#!/usr/bin/env python3
"""Plot STAR-Fusion heatmap and per-sample Circos plots.

Inputs are discovered under <results-dir>/fusion/<sample>/ and outputs are written to
<results-dir>/plots/fusion/.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path
from typing import Iterable

SKIP_DIRS = {"_starF_checkpoints", "tmp_chim_read_mappings_dir", "star-fusion.preliminary"}
PREDICTION_FILES = ["star-fusion.fusion_predictions.tsv", "star-fusion.fusion_predictions.abridged.tsv"]
LINK_HEADER = [
    "Sample", "FusionName", "LeftGene", "RightGene", "LeftBreakpoint", "RightBreakpoint",
    "left_chr", "left_pos", "right_chr", "right_pos", "JunctionReadCount", "SpanningFragCount", "FFPM",
]

HUMAN_CHR_SIZES = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555, "chr5": 181538259,
    "chr6": 170805979, "chr7": 159345973, "chr8": 145138636, "chr9": 138394717, "chr10": 133797422,
    "chr11": 135086622, "chr12": 133275309, "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345, "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}
MOUSE_CHR_SIZES = {
    "chr1": 195154279, "chr2": 181755017, "chr3": 159745316, "chr4": 156860686, "chr5": 151758149,
    "chr6": 149588044, "chr7": 144995196, "chr8": 130127694, "chr9": 124359700, "chr10": 130530862,
    "chr11": 121973369, "chr12": 120092757, "chr13": 120883175, "chr14": 125139656, "chr15": 104073951,
    "chr16": 98008968, "chr17": 95294699, "chr18": 90720763, "chr19": 61420004, "chrX": 169476592,
    "chrY": 91455967,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot STAR-Fusion FFPM heatmap and Circos fusion links.")
    parser.add_argument("--results-dir", required=True, help="RNA-seq results directory containing fusion/<sample>/ outputs.")
    parser.add_argument("--organism", default="human", help="human/hsa or mouse/mmu; controls chromosome order and sizes.")
    parser.add_argument("--min-ffpm", type=float, default=0.1, help="Minimum FFPM for retaining fusion events.")
    parser.add_argument("--min-junction-reads", type=float, default=1.0, help="Minimum JunctionReadCount for retaining fusion events.")
    parser.add_argument("--link-width", default="junction", choices=["junction", "ffpm"], help="Value used for Circos link width.")
    return parser.parse_args()


def chr_sizes_for_organism(organism: str) -> dict[str, int]:
    org = organism.lower()
    return MOUSE_CHR_SIZES if org in {"mouse", "mmu", "mus_musculus"} else HUMAN_CHR_SIZES


def canonical_chr(chrom: str) -> str:
    chrom = str(chrom).strip()
    if not chrom:
        return ""
    chrom = chrom.split()[0]
    return chrom if chrom.lower().startswith("chr") else f"chr{chrom}"


def parse_breakpoint(value: str) -> tuple[str, int]:
    parts = str(value).strip().split(":")
    if len(parts) < 2:
        return "", 0
    chrom = canonical_chr(parts[0])
    try:
        pos = int(float(parts[1]))
    except ValueError:
        pos = 0
    return chrom, pos


def to_float(value: object, default: float = 0.0) -> float:
    try:
        if value is None or value == "":
            return default
        out = float(value)
        return out if math.isfinite(out) else default
    except (TypeError, ValueError):
        return default


def find_prediction_file(sample_dir: Path) -> Path | None:
    for name in PREDICTION_FILES:
        path = sample_dir / name
        if path.is_file():
            return path
    return None


def sample_dirs(fusion_dir: Path) -> list[Path]:
    if not fusion_dir.exists():
        return []
    return sorted([p for p in fusion_dir.iterdir() if p.is_dir() and p.name not in SKIP_DIRS], key=lambda p: p.name)


def split_fusion_genes(fusion_name: str) -> tuple[str, str]:
    if "--" in fusion_name:
        left, right = fusion_name.split("--", 1)
        return left, right
    return fusion_name, ""


def read_fusions(sample: str, prediction_file: Path | None, min_ffpm: float, min_junction_reads: float) -> list[dict[str, object]]:
    if prediction_file is None:
        return []
    events: list[dict[str, object]] = []
    with prediction_file.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            fusion_name = row.get("FusionName", "")
            left_bp = row.get("LeftBreakpoint", "")
            right_bp = row.get("RightBreakpoint", "")
            left_chr, left_pos = parse_breakpoint(left_bp)
            right_chr, right_pos = parse_breakpoint(right_bp)
            ffpm = to_float(row.get("FFPM"))
            junction = to_float(row.get("JunctionReadCount"))
            spanning = to_float(row.get("SpanningFragCount"))
            if ffpm < min_ffpm and junction < min_junction_reads:
                continue
            if not left_chr or not right_chr or left_pos <= 0 or right_pos <= 0:
                continue
            left_gene = row.get("LeftGene") or split_fusion_genes(fusion_name)[0]
            right_gene = row.get("RightGene") or split_fusion_genes(fusion_name)[1]
            events.append({
                "Sample": sample, "FusionName": fusion_name, "LeftGene": left_gene, "RightGene": right_gene,
                "LeftBreakpoint": left_bp, "RightBreakpoint": right_bp, "left_chr": left_chr, "left_pos": left_pos,
                "right_chr": right_chr, "right_pos": right_pos, "JunctionReadCount": junction,
                "SpanningFragCount": spanning, "FFPM": ffpm,
            })
    return events


def write_links(path: Path, events: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LINK_HEADER, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for event in events:
            writer.writerow(event)


def write_heatmap_matrix(path: Path, all_events: list[dict[str, object]], samples: list[str], min_ffpm: float) -> list[str]:
    matrix: dict[str, dict[str, float]] = {}
    for event in all_events:
        fusion = str(event["FusionName"])
        sample = str(event["Sample"])
        matrix.setdefault(fusion, {})[sample] = max(matrix.setdefault(fusion, {}).get(sample, 0.0), float(event["FFPM"]))
    kept = sorted([fusion for fusion, values in matrix.items() if max(values.values() or [0.0]) > min_ffpm])
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Fusion", *samples])
        for fusion in kept:
            writer.writerow([fusion, *[matrix.get(fusion, {}).get(sample, 0.0) for sample in samples]])
    return kept


def plot_heatmap(matrix_path: Path, pdf_path: Path) -> None:
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        plot_heatmap_basic_pdf(matrix_path, pdf_path)
        return

    with matrix_path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        rows = list(reader)
    samples = header[1:]
    fusions = [row[0] for row in rows]
    values = np.array([[to_float(v) for v in row[1:]] for row in rows], dtype=float) if rows else np.zeros((0, len(samples)))
    height = max(4.0, min(18.0, 0.28 * max(1, len(fusions)) + 2.0))
    width = max(6.0, min(18.0, 0.45 * max(1, len(samples)) + 4.0))
    fig, ax = plt.subplots(figsize=(width, height))
    if values.size:
        im = ax.imshow(values, aspect="auto", cmap="Reds")
        fig.colorbar(im, ax=ax, label="FFPM")
    else:
        ax.text(0.5, 0.5, "No fusion with FFPM > threshold", ha="center", va="center", transform=ax.transAxes)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=45, ha="right")
    ax.set_yticks(range(len(fusions)))
    ax.set_yticklabels(fusions)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Fusion")
    ax.set_title("Fusion × Sample FFPM heatmap")
    fig.tight_layout()
    fig.savefig(pdf_path)
    plt.close(fig)


def pdf_escape(text: object) -> str:
    return str(text).replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def plot_heatmap_basic_pdf(matrix_path: Path, pdf_path: Path) -> None:
    """Dependency-free fallback heatmap PDF used when matplotlib is unavailable."""
    with matrix_path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        rows = list(reader)
    samples = header[1:]
    fusions = [row[0] for row in rows]
    values = [[to_float(v) for v in row[1:]] for row in rows]
    max_value = max([max(row) for row in values if row] or [1.0]) or 1.0
    width = 900
    height = max(360, min(1600, 120 + 20 * max(1, len(fusions))))
    left = 180
    top = height - 90
    cell_w = max(18, min(80, int((width - left - 60) / max(1, len(samples)))))
    cell_h = 16
    commands = ["1 1 1 rg 0 0 %d %d re f" % (width, height), "0 0 0 rg /F1 14 Tf 40 %d Td (Fusion x Sample FFPM heatmap) Tj" % (height - 35)]
    commands.append("/F1 8 Tf")
    for j, sample in enumerate(samples):
        x = left + j * cell_w + 2
        commands.append("0 0 0 rg 1 0 0 1 %d %d Tm (%s) Tj" % (x, top + 12, pdf_escape(sample[:16])))
    for i, fusion in enumerate(fusions[:80]):
        y = top - (i + 1) * cell_h
        commands.append("0 0 0 rg 1 0 0 1 40 %d Tm (%s) Tj" % (y + 4, pdf_escape(fusion[:35])))
        for j, value in enumerate(values[i][:len(samples)]):
            intensity = min(1.0, value / max_value)
            red = 1.0
            green = 1.0 - 0.75 * intensity
            blue = 1.0 - 0.75 * intensity
            x = left + j * cell_w
            commands.append("%.3f %.3f %.3f rg %d %d %d %d re f" % (red, green, blue, x, y, cell_w - 1, cell_h - 1))
    if not fusions:
        commands.append("0 0 0 rg /F1 12 Tf 1 0 0 1 260 %d Tm (No fusion with FFPM > threshold) Tj" % (height // 2))
    stream = "\n".join(commands).encode()
    objects = []
    objects.append(b"<< /Type /Catalog /Pages 2 0 R >>")
    objects.append(b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>")
    objects.append(f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {width} {height}] /Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>".encode())
    objects.append(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>")
    objects.append(b"<< /Length " + str(len(stream)).encode() + b" >>\nstream\n" + stream + b"\nendstream")
    out = bytearray(b"%PDF-1.4\n")
    offsets = [0]
    for idx, obj in enumerate(objects, 1):
        offsets.append(len(out))
        out.extend(f"{idx} 0 obj\n".encode() + obj + b"\nendobj\n")
    xref = len(out)
    out.extend(f"xref\n0 {len(objects)+1}\n0000000000 65535 f \n".encode())
    for offset in offsets[1:]:
        out.extend(f"{offset:010d} 00000 n \n".encode())
    out.extend(f"trailer << /Size {len(objects)+1} /Root 1 0 R >>\nstartxref\n{xref}\n%%EOF\n".encode())
    pdf_path.write_bytes(out)


def plot_circos(sample: str, events: list[dict[str, object]], chr_sizes: dict[str, int], out_pdf: Path, out_png: Path, link_width: str) -> None:
    if not events:
        print(f"[{sample}] No fusion events found for circos plot.")
        return
    try:
        from pycirclize import Circos
    except ImportError as exc:
        raise SystemExit("Please install pyCirclize: pip install pycirclize") from exc

    circos = Circos(chr_sizes, space=2)
    for sector in circos.sectors:
        track = sector.add_track((95, 100))
        track.axis(fc="#f5f5f5", ec="#555555", lw=0.5)
        track.text(sector.name.replace("chr", ""), size=7)
    for event in events:
        left_chr = str(event["left_chr"])
        right_chr = str(event["right_chr"])
        if left_chr not in chr_sizes or right_chr not in chr_sizes:
            continue
        left_pos = min(max(int(event["left_pos"]), 1), chr_sizes[left_chr])
        right_pos = min(max(int(event["right_pos"]), 1), chr_sizes[right_chr])
        weight = float(event["JunctionReadCount"] if link_width == "junction" else event["FFPM"])
        if weight <= 0:
            weight = float(event["FFPM"])
        lw = max(0.5, min(6.0, 0.5 + math.log1p(weight)))
        circos.link((left_chr, left_pos, min(left_pos + 1, chr_sizes[left_chr])), (right_chr, right_pos, min(right_pos + 1, chr_sizes[right_chr])), color="#c06c84aa", lw=lw)
    fig = circos.plotfig()
    fig.suptitle(f"{sample} Fusion Circos Plot", fontsize=12)
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=180)


def main() -> int:
    args = parse_args()
    results_dir = Path(args.results_dir)
    fusion_dir = results_dir / "fusion"
    plot_dir = results_dir / "plots" / "fusion"
    plot_dir.mkdir(parents=True, exist_ok=True)
    chr_sizes = chr_sizes_for_organism(args.organism)
    samples = [p.name for p in sample_dirs(fusion_dir)]
    all_events: list[dict[str, object]] = []
    for sample_dir in sample_dirs(fusion_dir):
        sample = sample_dir.name
        events = read_fusions(sample, find_prediction_file(sample_dir), args.min_ffpm, args.min_junction_reads)
        all_events.extend(events)
        write_links(plot_dir / f"{sample}.fusion_circos_links.tsv", events)
        plot_circos(sample, events, chr_sizes, plot_dir / f"{sample}.fusion_circos.pdf", plot_dir / f"{sample}.fusion_circos.png", args.link_width)
    matrix_path = plot_dir / "fusion_heatmap.tsv"
    write_heatmap_matrix(matrix_path, all_events, samples, args.min_ffpm)
    plot_heatmap(matrix_path, plot_dir / "fusion_heatmap.pdf")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

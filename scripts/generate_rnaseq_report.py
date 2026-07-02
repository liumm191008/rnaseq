#!/usr/bin/env python3
"""Generate an offline Chinese static HTML RNA-seq analysis report.

The report is intentionally dependency-free: it uses only the Python standard
library and embeds all CSS/JavaScript in a single HTML file. Data files remain
next to the analysis results and are linked by relative paths so the whole
results directory can be delivered to customers and browsed offline.
"""

from __future__ import annotations

import argparse
import base64
import csv
import html
import json
import os
import re
from pathlib import Path
from typing import Iterable, Sequence

IMAGE_EXTENSIONS = {".png", ".jpg", ".jpeg", ".gif", ".svg", ".pdf"}
TABLE_EXTENSIONS = {".csv", ".tsv", ".txt"}
DEFAULT_SMALL_TABLE_ROWS = 200
DEFAULT_BIG_TABLE_PREVIEW_ROWS = 30
REPORT_INTRO = "RNA测序（RNA sequencing，RNA-seq）是一种基于高通量测序技术的转录组研究方法，可对样本中的RNA进行全面检测和定量分析。通过对测序数据进行质量控制、序列比对、基因表达定量、差异表达分析及功能富集分析，可以系统揭示不同实验条件下基因表达变化及其潜在生物学机制。本报告展示了RNA-seq数据的主要分析结果，包括测序质量评估、表达定量、差异表达分析、功能富集分析等内容，为相关生物学研究和机制探索提供参考依据。"
STRING_ANALYSIS_TEXT = "STRING（Search Tool for the Retrieval of Interacting Genes/Proteins）是目前应用最广泛的蛋白互作数据库之一，整合了实验验证、数据库注释、文献挖掘、共表达分析及生物信息学预测等多种来源的蛋白相互作用信息。建议用户通过将差异表达基因、WGCNA关键模块基因或Hub Gene导入STRING数据库，构建蛋白-蛋白相互作用网络（Protein-Protein Interaction Network，PPI），从网络层面分析基因间的功能关联关系，识别网络中的关键节点和潜在核心调控因子。"
STRING_URL = "https://string-db.org/"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="生成可离线查阅的中文 RNA-seq 静态 HTML 分析报告。",
    )
    parser.add_argument("--input-json", required=True, help="WDL/Cromwell input JSON，例如 input.json。")
    parser.add_argument("--results-dir", required=True, help="RNA-seq 分析结果目录，例如 /path/rnaseq_results。")
    parser.add_argument("--report-dir", required=True, help="报告输出目录，建议放在 results-dir/report。")
    parser.add_argument("--project-name", default="RNA-seq 数据分析项目", help="项目名称。")
    parser.add_argument("--client-name", default="", help="客户名称，可选。")
    parser.add_argument("--report-title", default="RNA-seq 数据分析报告", help="报告标题。")
    parser.add_argument("--small-table-rows", type=int, default=DEFAULT_SMALL_TABLE_ROWS, help="小表格最多内嵌展示行数。")
    parser.add_argument("--big-table-preview-rows", type=int, default=DEFAULT_BIG_TABLE_PREVIEW_ROWS, help="大文件预览行数。")
    return parser.parse_args()


def load_input_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def esc(value: object) -> str:
    return html.escape("" if value is None else str(value), quote=True)


def rel_link(target: Path, report_dir: Path) -> str:
    return Path(os.path.relpath(target, report_dir)).as_posix()


def file_size(path: Path) -> str:
    size = path.stat().st_size
    units = ["B", "KB", "MB", "GB", "TB"]
    value = float(size)
    for unit in units:
        if value < 1024 or unit == units[-1]:
            return f"{value:.1f} {unit}" if unit != "B" else f"{int(value)} {unit}"
        value /= 1024
    return f"{size} B"


def find_existing(paths: Iterable[Path]) -> list[Path]:
    return [path for path in paths if path.exists()]


def collect_files(root: Path, extensions: set[str] | None = None) -> list[Path]:
    if not root.exists():
        return []
    files: list[Path] = []
    for path in root.rglob("*"):
        if path.is_file() and (extensions is None or path.suffix.lower() in extensions):
            files.append(path)
    return sorted(files)


def detect_delimiter(path: Path, sample_line: str | None = None) -> str:
    if path.suffix.lower() == ".csv":
        return ","
    if sample_line and "," in sample_line and "\t" not in sample_line:
        return ","
    return "\t"


def read_table(path: Path, small_limit: int, preview_limit: int) -> tuple[list[str], list[list[str]], bool, int | None]:
    """Read a CSV/TSV-like table.

    Returns (header, rows, truncated, total_rows). For big files, only preview
    rows are returned and total_rows may be None if the file is too large to
    count cheaply during generation.
    """
    raw_lines: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            raw_lines.append(line.rstrip("\n\r"))
            if len(raw_lines) > small_limit:
                break

    if not raw_lines:
        return ["提示"], [["文件为空或仅包含注释行"]], False, 0

    delimiter = detect_delimiter(path, raw_lines[0])
    parsed = list(csv.reader(raw_lines, delimiter=delimiter))
    header = parsed[0]
    rows = parsed[1:]
    truncated = len(raw_lines) > small_limit
    if truncated:
        rows = rows[:preview_limit]
        total_rows = None
    else:
        total_rows = len(rows)
    return header, rows, truncated, total_rows


def render_download(path: Path, report_dir: Path, label: str = "下载原始数据") -> str:
    link_text = f"{label}：{path.name}"
    return f'<a class="download" href="{esc(rel_link(path, report_dir))}" download>{esc(link_text)}（{esc(file_size(path))}）</a>'


def render_file_link(path: Path, report_dir: Path, text: str | None = None, class_name: str = "cell-link") -> str:
    label = text or path.name
    return f'<a class="{esc(class_name)}" href="{esc(rel_link(path, report_dir))}" download title="{esc(label)}">{esc(label)}</a>'


def basename_or_value(value: object) -> str:
    if not value or value == "-":
        return "-"
    return Path(str(value)).name


def cell_width_style(width: int | None) -> str:
    if width is None:
        return ""
    return f' style="width:{width}ch;min-width:{width}ch;max-width:{width}ch"'


def render_cell(value: object, tag: str = "td", href: str = "", width: int | None = None) -> str:
    text = esc(value)
    if href and tag == "td":
        content = f'<a href="{esc(href)}" target="_blank" rel="noopener" class="cell-link"><span class="cell-text">{text}</span></a>'
    else:
        content = f'<span class="cell-text">{text}</span>'
    return f'<{tag}{cell_width_style(width)} data-full-text="{text}">{content}</{tag}>'


def one_line_text(text: str, class_name: str = "muted") -> str:
    return f'<p class="{class_name} text-ellipsis" title="{esc(text)}">{esc(text)}</p>'


def describe_table_file(path: Path) -> str:
    name = path.name
    lower = name.lower()
    if name == "gene_counts.tsv":
        return "featureCounts 基因表达定量矩阵：每一行代表一个基因；Geneid 为基因 ID，Chr/Start/End/Strand/Length 为注释坐标与基因长度，后续样本列为该基因分配到的 paired-end fragment/read count，可用于差异表达、富集和共表达分析。"
    if name.endswith(".summary"):
        return "featureCounts 计数汇总表：每一行代表一种计数状态；Assigned 为成功分配到基因的 reads/fragments，Unassigned_* 表示因多重比对、无特征、嵌合、低质量或歧义等原因未计数的 reads/fragments；每一列对应一个样本 BAM。"
    if name == "pairwise_de_summary.tsv":
        return "差异表达比较汇总表：每一行代表一个组间比较；comparison 为比较名称，group_a/group_b 为比较双方，total_tested 为参与 DESeq2 检验的基因数，deg_total 为显著差异基因总数，deg_up/deg_down 分别为相对 group_a 上调/下调的基因数。"
    if name == "filtered_count_matrix.tsv":
        return "过滤后表达矩阵：每一行代表通过低表达过滤的基因，每一列代表样本表达计数；过滤标准为所有样本 count 总和 rowSums(count_data) >= min_count（流程默认 min_count=10，可由输入参数调整），该矩阵是 DESeq2 差异表达建模使用的输入矩阵。"
    if name == "vst_count_matrix.tsv":
        return "VST 矩阵：每一行代表通过过滤并完成 DESeq2 variance stabilizing transformation 的基因，每一列代表样本；该文件用于 PCA、样本距离和 WGCNA，不在报告正文展示。"
    if name == "normalized_count_matrix.tsv":
        return "DESeq2 标准化计数矩阵：每一行代表一个通过低表达过滤的基因，每一列代表一个样本；数值为使用 DESeq2 size factor 校正文库测序深度后的计数，适合样本间表达比较和差异表达建模。"
    if name == "fpkm_count_matrix.tsv":
        return "标准化表达矩阵：正文展示 FPKM 数据，每一行代表一个基因，每一列代表一个样本；数值同时按基因长度和样本文库大小校正，适合查看同一样本中的相对表达水平。下方同时提供 FPKM、TPM 和 DESeq2 标准化矩阵下载。"
    if name == "tpm_count_matrix.tsv":
        return "TPM 表达矩阵：每一行代表一个基因，每一列代表一个样本；数值先按基因长度校正，再将每个样本缩放到总和一百万，便于比较不同样本中的相对表达组成。"
    if lower.endswith(".deg.tsv"):
        return "显著差异表达基因表：每一行代表一个通过阈值筛选的差异基因；gene_symbol 为基因符号，gene_id 为基因编号，baseMean 为平均标准化表达量，log2FoldChange 为处理组相对对照组的 log2 倍数变化，lfcSE/stat 为模型统计量，pvalue 为原始显著性，padj 为多重检验校正后的 FDR，末尾 FPKM_<样本名> 列为各样本的 FPKM 表达量；下方提供完整差异表达表下载。"
    if lower.endswith(".all.tsv"):
        return "完整 DESeq2 差异表达结果表：每一行代表一个参与检验的基因；首列 gene_symbol 为基因符号，其后 gene_id 为基因编号，中间列为表达变化倍数及显著性指标，末尾 FPKM_<样本名> 列为各样本的 FPKM 表达量。"
    if lower.endswith("_go.csv"):
        return "GO 富集结果表：每一行代表一个 GO Biological Process 功能条目；ID 为 GO 编号，Description 为功能描述，Count/Overlap 为命中差异基因数，BgRatio/Background 为背景比例，pvalue/p.adjust/qvalue 为富集显著性，RichFactor 为命中基因数占该条目背景基因数的比例，geneID/SYMBOL 展示贡献该条目的基因。"
    if lower.endswith("_kegg.csv"):
        return "KEGG 富集结果表：每一行代表一条 KEGG pathway；ID 为通路编号且可点击跳转到 Pathway_Link，Description 为通路名称，Class/Subclass 为通路分类，Count/geneID 为富集到的差异基因，pvalue/p.adjust/qvalue 为显著性，RichFactor 表示差异基因在该通路背景中的占比。"
    if lower.endswith("_reactome.csv"):
        return "Reactome 富集结果表：每一行代表一个 Reactome pathway；ID 为通路编号，Description 为通路描述，Count/geneID 为命中基因，pvalue/p.adjust/qvalue 为富集显著性，RichFactor 为命中比例，用于判断差异基因集中涉及的反应和信号通路。"
    if "rmats" in lower or lower.endswith((".mats.jc.txt", ".mats.jcec.txt")):
        return "rMATS 可变剪切事件表：每一行代表一个剪切事件；ID/GeneID/geneSymbol 为事件和基因信息，chr/strand/exonStart_0base/exonEnd/upstreamES/downstreamEE 等为基因组坐标，IJC/SJC 为 inclusion/skipping junction counts，IncLevel 为各样本包含比例，IncLevelDifference 为两组剪切比例差，PValue/FDR 为显著性。"
    if "fusion" in lower:
        return "融合基因候选结果表：每一行代表一个候选融合事件，列标题展示融合基因、断点、支持 reads 和注释信息。"
    if name == "gene_modules.tsv":
        return "基因模块对应表：每一行代表一个基因，列标题展示基因 ID、模块颜色/编号等共表达模块归属信息，用于查看基因被分配到哪个 WGCNA 模块。"
    if name == "module_eigengenes.tsv":
        return "WGCNA 模块特征向量表：每一行代表一个样本，列标题为不同模块的 eigengene，用于比较模块表达模式。"
    if name == "important_modules.tsv":
        return "关键共表达模块列表：每一行代表一个与实验分组或表型显著相关的共表达模块；列标题通常包含模块名称、相关表型、相关系数、P 值或校正后显著性，用于筛选与实验设计关系最强、值得重点解读的模块。"
    if name == "hub_genes.tsv":
        return "WGCNA 核心调控基因表：每一行代表一个模块内候选 hub gene；列标题通常包含基因 ID、模块、模块成员关系 kME、连接度或排序信息，用于定位模块内代表性/关键调控基因。"
    return "分析结果表：每一行代表一个结果记录，列标题解释对应统计字段；可在页面内搜索、筛选、分页查看，完整文件可下载。"


def describe_image_file(path: Path) -> str:
    lower = path.name.lower()
    if "pca" in lower:
        return "PCA 主成分分析图：横轴为 PC1，纵轴为 PC2，括号中的百分比表示该主成分解释的表达变异比例；每个点代表一个样本，点旁文字为样本名，颜色/形状通常代表分组，用于观察组内重复和组间分离。"
    if "sample_distance" in lower:
        return "样本距离热图：横轴和纵轴均为样本，颜色表示样本间表达距离或相似度；距离越小表示表达谱越相近，聚类分支用于判断样本整体一致性。"
    if "heatmap" in lower and "wgcna" not in lower:
        return "表达热图：横轴通常为样本，纵轴为基因或功能条目，颜色表示标准化表达量或统计值高低，用于查看样本聚类和基因表达模式。"
    if "volcano" in lower:
        return "火山图：横轴为 log2FoldChange，表示处理组相对对照组的表达变化方向和倍数；纵轴通常为 -log10(padj 或 pvalue)，越高表示差异越显著，点代表基因。"
    if lower.endswith(".ma.pdf") or ".ma." in lower:
        return "MA 图：横轴为基因平均表达量或标准化平均 count，纵轴为 log2FoldChange；每个点代表一个基因，用于观察不同表达丰度下的差异变化。"
    if "go" in lower and "dag" in lower:
        return "GO DAG 有向无环图：节点代表 GO 功能条目，连线表示 GO 术语之间的层级/从属关系；颜色通常表示富集显著性或富集程度，用于查看显著 GO 条目在功能层级中的位置和上下游关系。"
    if "go" in lower and "dotplot" in lower:
        return "GO 富集气泡图：横轴为 RichFactor 或富集比例，纵轴为 GO 功能条目；气泡大小通常表示富集基因数，颜色表示校正后显著性。"
    if "kegg" in lower and "dotplot" in lower:
        return "KEGG 富集气泡图：横轴为 RichFactor 或富集比例，纵轴为 KEGG 通路名称；气泡大小表示通路中的差异基因数，颜色表示显著性。"
    if "barplot" in lower and "go" in lower:
        return "GO 分类/富集柱状图：横轴通常为基因数量或条目数量，纵轴为 GO 分类或功能条目，用于概览差异基因涉及的功能类别。"
    if "barplot" in lower and "kegg" in lower:
        return "KEGG 富集柱状图：横轴通常为基因数量、RichFactor 或 -log10 显著性，纵轴为通路名称，用于展示主要富集通路。"
    if "reactome" in lower and "dotplot" in lower:
        return "Reactome 富集气泡图：横轴通常为 RichFactor 或富集比例，纵轴为 Reactome 通路名称；气泡大小表示富集基因数，颜色表示显著性。"
    if "reactome" in lower and "barplot" in lower:
        return "Reactome 富集柱状图：横轴通常为基因数量、RichFactor 或 -log10 显著性，纵轴为 Reactome 通路名称，用于展示主要富集通路。"
    if "fusion_heatmap" in lower:
        return "Fusion Heatmap：横轴为样本，纵轴为融合基因事件，颜色表示 FFPM 值；仅保留至少一个样本 FFPM > 0.1 的 Fusion，用于查看融合事件在不同样本中的表达支持强度。"
    if "fusion_circos" in lower:
        return "Circos plot shows genomic connections between fusion gene breakpoints. Each arc represents one fusion event. Thicker arcs indicate stronger fusion support based on junction read count or FFPM."
    if "wgcna_module_trait_heatmap" in lower:
        return "WGCNA 模块-表型相关性热图：横轴为实验分组或表型变量，纵轴为共表达模块；每个色块表示模块特征向量与表型的相关系数，颜色深浅和方向表示相关强弱及正负关系，色块中的数字通常为相关系数和显著性 P 值，用于筛选与实验分组最相关的模块。"
    if "wgcna_sample_clustering" in lower:
        return "WGCNA 样本聚类图：横轴为样本名称或样本聚类树叶节点，纵轴为样本间表达差异对应的聚类高度；分支越近表示样本整体表达模式越相似，图中颜色注释通常表示分组或其他表型，用于发现离群样本和检查分组一致性。"
    if "wgcna_module_size_barplot" in lower:
        return "WGCNA 模块大小柱状图：横轴为共表达模块名称或模块颜色，纵轴为该模块包含的基因数量；柱子越高表示模块规模越大，用于了解模块划分结构，但模块大小不等同于生物学重要性，需要结合模块-表型相关性和 hub gene 一起解读。"
    if "wgcna_module_eigengene_heatmap" in lower:
        return "WGCNA 模块特征向量热图：横轴为样本，纵轴为共表达模块 eigengene；颜色表示模块整体表达趋势在不同样本中的高低，用于比较模块在样本或分组间的表达模式。"
    if "wgcna_gene_dendrogram" in lower:
        return "WGCNA 基因聚类树：横轴为参与网络分析的基因，纵轴为基于拓扑重叠或表达相似性计算的聚类高度；下方颜色条通常代表模块归属，用于展示基因如何被切分到不同共表达模块。"
    if "wgcna_soft_threshold" in lower:
        return "WGCNA 软阈值选择图：横轴为候选 soft-threshold power，纵轴通常展示无尺度网络拟合度和平均连接度；用于判断构建共表达网络时选用的软阈值是否兼顾无尺度特征和网络连接性。"
    if "wgcna" in lower or "dendrogram" in lower:
        return "WGCNA 共表达网络图：横轴、纵轴和颜色含义以图内标签为准，通常用于展示基因共表达模块、模块间关系或模块与样本表型的联系；解读时应结合样本量、分组信息和 hub gene 表格。"
    return "分析图表：横轴和纵轴含义见图内坐标轴标题或图例；图中点、线、颜色或分面表示对应分析对象及统计值，可下载离线查看。"


def kegg_link_for_cell(path: Path, header: str, row_map: dict[str, str]) -> str:
    if path.name.lower().endswith("_kegg.csv") and header == "ID":
        return row_map.get("Pathway_Link", "")
    return ""


def chinese_table_title(path: Path) -> str:
    name = path.name
    lower = name.lower()
    comparison = de_comparison_name(path)
    if name == "gene_counts.tsv":
        return "基因表达定量矩阵"
    if name.endswith(".summary"):
        return "基因计数汇总统计"
    if name == "pairwise_de_summary.tsv":
        return "差异表达比较汇总"
    if name == "filtered_count_matrix.tsv":
        return "过滤后表达矩阵"
    if name == "vst_count_matrix.tsv":
        return "VST 矩阵"
    if name == "normalized_count_matrix.tsv":
        return "DESeq2 标准化计数矩阵"
    if name == "fpkm_count_matrix.tsv":
        return "标准化表达矩阵"
    if name == "tpm_count_matrix.tsv":
        return "TPM 表达矩阵"
    if name == "sample_metadata.tsv":
        return "样本分组信息"
    if comparison and lower.endswith(".deg.tsv"):
        return f"{comparison} 显著差异基因表"
    if comparison and lower.endswith(".all.tsv"):
        return f"{comparison} 完整差异表达表"
    enrich = enrichment_comparison_name(path)
    if enrich:
        if lower.endswith("_go.csv"):
            return f"GO 富集结果（{enrich} ）"
        if lower.endswith("_go_classification.csv"):
            return f"GO 二级分类（{enrich} ）"
        if lower.endswith("_kegg.csv"):
            return f"KEGG 富集结果（{enrich} ）"
        if lower.endswith("_reactome.csv"):
            return f"Reactome 富集结果（{enrich} ）"
    if "rmats" in lower or lower.endswith((".mats.jc.txt", ".mats.jcec.txt")):
        return "rMATS 可变剪切结果"
    if "fusion" in lower:
        if name == "fusion_heatmap.tsv":
            return "Fusion × Sample FFPM矩阵"
        if name.endswith(".fusion_circos_links.tsv"):
            return f"{name.replace('.fusion_circos_links.tsv', '')} Fusion Circos绘图数据"
        return "融合基因候选结果"
    if name == "gene_modules.tsv":
        return "基因模块对应表"
    if name == "module_eigengenes.tsv":
        return "WGCNA 模块特征向量表"
    if name == "important_modules.tsv":
        return "关键共表达模块列表"
    if name == "hub_genes.tsv":
        return "核心调控基因（Hub Gene）"
    return "分析结果表"


def display_text_for_width(value: object) -> str:
    return re.sub(r"<[^>]+>", "", str(value))


def column_widths(rows: Sequence[Sequence[object]], column_count: int) -> list[int]:
    widths = []
    for idx in range(column_count):
        values = [
            display_text_for_width(row[idx])
            for row in rows
            if idx < len(row) and row[idx] is not None
        ]
        max_len = max((len(value) for value in values), default=8)
        widths.append(max(8, min(max_len, 42)))
    return widths


def apply_table_width_overrides(path: Path, header: Sequence[str], widths: list[int]) -> list[int]:
    adjusted = list(widths)
    if path.name == "gene_counts.tsv":
        compact_columns = {"chr", "start", "end", "strand"}
        for idx, column in enumerate(header):
            if column.strip().lower() in compact_columns and idx < len(adjusted):
                adjusted[idx] = 8
    return adjusted


def render_colgroup(widths: Sequence[int]) -> str:
    return "<colgroup>" + "".join(f'<col style="width:{width}ch;min-width:{width}ch;max-width:{width}ch">' for width in widths) + "</colgroup>"


def priority_de_plot(path: Path) -> int:
    name = path.name.lower()
    if name == "pca.pdf":
        return 0
    if name == "sample_distance_heatmap.pdf":
        return 1
    return 2


def enrichment_comparison_name(path: Path) -> str | None:
    name = path.name
    patterns = [
        r"(.+)_GO\.csv$",
        r"(.+)_GO_classification\.csv$",
        r"(.+)_KEGG\.csv$",
        r"(.+)_Reactome\.csv$",
        r"(.+)_KEGG\.html$",
        r"(.+)_Reactome\.html$",
        r"(.+)_(?:GO|KEGG|Reactome)_(?:dotplot|barplot|DAG)\.pdf$",
    ]
    for pattern in patterns:
        match = re.match(pattern, name)
        if match:
            return match.group(1)
    return None


def render_table(
    title: str,
    path: Path,
    report_dir: Path,
    small_limit: int,
    preview_limit: int,
    description: str = "",
    extra_downloads: Sequence[tuple[Path, str]] | None = None,
    primary_download_label: str | None = None,
) -> str:
    header, rows, truncated, total_rows = read_table(path, small_limit, preview_limit)
    note = ""
    if truncated:
        note = f"<p class='note'>该文件较大，仅展示前 {preview_limit} 行示例；请点击下载查看完整文件。</p>"
    elif total_rows is not None:
        note = f"<p class='note'>共展示 {total_rows} 行，可搜索、筛选和分页。</p>"

    display_title = title or chinese_table_title(path)
    widths = apply_table_width_overrides(path, header, column_widths(rows, len(header)))
    colgroup = render_colgroup(widths)
    thead = "".join(render_cell(col, "th", width=width) for col, width in zip(header, widths))
    tbody_rows = []
    for row in rows:
        cells = list(row) + [""] * max(0, len(header) - len(row))
        row_values = cells[: len(header)]
        row_map = dict(zip(header, row_values))
        rendered_cells = [
            render_cell(cell, href=kegg_link_for_cell(path, col, row_map), width=width)
            for col, cell, width in zip(header, row_values, widths)
        ]
        tbody_rows.append("<tr>" + "".join(rendered_cells) + "</tr>")
    tbody = "\n".join(tbody_rows)
    primary_label = primary_download_label or ("下载显著差异表达表" if path.name.lower().endswith(".deg.tsv") else "下载原始数据")
    downloads = [render_download(path, report_dir, primary_label)]
    for extra_path, label in extra_downloads or []:
        downloads.append(render_download(extra_path, report_dir, label))
    return f"""
    <article class="card table-card">
      <div class="card-header">
        <div>
          <h4>{esc(display_title)}</h4>
          {one_line_text(description or describe_table_file(path))}
          <div class="download-group">{"".join(downloads)}</div>
        </div>
      </div>
      {note}
      <div class="table-tools">
        <label>搜索/筛选：<input type="search" class="table-search" placeholder="输入关键词"></label>
        <label>每页：<select class="page-size"><option selected>10</option><option>20</option><option>50</option><option>100</option></select></label>
        <span class="page-info"></span>
        <button type="button" class="prev-page">上一页</button>
        <button type="button" class="next-page">下一页</button>
      </div>
      <div class="table-wrap">
        <table class="report-table">
          {colgroup}
          <thead><tr>{thead}</tr></thead>
          <tbody>{tbody}</tbody>
        </table>
      </div>
    </article>
    """


def render_simple_table_cell(cell: object, width: int | None = None) -> str:
    cell_text = str(cell)
    if cell_text.startswith("<a "):
        plain_text = re.sub(r"<[^>]+>", "", cell_text)
        return f'<td{cell_width_style(width)} data-full-text="{esc(plain_text)}"><span class="cell-text">{cell_text}</span></td>'
    return render_cell(cell, width=width)


def render_simple_table(title: str, headers: Sequence[str], rows: Sequence[Sequence[object]], description: str = "") -> str:
    widths = column_widths(rows, len(headers))
    colgroup = render_colgroup(widths)
    thead = "".join(render_cell(col, "th", width=width) for col, width in zip(headers, widths))
    tbody = "\n".join(
        "<tr>" + "".join(render_simple_table_cell(cell, width) for cell, width in zip(row, widths)) + "</tr>" for row in rows
    )
    return f"""
    <article class="card table-card">
      <h4>{esc(title)}</h4>
      {one_line_text(description or "表格支持搜索、筛选和分页；长字段悬停可查看完整内容。")}
      <div class="table-tools">
        <label>搜索/筛选：<input type="search" class="table-search" placeholder="输入关键词"></label>
        <label>每页：<select class="page-size"><option selected>10</option><option>20</option><option>50</option><option>100</option></select></label>
        <span class="page-info"></span>
        <button type="button" class="prev-page">上一页</button>
        <button type="button" class="next-page">下一页</button>
      </div>
      <div class="table-wrap"><table class="report-table">{colgroup}<thead><tr>{thead}</tr></thead><tbody>{tbody}</tbody></table></div>
    </article>
    """


def render_image(path: Path, report_dir: Path) -> str:
    href = esc(rel_link(path, report_dir))
    suffix = path.suffix.lower()
    caption = esc(path.name)
    if suffix == ".pdf":
        pdf_href = href + "#toolbar=0&navpanes=0&scrollbar=1"
        body = f'<object data="{pdf_href}" type="application/pdf"><p>浏览器无法内嵌 PDF，<a href="{href}">点击打开 {caption}</a></p></object>'
    else:
        body = f'<img src="{href}" alt="{caption}" loading="lazy">'
    description = describe_image_file(path)
    return f'<figure class="figure">{body}<figcaption title="{esc(path.name)}">{caption} · <a href="{href}" download>下载</a></figcaption>{one_line_text(description, "figure-desc")}</figure>'


def render_gallery(title: str, files: Sequence[Path], report_dir: Path, description: str = "") -> str:
    if not files:
        return render_missing(title, "未发现图片或 PDF 图表。")
    figures = "".join(f'<div class="gallery-item">{render_image(path, report_dir)}</div>' for path in files)
    return f"""
    <article class="card gallery-card">
      <h4>{esc(title)}</h4>
      {one_line_text(description or "图表按两列平铺展示；图片和 PDF 均使用相对路径，便于离线查阅。")}
      <div class="gallery">{figures}</div>
    </article>
    """


def render_embedded_static_image(filename: str, title: str, description: str) -> str:
    path = Path(__file__).with_name(filename)
    if not path.exists():
        return render_missing(title, f"未发现内嵌图片资源 {filename}。")
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    src = f"data:image/png;base64,{encoded}"
    caption = esc(filename)
    return f"""
    <figure class="figure embedded-static-figure">
      <img src="{src}" alt="{esc(title)}" loading="lazy">
      <figcaption title="{caption}">{caption}</figcaption>
      {one_line_text(description, "figure-desc")}
    </figure>
    """


def render_ppi_section() -> str:
    figures = "".join([
        '<div class="gallery-item">' + render_embedded_static_image(
            "string.input.png",
            "STRING 操作方法",
            "STRING操作方法：录入基因列表和物种类型，进行在线快速分析",
        ) + '</div>',
        '<div class="gallery-item">' + render_embedded_static_image(
            "string.demo.png",
            "PPI 网络示例图",
            "PPI网络示例图：节点（Node）代表蛋白或基因，连线（Edge）代表蛋白间已知或预测的相互作用。连接度较高的节点通常被认为是网络中的核心基因（Hub Gene），可能在相关生物学过程或疾病机制中发挥重要调控作用。",
        ) + '</div>',
    ])
    return f"""
    <article class="card ppi-section">
      <h4>在线STRING分析（用户自操作）</h4>
      <p class="note">{esc(STRING_ANALYSIS_TEXT)}</p>
      <p><a class="cell-link" href="{esc(STRING_URL)}" target="_blank" rel="noopener">在线STRING分析: {esc(STRING_URL)}</a></p>
      <div class="gallery">{figures}</div>
    </article>
    """


def render_file_links(title: str, files: Sequence[Path], report_dir: Path, description: str = "") -> str:
    if not files:
        return render_missing(title, "未发现相关文件。")
    items = "".join(
        f'<li><a href="{esc(rel_link(path, report_dir))}">{esc(path.name)}</a> <span class="muted">{esc(file_size(path))}</span></li>'
        for path in files
    )
    return f"""
    <article class="card">
      <h4>{esc(title)}</h4>
      {one_line_text(description or "相关结果文件可点击打开或下载，用于离线交付。")}
      <ul class="file-list">{items}</ul>
    </article>
    """


def file_by_name(files: Sequence[Path], name: str) -> Path | None:
    for path in files:
        if path.name == name:
            return path
    return None


def render_single_figure_card(title: str, path: Path | None, report_dir: Path, description: str = "") -> str:
    if path is None:
        return render_missing(title, "未发现指定图表。")
    extra_class = " module-size-figure-card" if path.name == "wgcna_module_size_barplot.pdf" else ""
    return f"""
    <article class="card gallery-card single-figure-card{extra_class}">
      <h4>{esc(title)}</h4>
      {one_line_text(description or describe_image_file(path))}
      <div class="single-figure">{render_image(path, report_dir)}</div>
    </article>
    """


def fastqc_sample_name(path: Path, fastqc_root: Path) -> str:
    try:
        relative = path.relative_to(fastqc_root)
    except ValueError:
        return path.parent.name
    if len(relative.parts) > 1:
        return relative.parts[0]
    return path.parent.name


def fastqc_read_label(name: str) -> str:
    lower = name.lower()
    if re.search(r"(?:^|[._-])r?1(?:[._-]|$)", lower) or "_val_1" in lower:
        return "Read1"
    if re.search(r"(?:^|[._-])r?2(?:[._-]|$)", lower) or "_val_2" in lower:
        return "Read2"
    return "FASTQ"


def collect_fastqc_entries(fastqc_root: Path) -> dict[str, list[dict[str, Path | str]]]:
    if not fastqc_root.exists():
        return {}
    grouped: dict[str, list[dict[str, Path | str]]] = {}
    for extracted_dir in sorted(path for path in fastqc_root.rglob("*_fastqc") if path.is_dir()):
        sample = fastqc_sample_name(extracted_dir, fastqc_root)
        if ".trimmed" in sample.lower() or ".trimmed" in extracted_dir.name.lower():
            continue
        images_dir = extracted_dir / "Images"
        quality_plot = images_dir / "per_base_quality.png"
        content_plot = images_dir / "per_base_sequence_content.png"
        if not (quality_plot.exists() or content_plot.exists()):
            continue
        html_path = extracted_dir.parent / f"{extracted_dir.name}.html"
        entry: dict[str, Path | str] = {
            "name": extracted_dir.name,
            "read_label": fastqc_read_label(extracted_dir.name),
            "quality_plot": quality_plot,
            "content_plot": content_plot,
        }
        if html_path.exists():
            entry["html"] = html_path
        grouped.setdefault(sample, []).append(entry)
    for entries in grouped.values():
        entries.sort(key=lambda item: (str(item["read_label"]), str(item["name"])))
    return dict(sorted(grouped.items()))


def parse_fastqc_basic_stats(extracted_dir: Path) -> dict[str, str]:
    data_path = extracted_dir / "fastqc_data.txt"
    if not data_path.exists():
        return {}
    stats: dict[str, str] = {}
    in_basic_stats = False
    for line in data_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">>Basic Statistics"):
            in_basic_stats = True
            continue
        if in_basic_stats and line.startswith(">>END_MODULE"):
            break
        if not in_basic_stats or not line or line.startswith("#"):
            continue
        if "\t" in line:
            key, value = line.split("\t", 1)
            stats[key.strip()] = value.strip()
    return stats


def collect_fastqc_stats(results_dir: Path) -> dict[str, list[dict[str, str]]]:
    fastqc_root = results_dir / "fastqc"
    if not fastqc_root.exists():
        return {}
    grouped: dict[str, list[dict[str, str]]] = {}
    for extracted_dir in sorted(path for path in fastqc_root.rglob("*_fastqc") if path.is_dir()):
        sample = fastqc_sample_name(extracted_dir, fastqc_root)
        if ".trimmed" in sample.lower() or ".trimmed" in extracted_dir.name.lower():
            continue
        stats = parse_fastqc_basic_stats(extracted_dir)
        if not stats:
            continue
        entry = {
            "read_label": fastqc_read_label(extracted_dir.name),
            "name": extracted_dir.name,
            **stats,
        }
        grouped.setdefault(sample, []).append(entry)
    for entries in grouped.values():
        entries.sort(key=lambda item: (item.get("read_label", ""), item.get("name", "")))
    return grouped


def format_fastqc_stat(entries: Sequence[dict[str, str]], metric: str) -> str:
    values = []
    for entry in entries:
        label = entry.get("read_label", "FASTQ")
        value = entry.get(metric, "-") or "-"
        values.append(f"{label}: {value}")
    return "; ".join(values) if values else "-"


def render_fastqc_plot(path: Path, report_dir: Path, title: str, description: str) -> str:
    href = esc(rel_link(path, report_dir))
    return f"""
    <figure class="figure fastqc-figure">
      <img src="{href}" alt="{esc(title)}" loading="lazy">
      <figcaption title="{esc(path.name)}">{esc(title)} · <a href="{href}" download>下载</a></figcaption>
      {one_line_text(description, "figure-desc")}
    </figure>
    """


def render_fastqc_read(entry: dict[str, Path | str], report_dir: Path) -> str:
    name = str(entry["name"])
    read_label = str(entry["read_label"])
    html_link = ""
    html_path = entry.get("html")
    if isinstance(html_path, Path):
        html_link = render_download(html_path, report_dir, "下载 FastQC HTML")
    plots = []
    quality_plot = entry.get("quality_plot")
    if isinstance(quality_plot, Path) and quality_plot.exists():
        plots.append(render_fastqc_plot(
            quality_plot,
            report_dir,
            "Per base sequence quality",
            "Per base sequence quality：横轴为 reads 从 5' 到 3' 的碱基位置，纵轴为 Phred 质量分数；箱线图展示每个位置的质量分布，用于判断测序质量是否随读长下降。",
        ))
    content_plot = entry.get("content_plot")
    if isinstance(content_plot, Path) and content_plot.exists():
        plots.append(render_fastqc_plot(
            content_plot,
            report_dir,
            "Per base sequence content",
            "Per base sequence content：横轴为 reads 碱基位置，纵轴为 A/C/G/T 各碱基百分比；不同颜色曲线表示碱基组成，用于识别建库偏倚或接头/低复杂度影响。",
        ))
    return f"""
    <div class="fastqc-read-card">
      <div class="fastqc-read-header">
        <strong>{esc(read_label)}</strong>
        <span class="muted text-ellipsis" title="{esc(name)}">{esc(name)}</span>
      </div>
      <div class="download-group">{html_link}</div>
      <div class="fastqc-plot-grid">{''.join(plots)}</div>
    </div>
    """


def render_fastqc_section(results_dir: Path, report_dir: Path) -> str:
    grouped = collect_fastqc_entries(results_dir / "fastqc")
    if not grouped:
        return render_missing("FastQC 质控结果", "未发现已解压的 FastQC 结果目录（*_fastqc/Images）。请先运行 FastQC 并解压 *_fastqc.zip。")
    panels = []
    for sample, entries in grouped.items():
        read_cards = "".join(render_fastqc_read(entry, report_dir) for entry in entries)
        panel = f"""
        <div class="fastqc-sample-panel">
          {one_line_text("每个样本仅展示原始 FASTQ 的 FastQC Per base sequence quality 和 Per base sequence content；同一样本的 Read1/Read2 质控结果并排显示，并提供 FastQC HTML 下载；trim 后 FastQC 结果不在报告中展示。")}
          <div class="fastqc-read-grid">{read_cards}</div>
        </div>
        """
        panels.append((sample, panel))
    return render_tabs(
        "FastQC 质控结果（按样本）",
        "点击标签页切换不同样本；本节仅展示原始数据 FastQC 解压结果中的关键质量图和 HTML 下载。",
        panels,
    )


def render_missing(title: str, message: str) -> str:
    return f"<article class='card missing'><h4>{esc(title)}</h4><p>{esc(message)}</p></article>"


def render_method_note(title: str, software: str, notes: str) -> str:
    text = f"使用软件：{software}；结果解读注意事项：{notes}"
    return f"""
    <article class="card method-card">
      <h4>{esc(title)}</h4>
      {one_line_text(text)}
    </article>
    """


def parse_star_log(path: Path) -> dict[str, str]:
    metrics: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "|" not in line:
            continue
        key, value = [part.strip() for part in line.split("|", 1)]
        if key and value:
            metrics[key] = value
    return metrics


def parse_trim_galore_report(path: Path) -> dict[str, str]:
    metrics = {
        "total_reads": "-",
        "reads_with_adapters": "-",
        "reads_written": "-",
        "total_bases": "-",
        "quality_trimmed": "-",
        "bases_written": "-",
    }
    key_map = {
        "total reads processed": "total_reads",
        "reads with adapters": "reads_with_adapters",
        "reads written (passing filters)": "reads_written",
        "total basepairs processed": "total_bases",
        "quality-trimmed": "quality_trimmed",
        "total written (filtered)": "bases_written",
    }
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if ":" not in line:
            continue
        key, value = [part.strip() for part in line.split(":", 1)]
        metric = key_map.get(key.lower())
        if metric:
            metrics[metric] = re.sub(r"\s+", " ", value)
    return metrics


def trim_report_read_label(path: Path) -> str:
    lower = path.name.lower()
    if "_val_1" in lower or re.search(r"(?:^|[._-])r?1(?:[._-]|$)", lower):
        return "Read1"
    if "_val_2" in lower or re.search(r"(?:^|[._-])r?2(?:[._-]|$)", lower):
        return "Read2"
    return "FASTQ"


def render_trim_galore_section(results_dir: Path, report_dir: Path) -> str:
    reports = sorted((results_dir / "trimmed").rglob("*_trimming_report.txt")) if (results_dir / "trimmed").exists() else []
    if not reports:
        return render_missing("Trim Galore 过滤统计", "未发现 trimmed/*/*_trimming_report.txt，无法展示接头和低质量过滤统计。")
    rows = []
    for report in reports:
        metrics = parse_trim_galore_report(report)
        rows.append([
            report.parent.name,
            trim_report_read_label(report),
            metrics["total_reads"],
            metrics["reads_with_adapters"],
            metrics["reads_written"],
            metrics["total_bases"],
            metrics["quality_trimmed"],
            metrics["bases_written"],
        ])
    return render_simple_table(
        "Trim Galore 过滤统计",
        ["样本", "端", "Total reads processed", "Reads with adapters", "Reads written", "Total basepairs processed", "Quality-trimmed", "Total written"],
        rows,
        "Trim Galore 过滤统计：每一行代表一个样本的一端 FASTQ；Total reads processed 为输入 reads 数，Reads with adapters 为检出接头的 reads 数及比例，Reads written 为通过过滤后保留的 reads 数及比例，Total basepairs processed 为输入碱基数，Quality-trimmed 为因质量剪切去除的碱基数，Total written 为过滤后保留的碱基数及比例。",
    )


def first_existing_file(paths: Sequence[Path]) -> Path | None:
    for path in paths:
        if path.exists() and path.is_file():
            return path
    return None


def star_download_link(path: Path | None, report_dir: Path, label: str) -> str:
    if path is None:
        return "-"
    return render_file_link(path, report_dir, label)


def render_alignment_section(results_dir: Path, report_dir: Path) -> str:
    logs = sorted((results_dir / "star").glob("*/Log.final.out")) if (results_dir / "star").exists() else []
    if not logs:
        return render_missing("STAR 比对结果", "未发现 star/*/Log.final.out。")
    wanted = [
        "Number of input reads",
        "Uniquely mapped reads number",
        "Uniquely mapped reads %",
        "% of reads mapped to multiple loci",
        "% of reads unmapped: too short",
    ]
    rows = []
    for log in logs:
        sample = log.parent.name
        metrics = parse_star_log(log)
        bam = first_existing_file([
            log.parent / f"{sample}.sorted.bam",
            log.parent / "Aligned.sortedByCoord.out.bam",
        ])
        rows.append(
            [sample]
            + [metrics.get(key, "-") for key in wanted]
            + [star_download_link(bam, report_dir, "下载 BAM")]
        )
    description = (
        "STAR 比对结果关键指标：样本为样本名称；Number of input reads 为进入 STAR 的 reads 数；"
        "Uniquely mapped reads number 为唯一比对 reads 数；Uniquely mapped reads % 为唯一比对比例；"
        "% of reads mapped to multiple loci 为多位置比对比例；% of reads unmapped: too short 为因序列过短未比对上的比例；"
        "BAM 下载为坐标排序后的比对结果文件，可用于 IGV 浏览、定量复核和下游分析。"
    )
    headers = ["样本"] + wanted + ["BAM 下载"]
    return render_simple_table("STAR 比对结果", headers, rows, description)


def raw_data_index(results_dir: Path) -> dict[str, Path]:
    return {path.name: path for path in collect_files(results_dir.parent / "seq_data", None)}


def sample_rows(input_data: dict) -> list[list[str]]:
    samples = input_data.get("RnaSeq.samples", [])
    rows = []
    for sample in samples:
        rows.append([
            sample.get("sample_id", ""),
            sample.get("group", ""),
            sample.get("read1_path", ""),
            sample.get("read2_path", ""),
        ])
    return rows


def render_sample_table(input_data: dict, results_dir: Path, report_dir: Path) -> str:
    raw_files = raw_data_index(results_dir)
    fastqc_stats = collect_fastqc_stats(results_dir)
    headers = ["样本编号", "分组", "Read1", "Read2", "Total Sequences", "Total Bases", "Sequence length", "%GC"]
    rows = []
    width_rows = []
    for sample in input_data.get("RnaSeq.samples", []):
        sample_id = str(sample.get("sample_id", ""))
        row = [sample_id, sample.get("group", "")]
        width_row = [sample_id, sample.get("group", "")]
        for key in ("read1_path", "read2_path"):
            value = sample.get(key, "")
            filename = Path(str(value)).name
            raw_path = raw_files.get(filename)
            if raw_path:
                row.append(render_file_link(raw_path, report_dir, filename))
            else:
                row.append(esc(value))
            width_row.append(filename or str(value))
        sample_stats = fastqc_stats.get(sample_id, [])
        for metric in ("Total Sequences", "Total Bases", "Sequence length", "%GC"):
            formatted = format_fastqc_stat(sample_stats, metric)
            row.append(formatted)
            width_row.append(formatted)
        rows.append(row)
        width_rows.append(width_row)
    widths = column_widths(width_rows, len(headers))
    colgroup = render_colgroup(widths)
    thead = "".join(render_cell(col, "th", width=width) for col, width in zip(headers, widths))
    body_rows = []
    for row in rows:
        cells = []
        for idx, (cell, width) in enumerate(zip(row, widths)):
            if idx in (2, 3):
                plain_text = re.sub(r"<[^>]+>", "", str(cell))
                cells.append(f'<td{cell_width_style(width)} data-full-text="{esc(plain_text)}"><span class="cell-text">{cell}</span></td>')
            else:
                cells.append(render_cell(cell, width=width))
        body_rows.append("<tr>" + "".join(cells) + "</tr>")
    tbody = "\n".join(body_rows)
    return f"""
    <article class="card table-card">
      <h4>样本信息</h4>
      {one_line_text("每一行代表一个样本；Read1/Read2 列为原始测序文件下载链接；Total Sequences、Total Bases、Sequence length 和 %GC 来自 FastQC Basic Statistics，用于快速查看每个样本双端 FASTQ 的数据量、读长和 GC 含量。")}
      <div class="table-tools">
        <label>搜索/筛选：<input type="search" class="table-search" placeholder="输入关键词"></label>
        <label>每页：<select class="page-size"><option selected>10</option><option>20</option><option>50</option><option>100</option></select></label>
        <span class="page-info"></span>
        <button type="button" class="prev-page">上一页</button>
        <button type="button" class="next-page">下一页</button>
      </div>
      <div class="table-wrap"><table class="report-table">{colgroup}<thead><tr>{thead}</tr></thead><tbody>{tbody}</tbody></table></div>
    </article>
    """


def project_rows(input_data: dict, args: argparse.Namespace) -> list[list[str]]:
    return [
        ["项目名称", args.project_name],
        ["客户名称", args.client_name or "-"],
        ["基因组", basename_or_value(input_data.get("RnaSeq.genome_fasta_path", "-"))],
        ["注释GTF", basename_or_value(input_data.get("RnaSeq.annotation_gtf_path", "-"))],
        ["物种", input_data.get("RnaSeq.enrichment_organism", "-")],
    ]


def section(title: str, section_id: str, body: str) -> str:
    return f'<section id="{esc(section_id)}"><h2>{esc(title)}</h2>{body}</section>'


def render_tables_for_files(title: str, files: Sequence[Path], report_dir: Path, small_limit: int, preview_limit: int) -> str:
    files = [path for path in files if path.name != "sample_metadata.tsv"]
    if not files:
        return render_missing(title, "未发现相关表格文件。")
    cards = []
    for path in files:
        cards.append(render_table(chinese_table_title(path), path, report_dir, small_limit, preview_limit, describe_table_file(path)))
    return f"<div class='stack'>{''.join(cards)}</div>"


def splicing_comparison_name(path: Path, splicing_root: Path) -> str:
    try:
        relative = path.relative_to(splicing_root)
    except ValueError:
        return path.parent.name
    if len(relative.parts) > 1:
        return relative.parts[0]
    return path.parent.name if path.parent != splicing_root else "可变剪切结果"


def is_jcec_file(path: Path) -> bool:
    return path.name.lower().endswith("jcec.txt")


def is_jc_file(path: Path) -> bool:
    lower = path.name.lower()
    return lower.endswith("jc.txt") and not lower.endswith("jcec.txt")


def is_splicing_summary_file(path: Path) -> bool:
    return path.name.lower() == "summary.txt"


def splicing_event_sort_key(path: Path) -> tuple[int, str]:
    name = path.name.lower()
    event_order = ["se", "a5ss", "a3ss", "mxe", "ri"]
    for idx, event in enumerate(event_order):
        if name.startswith(event + ".") or f".{event}." in name:
            return (idx, path.name)
    return (len(event_order), path.name)


def splicing_result_description(result_type: str) -> str:
    source = "JCEC 同时使用 junction counts 和 exon body counts" if result_type == "JCEC" else "JC 仅使用 junction counts"
    return (
        f"差异剪接结果（{result_type}）：{source} 估计剪切差异；每一行代表一个 rMATS 可变剪切事件。"
        "常见字段含义：ID 为事件编号；GeneID/geneSymbol 为基因 ID 和基因名称；chr/strand 为染色体和链方向；"
        "exonStart_0base/exonEnd/upstreamES/upstreamEE/downstreamES/downstreamEE 等为事件相关外显子或上下游外显子坐标；"
        "IJC_SAMPLE_*/SJC_SAMPLE_* 分别为 inclusion junction counts 和 skipping junction counts；"
        "IncFormLen/SkipFormLen 为包含/跳跃 isoform 的有效长度；PValue 为原始显著性，FDR 为多重校正显著性；"
        "IncLevel_SAMPLE_* 为各样本或各组的包含比例，IncLevelDifference 为两组包含比例差值，正负方向需结合比较顺序解读。"
    )


def render_splicing_file_tabs(title: str, files: Sequence[Path], result_type: str, report_dir: Path, preview_limit: int) -> str:
    if not files:
        return render_missing(title, f"未发现 {result_type}.txt 结果文件。")
    panels = []
    description = splicing_result_description(result_type)
    for path in sorted(files, key=splicing_event_sort_key):
        panels.append((path.name, render_table(f"{title}：{path.name}", path, report_dir, preview_limit, preview_limit, description)))
    return render_tabs(title, f"点击文件标签切换不同 {result_type}.txt 事件表；差异剪接事件表通常较大，页面仅加载前 {preview_limit} 行样例数据，完整结果请通过表格下载链接获取。", panels)


def render_splicing_section(splicing_tables: Sequence[Path], report_dir: Path, small_limit: int, preview_limit: int) -> str:
    root = None
    for path in splicing_tables:
        for parent in [path.parent, *path.parents]:
            if parent.name == "alternative_splicing":
                root = parent
                break
        if root is not None:
            break
    if root is None:
        return render_missing("rMATS 可变剪切分析", "未发现 alternative_splicing 结果目录。")

    grouped: dict[str, dict[str, list[Path]]] = {}
    for path in splicing_tables:
        if not (is_splicing_summary_file(path) or is_jcec_file(path) or is_jc_file(path)):
            continue
        comparison = splicing_comparison_name(path, root)
        bucket = grouped.setdefault(comparison, {"summary": [], "jcec": [], "jc": []})
        if is_splicing_summary_file(path):
            bucket["summary"].append(path)
        elif is_jcec_file(path):
            bucket["jcec"].append(path)
        elif is_jc_file(path):
            bucket["jc"].append(path)

    panels = []
    for comparison in sorted(grouped):
        items = grouped[comparison]
        blocks = []
        for summary_path in sorted(items["summary"]):
            blocks.append(render_table(
                f"可变剪切统计汇总（{comparison} ）",
                summary_path,
                report_dir,
                small_limit,
                preview_limit,
                "rMATS summary.txt 统计表：每一行代表一种可变剪切事件类型，列标题说明该组间比较中检测到的事件数量、显著事件数量或 FDR 阈值统计。",
            ))
        blocks.append(render_splicing_file_tabs("差异剪接结果（JCEC）", items["jcec"], "JCEC", report_dir, preview_limit))
        blocks.append(render_splicing_file_tabs("差异剪接结果（JC）", items["jc"], "JC", report_dir, preview_limit))
        panels.append((comparison, f"<div class='stack'>{''.join(blocks)}</div>"))

    return render_tabs(
        "可变剪切分析结果（组间两两比较）",
        "点击标签页切换每个组间比较；每个比较展示 summary.txt 统计表，并分别通过文件标签切换 JCEC.txt 和 JC.txt 差异剪接事件样例结果，完整事件表请下载原始文件。",
        panels,
    )


def de_comparison_name(path: Path) -> str | None:
    match = re.match(r"(.+)\.(all|deg)\.(tsv|csv|txt)$", path.name)
    return match.group(1) if match else None


def render_tabs(title: str, intro: str, panels: Sequence[tuple[str, str]]) -> str:
    if not panels:
        return render_missing(title, "未发现可切换的组间比较结果。")
    buttons = []
    panel_html = []
    for idx, (label, body) in enumerate(panels):
        active = " active" if idx == 0 else ""
        safe_label = esc(label)
        buttons.append(f'<button type="button" class="tab-button{active}" data-tab="{safe_label}">{safe_label}</button>')
        panel_html.append(f'<div class="tab-panel{active}" data-tab="{safe_label}">{body}</div>')
    return f"""
    <article class="card tab-card">
      <h4>{esc(title)}</h4>
      {one_line_text(intro)}
      <div class="tab-buttons">{"".join(buttons)}</div>
      <div class="tab-panels">{"".join(panel_html)}</div>
    </article>
    """


def overview_table_sort_key(path: Path) -> tuple[int, str]:
    if path.name == "pairwise_de_summary.tsv":
        return (0, path.name)
    if path.name == "filtered_count_matrix.tsv":
        return (1, path.name)
    if path.name == "vst_count_matrix.tsv":
        return (2, path.name)
    if path.name == "normalized_count_matrix.tsv":
        return (3, path.name)
    if path.name == "fpkm_count_matrix.tsv":
        return (4, path.name)
    if path.name == "tpm_count_matrix.tsv":
        return (5, path.name)
    return (6, path.name)


def comparison_plot_sort_key(path: Path) -> tuple[int, str]:
    name = path.name.lower()
    if "volcano" in name:
        return (0, path.name)
    if name.endswith(".ma.pdf") or ".ma." in name:
        return (1, path.name)
    return (2, path.name)


def render_de_section(
    de_tables: Sequence[Path],
    plot_files: Sequence[Path],
    report_dir: Path,
    small_limit: int,
    preview_limit: int,
) -> str:
    de_tables = [path for path in de_tables if path.name != "sample_metadata.tsv"]
    matrix_paths = {path.name: path for path in de_tables if path.name in {
        "vst_count_matrix.tsv",
        "normalized_count_matrix.tsv",
        "fpkm_count_matrix.tsv",
        "tpm_count_matrix.tsv",
    }}
    hidden_matrix_names = {"vst_count_matrix.tsv", "normalized_count_matrix.tsv", "fpkm_count_matrix.tsv", "tpm_count_matrix.tsv"}
    grouped_deg: dict[str, list[Path]] = {}
    grouped_all: dict[str, Path] = {}
    overview: list[Path] = []
    for path in de_tables:
        comparison = de_comparison_name(path)
        if comparison and path.name.lower().endswith(".deg.tsv"):
            grouped_deg.setdefault(comparison, []).append(path)
        elif comparison and path.name.lower().endswith(".all.tsv"):
            grouped_all[comparison] = path
        elif path.name not in hidden_matrix_names:
            overview.append(path)

    comparisons = sorted(set(grouped_deg) | set(grouped_all))
    comparison_plots = {
        comparison: sorted([p for p in plot_files if p.name.startswith(comparison + ".") or p.name.startswith(comparison + "_")], key=comparison_plot_sort_key)
        for comparison in comparisons
    }
    used_plot_names = {p.name for files in comparison_plots.values() for p in files}
    overview_plots = sorted([p for p in plot_files if p.name not in used_plot_names], key=priority_de_plot)
    priority_plots = [p for p in overview_plots if priority_de_plot(p) < 2]
    other_overview_plots = [p for p in overview_plots if priority_de_plot(p) >= 2]

    body = ""
    if priority_plots:
        body += render_gallery("样本总体关系图", priority_plots, report_dir, "PCA 图横轴为 PC1、纵轴为 PC2；样本距离热图横轴和纵轴均为样本，优先用于查看样本整体关系。")
    ordered_overview = sorted(overview, key=overview_table_sort_key)
    body += render_tables_for_files("差异表达汇总/辅助表", ordered_overview, report_dir, small_limit, preview_limit)
    fpkm_path = matrix_paths.get("fpkm_count_matrix.tsv")
    if fpkm_path:
        matrix_downloads = []
        if "tpm_count_matrix.tsv" in matrix_paths:
            matrix_downloads.append((matrix_paths["tpm_count_matrix.tsv"], "下载TPM矩阵"))
        if "normalized_count_matrix.tsv" in matrix_paths:
            matrix_downloads.append((matrix_paths["normalized_count_matrix.tsv"], "下载标准化矩阵"))
        body += render_table(
            "标准化表达矩阵",
            fpkm_path,
            report_dir,
            small_limit,
            preview_limit,
            describe_table_file(fpkm_path),
            matrix_downloads,
            "下载FPKM矩阵",
        )

    panels = []
    for comparison in comparisons:
        panel = ""
        if comparison_plots[comparison]:
            panel += render_gallery(f"{comparison} 差异表达图", comparison_plots[comparison], report_dir, "火山图横轴为 log2FoldChange、纵轴为 -log10 显著性；MA 图横轴为平均表达量、纵轴为 log2FoldChange，二者均先于显著差异基因表展示。")
        table_cards = []
        for deg_path in sorted(grouped_deg.get(comparison, [])):
            extras = []
            if comparison in grouped_all:
                extras.append((grouped_all[comparison], "下载完整差异表达表"))
            table_cards.append(render_table(chinese_table_title(deg_path), deg_path, report_dir, 10**12, preview_limit, describe_table_file(deg_path), extras))
        if table_cards:
            panel += f"<div class='stack'>{''.join(table_cards)}</div>"
        elif comparison in grouped_all:
            panel += render_missing(f"{comparison} 显著差异基因表", "未发现 deg.tsv，仅提供完整差异表达表下载。")
            panel += render_download(grouped_all[comparison], report_dir, "下载完整差异表达表")
        panels.append((comparison, panel))
    body += render_tabs("差异表达分析结果（组间两两比较）", "点击标签页切换不同组间比较；每个标签中先展示 MA/火山图，再展示显著差异基因表，完整差异表达表仅作为下载链接提供。", panels)

    return body


def enrichment_plot_sort_key(path: Path) -> tuple[int, str]:
    name = path.name.lower()
    if "go_barplot" in name:
        return (0, path.name)
    if "go_dotplot" in name:
        return (1, path.name)
    if "go_dag" in name:
        return (2, path.name)
    if "kegg_barplot" in name:
        return (3, path.name)
    if "kegg_dotplot" in name:
        return (4, path.name)
    if "reactome_barplot" in name:
        return (5, path.name)
    if "reactome_dotplot" in name:
        return (6, path.name)
    return (7, path.name)


def enrichment_table_sort_key(path: Path) -> tuple[int, str]:
    lower = path.name.lower()
    if lower.endswith("_go.csv"):
        return (0, path.name)
    if lower.endswith("_kegg.csv"):
        return (1, path.name)
    if lower.endswith("_reactome.csv"):
        return (2, path.name)
    return (3, path.name)



def enrichment_kind(path: Path) -> str | None:
    lower = path.name.lower()
    if "_go" in lower:
        return "go"
    if "_kegg" in lower:
        return "kegg"
    if "_reactome" in lower:
        return "reactome"
    return None


def matching_enrichment_plots(table_path: Path, plots: Sequence[Path]) -> list[Path]:
    kind = enrichment_kind(table_path)
    if not kind:
        return []
    comparison = enrichment_comparison_name(table_path)
    if not comparison:
        return []
    prefix = f"{comparison}_{kind}".lower()
    return sorted([path for path in plots if path.name.lower().startswith(prefix)], key=enrichment_plot_sort_key)


def render_inline_enrichment_plots(table_path: Path, plots: Sequence[Path], report_dir: Path) -> str:
    matched = matching_enrichment_plots(table_path, plots)
    if not matched:
        return ""
    figures = "".join(f'<div class="gallery-item">{render_image(path, report_dir)}</div>' for path in matched)
    description = describe_image_file(matched[0]) if len(matched) == 1 else "富集 barplot/dotplot 紧跟对应结果表展示；横轴通常为 RichFactor、富集比例或基因数量，纵轴为功能条目或通路名称。"
    return f"""
    <div class="enrichment-inline-plots">
      {one_line_text(description, "figure-desc")}
      <div class="gallery">{figures}</div>
    </div>
    """

def should_display_enrichment_table(path: Path) -> bool:
    lower = path.name.lower()
    return lower.endswith("_go.csv") or lower.endswith("_kegg.csv") or lower.endswith("_reactome.csv")


def render_enrichment_section(
    enrichment_tables: Sequence[Path],
    enrichment_html: Sequence[Path],
    enrichment_plots: Sequence[Path],
    report_dir: Path,
    small_limit: int,
    preview_limit: int,
) -> str:
    grouped: dict[str, dict[str, list[Path]]] = {}
    overview_tables: list[Path] = []
    for path in enrichment_tables:
        if not should_display_enrichment_table(path):
            continue
        comparison = enrichment_comparison_name(path)
        if comparison:
            grouped.setdefault(comparison, {"tables": [], "html": [], "plots": []})["tables"].append(path)
        else:
            overview_tables.append(path)
    for path in enrichment_plots:
        comparison = enrichment_comparison_name(path)
        if comparison:
            grouped.setdefault(comparison, {"tables": [], "html": [], "plots": []})["plots"].append(path)

    panels = []
    for comparison in sorted(grouped):
        items = grouped[comparison]
        table_blocks = []
        for table_path in sorted(items["tables"], key=enrichment_table_sort_key):
            table_blocks.append(render_table(chinese_table_title(table_path), table_path, report_dir, 10**12, preview_limit, describe_table_file(table_path)))
            table_blocks.append(render_inline_enrichment_plots(table_path, items["plots"], report_dir))
        panel = f"<div class='stack'>{''.join(table_blocks)}</div>" if table_blocks else render_missing(f"{comparison} 富集结果表", "未发现 GO/KEGG/Reactome 富集结果表。")
        panels.append((comparison, panel))

    body = render_tabs("富集分析结果（组间两两比较）", "点击标签页切换不同组间比较；仅展示 GO、KEGG 和 Reactome 富集结果，表格不截断，全部行都会加载。", panels)
    if overview_tables:
        body += render_tables_for_files("富集分析其他结果表", overview_tables, report_dir, 10**12, preview_limit)
    return body


def render_report_tabs(panels: Sequence[tuple[str, str, str]]) -> str:
    if not panels:
        return ""
    links = []
    sections = []
    for idx, (title, section_id, body) in enumerate(panels):
        active = " active" if idx == 0 else ""
        safe_id = esc(section_id)
        links.append(f'<a class="tab-button report-tab-button{active}" href="#{safe_id}" data-target="{safe_id}">{esc(title)}</a>')
        sections.append(f'<section id="{safe_id}" class="report-tab-panel"><h2>{esc(title)}</h2>{body}</section>')
    return f"""
    <div class="report-layout">
      <aside class="report-tabs" aria-label="报告内容导航">
        <div class="tab-buttons report-tab-buttons">{"".join(links)}</div>
      </aside>
      <div class="report-sections">{"".join(sections)}</div>
    </div>
    """


def render_coexpression_section(
    coexpression_tables: Sequence[Path],
    plot_files: Sequence[Path],
    report_dir: Path,
    small_limit: int,
    preview_limit: int,
) -> str:
    table_map = {path.name: path for path in coexpression_tables}
    plot_map = {path.name: path for path in plot_files}

    trait_plots = [
        plot_map[name]
        for name in ["wgcna_module_trait_heatmap.pdf", "wgcna_sample_clustering.pdf"]
        if name in plot_map
    ]
    phenotype_block = render_gallery(
        "模块-表型热图和样本聚类图",
        trait_plots,
        report_dir,
        "模块-表型热图横轴为实验分组/表型，纵轴为共表达模块，颜色表示相关系数并用于定位与分组最相关的模块；样本聚类图横轴为样本、纵轴为聚类高度，用于检查样本整体表达相似性和离群情况。",
    )
    important_modules = table_map.get("important_modules.tsv")
    phenotype_block += render_table(
        "关键共表达模块列表",
        important_modules,
        report_dir,
        small_limit,
        preview_limit,
        describe_table_file(important_modules),
    ) if important_modules else render_missing("关键共表达模块列表", "未发现 important_modules.tsv。")

    hub_genes = table_map.get("hub_genes.tsv")
    hub_block = render_table(
        "核心调控基因（Hub Gene）",
        hub_genes,
        report_dir,
        small_limit,
        preview_limit,
        describe_table_file(hub_genes),
    ) if hub_genes else render_missing("核心调控基因（Hub Gene）", "未发现 hub_genes.tsv。")
    hub_block += render_single_figure_card(
        "模块大小统计图",
        plot_map.get("wgcna_module_size_barplot.pdf"),
        report_dir,
        "模块大小柱状图横轴为模块，纵轴为模块包含基因数量，用于判断各共表达模块规模和后续重点解读模块。",
    )

    more_plots = [
        plot_map[name]
        for name in ["wgcna_gene_dendrogram.pdf", "wgcna_module_eigengene_heatmap.pdf"]
        if name in plot_map
    ]
    more_block = render_gallery(
        "更多结果图",
        more_plots,
        report_dir,
        "基因聚类树横轴为参与网络分析的基因、纵轴为聚类高度，用于查看基因模块划分；模块特征向量热图横轴为样本，纵轴为模块 eigengene，颜色表示模块整体表达趋势。",
    )
    gene_modules = table_map.get("gene_modules.tsv")
    more_block += render_table(
        "基因模块对应表",
        gene_modules,
        report_dir,
        small_limit,
        preview_limit,
        describe_table_file(gene_modules),
    ) if gene_modules else render_missing("基因模块对应表", "未发现 gene_modules.tsv。")

    sections = [
        ("模块与实验分组相关性分析", phenotype_block),
        ("核心调控基因（Hub Gene）", hub_block),
        ("更多结果", more_block),
    ]
    return "".join(f"<article class='card'><h4>{esc(title)}</h4>{content}</article>" for title, content in sections)


def render_fusion_section(
    fusion_tables: Sequence[Path],
    fusion_plots: Sequence[Path],
    report_dir: Path,
    small_limit: int,
    preview_limit: int,
) -> str:
    table_map = {path.name: path for path in fusion_tables}
    plot_map = {path.name: path for path in fusion_plots}
    blocks: list[str] = []

    heatmap_pdf = plot_map.get("fusion_heatmap.pdf")
    if heatmap_pdf:
        blocks.append(render_gallery(
            "Fusion Heatmap",
            [heatmap_pdf],
            report_dir,
            "Fusion × Sample 热图以 FFPM 为数值，展示各融合事件在不同样本中的表达支持强度；矩阵保留至少一个样本 FFPM > 0.1 的 Fusion。",
        ))
    fusion_heatmap = table_map.get("fusion_heatmap.tsv")
    if fusion_heatmap:
        blocks.append(render_table(
            "Fusion × Sample FFPM矩阵",
            fusion_heatmap,
            report_dir,
            small_limit,
            preview_limit,
            "Fusion × Sample 矩阵：每一行代表一个 Fusion，列为样本，数值为 FFPM；仅保留至少一个样本 FFPM > 0.1 的 Fusion。",
        ))

    circos_pdfs = sorted([p for p in fusion_plots if p.name.endswith(".fusion_circos.pdf")], key=lambda p: p.name)
    if circos_pdfs:
        blocks.append(render_gallery(
            "Fusion Circos Plot",
            circos_pdfs,
            report_dir,
            "Circos plot shows genomic connections between fusion gene breakpoints. Each arc represents one fusion event. Thicker arcs indicate stronger fusion support based on junction read count or FFPM.",
        ))
    circos_pngs = sorted([p for p in fusion_plots if p.name.endswith(".fusion_circos.png")], key=lambda p: p.name)
    if circos_pngs:
        blocks.append(render_file_links("Fusion Circos PNG下载", circos_pngs, report_dir, "每个样本的 Fusion Circos PNG 图片，可用于报告复用或快速预览。"))
    circos_links = sorted([p for p in fusion_tables if p.name.endswith(".fusion_circos_links.tsv")], key=lambda p: p.name)
    if circos_links:
        blocks.append(render_file_links("Fusion Circos绘图数据", circos_links, report_dir, "每个样本用于绘制 Circos 连接线的数据表，包含断点、染色体位置、JunctionReadCount、SpanningFragCount 和 FFPM。"))

    return "".join(blocks) if blocks else render_missing("Fusion Gene Analysis", "未发现 Fusion Heatmap 或 Fusion Circos 结果。")


def build_report(args: argparse.Namespace) -> str:
    input_path = Path(args.input_json).resolve()
    results_dir = Path(args.results_dir).resolve()
    report_dir = Path(args.report_dir).resolve()
    input_data = load_input_json(input_path)

    sample_table = render_sample_table(input_data, results_dir, report_dir)
    project_table = render_simple_table("项目信息与分析参数", ["项目", "内容"], project_rows(input_data, args), "项目基础信息、参考基因组和主要分析参数。")

    count_tables = collect_files(results_dir / "counts", TABLE_EXTENSIONS)
    de_tables = [path for path in collect_files(results_dir / "differential_expression", TABLE_EXTENSIONS) if path.name not in {"sample_metadata.tsv", "filtered_count_matrix.tsv"}]
    enrichment_tables = collect_files(results_dir / "enrichment", {".csv", ".tsv", ".txt"})
    enrichment_html = collect_files(results_dir / "enrichment", {".html"})
    splicing_tables = collect_files(results_dir / "alternative_splicing", TABLE_EXTENSIONS)
    coexpression_tables = collect_files(results_dir / "coexpression", TABLE_EXTENSIONS)
    fusion_plot_tables = collect_files(results_dir / "plots" / "fusion", TABLE_EXTENSIONS)
    plot_files = collect_files(results_dir / "plots", IMAGE_EXTENSIONS)
    coexpression_plots = collect_files(results_dir / "coexpression" / "plots", IMAGE_EXTENSIONS)
    enrichment_plots = collect_files(results_dir / "enrichment" / "plots", IMAGE_EXTENSIONS)
    fusion_plots = collect_files(results_dir / "plots" / "fusion", IMAGE_EXTENSIONS)

    report_panels = [
        ("1. 项目信息", "project", project_table),
        ("2. 样本信息", "samples", sample_table),
        ("3. 质控结果", "qc", render_method_note("质控与剪切方法", "FastQC / Trim Galore", "重点关注原始数据每碱基质量是否随读长下降、碱基组成是否异常，以及 Trim Galore 后 reads/bases 保留比例和接头污染比例；低保留率或高接头比例提示样本或建库质量需重点排查。") + render_fastqc_section(results_dir, report_dir) + render_trim_galore_section(results_dir, report_dir)),
        ("4. 比对结果", "alignment", render_method_note("比对方法", "STAR", "重点关注唯一比对比例、多重比对比例和过短未比对比例；唯一比对比例偏低可能与样本污染、物种/参考基因组不匹配、reads 质量或 rRNA/重复序列比例偏高有关。") + render_alignment_section(results_dir, report_dir)),
        ("5. 定量结果", "counts", render_method_note("基因定量方法", "featureCounts", "重点关注 Assigned 比例和各类 Unassigned 原因；Assigned 偏低时需结合比对结果、GTF 注释版本和文库类型判断是否存在注释不匹配或比对质量问题。") + render_tables_for_files("featureCounts 结果", count_tables, report_dir, args.small_table_rows, args.big_table_preview_rows)),
        ("6. 差异表达分析", "de", render_method_note("差异表达分析方法", "DESeq2 / ggplot2 / pheatmap", "重点结合 PCA/样本距离判断重复一致性，再解读火山图、MA 图和 DEG 表；显著性以 padj/FDR 和表达变化倍数共同判断，低表达基因的倍数变化需谨慎解释。") + render_de_section(de_tables, plot_files, report_dir, args.small_table_rows, args.big_table_preview_rows)),
        ("7. GO / KEGG / Reactome 富集分析", "enrichment", render_method_note("功能富集分析方法", "bioconductor 镜像内 /kegg_data/gene_cluster_enrich.R", "富集结果反映差异基因集合中功能条目/通路的过度代表性；解读时应同时关注校正后显著性、RichFactor、命中基因数和通路生物学相关性，避免只按 p 值排序。") + render_enrichment_section(enrichment_tables, enrichment_html, enrichment_plots, report_dir, args.small_table_rows, args.big_table_preview_rows)),
        ("8. 可变剪切分析", "splicing", render_method_note("可变剪切分析方法", "rMATS", "重点关注 FDR、IncLevelDifference 和支持剪切事件的 reads 数；低 reads 支持或重复间差异较大的事件应谨慎解读，必要时结合基因浏览器验证。") + render_splicing_section(splicing_tables, report_dir, args.small_table_rows, args.big_table_preview_rows)),
        ("9. Fusion Gene Analysis", "fusion", render_method_note("Fusion Gene Analysis", "STAR-Fusion / pyCirclize", "重点关注 FFPM、JunctionReadCount 和 SpanningFragCount 支持强度；Circos 图用于观察融合断点在染色体之间的连接关系，低支持度事件需谨慎解释。") + render_fusion_section(fusion_plot_tables, fusion_plots, report_dir, args.small_table_rows, args.big_table_preview_rows)),
        ("10. 基因共表达网络分析", "coexpression", render_method_note("共表达网络分析方法", "WGCNA", "重点关注模块内基因共表达一致性、模块特征向量与样本分组/表型的关系以及 hub genes；样本量较少时模块稳定性有限，应结合生物学背景解释。") + render_coexpression_section(coexpression_tables, coexpression_plots, report_dir, args.small_table_rows, args.big_table_preview_rows)),
        ("11. 蛋白互作网络分析", "ppi", render_ppi_section()),
    ]
    sections = render_report_tabs(report_panels)
    nav = ""

    return HTML_TEMPLATE.format(
        title=esc(args.report_title),
        project=esc(args.project_name),
        intro=esc(REPORT_INTRO),
        nav=nav,
        sections=sections,
    )


HTML_TEMPLATE = r"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{title}</title>
  <style>
    :root {{ --bg:#fff7f8; --card:#fffdfc; --text:#3f2d32; --muted:#8a6b72; --brand:#c06c84; --brand-dark:#8d4054; --line:#ead5db; --soft:#f8e8ec; }}
    * {{ box-sizing: border-box; }}
    body {{ margin:0; font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","Microsoft YaHei",Arial,sans-serif; color:var(--text); background:var(--bg); }}
    header {{ padding:32px 40px; background:linear-gradient(135deg,#9f4f63,#d89aaa); color:white; }}
    header h1 {{ margin:0 0 8px; font-size:32px; }}
    header p {{ margin:4px 0; opacity:.95; }}
    nav {{ display:none; }}
    main {{ max-width:1440px; margin:0 auto; padding:24px; }}
    section {{ margin:0 0 28px; }}
    h2 {{ margin:24px 0 12px; padding-left:10px; border-left:5px solid var(--brand); }}
    h4 {{ margin:0 0 8px; }}
    .card {{ background:var(--card); border:1px solid var(--line); border-radius:14px; padding:18px; margin:14px 0; box-shadow:0 6px 20px rgba(15,23,42,.05); }}
    .card-header {{ display:flex; justify-content:space-between; gap:12px; align-items:flex-start; min-width:0; }}
    .card-header > div {{ min-width:0; max-width:100%; }}
    .table-card .text-ellipsis,.gallery-card .text-ellipsis,.tab-card .text-ellipsis,.method-card .text-ellipsis {{ display:block; max-width:100%; min-width:0; }}
    .muted,.note {{ color:var(--muted); font-size:14px; }}
    .download {{ color:var(--brand); white-space:nowrap; }}
    .download-group {{ display:flex; gap:10px; flex-wrap:wrap; justify-content:flex-start; margin-top:8px; }}
    .table-tools {{ display:flex; gap:12px; align-items:center; flex-wrap:wrap; margin:10px 0; }}
    input,select,button {{ padding:6px 8px; border:1px solid var(--line); border-radius:8px; background:white; }}
    button {{ cursor:pointer; }}
    .table-wrap {{ width:100%; overflow:auto; border:1px solid var(--line); border-radius:10px; }}
    table {{ width:max-content; min-width:100%; table-layout:fixed; border-collapse:collapse; font-size:13px; }}
    th,td {{ padding:8px 10px; border-bottom:1px solid var(--line); vertical-align:top; max-width:280px; min-width:0; overflow:hidden; }}
    th {{ background:#fff1df; position:sticky; top:0; z-index:1; text-align:left; }}
    .cell-text {{ display:block; max-width:min(260px, 100%); overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }}
    .cell-popover {{ position:fixed; z-index:9999; max-width:min(720px, calc(100vw - 32px)); max-height:260px; overflow:auto; padding:10px 12px; border:1px solid var(--line); border-radius:10px; background:#fffdfc; color:var(--text); box-shadow:0 12px 32px rgba(63,45,50,.18); white-space:pre-wrap; word-break:break-word; user-select:text; font-size:13px; line-height:1.55; }}
    .text-ellipsis {{ max-width:100%; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }}
    .cell-link {{ color:var(--brand-dark); text-decoration:underline; }}
    tr:nth-child(even) td {{ background:#fffdf9; }}
    .tab-buttons {{ display:flex; gap:18px; flex-wrap:wrap; margin:12px 0; border-bottom:1px solid var(--line); }}
    .tab-card > .tab-buttons {{ position:sticky; top:0; z-index:4; background:rgba(255,253,252,.97); padding-top:8px; }}
    .tab-button {{ border:0; border-bottom:3px solid transparent; border-radius:0; background:transparent; color:var(--muted); padding:10px 2px 9px; }}
    .tab-button.active {{ color:var(--brand-dark); border-bottom-color:var(--brand); font-weight:700; }}
    .report-layout {{ display:grid; grid-template-columns:220px minmax(0, 1fr); gap:24px; align-items:start; }}
    .report-tabs {{ position:sticky; top:16px; background:var(--card); border:1px solid var(--line); border-radius:14px; padding:12px; box-shadow:0 6px 20px rgba(15,23,42,.05); }}
    .report-tab-buttons {{ display:flex; flex-direction:column; gap:4px; margin:0; border-bottom:0; }}
    .report-tab-button {{ text-decoration:none; flex:0 0 auto; border-left:3px solid transparent; border-bottom:0; padding:9px 10px; border-radius:8px; }}
    .report-tab-button.active {{ border-left-color:var(--brand); border-bottom-color:transparent; background:var(--soft); }}
    .report-tab-panel {{ margin:0 0 28px; scroll-margin-top:24px; }}
    html {{ scroll-behavior:smooth; }}
    .tab-panel {{ display:none; margin-top:14px; }}
    .tab-panel.active {{ display:block; }}
    .figure {{ margin:0; }}
    .figure img,.figure object {{ max-width:100%; width:100%; border:1px solid var(--line); border-radius:12px; background:white; }}
    .figure img {{ min-height:360px; object-fit:contain; }}
    .embedded-static-figure img {{ min-height:0; }}
    .figure object {{ min-height:518px; }}
    .module-size-figure-card .figure object {{ min-height:259px; }}
    .module-size-figure-card .figure img {{ min-height:180px; }}
    figcaption {{ color:var(--muted); font-size:13px; margin-top:6px; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }}
    .figure-desc {{ color:var(--muted); font-size:13px; margin:4px 0 0; }}
    .gallery {{ display:grid; grid-template-columns:repeat(2, minmax(0, 1fr)); gap:16px; }}
    .gallery-item {{ min-width:0; }}
    .fastqc-read-grid {{ display:grid; grid-template-columns:repeat(2, minmax(0, 1fr)); gap:16px; }}
    .fastqc-read-card {{ border:1px solid var(--line); border-radius:12px; padding:12px; background:#fffdfc; min-width:0; }}
    .fastqc-read-header {{ display:flex; gap:10px; justify-content:space-between; align-items:center; margin-bottom:8px; }}
    .fastqc-plot-grid {{ display:grid; grid-template-columns:repeat(2, minmax(0, 1fr)); gap:12px; margin-top:10px; }}
    .fastqc-figure img {{ min-height:260px; }}
    .enrichment-inline-plots {{ margin:-4px 0 18px; padding:0 18px 12px; }}
    .file-list {{ columns:2; }}
    .file-list li {{ margin:0 0 8px; break-inside:avoid; }}
    .missing {{ border-style:dashed; color:var(--muted); }}
    footer {{ text-align:center; color:var(--muted); padding:24px; }}
    @media (max-width: 760px) {{ header {{ padding:24px; }} main {{ padding:14px; }} .report-layout {{ grid-template-columns:1fr; }} .report-tabs {{ position:static; }} .report-tab-buttons {{ flex-direction:row; overflow:auto; }} .file-list {{ columns:1; }} .card-header {{ display:block; }} .gallery,.fastqc-read-grid,.fastqc-plot-grid {{ grid-template-columns:1fr; }} }}
  </style>
</head>
<body>
  <header>
    <h1>{title}</h1>
    <p>项目：{project}</p>
    <p>{intro}</p>
  </header>
  <nav>{nav}</nav>
  <main>{sections}</main>
  <footer>RNA-seq 静态分析报告</footer>
  <script>
    function setupTables() {{
      document.querySelectorAll('.table-card').forEach(function(card) {{
        const table = card.querySelector('table');
        if (!table) return;
        const rows = Array.from(table.querySelectorAll('tbody tr'));
        const search = card.querySelector('.table-search');
        const pageSize = card.querySelector('.page-size');
        const prev = card.querySelector('.prev-page');
        const next = card.querySelector('.next-page');
        const info = card.querySelector('.page-info');
        let page = 1;
        function render() {{
          const q = (search.value || '').toLowerCase();
          const filtered = rows.filter(row => row.textContent.toLowerCase().includes(q));
          const size = parseInt(pageSize.value, 10) || 10;
          const pages = Math.max(1, Math.ceil(filtered.length / size));
          page = Math.min(page, pages);
          rows.forEach(row => row.style.display = 'none');
          filtered.slice((page - 1) * size, page * size).forEach(row => row.style.display = '');
          info.textContent = `第 ${{page}} / ${{pages}} 页，共 ${{filtered.length}} 行`;
        }}
        search.addEventListener('input', () => {{ page = 1; render(); }});
        pageSize.addEventListener('change', () => {{ page = 1; render(); }});
        prev.addEventListener('click', () => {{ page = Math.max(1, page - 1); render(); }});
        next.addEventListener('click', () => {{ page += 1; render(); }});
        render();
      }});
    }}
    function setupReportNav() {{
      const links = Array.from(document.querySelectorAll('.report-tab-button'));
      const sections = links.map(link => document.querySelector(link.getAttribute('href'))).filter(Boolean);
      function activate(id) {{
        links.forEach(link => link.classList.toggle('active', link.getAttribute('href') === '#' + id));
      }}
      function updateActiveFromScroll() {{
        let current = sections[0];
        sections.forEach(section => {{
          if (section.getBoundingClientRect().top <= 120) current = section;
        }});
        if (current) activate(current.id);
      }}
      links.forEach(link => link.addEventListener('click', () => activate(link.dataset.target)));
      window.addEventListener('scroll', updateActiveFromScroll, {{ passive: true }});
      window.addEventListener('resize', updateActiveFromScroll);
      updateActiveFromScroll();
    }}
    function setupCellPopovers() {{
      const popover = document.createElement('div');
      popover.className = 'cell-popover';
      popover.hidden = true;
      document.body.appendChild(popover);
      let hideTimer = null;
      function clearHideTimer() {{
        if (hideTimer) window.clearTimeout(hideTimer);
        hideTimer = null;
      }}
      function scheduleHide() {{
        clearHideTimer();
        hideTimer = window.setTimeout(() => {{ popover.hidden = true; }}, 180);
      }}
      function showForCell(cell) {{
        const text = cell.dataset.fullText || cell.textContent || '';
        const textSpan = cell.querySelector('.cell-text');
        if (!text.trim() || !textSpan || textSpan.scrollWidth <= textSpan.clientWidth + 2) return;
        clearHideTimer();
        popover.textContent = text;
        popover.hidden = false;
        const rect = cell.getBoundingClientRect();
        const top = Math.min(window.innerHeight - popover.offsetHeight - 12, rect.bottom + 6);
        const left = Math.min(window.innerWidth - popover.offsetWidth - 12, Math.max(12, rect.left));
        popover.style.top = Math.max(12, top) + 'px';
        popover.style.left = left + 'px';
      }}
      document.querySelectorAll('td[data-full-text], th[data-full-text]').forEach(cell => {{
        cell.addEventListener('mouseenter', () => showForCell(cell));
        cell.addEventListener('focusin', () => showForCell(cell));
        cell.addEventListener('mouseleave', scheduleHide);
        cell.addEventListener('focusout', scheduleHide);
      }});
      popover.addEventListener('mouseenter', clearHideTimer);
      popover.addEventListener('mouseleave', scheduleHide);
      window.addEventListener('scroll', () => {{ popover.hidden = true; }}, {{ passive: true }});
    }}
    function setupTabs() {{
      document.querySelectorAll('.tab-card').forEach(function(card) {{
        const buttonBar = card.querySelector(':scope > .tab-buttons');
        const panelWrap = card.querySelector(':scope > .tab-panels');
        if (!buttonBar || !panelWrap) return;
        const buttons = Array.from(buttonBar.querySelectorAll(':scope > .tab-button'));
        const panels = Array.from(panelWrap.querySelectorAll(':scope > .tab-panel'));
        function show(value) {{
          buttons.forEach(button => button.classList.toggle('active', button.dataset.tab === value));
          panels.forEach(panel => panel.classList.toggle('active', panel.dataset.tab === value));
        }}
        buttons.forEach(button => button.addEventListener('click', () => show(button.dataset.tab)));
        if (buttons.length) show(buttons[0].dataset.tab);
      }});
    }}
    setupTables();
    setupReportNav();
    setupCellPopovers();
    setupTabs();
  </script>
</body>
</html>
"""


def main() -> None:
    args = parse_args()
    report_dir = Path(args.report_dir).resolve()
    report_dir.mkdir(parents=True, exist_ok=True)
    html_text = build_report(args)
    output_path = report_dir / "index.html"
    output_path.write_text(html_text, encoding="utf-8")
    print(f"报告已生成: {output_path}")


if __name__ == "__main__":
    main()

# RNA-seq WDL workflow

This repository contains a WDL 1.0 paired-end RNA-seq workflow in `workflows/rna_seq.wdl`.

The workflow performs:

1. STAR genome-index generation, or reuse of a supplied STAR index directory.
2. Raw-read FastQC.
3. Adapter/quality trimming with Trim Galore.
4. Trimmed-read FastQC.
5. STAR alignment to coordinate-sorted BAM and per-gene STAR counts.
6. Gene-level summarization with featureCounts.
7. Pairwise differential-expression analysis across all sample groups with DESeq2.
8. GO, KEGG, and Reactome pathway enrichment display for each pairwise DEG set.
9. Alternative splicing analysis with rMATS across all pairwise sample groups.
10. Fusion-gene analysis with STAR-Fusion when a CTAT genome library is supplied, plus Fusion × Sample FFPM heatmap and per-sample Fusion Circos plots.
11. Gene co-expression network analysis with WGCNA.
12. PCA, sample-distance heatmap, volcano, MA, DEG heatmap, enrichment bar/dot-plot, and curated WGCNA visualization.
13. Aggregated QC reporting with MultiQC.

## Docker execution model

The Docker runner command is configured once with `docker_run`:

```wdl
String docker_run = "docker run --rm --security-opt seccomp=unconfined  -v /home/data/vip01/work:/home/data/vip01/work"
```

Each tool image parameter contains only the image name, for example:

```wdl
String fastqc_image = "registry.cn-guangzhou.aliyuncs.com/origen/fastqc"
```

Tasks combine these two values directly:

```bash
~{docker_run} ~{image} bash -c "fastqc --extract --threads ~{threads} --outdir '~{output_dir}/fastqc/~{sample_id}' '~{read1_path}' '~{read2_path}'"
~{docker_run} ~{image} bash -c 'fastqc_dir="~{output_dir}/fastqc/~{sample_id}"; for zip_file in "${fastqc_dir}"/*_fastqc.zip; do [ -e "${zip_file}" ] || continue; unzip -o "${zip_file}" -d "${fastqc_dir}"; done'
```

Because the default Docker runner mounts `/home/data/vip01/work`, all default input and output paths should be absolute paths under `/home/data/vip01/work` unless you override `docker_run` with a different mount. The workflow defaults to the mm39 FASTA `/home/data/vip01/work/pipeline/database/mm39/Mus_musculus.GRCm39.dna.toplevel.fa` and GTF `/home/data/vip01/work/pipeline/database/mm39/Mus_musculus.GRCm39.115.gtf`.

## STAR index reuse

Set `RnaSeq.star_index_path` to an existing STAR index directory to skip `BuildStarIndex` and use that index directly. For the provided mm39 setup, use:

```json
"RnaSeq.star_index_path": "/home/data/vip01/work/pipeline/database/mm39/star_index_seqlen150"
```

If `RnaSeq.star_index_path` is empty, the workflow builds a new STAR index under `output_dir/star_index` using `RnaSeq.genome_fasta_path` and `RnaSeq.annotation_gtf_path`.

The mouse and human reference FASTA files used here include primary chromosomes plus additional contigs, and chromosome names do not have a `chr` prefix. Keep the FASTA, GTF, and STAR index naming conventions consistent within each species-specific run.


## Sample groups and downstream analysis

Each sample object must include `sample_id`, `group`, `read1_path`, and `read2_path`. If two or more groups are present, the workflow runs all pairwise group comparisons and writes:

- `differential_expression/*.all.tsv`: full DESeq2 result table for each comparison; the first columns are `gene_symbol` and `gene_id`, and the final columns are per-sample `FPKM_<sample_id>` values.
- `differential_expression/*.deg.tsv`: filtered DEG table for each comparison using `padj_cutoff` and `log2fc_cutoff`, with the same gene-symbol and per-sample FPKM columns as the full table.
- `differential_expression/pairwise_de_summary.tsv`: DEG counts for all pairwise comparisons.
- `differential_expression/normalized_count_matrix.tsv`: DESeq2 size-factor-normalized count matrix.
- `differential_expression/fpkm_count_matrix.tsv`: FPKM matrix calculated from filtered raw counts, featureCounts gene lengths, and library sizes.
- `differential_expression/tpm_count_matrix.tsv`: TPM matrix calculated from filtered raw counts and featureCounts gene lengths.
- `differential_expression/vst_count_matrix.tsv`: DESeq2 variance-stabilized matrix used by PCA, sample-distance plots, and WGCNA.
- `enrichment/<comparison>_GO.csv`, `<comparison>_KEGG.csv`, and `<comparison>_Reactome.csv`: GO BP, KEGG, and Reactome enrichment tables generated from each pairwise DEG list when available.
- `enrichment/plots/*.pdf`: GO/KEGG/Reactome enrichment plots produced by `/kegg_data/gene_cluster_enrich.R`; the report focuses on `<comparison>_GO.csv`, `<comparison>_KEGG.csv`, and `<comparison>_Reactome.csv`, while mapping files, gene-list files, GO classification tables, and KEGG HTML files remain downloadable in the results directory but are not embedded in the report.
- `plots/*.pdf`: PCA, sample distance heatmap, volcano plots, MA plots, top-DEG heatmaps, and enrichment plots.

Functional enrichment no longer requires a custom annotation table. The `FunctionalEnrichment` task extracts the `gene_id` column from every `*.deg.tsv` file, writes `enrichment/gene_lists/<comparison>.genes.txt`, and runs the bundled bioconductor-image script:

```bash
Rscript /kegg_data/gene_cluster_enrich.R gene_list.txt gene_id_type organism output_dir prefix
```

Set `RnaSeq.enrichment_gene_id_type` to the DEG gene identifier type accepted by the selected OrgDb, for example `ENSEMBL`, `SYMBOL`, `ENTREZID`, or `REFSEQ`. Set `RnaSeq.enrichment_organism` to `mouse`/`mmu` or `human`/`hsa`; the provided mouse and human inputs default to `ENSEMBL` plus `mouse` or `human` respectively.

Additional downstream tasks are included: `AlternativeSplicing` runs rMATS pairwise for every group pair from the coordinate-sorted BAMs; `FusionGeneAnalysis` runs STAR-Fusion from trimmed FASTQs when `fusion_ctat_lib_path` points to a CTAT genome library and `FusionVisualization` runs `scripts/plot_fusion_circos.py` to create `plots/fusion/fusion_heatmap.tsv`, `plots/fusion/fusion_heatmap.pdf`, per-sample `*.fusion_circos_links.tsv`, `*.fusion_circos.pdf`, and `*.fusion_circos.png` (mouse default `/home/data/vip01/work/pipeline/database/mm39/ctat_mm39_lib`; human input default `/home/data/vip01/work/pipeline/database/hg38/ctat_hg38_lib`); `CoexpressionNetwork` runs `/kegg_data/wgcna_analysis.R` from the bioconductor image on the DESeq2 VST matrix and sample metadata, then writes gene modules, module eigengenes, module-trait relationships, hub genes, and WGCNA plots. rMATS uses `registry.cn-guangzhou.aliyuncs.com/origen/rmats` by default, while WGCNA uses `registry.cn-guangzhou.aliyuncs.com/origen/bioconductor`. STAR alignment now emits chimeric junctions (`Chimeric.out.junction`) for fusion-aware outputs, and `pca.pdf` labels each point with the sample name, expands the PC2 coordinate span, and uses a taller aspect ratio so the vertical axis is visibly longer rather than only increasing the PDF canvas size.

## Example inputs

See `docs/rna_seq.inputs.json` for a minimal paired-end example. Update the sample FASTQ paths, sample groups, reference FASTA, annotation GTF, STAR index path, enrichment gene ID type/organism, thread count, output directory, thresholds, `docker_run`, and image names before running. The root-level `input.json` is a larger project-specific example with multiple groups.

## Example Cromwell run

```bash
java -jar cromwell.jar run workflows/rna_seq.wdl --inputs docs/rna_seq.inputs.json
```

The default output directory is `/home/data/vip01/work/rnaseq_results`.

## Standalone TrimReads script

For temporary direct execution outside WDL, run the generated Trim Galore script:

```bash
bash scripts/run_trimreads_MJ20260515137.sh
```

The script uses the same `/home/data/vip01/work` Docker mount, Trim Galore image, project FASTQ paths, output directory, and thread count as `input.json`.

## Standalone StarAlign script

After trimming has completed, run STAR alignment directly outside WDL with:

```bash
bash scripts/run_staralign_MJ20260515137.sh
```

The script reads `${OUTPUT_DIR}/trimmed/<sample>/<sample>_val_1.fq.gz` and `${OUTPUT_DIR}/trimmed/<sample>/<sample>_val_2.fq.gz`, uses the configured mm39 STAR index, and writes each sorted BAM to `${OUTPUT_DIR}/star/<sample>/<sample>.sorted.bam`.

## Standalone scripts for every task

For ad-hoc execution without Cromwell, this repository includes one shell helper per WDL task:

1. Optional STAR index build: `bash scripts/run_buildstarindex_MJ20260515137.sh`
2. Raw FastQC: `bash scripts/run_fastqc_MJ20260515137.sh raw`
3. TrimReads / Trim Galore: `bash scripts/run_trimreads_MJ20260515137.sh`
4. Trimmed FastQC: `bash scripts/run_fastqc_MJ20260515137.sh trimmed`
5. STAR alignment: `bash scripts/run_staralign_MJ20260515137.sh`
6. featureCounts: `bash scripts/run_featurecounts_MJ20260515137.sh` (the count matrix header is rewritten from BAM paths to sample IDs)
7. Differential expression, VST matrix, and DE plots: `bash scripts/run_differential_expression_MJ20260515137.sh`
8. GO/KEGG/Reactome enrichment with gene_cluster_enrich.R: `bash scripts/run_functional_enrichment_MJ20260515137.sh`
9. Alternative splicing with rMATS: `bash scripts/run_alternative_splicing_MJ20260515137.sh`
10. Fusion genes with STAR-Fusion plus heatmap/circos plots (mouse default CTAT library: `/home/data/vip01/work/pipeline/database/mm39/ctat_mm39_lib`): `bash scripts/run_fusion_gene_MJ20260515137.sh`
11. Gene co-expression network with WGCNA from DESeq2 VST data: `bash scripts/run_coexpression_network_MJ20260515137.sh`
12. MultiQC: `bash scripts/run_multiqc_MJ20260515137.sh`

`run_fastqc_MJ20260515137.sh all` runs both raw and trimmed FastQC. FastQC runs with `--extract`, and `.zip` outputs are also explicitly decompressed in each sample FastQC directory so the extracted `*_fastqc/` result folders are available for direct inspection and downstream report delivery. Run the trimmed FastQC and downstream scripts only after the preceding scripts have produced their expected files.

The differential-expression helpers create and append `${OUTPUT_DIR}/logs/differential_expression.<timestamp>.log` from inside the Docker container to avoid host-side permission problems on mounted project directories. R-based helpers run Docker with `-i` for heredoc stdin, otherwise `Rscript -` receives no script and exits immediately. If `counts/gene_counts.tsv` is missing or empty, the script stops immediately with a clear error in that log; run the matching `run_featurecounts_*` script first and then rerun differential expression.


## Static Chinese HTML report

After the analysis finishes, generate an offline deliverable report with the dependency-free Python helper:

```bash
python scripts/generate_rnaseq_report.py \
  --input-json input.json \
  --results-dir /home/data/vip01/work/bioproject/MJ20260515137/rnaseq_results \
  --report-dir /home/data/vip01/work/bioproject/MJ20260515137/rnaseq_results/report \
  --project-name MJ20260515137
```

The generated `report/index.html` is a standalone Chinese report using a Nature Rose palette and opens with an RNA-seq background/introduction paragraph instead of a report-generation timestamp. It summarizes project information with Chinese labels (`基因组`, `注释GTF`, `物种`) and only displays reference file basenames, embeds raw sequencing download links directly in the sample table `Read1`/`Read2` columns by scanning the sibling `seq_data` directory next to `--results-dir`, and no longer creates a separate raw-data section. The sample table also reads FastQC `Basic Statistics` and displays `Total Sequences`, `Total Bases`, `Sequence length`, and `%GC` for Read1/Read2. The QC section switches by sample and shows FastQC extracted `Per base sequence quality` and `Per base sequence content` plots for each FASTQ, with Read1/Read2 displayed in one row and FastQC HTML downloads available; MultiQC results are not embedded in that section. Plots/PDFs use relative paths for offline viewing, PDF sidebars are closed by default, PDF display boxes are 20% taller again, and figure descriptions include axis meanings such as PC1/PC2, log2FoldChange, -log10 significance, RichFactor, sample axes, and heatmap color meanings; the STAR alignment result table explains each displayed metric and provides per-sample coordinate-sorted BAM download links when those files are present. Report sections use a sticky left-side navigation bar that jumps to each section and automatically highlights the section currently visible while keeping all sections visible during scrolling; sample/group comparisons and enrichment comparisons still use sticky underline-style switching tabs so comparison controls remain accessible while scrolling within long sections. Tables default to 10 rows per page, left-align headers, truncate long headers and cells with selectable floating hover previews instead of expanding table row height, size columns from data values, keep featureCounts `Chr`/`Start`/`End`/`Strand` columns compact, use Chinese table titles, place download links below the table explanation text, and include detailed row/column explanation text for interpreting each chart and table. Trim Galore filtering statistics from `*_trimming_report.txt` are shown in the QC section. Each analysis section states the software used and result-interpretation notes, and differential-expression tabs show MA/volcano plots before the `deg.tsv` table, load significant DEG tables fully regardless of size, label DEG downloads as `下载显著差异表达表`, hide `*.all.tsv` from embedded display while adding a `下载完整差异表达表` link beside the DEG download, hide `filtered_count_matrix.tsv` and `vst_count_matrix.tsv` from the report body, display `fpkm_count_matrix.tsv` as `标准化表达矩阵`, and provide dedicated FPKM, TPM, and DESeq2 normalized-matrix download links. Enrichment tabs embed GO/KEGG/Reactome CSV tables and plots, place GO/KEGG/Reactome barplot and dotplot PDFs plus GO DAG PDFs directly below their corresponding enrichment result tables, load enrichment CSV rows without truncation, and omit `sample_metadata.tsv`, GO classification tables, KEGG HTML files, gene ID mapping files, and gene-list text files from the report body. Fusion Gene Analysis is embedded in the report and displays Fusion Heatmap plus Fusion Circos Plot outputs from `plots/fusion`. WGCNA plots are read from `coexpression/plots`. The co-expression section is limited to three customer-facing blocks: `模块与实验分组相关性分析` (`模块-表型热图和样本聚类图`, plus `关键共表达模块列表` from `important_modules.tsv`), `核心调控基因（Hub Gene）` (`hub_genes.tsv` plus a half-height `wgcna_module_size_barplot.pdf` preview), and `更多结果` (`wgcna_gene_dendrogram.pdf`, `wgcna_module_eigengene_heatmap.pdf`, plus `基因模块对应表` from `gene_modules.tsv`). WGCNA PDF descriptions now explain each plot's axes, colors, and interpretation separately. Alternative-splicing tabs switch by pairwise comparison; each comparison shows the rMATS `summary.txt` statistics table and nested file tabs for all `JCEC.txt` and `JC.txt` differential splicing event tables, loading only preview rows for these large event files while providing detailed field explanations and original-file downloads. The final `蛋白互作网络分析` section provides an online STRING analysis link and embeds `scripts/string.input.png` plus `scripts/string.demo.png` as Base64 data URIs so the report does not depend on external image files.

## Separate human-side RNA-seq subset

`N1_mix`, `N2`, `N3`, `E1`, `E2`, and `E3` are separated into `input_human.json` for an independent human-side analysis. The original `input.json` now keeps only the mouse-side `LPS`, `exo`, and `con` samples.

The human-side input uses the requested hg38 paths:

```json
"RnaSeq.genome_fasta_path": "/home/data/vip01/work/pipeline/database/hg38/Homo_sapiens.GRCh38.dna.toplevel.fa",
"RnaSeq.annotation_gtf_path": "/home/data/vip01/work/pipeline/database/hg38/Homo_sapiens.GRCh38.113.gtf",
"RnaSeq.star_index_path": "/home/data/vip01/work/pipeline/database/hg38/star_index_seqlen150"
```

Run the human-side standalone scripts with the `_human` suffix, for example:

```bash
bash scripts/run_trimreads_MJ20260515137_human.sh
bash scripts/run_staralign_MJ20260515137_human.sh
bash scripts/run_featurecounts_MJ20260515137_human.sh
bash scripts/run_differential_expression_MJ20260515137_human.sh
bash scripts/run_functional_enrichment_MJ20260515137_human.sh
bash scripts/run_alternative_splicing_MJ20260515137_human.sh
bash scripts/run_fusion_gene_MJ20260515137_human.sh  # defaults to /home/data/vip01/work/pipeline/database/hg38/ctat_hg38_lib
bash scripts/run_coexpression_network_MJ20260515137_human.sh
```

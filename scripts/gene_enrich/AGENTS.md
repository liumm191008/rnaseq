# Gene Cluster Enrichment Tool

## 工具介绍

本仓库提供 `gene_cluster_enrich.R`，用于对输入基因列表进行聚类/功能富集分析。该脚本支持：

- GO Biological Process 富集分析；
- GO 二级分类统计；
- 离线 KEGG 富集分析，使用镜像内 `/kegg_data` 下的 `gson` 与 KEGG class TSV 文件；
- Reactome pathway 富集分析；
- 可选背景基因列表；
- 富集结果 CSV、KEGG HTML、GO/KEGG/Reactome PDF 图输出。

Docker 镜像构建后，脚本会被放置在容器内：

```text
/kegg_data/gene_cluster_enrich.R
```

## 使用方法

```bash
Rscript /kegg_data/gene_cluster_enrich.R \
  gene_list.txt \
  gene_id_type \
  organism \
  output_dir \
  prefix \
  [background_gene_list.txt]
```

参数说明：

- `gene_list.txt`：输入基因列表文件；支持一行一个基因，也支持逗号、空格或制表符分隔；
- `gene_id_type`：输入基因 ID 类型，需为对应 OrgDb 支持的 keytype，例如 `ENTREZID`、`SYMBOL`、`ENSEMBL`、`REFSEQ`；
- `organism`：物种，支持 `mouse`/`mmu` 或 `human`/`hsa`；
- `output_dir`：输出目录；
- `prefix`：输出文件前缀；
- `background_gene_list.txt`：可选背景基因列表，ID 类型需与 `gene_id_type` 一致；不填写时使用默认富集背景。

示例：

```bash
Rscript /kegg_data/gene_cluster_enrich.R \
  cluster_genes.txt \
  SYMBOL \
  mouse \
  results \
  cluster1
```

带背景基因列表示例：

```bash
Rscript /kegg_data/gene_cluster_enrich.R \
  cluster_genes.txt \
  SYMBOL \
  mouse \
  results \
  cluster1 \
  background_genes.txt
```

主要输出：

- `<prefix>_gene_id_mapping.csv`：输入基因 ID 映射结果；
- `<prefix>_background_gene_id_mapping.csv`：可选背景基因 ID 映射结果；
- `<prefix>_GO.csv`：GO BP 富集结果；
- `<prefix>_GO_classification.csv`：GO 分类统计结果；
- `<prefix>_KEGG.csv`：KEGG 富集结果；
- `<prefix>_KEGG.html`：带 KEGG pathway 链接的 HTML 结果；
- `<prefix>_Reactome.csv`：Reactome 富集结果；
- `plots/`：GO、KEGG、Reactome 的 dotplot/barplot 以及 GO DAG PDF 图。

## 维护约定

# User-provided custom instructions

所有代码作者写：liumm

FROM registry.cn-guangzhou.aliyuncs.com/origen/bioconductor

LABEL maintainer="liumm" \
      org.opencontainers.image.authors="liumm" \
      org.opencontainers.image.description="Gene-list enrichment pipeline with offline KEGG support"

COPY gene_cluster_enrich.R /kegg_data/gene_cluster_enrich.R
RUN chmod 0755 /kegg_data/gene_cluster_enrich.R

WORKDIR /kegg_data

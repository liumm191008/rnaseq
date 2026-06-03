version 1.0


## Paired-end RNA-seq sample definition.
## All paths must be absolute host paths available through the docker_run mount.
struct RnaSeqSample {
  String sample_id
  String group
  String read1_path
  String read2_path
}

workflow RnaSeq {
  input {
    Array[RnaSeqSample] samples
    String test1
    String test2
    String genome_fasta_path = "/home/data/vip01/work/pipeline/database/mm39/Mus_musculus.GRCm39.dna.toplevel.fa"
    String annotation_gtf_path = "/home/data/vip01/work/pipeline/database/mm39/Mus_musculus.GRCm39.115.gtf"
    String star_index_path = ""
    String output_dir = "/home/data/vip01/work/rnaseq_results"
    String enrichment_gene_id_type = "ENSEMBL"
    String enrichment_organism = "mouse"
    String fusion_ctat_lib_path = ""
    Int threads = 8
    Int sjdb_overhang = 149
    Int min_count = 10
    Float padj_cutoff = 0.05
    Float log2fc_cutoff = 1.0

    String docker_run = "docker run --rm --security-opt seccomp=unconfined  -v /home/data/vip01/work:/home/data/vip01/work"
    String fastqc_image = "registry.cn-guangzhou.aliyuncs.com/origen/fastqc"
    String trim_galore_image = "registry.cn-guangzhou.aliyuncs.com/origen/trim-galore"
    String star_image = "registry.cn-guangzhou.aliyuncs.com/origen/star"
    String subread_image = "registry.cn-guangzhou.aliyuncs.com/origen/subread"
    String multiqc_image = "registry.cn-guangzhou.aliyuncs.com/origen/multiqc"
    String differential_analysis_image = "registry.cn-guangzhou.aliyuncs.com/origen/bioconductor"
    String enrichment_image = "registry.cn-guangzhou.aliyuncs.com/origen/bioconductor"
    String rmats_image = "registry.cn-guangzhou.aliyuncs.com/origen/rmats"
    String star_fusion_image = "registry.cn-guangzhou.aliyuncs.com/origen/star-fusion"
    String coexpression_image = "registry.cn-guangzhou.aliyuncs.com/origen/bioconductor"
  }

  if (star_index_path == "") {
    call BuildStarIndex {
      input:
        genome_fasta_path = genome_fasta_path,
        annotation_gtf_path = annotation_gtf_path,
        output_dir = output_dir,
        threads = threads,
        sjdb_overhang = sjdb_overhang,
        docker_run = docker_run,
        image = star_image
    }
  }

  String resolved_star_index_path = if star_index_path != "" then star_index_path else select_first([BuildStarIndex.genome_dir])

  scatter (sample in samples) {
    String current_sample_id = sample.sample_id
    String current_sample_group = sample.group

    call FastQc as RawFastQc {
      input:
        sample_id = sample.sample_id,
        read1_path = sample.read1_path,
        read2_path = sample.read2_path,
        output_dir = output_dir,
        threads = threads,
        docker_run = docker_run,
        image = fastqc_image
    }

    call TrimReads {
      input:
        sample_id = sample.sample_id,
        read1_path = sample.read1_path,
        read2_path = sample.read2_path,
        output_dir = output_dir,
        threads = threads,
        docker_run = docker_run,
        image = trim_galore_image
    }

    call FastQc as TrimmedFastQc {
      input:
        sample_id = sample.sample_id + ".trimmed",
        read1_path = TrimReads.trimmed_read1_path,
        read2_path = TrimReads.trimmed_read2_path,
        output_dir = output_dir,
        threads = threads,
        docker_run = docker_run,
        image = fastqc_image
    }

    call StarAlign {
      input:
        sample_id = sample.sample_id,
        read1_path = TrimReads.trimmed_read1_path,
        read2_path = TrimReads.trimmed_read2_path,
        genome_dir = resolved_star_index_path,
        output_dir = output_dir,
        threads = threads,
        docker_run = docker_run,
        image = star_image
    }
  }

  call FeatureCounts {
    input:
      bam_paths = StarAlign.sorted_bam_path,
      sample_ids = current_sample_id,
      annotation_gtf_path = annotation_gtf_path,
      output_dir = output_dir,
      threads = threads,
      docker_run = docker_run,
      image = subread_image
  }

  call DifferentialExpression {
    input:
      count_matrix_path = FeatureCounts.count_matrix_path,
      sample_ids = current_sample_id,
      sample_groups = current_sample_group,
      output_dir = output_dir,
      min_count = min_count,
      padj_cutoff = padj_cutoff,
      log2fc_cutoff = log2fc_cutoff,
      docker_run = docker_run,
      image = differential_analysis_image
  }

  call FunctionalEnrichment {
    input:
      de_summary_path = DifferentialExpression.de_summary_path,
      de_results_dir = DifferentialExpression.de_results_dir,
      gene_id_type = enrichment_gene_id_type,
      organism = enrichment_organism,
      output_dir = output_dir,
      docker_run = docker_run,
      image = enrichment_image
  }

  call AlternativeSplicing {
    input:
      bam_paths = StarAlign.sorted_bam_path,
      sample_ids = current_sample_id,
      sample_groups = current_sample_group,
      annotation_gtf_path = annotation_gtf_path,
      output_dir = output_dir,
      read_length = sjdb_overhang + 1,
      threads = threads,
      docker_run = docker_run,
      image = rmats_image
  }

  call FusionGeneAnalysis {
    input:
      sample_ids = current_sample_id,
      read1_paths = TrimReads.trimmed_read1_path,
      read2_paths = TrimReads.trimmed_read2_path,
      output_dir = output_dir,
      fusion_ctat_lib_path = fusion_ctat_lib_path,
      threads = threads,
      docker_run = docker_run,
      image = star_fusion_image
  }

  call CoexpressionNetwork {
    input:
      vst_count_matrix_path = DifferentialExpression.vst_count_matrix_path,
      sample_metadata_path = DifferentialExpression.sample_metadata_path,
      output_dir = output_dir,
      docker_run = docker_run,
      image = coexpression_image
  }

  call MultiQc {
    input:
      output_dir = output_dir,
      count_matrix_path = FeatureCounts.count_matrix_path,
      docker_run = docker_run,
      image = multiqc_image
  }

  output {
    String star_index_dir = resolved_star_index_path
    Array[String] raw_fastqc_dirs = RawFastQc.fastqc_dir_path
    Array[String] trimmed_fastqc_dirs = TrimmedFastQc.fastqc_dir_path
    Array[String] trimmed_read1 = TrimReads.trimmed_read1_path
    Array[String] trimmed_read2 = TrimReads.trimmed_read2_path
    Array[String] sorted_bams = StarAlign.sorted_bam_path
    Array[String] alignment_logs = StarAlign.final_log_path
    String count_matrix = FeatureCounts.count_matrix_path
    String count_summary = FeatureCounts.summary_path
    String sample_metadata = DifferentialExpression.sample_metadata_path
    String de_summary = DifferentialExpression.de_summary_path
    String de_results_dir = DifferentialExpression.de_results_dir
    String vst_count_matrix = DifferentialExpression.vst_count_matrix_path
    String enrichment_results_dir = FunctionalEnrichment.enrichment_results_dir
    String alternative_splicing_dir = AlternativeSplicing.alternative_splicing_dir
    String fusion_results_dir = FusionGeneAnalysis.fusion_results_dir
    String coexpression_results_dir = CoexpressionNetwork.coexpression_results_dir
    String coexpression_modules = CoexpressionNetwork.gene_module_path
    String coexpression_module_sizes = CoexpressionNetwork.module_sizes_path
    String coexpression_module_eigengenes = CoexpressionNetwork.module_eigengenes_path
    String coexpression_hub_genes = CoexpressionNetwork.hub_genes_path
    String coexpression_input_type = CoexpressionNetwork.input_type_path
    String visualization_dir = DifferentialExpression.visualization_dir
    String pca_plot = DifferentialExpression.pca_plot_path
    String sample_distance_heatmap = DifferentialExpression.sample_distance_heatmap_path
    String multiqc_report = MultiQc.report_path
  }
}

task BuildStarIndex {
  input {
    String genome_fasta_path
    String annotation_gtf_path
    String output_dir
    Int threads
    Int sjdb_overhang
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c "mkdir -p '~{output_dir}/star_index'"
    ~{docker_run} ~{image} bash -c "STAR --runThreadN ~{threads} --runMode genomeGenerate --genomeDir '~{output_dir}/star_index' --genomeFastaFiles '~{genome_fasta_path}' --sjdbGTFfile '~{annotation_gtf_path}' --sjdbOverhang ~{sjdb_overhang}"
  >>>

  output {
    String genome_dir = output_dir + "/star_index"
  }
}

task FastQc {
  input {
    String sample_id
    String read1_path
    String read2_path
    String output_dir
    Int threads
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c "mkdir -p '~{output_dir}/fastqc/~{sample_id}'"
    ~{docker_run} ~{image} bash -c "fastqc --extract --threads ~{threads} --outdir '~{output_dir}/fastqc/~{sample_id}' '~{read1_path}' '~{read2_path}'"
    ~{docker_run} ~{image} bash -c 'fastqc_dir="~{output_dir}/fastqc/~{sample_id}"; for zip_file in "${fastqc_dir}"/*_fastqc.zip; do [ -e "${zip_file}" ] || continue; unzip -o "${zip_file}" -d "${fastqc_dir}"; done'
  >>>

  output {
    String fastqc_dir_path = output_dir + "/fastqc/" + sample_id
  }
}

task TrimReads {
  input {
    String sample_id
    String read1_path
    String read2_path
    String output_dir
    Int threads
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c "mkdir -p '~{output_dir}/trimmed/~{sample_id}'"
    ~{docker_run} ~{image} bash -c "trim_galore --paired --cores ~{threads} --gzip --basename '~{sample_id}' --output_dir '~{output_dir}/trimmed/~{sample_id}' '~{read1_path}' '~{read2_path}'"
    ~{docker_run} ~{image} bash -c "test -s '~{output_dir}/trimmed/~{sample_id}/~{sample_id}_val_1.fq.gz'"
    ~{docker_run} ~{image} bash -c "test -s '~{output_dir}/trimmed/~{sample_id}/~{sample_id}_val_2.fq.gz'"
  >>>

  output {
    String trimmed_read1_path = output_dir + "/trimmed/" + sample_id + "/" + sample_id + "_val_1.fq.gz"
    String trimmed_read2_path = output_dir + "/trimmed/" + sample_id + "/" + sample_id + "_val_2.fq.gz"
    String trim_dir_path = output_dir + "/trimmed/" + sample_id
  }
}

task StarAlign {
  input {
    String sample_id
    String read1_path
    String read2_path
    String genome_dir
    String output_dir
    Int threads
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c "mkdir -p '~{output_dir}/star/~{sample_id}'"
    ~{docker_run} ~{image} bash -c "STAR --runThreadN ~{threads} --genomeDir '~{genome_dir}' --readFilesIn '~{read1_path}' '~{read2_path}' --readFilesCommand zcat --outFileNamePrefix '~{output_dir}/star/~{sample_id}/' --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimOutType Junctions"
    ~{docker_run} ~{image} bash -c "mv '~{output_dir}/star/~{sample_id}/Aligned.sortedByCoord.out.bam' '~{output_dir}/star/~{sample_id}/~{sample_id}.sorted.bam'"
  >>>

  output {
    String sorted_bam_path = output_dir + "/star/" + sample_id + "/" + sample_id + ".sorted.bam"
    String final_log_path = output_dir + "/star/" + sample_id + "/Log.final.out"
    String gene_counts_path = output_dir + "/star/" + sample_id + "/ReadsPerGene.out.tab"
    String chimeric_junction_path = output_dir + "/star/" + sample_id + "/Chimeric.out.junction"
  }
}

task FeatureCounts {
  input {
    Array[String] bam_paths
    Array[String] sample_ids
    String annotation_gtf_path
    String output_dir
    Int threads
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c "mkdir -p '~{output_dir}/counts'"
    ~{docker_run} ~{image} bash -c "featureCounts -T ~{threads} -p -B -C -a '~{annotation_gtf_path}' -o '~{output_dir}/counts/gene_counts.tsv' ~{sep=" " bam_paths}"
    ~{docker_run} ~{image} bash -c "awk -v FS='\t' -v OFS='\t' -v samples='~{sep="," sample_ids}' 'BEGIN { split(samples, sample_names, \",\") } /^#/ { print; next } \$1 == \"Geneid\" { for (i = 7; i <= NF; i++) { idx = i - 6; if (idx in sample_names) \$i = sample_names[idx] } } { print }' '~{output_dir}/counts/gene_counts.tsv' > '~{output_dir}/counts/gene_counts.tsv.tmp' && mv '~{output_dir}/counts/gene_counts.tsv.tmp' '~{output_dir}/counts/gene_counts.tsv'"
  >>>

  output {
    String count_matrix_path = output_dir + "/counts/gene_counts.tsv"
    String summary_path = output_dir + "/counts/gene_counts.tsv.summary"
  }
}

task DifferentialExpression {
  input {
    String count_matrix_path
    Array[String] sample_ids
    Array[String] sample_groups
    String output_dir
    Int min_count
    Float padj_cutoff
    Float log2fc_cutoff
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c "mkdir -p '~{output_dir}/differential_expression' '~{output_dir}/plots'"
    ~{docker_run} -i ~{image} Rscript - <<'RSCRIPT'
    suppressPackageStartupMessages(library(DESeq2))
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(pheatmap))

    count_file <- "~{count_matrix_path}"
    output_dir <- "~{output_dir}"
    min_count <- ~{min_count}
    padj_cutoff <- ~{padj_cutoff}
    log2fc_cutoff <- ~{log2fc_cutoff}
    sample_ids <- c("~{sep='","' sample_ids}")
    sample_groups_raw <- c("~{sep='","' sample_groups}")
    sample_groups <- make.names(sample_groups_raw)

    de_dir <- file.path(output_dir, "differential_expression")
    plot_dir <- file.path(output_dir, "plots")
    dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

    metadata <- data.frame(
      sample_id = sample_ids,
      group = factor(sample_groups),
      original_group = sample_groups_raw,
      stringsAsFactors = FALSE
    )
    rownames(metadata) <- metadata$sample_id
    write.table(metadata, file.path(de_dir, "sample_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    feature_counts <- read.delim(count_file, comment.char = "#", check.names = FALSE)
    if (!"Geneid" %in% colnames(feature_counts)) {
      stop("featureCounts output must contain a Geneid column")
    }
    if (ncol(feature_counts) < 7) {
      stop("featureCounts output must contain annotation columns plus at least one count column")
    }

    gene_ids <- feature_counts$Geneid
    count_data <- feature_counts[, 7:ncol(feature_counts), drop = FALSE]
    if (ncol(count_data) != length(sample_ids)) {
      stop(sprintf("count column number (%s) does not match sample number (%s)", ncol(count_data), length(sample_ids)))
    }
    rownames(count_data) <- gene_ids
    colnames(count_data) <- sample_ids
    count_data <- round(as.matrix(count_data))
    keep <- rowSums(count_data) >= min_count
    count_data <- count_data[keep, , drop = FALSE]
    write.table(count_data, file.path(de_dir, "filtered_count_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)

    if (length(unique(metadata$group)) < 2) {
      write.table(data.frame(), file.path(de_dir, "pairwise_de_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
      writeLines("Only one group was provided; pairwise differential expression was skipped.", file.path(de_dir, "differential_expression.skipped.txt"))
      quit(save = "no", status = 0)
    }

    dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~ group)
    dds <- DESeq(dds)
    vst_counts <- varianceStabilizingTransformation(dds, blind = FALSE)
    write.table(assay(vst_counts), file.path(de_dir, "vst_count_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)

    pca_data <- plotPCA(vst_counts, intgroup = "group", returnData = TRUE)
    if (!"name" %in% colnames(pca_data)) pca_data$name <- rownames(pca_data)
    percent_var <- round(100 * attr(pca_data, "percentVar"))
    x_span <- diff(range(pca_data$PC1, na.rm = TRUE))
    y_span <- diff(range(pca_data$PC2, na.rm = TRUE))
    if (!is.finite(x_span) || x_span == 0) x_span <- 1
    if (!is.finite(y_span) || y_span == 0) y_span <- 1
    target_y_span <- max(y_span, x_span * 0.75)
    y_mid <- mean(range(pca_data$PC2, na.rm = TRUE))
    pca_y_limits <- c(y_mid - target_y_span / 2, y_mid + target_y_span / 2)
    pdf(file.path(plot_dir, "pca.pdf"), width = 7, height = 8)
    print(
      ggplot(pca_data, aes(PC1, PC2, color = group)) +
        geom_point(size = 3) +
        geom_text(aes(label = name), vjust = -0.8, size = 3, show.legend = FALSE) +
        xlab(paste0("PC1: ", percent_var[1], "% variance")) +
        ylab(paste0("PC2: ", percent_var[2], "% variance")) +
        coord_cartesian(ylim = pca_y_limits) +
        theme_bw() +
        theme(aspect.ratio = 1.2)
    )
    dev.off()

    sample_dist <- dist(t(assay(vst_counts)))
    pdf(file.path(plot_dir, "sample_distance_heatmap.pdf"), width = 8, height = 7)
    pheatmap(as.matrix(sample_dist), clustering_distance_rows = sample_dist, clustering_distance_cols = sample_dist)
    dev.off()

    groups <- levels(metadata$group)
    comparisons <- combn(groups, 2, simplify = FALSE)
    summary_list <- list()

    for (comparison in comparisons) {
      group_a <- comparison[1]
      group_b <- comparison[2]
      comparison_name <- paste0(group_b, "_vs_", group_a)
      res <- as.data.frame(results(dds, contrast = c("group", group_b, group_a)))
      res$gene_id <- rownames(res)
      res <- res[, c("gene_id", setdiff(colnames(res), "gene_id"))]
      res <- res[order(res$padj), ]
      all_path <- file.path(de_dir, paste0(comparison_name, ".all.tsv"))
      deg_path <- file.path(de_dir, paste0(comparison_name, ".deg.tsv"))
      write.table(res, all_path, sep = "\t", quote = FALSE, row.names = FALSE)

      deg <- res[!is.na(res$padj) & res$padj <= padj_cutoff & abs(res$log2FoldChange) >= log2fc_cutoff, ]
      write.table(deg, deg_path, sep = "\t", quote = FALSE, row.names = FALSE)
      summary_list[[comparison_name]] <- data.frame(
        comparison = comparison_name,
        group_a = group_a,
        group_b = group_b,
        total_tested = nrow(res),
        deg_total = nrow(deg),
        deg_up = sum(deg$log2FoldChange > 0),
        deg_down = sum(deg$log2FoldChange < 0),
        stringsAsFactors = FALSE
      )

      plot_df <- res
      plot_df$significant <- !is.na(plot_df$padj) & plot_df$padj <= padj_cutoff & abs(plot_df$log2FoldChange) >= log2fc_cutoff
      volcano <- ggplot(plot_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) + geom_point(alpha = 0.6, size = 1) + theme_bw() + ggtitle(comparison_name) + xlab("log2 fold change") + ylab("-log10 adjusted P value")
      ggsave(file.path(plot_dir, paste0(comparison_name, ".volcano.pdf")), volcano, width = 7, height = 6)

      ma <- ggplot(plot_df, aes(x = baseMean, y = log2FoldChange, color = significant)) + geom_point(alpha = 0.6, size = 1) + scale_x_log10() + theme_bw() + ggtitle(comparison_name) + xlab("mean normalized count") + ylab("log2 fold change")
      ggsave(file.path(plot_dir, paste0(comparison_name, ".ma.pdf")), ma, width = 7, height = 6)

      top_genes <- head(deg$gene_id, 50)
      if (length(top_genes) > 1) {
        pdf(file.path(plot_dir, paste0(comparison_name, ".top_deg_heatmap.pdf")), width = 9, height = 9)
        pheatmap(assay(vst_counts)[top_genes, , drop = FALSE], scale = "row", annotation_col = metadata[, "group", drop = FALSE])
        dev.off()
      }
    }

    summary_table <- do.call(rbind, summary_list)
    write.table(summary_table, file.path(de_dir, "pairwise_de_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    RSCRIPT
  >>>

  output {
    String sample_metadata_path = output_dir + "/differential_expression/sample_metadata.tsv"
    String de_summary_path = output_dir + "/differential_expression/pairwise_de_summary.tsv"
    String de_results_dir = output_dir + "/differential_expression"
    String vst_count_matrix_path = output_dir + "/differential_expression/vst_count_matrix.tsv"
    String visualization_dir = output_dir + "/plots"
    String pca_plot_path = output_dir + "/plots/pca.pdf"
    String sample_distance_heatmap_path = output_dir + "/plots/sample_distance_heatmap.pdf"
  }
}

task FunctionalEnrichment {
  input {
    String de_summary_path
    String de_results_dir
    String gene_id_type
    String organism
    String output_dir
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} -i ~{image} bash <<'BASH'
    set -euo pipefail
    de_summary_path="~{de_summary_path}"
    de_results_dir="~{de_results_dir}"
    gene_id_type="~{gene_id_type}"
    organism="~{organism}"
    enrich_dir="~{output_dir}/enrichment"
    mkdir -p "${enrich_dir}" "${enrich_dir}/gene_lists" "${enrich_dir}/plots"

    if [ ! -s "${de_summary_path}" ]; then
      echo "Differential-expression summary is missing or empty; enrichment was skipped." > "${enrich_dir}/enrichment.skipped.txt"
      exit 0
    fi

    comparisons=$(awk -F '\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i=="comparison") c=i; next} c && $c!="" {print $c}' "${de_summary_path}")
    if [ -z "${comparisons}" ]; then
      echo "No pairwise comparisons were available; enrichment was skipped." > "${enrich_dir}/enrichment.skipped.txt"
      exit 0
    fi

    while IFS= read -r comparison; do
      [ -n "${comparison}" ] || continue
      deg_path="${de_results_dir}/${comparison}.deg.tsv"
      gene_list="${enrich_dir}/gene_lists/${comparison}.genes.txt"
      if [ ! -s "${deg_path}" ]; then
        echo "DEG file is missing or empty: ${deg_path}" > "${enrich_dir}/${comparison}.skipped.txt"
        continue
      fi
      awk -F '\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i=="gene_id") c=i; next} c && $c!="" {print $c}' "${deg_path}" | sort -u > "${gene_list}"
      if [ ! -s "${gene_list}" ]; then
        echo "No gene_id values were found in ${deg_path}." > "${enrich_dir}/${comparison}.skipped.txt"
        continue
      fi
      Rscript /kegg_data/gene_cluster_enrich.R "${gene_list}" "${gene_id_type}" "${organism}" "${enrich_dir}" "${comparison}"
    done <<< "${comparisons}"
    BASH
  >>>

  output {
    String enrichment_results_dir = output_dir + "/enrichment"
  }
}


task AlternativeSplicing {
  input {
    Array[String] bam_paths
    Array[String] sample_ids
    Array[String] sample_groups
    String annotation_gtf_path
    String output_dir
    Int read_length
    Int threads
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} -i ~{image} bash <<'BASH'
    set -euo pipefail
    output_dir="~{output_dir}"
    gtf="~{annotation_gtf_path}"
    read_length="~{read_length}"
    threads="~{threads}"
    sample_ids=("~{sep='" "' sample_ids}")
    sample_groups=("~{sep='" "' sample_groups}")
    bam_paths=("~{sep='" "' bam_paths}")
    as_dir="${output_dir}/alternative_splicing"
    mkdir -p "${as_dir}"

    mapfile -t groups < <(printf '%s\n' "${sample_groups[@]}" | sort -u)
    if [ "${#groups[@]}" -lt 2 ]; then
      echo "Only one group was provided; rMATS was skipped." > "${as_dir}/alternative_splicing.skipped.txt"
      exit 0
    fi

    for ((i=0; i<${#groups[@]}; i++)); do
      for ((j=i+1; j<${#groups[@]}; j++)); do
        group_a="${groups[$i]}"
        group_b="${groups[$j]}"
        comparison="${group_b}_vs_${group_a}"
        comparison_dir="${as_dir}/${comparison}"
        tmp_dir="${comparison_dir}/tmp"
        mkdir -p "${comparison_dir}" "${tmp_dir}"
        b1=()
        b2=()
        for ((k=0; k<${#sample_ids[@]}; k++)); do
          if [ "${sample_groups[$k]}" = "${group_a}" ]; then
            b1+=("${bam_paths[$k]}")
          elif [ "${sample_groups[$k]}" = "${group_b}" ]; then
            b2+=("${bam_paths[$k]}")
          fi
        done
        (IFS=,; printf '%s\n' "${b1[*]}") > "${comparison_dir}/${group_a}.bam.list"
        (IFS=,; printf '%s\n' "${b2[*]}") > "${comparison_dir}/${group_b}.bam.list"
        rmats.py --b1 "${comparison_dir}/${group_a}.bam.list" --b2 "${comparison_dir}/${group_b}.bam.list" --gtf "${gtf}" --od "${comparison_dir}" --tmp "${tmp_dir}" -t paired --readLength "${read_length}" --nthread "${threads}"
      done
    done
    BASH
  >>>

  output {
    String alternative_splicing_dir = output_dir + "/alternative_splicing"
  }
}

task FusionGeneAnalysis {
  input {
    Array[String] sample_ids
    Array[String] read1_paths
    Array[String] read2_paths
    String output_dir
    String fusion_ctat_lib_path
    Int threads
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} -i ~{image} bash <<'BASH'
    set -euo pipefail
    output_dir="~{output_dir}"
    ctat_lib="~{fusion_ctat_lib_path}"
    threads="~{threads}"
    sample_ids=("~{sep='" "' sample_ids}")
    read1_paths=("~{sep='" "' read1_paths}")
    read2_paths=("~{sep='" "' read2_paths}")
    fusion_dir="${output_dir}/fusion"
    mkdir -p "${fusion_dir}"
    if [ -z "${ctat_lib}" ] || [ ! -d "${ctat_lib}" ]; then
      echo "fusion_ctat_lib_path is empty or missing; STAR-Fusion was skipped." > "${fusion_dir}/fusion.skipped.txt"
      exit 0
    fi
    for ((i=0; i<${#sample_ids[@]}; i++)); do
      sample_id="${sample_ids[$i]}"
      sample_dir="${fusion_dir}/${sample_id}"
      mkdir -p "${sample_dir}"
      STAR-Fusion --genome_lib_dir "${ctat_lib}" --left_fq "${read1_paths[$i]}" --right_fq "${read2_paths[$i]}" --CPU "${threads}" --output_dir "${sample_dir}"
    done
    BASH
  >>>

  output {
    String fusion_results_dir = output_dir + "/fusion"
  }
}

task CoexpressionNetwork {
  input {
    String vst_count_matrix_path
    String sample_metadata_path
    String output_dir
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c "mkdir -p '~{output_dir}/coexpression' '~{output_dir}/plots' && test -s '~{vst_count_matrix_path}' && test -s '~{sample_metadata_path}'"
    ~{docker_run} ~{image} Rscript /kegg_data/wgcna_analysis.R --vst-file '~{vst_count_matrix_path}' --output-dir '~{output_dir}/coexpression' --metadata-file '~{sample_metadata_path}'
  >>>

  output {
    String coexpression_results_dir = output_dir + "/coexpression"
    String gene_module_path = output_dir + "/coexpression/gene_modules.tsv"
    String module_sizes_path = output_dir + "/coexpression/module_sizes.tsv"
    String module_eigengenes_path = output_dir + "/coexpression/module_eigengenes.tsv"
    String hub_genes_path = output_dir + "/coexpression/hub_genes.tsv"
    String input_type_path = output_dir + "/coexpression/input_type.txt"
  }
}
task MultiQc {
  input {
    String output_dir
    String count_matrix_path
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c "test -s '~{count_matrix_path}'"
    ~{docker_run} ~{image} bash -c "mkdir -p '~{output_dir}/multiqc'"
    ~{docker_run} ~{image} bash -c "multiqc --outdir '~{output_dir}/multiqc' --filename multiqc_report.html '~{output_dir}'"
  >>>

  output {
    String report_path = output_dir + "/multiqc/multiqc_report.html"
  }
}

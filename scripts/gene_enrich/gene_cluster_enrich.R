#!/usr/bin/env Rscript

############################################################
# Gene-list enrichment pipeline
# Author: ChatGPT
# Version: gene-list stable offline KEGG
############################################################

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(enrichplot)
    library(gson)
    library(ReactomePA)
    library(AnnotationDbi)
    library(ggplot2)
})

usage <- function() {
    cat("Usage:\n\n",
        "Rscript gene_cluster_enrich.R \\\n",
        "    gene_list.txt \\\n",
        "    gene_id_type \\\n",
        "    organism \\\n",
        "    output_dir \\\n",
        "    prefix\n\n",
        "Arguments:\n",
        "  gene_list.txt   Gene list file. Gene IDs may be one per line, or separated by comma/tab/space.\n",
        "                  Lines beginning with # are ignored. If a header named gene/gene_id/id is present,\n",
        "                  that column is used; otherwise all tokens are parsed as gene IDs.\n",
        "  gene_id_type    Input gene ID type in the selected OrgDb, e.g. ENTREZID, SYMBOL, ENSEMBL, REFSEQ.\n",
        "  organism        Species: mouse/mmu or human/hsa. This selects OrgDb and KEGG files automatically.\n",
        "  output_dir      Output directory.\n",
        "  prefix          Output file prefix.\n\n",
        "Outputs:\n",
        "  output_dir/<prefix>_gene_id_mapping.csv       Input gene IDs mapped to ENTREZID/SYMBOL.\n",
        "  output_dir/<prefix>_GO.csv                    GO BP enrichment table with RichFactor.\n",
        "  output_dir/<prefix>_GO_classification.csv     GO level-2 classification table.\n",
        "  output_dir/<prefix>_KEGG.csv                  KEGG enrichment table with class metadata, pathway link, RichFactor, and symbol geneID.\n",
        "  output_dir/<prefix>_KEGG.html                 KEGG enrichment HTML table with clickable pathway links.\n",
        "  output_dir/<prefix>_Reactome.csv              Reactome pathway enrichment table with RichFactor.\n",
        "  output_dir/plots/<prefix>_GO_dotplot.pdf      GO enrichment bubble plot; x-axis is RichFactor.\n",
        "  output_dir/plots/<prefix>_GO_DAG.pdf          GO enrichment DAG plot for top significant GO terms.\n",
        "  output_dir/plots/<prefix>_GO_barplot.pdf      GO classification bar plot.\n",
        "  output_dir/plots/<prefix>_KEGG_dotplot.pdf    KEGG enrichment bubble plot; x-axis is RichFactor.\n",
        "  output_dir/plots/<prefix>_KEGG_barplot.pdf    KEGG enrichment bar plot.\n",
        "  output_dir/plots/<prefix>_Reactome_dotplot.pdf Reactome enrichment bubble plot; x-axis is RichFactor.\n",
        "  output_dir/plots/<prefix>_Reactome_barplot.pdf Reactome enrichment bar plot.\n",
        sep = "")
}

log_step <- function(...) {
    message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
}

validate_file <- function(path, label) {
    if (!file.exists(path)) {
        stop(label, " not found: ", path, call. = FALSE)
    }
    if (file.info(path)$size == 0) {
        stop(label, " is empty: ", path, call. = FALSE)
    }
    invisible(path)
}

write_result <- function(result_df, path) {
    write.csv(result_df, path, row.names = FALSE)
    log_step("Wrote ", nrow(result_df), " rows: ", path)
}

empty_kegg_class_df <- function() {
    data.frame(
        ID = character(0),
        category = character(0),
        subcategory = character(0),
        Description = character(0),
        stringsAsFactors = FALSE
    )
}

normalize_gson_term2name <- function(kk_gson) {
    term2name <- kk_gson@gsid2name

    if (is.data.frame(term2name) && ncol(term2name) >= 2) {
        return(data.frame(
            ID = as.character(term2name[[1]]),
            Description = as.character(term2name[[2]]),
            stringsAsFactors = FALSE
        ))
    }

    if (is.list(term2name) && !is.null(names(term2name))) {
        return(data.frame(
            ID = names(term2name),
            Description = as.character(unlist(term2name, use.names = FALSE)),
            stringsAsFactors = FALSE
        ))
    }

    if (is.atomic(term2name) && !is.null(names(term2name))) {
        return(data.frame(
            ID = names(term2name),
            Description = as.character(term2name),
            stringsAsFactors = FALSE
        ))
    }

    empty_kegg_class_df()[, c("ID", "Description")]
}

normalize_gson_term_ids <- function(kk_gson) {
    term2gene <- kk_gson@gsid2gene

    if (is.data.frame(term2gene) && ncol(term2gene) >= 1) {
        return(unique(as.character(term2gene[[1]])))
    }

    if (is.list(term2gene) && !is.null(names(term2gene))) {
        return(names(term2gene))
    }

    character(0)
}

make_kegg_class_from_gson <- function(kk_gson) {
    term_ids <- normalize_gson_term_ids(kk_gson)
    term2name <- normalize_gson_term2name(kk_gson)

    if (length(term_ids) == 0 && nrow(term2name) > 0) {
        term_ids <- term2name$ID
    }

    if (length(term_ids) == 0) {
        stop("Unable to infer KEGG pathway IDs from gson; please provide a non-empty kegg_class.tsv", call. = FALSE)
    }

    class_df <- data.frame(
        ID = unique(term_ids),
        category = "Unclassified",
        subcategory = "Unclassified",
        stringsAsFactors = FALSE
    )
    class_df <- merge(class_df, term2name, by = "ID", all.x = TRUE, sort = FALSE)
    missing_description <- is.na(class_df$Description) | !nzchar(class_df$Description)
    class_df$Description[missing_description] <- class_df$ID[missing_description]
    class_df[, c("ID", "category", "subcategory", "Description")]
}

prepare_kegg_class_df <- function(class_file, kk_gson) {
    validate_file(class_file, "KEGG class file")
    class_df <- read.delim(
        class_file,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    if (!"ID" %in% colnames(class_df)) {
        stop("KEGG class file must contain an ID column", call. = FALSE)
    }

    if (nrow(class_df) == 0) {
        warning(
            "KEGG class file contains 0 rows; using gson-derived descriptions with Unclassified categories.",
            call. = FALSE
        )
        class_df <- make_kegg_class_from_gson(kk_gson)
    }

    class_df
}

html_escape <- function(x) {
    x <- as.character(x)
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x <- gsub('"', "&quot;", x, fixed = TRUE)
    x
}

extract_kegg_pathway_id <- function(pathway_id, organism = "mmu") {
    pathway_id <- sub("^path:", "", pathway_id)
    pathway_id <- sub("^map", organism, pathway_id)

    if (grepl("^[0-9]{5}$", pathway_id)) {
        pathway_id <- paste0(organism, pathway_id)
    }

    pathway_id
}

split_gene_ids <- function(gene_ids) {
    gene_ids <- as.character(gene_ids)
    if (is.na(gene_ids) || !nzchar(gene_ids)) {
        return(character(0))
    }
    unique(strsplit(gene_ids, "/", fixed = TRUE)[[1]])
}

make_kegg_pathway_url <- function(pathway_id, gene_ids, organism = "mmu",
                                  bg_color = "#FF9999", fg_color = "#000000",
                                  max_url_length = 1900) {
    pathway_id <- extract_kegg_pathway_id(pathway_id, organism)
    base_url <- paste0("https://www.kegg.jp/kegg-bin/show_pathway?map=", pathway_id)
    genes <- split_gene_ids(gene_ids)
    genes <- genes[!is.na(genes) & nzchar(genes)]

    if (length(genes) == 0) {
        return(base_url)
    }

    genes <- ifelse(grepl(":", genes, fixed = TRUE), genes, paste0(organism, ":", genes))
    color_lines <- paste0(genes, " ", bg_color, ",", fg_color)
    color_query <- utils::URLencode(paste(color_lines, collapse = "\r\n"), reserved = TRUE)
    color_query <- gsub("%20", "+", color_query, fixed = TRUE)

    highlighted_url <- paste0(
        "https://www.kegg.jp/kegg-bin/show_pathway?map=",
        pathway_id,
        "&multi_query=",
        color_query
    )

    if (nchar(highlighted_url, type = "bytes") > max_url_length) {
        return(base_url)
    }

    highlighted_url
}

make_kegg_pathway_link_note <- function(pathway_id, gene_ids, organism = "mmu",
                                        max_url_length = 1900) {
    pathway_id <- extract_kegg_pathway_id(pathway_id, organism)
    genes <- split_gene_ids(gene_ids)
    genes <- genes[!is.na(genes) & nzchar(genes)]

    if (length(genes) == 0) {
        return("No genes available for highlighting")
    }

    genes <- ifelse(grepl(":", genes, fixed = TRUE), genes, paste0(organism, ":", genes))
    color_lines <- paste0(genes, " #FF9999,#000000")
    color_query <- utils::URLencode(paste(color_lines, collapse = "\r\n"), reserved = TRUE)
    color_query <- gsub("%20", "+", color_query, fixed = TRUE)
    highlighted_url <- paste0(
        "https://www.kegg.jp/kegg-bin/show_pathway?map=",
        pathway_id,
        "&multi_query=",
        color_query
    )

    if (nchar(highlighted_url, type = "bytes") > max_url_length) {
        return(paste0(
            "Gene list too long for a stable KEGG URL; ",
            "Pathway_Link opens the pathway without URL-based highlighting (",
            length(genes),
            " genes)."
        ))
    }

    "Highlighted KEGG URL"
}

add_kegg_pathway_links <- function(kegg_df, organism = "mmu") {
    if (nrow(kegg_df) == 0) {
        kegg_df$Pathway_Link <- character(0)
        kegg_df$Pathway_Link_Note <- character(0)
        return(kegg_df)
    }

    if (!all(c("ID", "geneID") %in% colnames(kegg_df))) {
        stop("KEGG result must contain ID and geneID columns to create pathway links", call. = FALSE)
    }

    kegg_df$Pathway_Link <- mapply(
        make_kegg_pathway_url,
        pathway_id = kegg_df$ID,
        gene_ids = kegg_df$geneID,
        MoreArgs = list(organism = organism),
        USE.NAMES = FALSE
    )
    kegg_df$Pathway_Link_Note <- mapply(
        make_kegg_pathway_link_note,
        pathway_id = kegg_df$ID,
        gene_ids = kegg_df$geneID,
        MoreArgs = list(organism = organism),
        USE.NAMES = FALSE
    )

    kegg_df
}

make_entrez_symbol_lookup <- function(gene_id_values, orgdb) {
    genes <- unique(unlist(lapply(gene_id_values, split_gene_ids), use.names = FALSE))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    genes <- unique(sub("^[^:]+:", "", genes))

    if (length(genes) == 0) {
        return(character(0))
    }

    symbols <- suppressMessages(AnnotationDbi::mapIds(
        orgdb,
        keys = genes,
        keytype = "ENTREZID",
        column = "SYMBOL",
        multiVals = "first"
    ))
    symbols <- unname(as.character(symbols))
    missing_symbols <- is.na(symbols) | !nzchar(symbols)
    symbols[missing_symbols] <- genes[missing_symbols]
    names(symbols) <- genes

    symbols
}

map_entrez_list_to_symbols <- function(gene_ids, symbol_lookup) {
    genes <- split_gene_ids(gene_ids)
    genes <- genes[!is.na(genes) & nzchar(genes)]

    if (length(genes) == 0) {
        return("")
    }

    genes_without_prefix <- sub("^[^:]+:", "", genes)
    symbols <- unname(symbol_lookup[genes_without_prefix])
    missing_symbols <- is.na(symbols) | !nzchar(symbols)
    symbols[missing_symbols] <- genes_without_prefix[missing_symbols]

    paste(unique(symbols), collapse = "/")
}

replace_kegg_gene_ids_with_symbols <- function(kegg_df, orgdb) {
    if (!"geneID" %in% colnames(kegg_df)) {
        return(kegg_df)
    }

    symbol_lookup <- make_entrez_symbol_lookup(kegg_df$geneID, orgdb)
    kegg_df$geneID <- vapply(
        kegg_df$geneID,
        map_entrez_list_to_symbols,
        character(1),
        symbol_lookup = symbol_lookup
    )

    kegg_df[, c(setdiff(colnames(kegg_df), "geneID"), "geneID"), drop = FALSE]
}

coalesce_text <- function(primary, fallback) {
    primary <- as.character(primary)
    fallback <- as.character(fallback)
    missing_primary <- is.na(primary) | !nzchar(primary)
    primary[missing_primary] <- fallback[missing_primary]
    primary
}

normalize_kegg_description_columns <- function(kegg_df) {
    if ("Description.y" %in% colnames(kegg_df) && "Description.x" %in% colnames(kegg_df)) {
        kegg_df$Description <- coalesce_text(kegg_df$Description.y, kegg_df$Description.x)
        kegg_df$Description.x <- NULL
        kegg_df$Description.y <- NULL
    } else if ("Description.y" %in% colnames(kegg_df)) {
        kegg_df$Description <- kegg_df$Description.y
        kegg_df$Description.y <- NULL
    } else if ("Description.x" %in% colnames(kegg_df)) {
        kegg_df$Description <- kegg_df$Description.x
        kegg_df$Description.x <- NULL
    }

    if (!"Description" %in% colnames(kegg_df)) {
        kegg_df$Description <- kegg_df$ID
    }

    description_missing <- is.na(kegg_df$Description) | !nzchar(kegg_df$Description)
    kegg_df$Description[description_missing] <- kegg_df$ID[description_missing]

    leading_cols <- intersect(c("ID", "category", "subcategory", "Description"), colnames(kegg_df))
    kegg_df[, c(leading_cols, setdiff(colnames(kegg_df), leading_cols)), drop = FALSE]
}

write_kegg_html_result <- function(kegg_df, path) {
    if (nrow(kegg_df) == 0) {
        writeLines(
            c(
                "<!doctype html>",
                "<html><head><meta charset=\"UTF-8\"><title>KEGG enrichment</title></head>",
                "<body><p>No KEGG enrichment terms found.</p></body></html>"
            ),
            path
        )
        log_step("Wrote KEGG HTML result: ", path)
        return(invisible(path))
    }

    display_df <- kegg_df
    link_col <- display_df$Pathway_Link
    display_df$Pathway_Link <- paste0(
        "<a href=\"", html_escape(link_col),
        "\" target=\"_blank\" rel=\"noopener noreferrer\">View highlighted pathway</a>"
    )

    header <- paste0("<th>", html_escape(colnames(display_df)), "</th>", collapse = "")
    rows <- apply(display_df, 1, function(row) {
        cells <- mapply(function(value, name) {
            if (name == "Pathway_Link") {
                paste0("<td>", value, "</td>")
            } else {
                paste0("<td>", html_escape(value), "</td>")
            }
        }, row, names(row), USE.NAMES = FALSE)
        paste0("<tr>", paste0(cells, collapse = ""), "</tr>")
    })

    html <- c(
        "<!doctype html>",
        "<html>",
        "<head>",
        "<meta charset=\"UTF-8\">",
        "<title>KEGG enrichment</title>",
        "<style>",
        "body{font-family:Arial,sans-serif;margin:24px;}",
        "table{border-collapse:collapse;font-size:13px;}",
        "th,td{border:1px solid #ddd;padding:6px 8px;vertical-align:top;}",
        "th{background:#f2f2f2;position:sticky;top:0;}",
        "tr:nth-child(even){background:#fafafa;}",
        "</style>",
        "</head>",
        "<body>",
        "<h1>KEGG enrichment result</h1>",
        "<p>Click <strong>Pathway_Link</strong> to open the KEGG pathway map. Short gene lists are highlighted on the KEGG interactive map; long gene lists fall back to the plain pathway URL to avoid browser/KEGG URL-length errors. Check <strong>Pathway_Link_Note</strong> for details.</p>",
        "<table>",
        paste0("<thead><tr>", header, "</tr></thead>"),
        "<tbody>",
        rows,
        "</tbody>",
        "</table>",
        "</body></html>"
    )

    writeLines(html, path)
    log_step("Wrote KEGG HTML result: ", path)
    invisible(path)
}

parse_ratio_numerator <- function(ratio) {
    ratio <- as.character(ratio)
    suppressWarnings(as.numeric(sub("/.*$", "", ratio)))
}

add_rich_factor <- function(result_df) {
    if (nrow(result_df) == 0) {
        result_df$RichFactor <- numeric(0)
        return(result_df)
    }

    if (all(c("Count", "BgRatio") %in% colnames(result_df))) {
        bg_count <- parse_ratio_numerator(result_df$BgRatio)
        result_df$RichFactor <- as.numeric(result_df$Count) / bg_count
    } else if ("GeneRatio" %in% colnames(result_df)) {
        result_df$RichFactor <- parse_ratio_numerator(result_df$GeneRatio)
    } else {
        result_df$RichFactor <- NA_real_
    }

    result_df
}

plot_enrichment <- function(enrich_result, result_df, plot_dir, prefix, label,
                            show_category = 15, width = 10, height = 8,
                            draw_barplot = TRUE) {
    if (nrow(result_df) == 0) {
        log_step("No ", label, " enrichment terms found; skipping plots")
        return(invisible(FALSE))
    }

    dotplot_file <- file.path(plot_dir, paste0(prefix, "_", label, "_dotplot.pdf"))
    barplot_file <- file.path(plot_dir, paste0(prefix, "_", label, "_barplot.pdf"))
    plot_df <- head(add_rich_factor(result_df), show_category)
    plot_df$Count <- as.numeric(plot_df$Count)
    plot_df$p.adjust <- as.numeric(plot_df$p.adjust)
    if (!"Description" %in% colnames(plot_df)) {
        if ("Description.y" %in% colnames(plot_df)) {
            plot_df$Description <- plot_df$Description.y
        } else if ("Description.x" %in% colnames(plot_df)) {
            plot_df$Description <- plot_df$Description.x
        } else if ("ID" %in% colnames(plot_df)) {
            plot_df$Description <- plot_df$ID
        } else {
            plot_df$Description <- seq_len(nrow(plot_df))
        }
    }
    plot_df$Description <- as.character(plot_df$Description)
    missing_description <- is.na(plot_df$Description) | !nzchar(plot_df$Description)
    if ("ID" %in% colnames(plot_df)) {
        plot_df$Description[missing_description] <- plot_df$ID[missing_description]
    }
    plot_df$Description <- factor(
        plot_df$Description,
        levels = rev(unique(plot_df$Description))
    )

    log_step("Plotting ", label, " dotplot: ", dotplot_file)
    pdf(dotplot_file, width = width, height = height)
    print(
        ggplot(plot_df, aes(x = RichFactor, y = Description)) +
            geom_point(aes(size = Count, color = p.adjust)) +
            scale_color_continuous(low = "red", high = "blue", name = "p.adjust") +
            labs(x = "RichFactor", y = NULL, size = "Count") +
            theme_bw(base_size = 11)
    )
    dev.off()

    if (draw_barplot) {
        log_step("Plotting ", label, " barplot: ", barplot_file)
        pdf(barplot_file, width = width, height = height)
        print(barplot(enrich_result, showCategory = show_category))
        dev.off()
    }

    invisible(TRUE)
}

plot_go_dag <- function(go_enrich_result, result_df, plot_dir, prefix,
                        show_category = 10, width = 12, height = 10) {
    if (nrow(result_df) == 0) {
        log_step("No GO enrichment terms found; skipping GO DAG plot")
        return(invisible(FALSE))
    }

    dag_file <- file.path(plot_dir, paste0(prefix, "_GO_DAG.pdf"))
    show_category <- min(show_category, nrow(result_df))
    log_step("Plotting GO DAG top ", show_category, " terms: ", dag_file)

    pdf(dag_file, width = width, height = height)
    print(enrichplot::goplot(go_enrich_result, showCategory = show_category))
    dev.off()

    invisible(TRUE)
}

go_ontology_labels <- c(
    BP = "biological_process",
    CC = "cellular_component",
    MF = "molecular_function"
)

go_ontology_colors <- c(
    biological_process = "#F8766D",
    cellular_component = "#00BA38",
    molecular_function = "#619CFF"
)

read_gene_list <- function(gene_file) {
    lines <- readLines(gene_file, warn = FALSE)
    lines <- trimws(lines)
    lines <- lines[nzchar(lines) & !startsWith(lines, "#")]

    if (length(lines) == 0) {
        stop("Gene list file contains no usable gene IDs: ", gene_file, call. = FALSE)
    }

    table_candidate <- tryCatch(
        read.delim(
            gene_file,
            header = TRUE,
            sep = "\t",
            stringsAsFactors = FALSE,
            check.names = FALSE,
            comment.char = "#"
        ),
        error = function(e) NULL
    )

    if (!is.null(table_candidate) && nrow(table_candidate) > 0) {
        header_names <- tolower(colnames(table_candidate))
        gene_col <- match(TRUE, header_names %in% c("gene", "gene_id", "geneid", "id"))
        if (!is.na(gene_col)) {
            genes <- table_candidate[[gene_col]]
            genes <- trimws(as.character(genes))
            return(unique(genes[!is.na(genes) & nzchar(genes)]))
        }
    }

    genes <- unlist(strsplit(lines, "[,[:space:]]+", perl = TRUE), use.names = FALSE)
    genes <- trimws(genes)
    unique(genes[!is.na(genes) & nzchar(genes)])
}

get_organism_config <- function(organism) {
    organism <- tolower(organism)

    if (organism %in% c("mouse", "mmu", "mus_musculus", "mus musculus")) {
        return(list(
            label = "mouse",
            kegg_code = "mmu",
            orgdb_package = "org.Mm.eg.db",
            orgdb_object = "org.Mm.eg.db",
            reactome_organism = "mouse",
            gson_file = "/kegg_data/kegg_mmu.gson",
            class_file = "/kegg_data/kegg_class_mmu.tsv"
        ))
    }

    if (organism %in% c("human", "hsa", "homo_sapiens", "homo sapiens")) {
        return(list(
            label = "human",
            kegg_code = "hsa",
            orgdb_package = "org.Hs.eg.db",
            orgdb_object = "org.Hs.eg.db",
            reactome_organism = "human",
            gson_file = "/kegg_data/kegg_hsa.gson",
            class_file = "/kegg_data/kegg_class_hsa.tsv"
        ))
    }

    stop("Unsupported organism: ", organism, ". Use mouse/mmu or human/hsa.", call. = FALSE)
}

load_orgdb <- function(config) {
    if (!requireNamespace(config$orgdb_package, quietly = TRUE)) {
        stop("Required OrgDb package is not installed: ", config$orgdb_package, call. = FALSE)
    }

    get(config$orgdb_object, envir = asNamespace(config$orgdb_package))
}

validate_gene_id_type <- function(gene_id_type, orgdb, orgdb_package) {
    gene_id_type <- toupper(gene_id_type)
    valid_keytypes <- AnnotationDbi::keytypes(orgdb)

    if (!gene_id_type %in% valid_keytypes) {
        stop(
            "Unsupported gene_id_type: ", gene_id_type,
            ". Supported ", orgdb_package, " key types include: ",
            paste(valid_keytypes, collapse = ", "),
            call. = FALSE
        )
    }

    gene_id_type
}

convert_genes_to_entrez <- function(input_genes, gene_id_type, orgdb) {
    input_genes <- unique(as.character(input_genes))
    input_genes <- input_genes[!is.na(input_genes) & nzchar(input_genes)]

    if (gene_id_type == "ENTREZID") {
        return(data.frame(
            input_gene_id = input_genes,
            gene_id_type = gene_id_type,
            ENTREZID = input_genes,
            SYMBOL = suppressMessages(AnnotationDbi::mapIds(
                orgdb,
                keys = input_genes,
                keytype = "ENTREZID",
                column = "SYMBOL",
                multiVals = "first"
            )),
            stringsAsFactors = FALSE,
            row.names = NULL
        ))
    }

    mapping <- suppressMessages(AnnotationDbi::select(
        orgdb,
        keys = input_genes,
        keytype = gene_id_type,
        columns = c("ENTREZID", "SYMBOL")
    ))
    colnames(mapping)[colnames(mapping) == gene_id_type] <- "input_gene_id"
    mapping$gene_id_type <- gene_id_type
    mapping <- mapping[, c("input_gene_id", "gene_id_type", "ENTREZID", "SYMBOL")]
    mapping <- mapping[!is.na(mapping$ENTREZID) & nzchar(mapping$ENTREZID), ]
    mapping <- unique(mapping)
    row.names(mapping) <- NULL
    mapping
}

run_go_classification <- function(genes, orgdb, level = 2) {
    go_class_list <- lapply(names(go_ontology_labels), function(ontology) {
        go_class <- groupGO(
            gene = genes,
            OrgDb = orgdb,
            keyType = "ENTREZID",
            ont = ontology,
            level = level,
            readable = TRUE
        )
        go_class_df <- as.data.frame(go_class)
        go_class_df$ONTOLOGY <- go_ontology_labels[[ontology]]
        go_class_df
    })

    go_class_df <- do.call(rbind, go_class_list)
    go_class_df$Count <- as.numeric(go_class_df$Count)
    go_class_df <- go_class_df[!is.na(go_class_df$Count) & go_class_df$Count > 0, ]
    row.names(go_class_df) <- NULL

    go_class_df
}

plot_go_classification_bar <- function(go_class_df, plot_dir, prefix,
                                       width = 10, height = 9) {
    if (nrow(go_class_df) == 0) {
        log_step("No GO classification terms found; skipping GO classification barplot")
        return(invisible(FALSE))
    }

    go_class_df$ONTOLOGY <- factor(
        go_class_df$ONTOLOGY,
        levels = unname(go_ontology_labels)
    )
    go_class_df$Description <- factor(
        go_class_df$Description,
        levels = rev(unique(go_class_df$Description))
    )

    barplot_file <- file.path(plot_dir, paste0(prefix, "_GO_barplot.pdf"))
    log_step("Plotting GO classification barplot: ", barplot_file)

    pdf(barplot_file, width = width, height = height)
    print(
        ggplot(go_class_df, aes(x = Count, y = Description, fill = ONTOLOGY)) +
            geom_col(width = 0.55) +
            facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
            scale_fill_manual(values = go_ontology_colors, guide = "none") +
            labs(x = "Number", y = "Description", caption = "GO 分类注释结果") +
            theme_gray(base_size = 11) +
            theme(
                plot.caption = element_text(hjust = 0.5, size = 14, margin = margin(t = 18)),
                axis.text.y = element_text(size = 8),
                strip.text.y = element_text(size = 9),
                panel.spacing.y = grid::unit(0.15, "lines")
            )
    )
    dev.off()

    invisible(TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
    usage()
    quit(save = "no", status = 1)
}

gene_file <- args[1]
organism_config <- get_organism_config(args[3])
orgdb <- load_orgdb(organism_config)
gene_id_type <- validate_gene_id_type(args[2], orgdb, organism_config$orgdb_package)
gson_file <- organism_config$gson_file
class_file <- organism_config$class_file
outdir <- args[4]
prefix <- args[5]

validate_file(gene_file, "Gene list file")
validate_file(gson_file, "KEGG gson file")
validate_file(class_file, "KEGG class file")
log_step("Organism: ", organism_config$label, " (KEGG ", organism_config$kegg_code, ", OrgDb ", organism_config$orgdb_package, ")")
log_step("Using KEGG gson: ", gson_file)
log_step("Using KEGG class file: ", class_file)

if (!nzchar(prefix)) {
    stop("prefix must not be empty", call. = FALSE)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(outdir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

log_step("Reading genes: ", gene_file)
input_genes <- read_gene_list(gene_file)
log_step("Input genes parsed: ", length(input_genes))

if (length(input_genes) == 0) {
    stop("No genes found in gene list", call. = FALSE)
}

log_step("Converting ", gene_id_type, " genes to ENTREZID")
gene_mapping <- convert_genes_to_entrez(input_genes, gene_id_type, orgdb)
write_result(gene_mapping, file.path(outdir, paste0(prefix, "_gene_id_mapping.csv")))

genes <- unique(gene_mapping$ENTREZID)
genes <- genes[!is.na(genes) & nzchar(genes)]
log_step("Unique ENTREZID genes: ", length(genes))

if (length(genes) == 0) {
    stop("No input genes could be converted to ENTREZID", call. = FALSE)
}

log_step("Running GO BP enrichment")
ego <- enrichGO(
    gene = genes,
    OrgDb = orgdb,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
)
go_df <- add_rich_factor(as.data.frame(ego))
write_result(go_df, file.path(outdir, paste0(prefix, "_GO.csv")))

log_step("Running GO classification")
go_class_df <- run_go_classification(genes, orgdb, level = 2)
write_result(
    go_class_df,
    file.path(outdir, paste0(prefix, "_GO_classification.csv"))
)

log_step("Running offline KEGG enrichment")
kk_gson <- read.gson(gson_file)
kegg <- enricher(
    gene = genes,
    TERM2GENE = kk_gson@gsid2gene,
    TERM2NAME = kk_gson@gsid2name,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
)
kegg_df <- add_rich_factor(as.data.frame(kegg))

log_step("Adding KEGG categories")
class_df <- prepare_kegg_class_df(
    class_file = class_file,
    kk_gson = kk_gson
)

kegg_df$.enrich_order <- seq_len(nrow(kegg_df))
kegg_df <- merge(class_df, kegg_df, by = "ID", all.y = TRUE, sort = FALSE)
kegg_df <- kegg_df[order(kegg_df$.enrich_order), , drop = FALSE]
kegg_df$.enrich_order <- NULL
kegg_df <- normalize_kegg_description_columns(kegg_df)
kegg_df <- add_kegg_pathway_links(kegg_df, organism = organism_config$kegg_code)
kegg_df <- replace_kegg_gene_ids_with_symbols(kegg_df, orgdb)

write_result(kegg_df, file.path(outdir, paste0(prefix, "_KEGG.csv")))
write_kegg_html_result(
    kegg_df,
    file.path(outdir, paste0(prefix, "_KEGG.html"))
)

log_step("Running Reactome pathway enrichment")
reactome <- ReactomePA::enrichPathway(
    gene = genes,
    organism = organism_config$reactome_organism,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
)
reactome_df <- add_rich_factor(as.data.frame(reactome))
write_result(reactome_df, file.path(outdir, paste0(prefix, "_Reactome.csv")))

plot_enrichment(ego, go_df, plot_dir, prefix, "GO", draw_barplot = FALSE)
plot_go_dag(ego, go_df, plot_dir, prefix)
plot_go_classification_bar(go_class_df, plot_dir, prefix)
plot_enrichment(kegg, kegg_df, plot_dir, prefix, "KEGG")
plot_enrichment(reactome, reactome_df, plot_dir, prefix, "Reactome")

log_step("================================")
log_step("Analysis finished")
log_step("================================")

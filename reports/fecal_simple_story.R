library(tidyverse)
library(data.table)
library(seqinr)


motifs_scored <- fread("fecal_simple/binnary/motifs-scored-read-methylation.tsv")

include <- fread("fecal_simple/binnary/include_contigs.tsv") %>%
    select(bin, contig) %>%
    mutate(
        included = TRUE
    )

assembly <- read.fasta("data/real_communities/fecal_simple/mmlong2_lite/results/mmlong2_lite_assembly.fasta", seqtype = "DNA")

assembly_info <- data.frame(
    contig = names(assembly),
    length = sapply(assembly, function(seq) {
        length(seq)
    }),
    gc = sapply(assembly, function(seq) {
        seqinr::GC(seq)
    })
)
rm(assembly)

cov <- fread("data/real_communities/fecal_simple/mmlong2_lite/tmp/binning/mapping/cov_all.tsv") %>%
    rename(
        contig = 1,
        length = 2,
        read_coverage_mean = 3
    ) %>%
    select(contig, length, read_coverage_mean)

contig_bin <- fread("data/real_communities/fecal_simple/mmlong2_lite/tmp/binning/contig_bin.tsv") %>%
    rename(
        contig = 1,
        bin = 2
    ) %>%
    select(contig, bin) %>%
    mutate(
        included = FALSE
    ) %>%
    bind_rows(include)

genomad <- fread("data/real_communities/fecal_simple/genomad/mmlong2_lite_assembly_aggregated_classification/mmlong2_lite_assembly_aggregated_classification.tsv") %>%
    rename(
        contig = 1
    ) %>%
    pivot_longer(cols = -contig, names_to = "mge", values_to = "count") %>%
    group_by(contig) %>%
    filter(count == max(count)) %>%
    select(-count) %>%
    mutate(
        mge = str_remove(mge, "_score")
    )

motifs_scored_f <- motifs_scored %>%
    filter(N_motif_obs >= 3) %>%
    left_join(contig_bin) %>%
    left_join(assembly_info %>% select(contig, length)) %>%
    left_join(genomad) %>%
    mutate(
        motif_mod = paste0(motif, "_", mod_type, "_", mod_position),
        contig_l = paste0(contig, " - ", length, " - ", mge)
    ) %>%
    mutate(bin = ifelse(is.na(bin), "unbinned", bin)) %>%
    mutate(
        included = ifelse(is.na(included), FALSE, included)
    )


gc_cov <- assembly_info %>%
    left_join(cov) %>%
    left_join(contig_bin) %>%
    left_join(genomad) %>%
    mutate(
        bin = ifelse(is.na(bin), "unbinned", bin)
    ) %>%
    mutate(
        included = ifelse(is.na(included), FALSE, included)
    )

BIN <- "mmlong2_lite.bin.1.6"

bin_gc_cov <- gc_cov %>%
    filter(bin == BIN) %>%
    arrange(desc(read_coverage_mean))

contigs_of_interest <- c("contig_349", "contig_251", "contig_916")

contig_gc_cov <- gc_cov %>%
    filter(contig %in% contigs_of_interest)

gc_cov %>%
    filter(bin != BIN) %>%
    ggplot(aes(x = gc, y = read_coverage_mean, shape = mge, color = included, size = length)) +
    geom_point(color = "gray50", alpha = 0.6) +
    geom_point(data = bin_gc_cov, alpha = 0.8) +
    geom_point(data = contig_gc_cov, alpha = 1, color = "#e9c600") +
    theme_minimal() +
    labs(
        x = "GC content",
        y = "Read coverage mean",
        color = "Bin",
        shape = "Included"
    ) +
    scale_y_log10()


motifs_scored_f %>%
    filter(bin == BIN) %>%
    ggplot(aes(x = motif_mod, y = contig_l, fill = median)) +
    geom_tile(color = "gray50")


bin_included <- motifs_scored_f %>%
    filter(included, bin == BIN) %>%
    mutate(
        contig_l = case_when(
            contig %in% contigs_of_interest ~ paste0("<span style='color:#e9c600'>", contig_l, "</span>"),
            TRUE ~ contig_l
        )
    ) %>%
    select(-bin, median, contig_l, motif_mod, included) %>%
    rename(
        bin = contig_l,
        mean = median
    )

motifs_scored_f %>%
    filter(!included, bin != "unbinned") %>%
    group_by(bin, motif_mod, included) %>%
    summarise(
        mean = mean(median),
    ) %>%
    bind_rows(bin_included) %>%
    mutate(
        included = ifelse(included, paste0("Added to ", BIN), "Bin")
    ) %>%
    ggplot(aes(x = motif_mod, y = bin, fill = mean)) +
    geom_tile(color = "gray50") +
    scale_fill_gradient(low = "white", high = "blue") +
    facet_grid(included ~ ., scales = "free", space = "free") +
    labs(
        x = "",
        y = "",
        fill = "Methylation degree"
    ) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_markdown()
    )

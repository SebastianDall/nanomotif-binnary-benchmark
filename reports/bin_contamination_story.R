library(tidyverse)
library(data.table)
library(seqinr)
library(ggtext)


motifs_scored <- fread("data/real_communities/ZymoFecal/nanomotif_binnary/motifs-scored-read-methylation.tsv")

contamination <- fread("data/real_communities/ZymoFecal/nanomotif_binnary/bin_contamination.tsv")

checkm <- fread("data/real_communities/ZymoFecal/mmlong2_lite/results/mmlong2_lite_bins.tsv") %>%
    rename(
        bin = 1,
        checkm_comp = 2,
        checkm_cont = 3
    ) %>%
    select(bin, checkm_comp, checkm_cont) %>%
    mutate(
        bin = str_remove(bin, "mmlong2_lite.")
    ) %>%
    mutate(
        qc = case_when(
            checkm_comp >= 90 & checkm_cont <= 5 ~ "HQ",
            checkm_comp >= 50 & checkm_cont <= 10 ~ "MQ",
            TRUE ~ "Low"
        )
    )
# assembly <- read.fasta("data/real_communities/ZymoFecal/mmlong2_lite/results/mmlong2_lite_assembly.fasta", seqtype = "DNA")

# assembly_info <- data.frame(
#     contig = names(assembly),
#     length = sapply(assembly, function(seq) {
#         length(seq)
#     }),
#     gc = sapply(assembly, function(seq) {
#         seqinr::GC(seq)
#     })
# )
# rm(assembly)

# assembly_info <- assembly_info %>% as_tibble()
# write_delim(assembly_info, "reports/ZymoFecal_assembly_info.txt", delim = "\t")
assembly_info <- fread("reports/ZymoFecal_assembly_info.txt")

cov <- fread("data/real_communities/ZymoFecal/mmlong2_lite/tmp/binning/mapping/cov_all.tsv") %>%
    rename(
        contig = 1,
        length = 2,
        read_coverage_mean = 3
    ) %>%
    select(contig, length, read_coverage_mean)

contig_bin <- fread("data/real_communities/ZymoFecal/mmlong2_lite/tmp/binning/contig_bin.tsv", header = FALSE) %>%
    rename(
        contig = V1,
        bin = V2
    ) %>%
    mutate(
        contamination = case_when(
            contig %in% contamination$contig ~ TRUE,
            TRUE ~ FALSE
        ),
        bin = str_remove(bin, "mmlong2_lite.")
    ) %>%
    left_join(checkm, by = "bin")


gc_coverage <- cov %>%
    left_join(assembly_info %>% select(contig, gc), by = "contig") %>%
    left_join(contig_bin, by = "contig")


MOD_TYPE_PRETTY <- c(
    a = "m6",
    m = "m5",
    `21839` = "m4"
)

motifs_scored_f <- motifs_scored %>%
    filter(N_motif_obs * mean_read_cov >= 17) %>%
    left_join(contig_bin, by = "contig") %>%
    left_join(assembly_info, by = "contig") %>%
    mutate(
        bin = ifelse(is.na(bin), "Unbinned", bin),
        contamination = ifelse(is.na(contamination), FALSE, contamination)
    ) %>%
    mutate(
        contig_label = paste0(contig, " - ", length),
        contig_label = case_when(
            contamination ~ paste0("<span style='color:red'><strong>", contig_label, "</stron></span>"),
            TRUE ~ contig_label
        ),
        motif_mod = paste0(motif, "_", mod_type, "_", mod_position)
    )




contaminated_bins <- contamination %>%
    arrange(bin) %>%
    mutate(bin = str_remove(bin, "mmlong2_lite.")) %>%
    pull(bin) %>%
    unique()


plot_motifs <- function(motifs_scored_f, bin2see, sample_motifs = TRUE) {
    df <- motifs_scored_f %>%
        filter(bin == bin2see)

    if (sample_motifs) {
        motifs_high <- df %>%
            filter(median > 0.15) %>%
            pull(motif_mod) %>%
            unique()

        motifs_low <- df %>%
            filter(median < 0.15) %>%
            pull(motif_mod) %>%
            unique()

        # sample 10 motifs from motifs_low
        motifs_low <- sample(motifs_low, 10)

        motifs <- c(motifs_high, motifs_low)
    } else {
        motifs <- df$motif_mod %>% unique()
    }

    df_cross <- crossing(contig_label = df$contig_label, motif_mod = df$motif_mod %>% unique()) %>%
        left_join(df, by = c("contig_label", "motif_mod")) %>%
        mutate(
            txt = case_when(
                is.na(median) ~ "\u00B7",
                TRUE ~ NA
            ),
            motif = case_when(
                is.na(motif) ~ str_split(motif_mod, "_") %>% map_chr(1),
                TRUE ~ motif
            ),
            mod_type = case_when(
                is.na(mod_type) ~ str_split(motif_mod, "_") %>% map_chr(2),
                TRUE ~ mod_type
            ),
            mod_position = case_when(
                is.na(mod_position) ~ str_split(motif_mod, "_") %>% map_chr(3) %>% as.integer(),
                TRUE ~ mod_position
            ),
            median = ifelse(is.na(median), 0, median)
        )

    bin_stats <- tibble(
        bin = bin2see,
        comp = df$checkm_comp[1],
        cont = df$checkm_cont[1],
        qc = df$qc[1]
    )

    df_cross %>%
        filter(motif_mod %in% motifs) %>%
        mutate(
            ns = str_extract(motif, "NNN*"),
            ns = str_length(ns),
            motif = str_replace(motif, "NNN*", paste0("(N)<sub>", ns, "</sub>")),
            motif_axis = paste0(
                str_sub(motif, 1, mod_position),
                "<strong>",
                str_sub(motif, mod_position + 1, mod_position + 1),
                "<sub>",
                map(mod_type, function(x) MOD_TYPE_PRETTY[[x]]) %>% unlist(),
                "</sup>",
                "</strong>",
                str_sub(motif, mod_position + 2, -1)
            )
        ) %>%
        ggplot(aes(x = motif_axis, y = contig_label, fill = median, label = txt)) +
        geom_tile(color = "gray60", ) +
        geom_text(aes(label = txt), size = 15) +
        scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
        labs(
            x = "",
            y = "",
            title = paste0(bin_stats$qc, " ", bin2see, " - Comp. ", bin_stats$comp, " - Cont. ", bin_stats$cont)
        ) +
        coord_fixed() +
        theme_minimal() +
        theme(
            axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_markdown()
        )
}


plot_gc_coverage <- function(gc_coverage, bin2see) {
    bin_gc_coverage <- gc_coverage %>%
        filter(bin == bin2see)

    bin_stats <- tibble(
        bin = bin2see,
        comp = bin_gc_coverage$checkm_comp[1],
        cont = bin_gc_coverage$checkm_cont[1],
        qc = bin_gc_coverage$qc[1]
    )

    gc_coverage %>%
        filter(length > 20000) %>%
        filter(read_coverage_mean > 10) %>%
        ggplot(aes(x = gc, y = read_coverage_mean, size = length)) +
        geom_point(alpha = 0.5, color = "#6e6e6e") +
        geom_point(data = bin_gc_coverage, aes(color = contamination)) +
        labs(
            x = "GC content (%)",
            y = "coverage",
            title = paste0(bin_stats$qc, " ", bin2see, " - Comp. ", bin_stats$comp, " - Cont. ", bin_stats$cont)
        ) +
        scale_y_log10() +
        theme_minimal()
}

plot_motifs(motifs_scored_f, "bin.2.129")
plot_gc_coverage(gc_coverage, "bin.2.129")

for (BIN in contaminated_bins) {
    mag_qual <- checkm %>%
        filter(bin == BIN) %>%
        pull(qc)

    plot_motifs(motifs_scored_f, BIN)
    ggsave(paste0("analysis/contamination_story/all_bins/", mag_qual, "_", BIN, "_motifs.png"), width = 10, height = 10, dpi = 300, create.dir = TRUE)

    plot_gc_coverage(gc_coverage, BIN)
    ggsave(paste0("analysis/contamination_story/all_bins/", mag_qual, "_", BIN, "_gc_coverage.png"), width = 10, height = 10, dpi = 300)
}

motifs <- c(
    "AACCNNNNNCTC_a_1",
    "AAAANGC_a_3",
    "CCWGG_21839_1",
    "CCAGG_m_1",
    "CTCGAG_a_4",
    "GAGNNNNNGGTY_a_1"
)
plot_motifs(motifs_scored_f %>% filter(motif_mod %in% motifs), "bin.1.13", sample_motifs = FALSE) +
    labs(
        title = "ZymoFecal",
        subtitle = "HQ bin.1.13 - Comp. 98.62 - Cont 2.02"
    )
ggsave("analysis/contamination_story/motifs_selected.png", width = 5, height = 5, dpi = 300)

library(data.table)
library(tidyverse)
library(ggtext)
source("src/colors.R")



samples <- c("fecal_simple_recovered", "ZymoFecal", "fecal_inhouse", "anaerobic_digester", "mfd02199_backup")
labels <- c("Fecal Simple", "ZymoFecal", "Fecal Complex", "Anaerobic Digester", "Soil")

df <- tibble(
    sample = samples,
    label = labels
) %>%
    mutate(
        before_file = map(sample, ~ list.files(paste0("data/real_communities/", .x, "/mmlong2_lite/results"), full.names = TRUE, pattern = "*bins.tsv")),
    ) %>%
    mutate(
        contamination = map(.x = sample, ~ fread(paste0("data/real_communities/", .x, "/nanomotif_binnary/bin_contamination_std.tsv"))),
        before = map(before_file, ~ fread(.x)),
        after = map(sample, ~ fread(paste0("data/real_communities/", .x, "/nanomotif_binnary/checkm2_decon_std/quality_report.tsv"))),
        coverage = map(sample, ~ fread(paste0("data/real_communities/", .x, "/mmlong2_lite/tmp/binning/bin_cov.tsv"))),
    ) %>%
    select(-sample) %>%
    rename(
        sample = label
    )




df_cov_soil <- df %>%
    filter(sample == "Soil") %>%
    select(sample, coverage) %>%
    unnest(coverage) %>%
    rename(
        sample = sample,
        bin = Genome,
        coverage = 3
    ) %>%
    select(sample:coverage)
df_coverage <- df %>%
    filter(sample != "Soil") %>%
    select(sample, coverage) %>%
    unnest(coverage) %>%
    rename(
        sample = sample,
        bin = Genome,
        coverage = 3
    ) %>%
    select(sample:coverage) %>%
    bind_rows(df_cov_soil)



df_contamination <- df %>%
    # filter(sample != "Fecal Simple") %>%
    select(sample, contamination) %>%
    mutate(
        tib_len = map_int(contamination, ~ nrow(.x))
    ) %>%
    filter(tib_len > 0) %>%
    select(-tib_len) %>%
    unnest(contamination)


df_before <- df %>%
    filter(sample != "Soil") %>%
    unnest(before) %>%
    select(sample, bin, completeness_checkm2, contamination_checkm2, sample) %>%
    bind_rows(
        df %>%
            filter(sample == "Soil") %>%
            unnest(before) %>%
            select(sample, bin, Completeness, Contamination) %>%
            rename(
                completeness_checkm2 = Completeness,
                contamination_checkm2 = Contamination
            )
    ) %>%
    mutate(
        mag_qual = case_when(
            completeness_checkm2 >= 90 & contamination_checkm2 < 5 ~ "HQ",
            completeness_checkm2 >= 50 & contamination_checkm2 < 10 ~ "MQ",
            TRUE ~ "LQ"
        ),
        decon = "before"
    ) %>%
    rename(
        completeness = completeness_checkm2,
        contamination = contamination_checkm2
    )


df_after <- df %>%
    unnest(after) %>%
    select(sample, Name, Completeness, Contamination) %>%
    rename(
        bin = Name,
        completeness = Completeness,
        contamination = Contamination
    ) %>%
    mutate(
        mag_qual = case_when(
            completeness >= 90 & contamination < 5 ~ "HQ",
            completeness >= 50 & contamination < 10 ~ "MQ",
            TRUE ~ "LQ"
        ),
        decon = "after"
    )

df_com <- df_contamination %>%
    distinct(sample, bin) %>%
    left_join(df_coverage) %>%
    left_join(bind_rows(df_before, df_after))

n_contaminants <- df_contamination %>%
    group_by(sample) %>%
    summarise(
        n_bins = n_distinct(bin),
        n_contaminants = n_distinct(contig)
    ) %>%
    mutate(
        label = paste0("N bins: ", n_bins, "\nN contaminants: ", n_contaminants)
    )

df_com %>% filter(sample == "Soil")

df_com %>%
    group_by(sample, bin) %>%
    # filter(mag_qual != "LQ") %>%
    filter(n_distinct(decon) == 2) %>%
    ggplot(aes(x = completeness, y = contamination, color = coverage, shape = decon, group = bin)) +
    geom_rect(xmin = 90, xmax = 100, ymin = 0, ymax = 5, fill = NA, color = "gray30", inherit.aes = F, lty = "dashed") +
    geom_line(color = "gray60") +
    geom_point(size = 3) +
    geom_label(data = n_contaminants, aes(label = label), x = 52, y = 8.5, size = 3, hjust = 0, vjust = 0, inherit.aes = F) +
    facet_wrap(~sample) +
    guides(shape = guide_legend(title = "Decontamination"), color = guide_legend(title = "Coverage (log10)")) +
    scale_color_gradient(trans = "log10", breaks = c(10, 100, 300), labels = c("10", "100", "300")) +
    theme_minimal()

ggsave("analysis/contamination_story/decontamination_before_after_nostd_w_simple.png", width = 7, height = 7, dpi = 300)


df_com %>%
    group_by(sample, bin) %>%
    filter(mag_qual != "LQ") %>%
    filter(n_distinct(decon) == 2) %>%
    summarise(
        n_mag_qual = n_distinct(mag_qual)
    ) %>%
    filter(n_mag_qual == 2) %>%
    group_by(sample) %>%
    summarise(
        n_bins = n_distinct(bin)
    )

df_com %>% filter(sample == "Soil")

sample_label <- tibble(
    label = samples,
    sample = labels
)

changed_bin_status <- df_com %>%
    group_by(sample, bin) %>%
    filter(mag_qual != "LQ") %>%
    filter(n_distinct(decon) == 2) %>%
    summarise(
        n_mag_qual = n_distinct(mag_qual)
    ) %>%
    filter(n_mag_qual == 2) %>%
    left_join(sample_label) %>%
    ungroup() %>%
    nest(.by = label)

changed_bin_status <- df_com %>%
    group_by(sample, bin) %>%
    # filter(mag_qual != "LQ") %>%
    filter(n_distinct(decon) == 2) %>%
    summarise(
        n_mag_qual = n_distinct(mag_qual)
    ) %>%
    left_join(sample_label) %>%
    ungroup() %>%
    nest(.by = label)

for (s in changed_bin_status$label) {
    motifs_scored <- read_delim(file.path("data/real_communities/", s, "nanomotif_binnary/motifs-scored-read-methylation.tsv"), "\t")
    contig_bin <- read_delim(file.path("data/real_communities/", s, "mmlong2_lite/tmp/binning/contig_bin.tsv"), "\t", col_names = c("contig", "bin"))

    bins <- changed_bin_status %>%
        filter(label == s) %>%
        unnest(data) %>%
        pull(bin)

    sample_name <- changed_bin_status %>%
        filter(label == s) %>%
        unnest(data)
    sample_name <- sample_name$sample[1]

    contamination <- df_contamination %>%
        filter(sample == sample_name)

    for (b in bins) {
        contigs <- contig_bin %>%
            filter(bin == b) %>%
            pull(contig)

        contaminants <- contamination %>%
            filter(bin == b) %>%
            pull(contig)

        motifs_scored_f <- motifs_scored %>%
            filter(contig %in% contigs) %>%
            filter((N_motif_obs * mean_read_cov) >= 24) %>%
            mutate(
                ns = str_extract(motif, "NNN*"),
                ns = str_length(ns),
                motif_axis = paste0(
                    str_sub(motif, 1, mod_position),
                    "<strong>",
                    str_sub(motif, mod_position + 1, mod_position + 1),
                    "<sub>",
                    map(mod_type, function(x) MOD_TYPE_PRETTY[[x]]) %>% unlist(),
                    "</sup>",
                    "</strong>",
                    str_sub(motif, mod_position + 2, -1)
                ),
                motif_mod = str_replace(motif_axis, "NNN*", paste0("(N)<sub>", ns, "</sub>")),
                contig_l = case_when(
                    contig %in% contaminants ~ paste0("<span style='color:red'>", contig, "</span>"),
                    TRUE ~ contig
                )
            )

        motifs <- motifs_scored_f %>%
            filter(!is.na(median)) %>%
            filter(median >= 0.1) %>%
            pull(motif_mod) %>%
            unique()

        df_wide <- motifs_scored_f %>%
            filter(motif_mod %in% motifs) %>%
            select(contig_l, motif_mod, median) %>%
            pivot_wider(names_from = motif_mod, values_from = median, values_fill = list(median = 0))

        # Convert to matrix
        data_matrix <- as.matrix(df_wide[, -1])
        rownames(data_matrix) <- df_wide$contig_l

        # Perform clustering
        dist_rows <- dist(data_matrix, method = "euclidean")
        dist_cols <- dist(t(data_matrix), method = "euclidean")
        hc_rows <- hclust(dist_rows, method = "average")
        hc_cols <- hclust(dist_cols, method = "average")

        mod_cols <- colnames(data_matrix)[hc_cols$order]
        bin_rows <- rownames(data_matrix)[hc_rows$order]

        df_cross <- crossing(
            motif_mod = motifs,
            contig_l = bin_rows
        ) %>%
            left_join(motifs_scored_f) %>%
            mutate(
                contigl_l = factor(contig_l, levels = bin_rows),
                motif_mod = factor(motif_mod, levels = c(mod_cols)),
                txt = case_when(
                    is.na(median) ~ "\u00B7",
                    TRUE ~ NA
                ),
                median = ifelse(is.na(median), 0, median),
            )

        df_cross %>%
            ggplot(aes(x = motif_mod, y = contig_l, fill = median)) +
            geom_tile(color = "gray30") +
            geom_text(aes(label = txt), size = 10) +
            scale_fill_gradientn(
                labels = scales::percent,
                limits = c(0, 1),
                colours = c("white", PLOT_COLORS[[4]], PLOT_COLORS[[5]]), na.value = "white",
                values = c(0, 0.5, 1)
            ) +
            theme_minimal() +
            labs(
                x = "",
                y = ""
            ) +
            theme(
                axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_markdown(),
                # add bg to strio
                strip.background = element_rect(fill = "gray70"),
            )
        binstats <- df_com %>%
            filter(sample == sample_name, bin == b)
        before <- binstats %>%
            filter(decon == "before") %>%
            pull(mag_qual)
        after <- binstats %>%
            filter(decon == "after") %>%
            pull(mag_qual)
        ggsave(file.path("analysis/contamination_story/change_mag_qual_check_wstd_all/", paste0(s, "_", b, "_", before[1], "_", after[1], "_motif_score.png")), width = 30, height = 12)
    }
}



changed_bin_status %>%
    unnest(data) %>%
    left_join(df_com %>% filter(decon == "before") %>% select(bin, mag_qual))


write_delim(df_com, "analysis/contamination_story/supplmentary_data_3.csv", delim = ",")

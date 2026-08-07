# Script given to me by Kris Wang
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Checking bins against a dataset of normal cells
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
# use HBCA cells to define regular coverage per bin for BS conversion
# Rsubread was used to count reads into the bins


library(amethyst)
library(rhdf5)
library(data.table)
library(ggplot2)
library(patchwork)
library(tibble)
library(tidyr)
library(plyr); library(dplyr)
library(future)
library(furrr)
library(purrr)
library(cowplot)
library(pheatmap)
library(optparse,lib.loc="/home/rmulqueen/R/x86_64-conda-linux-gnu-library/4.4") #add this

#BiocManager::install("Rsubread")
library(Rsubread)
#devtools::install_github("navinlabcode/copykit")
library(copykit)

obj<-readRDS("batch1_homebrew_plate11.amethyst.rds")

# annotation for Rsubread with a pseudogeneID
ann <- as.data.frame(hg38_scaffold_200k_gr)

names(ann)[1] <- c("chr")

ann <- ann %>%
    mutate(geneID = 1:nrow(ann),
        strand = "*") %>%
    dplyr::rename(
        Chr = "chr",
        Start = "start",
        End = "end",
        Strand = "strand"
    )


files <- list.files("/mnt/lab/users/dminussi/projects/copy_number_pipeline_hg38/test_data/marked/", full.names = T, pattern = "*.bam$")
files_names <- list.files("/mnt/lab/users/dminussi/projects/copy_number_pipeline_hg38/test_data/marked/", full.names = F, pattern = "*.bam$")
files <- files[!str_detect(files, "bai")]
files_names <- files_names[!str_detect(files_names, "bai")]

read_counts <- parallel::mclapply(files, function(x) {
    fc <- featureCounts(x, annot.ext = ann, countMultiMappingReads = FALSE, ignoreDup = T, useMetaFeatures = FALSE)$counts
}, mc.cores = 70)

names(read_counts) <- files_names

read_counts_df <- bind_cols(read_counts)
read_counts_df <- read_counts_df %>%
    mutate(
        bin_index = 1:nrow(read_counts_df),
        Chr = ann$Chr,
        Start = ann$Start,
        End = ann$End
    )

read_counts_l <- gather(read_counts_df,
    key = "cell",
    value = "bin_count",
    -bin_index,
    -c(Chr:End)
)


# Cells are filtered according to bin count
# - Need > 30 mean bin counts
# - Need a coefficient of variation < than the median CV plus 3 times the median absolute deviation of the CV. This step also guarantees that only normal cells remain in the dataset
#
# ChromosomeY and chrX are excluded to avoid influencing the mean

stats <- read_counts_l %>%
    filter(Chr %!in% c("chrY", "chrX")) %>%
    group_by(cell) %>%
    summarise(
        mean_bc = mean(bin_count),
        median_bc = median(bin_count),
        sd_bc = sd(bin_count),
        cv = sd_bc / mean_bc,
        sd_up = mean_bc + 3 * (sd_bc),
        sd_down = mean_bc - 3 * (sd_bc)
    )

stats_hq <- stats %>%
    filter(
        mean_bc >= 30,
        median_bc >= 30,
        cv < (median(cv) + 2 * mad(cv))
    )

hist(stats_hq$cv, breaks = "fd")

read_counts_f <- read_counts_l %>%
    filter(cell %in% stats_hq$cell)

#
# Plotting some of the cells
#
# Mean is marked by the red line
# Blue lines are marking plus or minus 3 standard deviations from the mean

set.seed(2)

smp_cells <- sample(unique(read_counts_f$cell), 12)

read_counts_f %>%
    left_join(stats_hq) %>%
    filter(cell %in% smp_cells) %>%
    ggplot(aes(
        x = bin_index,
        y = bin_count
    )) +
    geom_point() +
    geom_hline(aes(yintercept = sd_up), color = "blue") +
    geom_hline(aes(yintercept = sd_down), color = "blue") +
    geom_hline(aes(yintercept = mean_bc), color = "red") +
    facet_wrap(vars(cell)) +
    cowplot::theme_cowplot()


# Marking bins that have either 3*sd up or down from the mean and counting those occurrences

read_counts_f <- read_counts_f %>%
    filter(Chr != "chrY") %>%
    group_by(cell) %>%
    mutate(
        mean_bc = mean(bin_count),
        sd_bc = sd(bin_count),
        cv = sd_bc / mean_bc,
        sd_up = mean_bc + 2.5 * (sd_bc),
        sd_down = mean_bc - 2.5 * (sd_bc)
    ) %>%
    ungroup() %>%
    mutate(
        bin_up = bin_count > sd_up,
        bin_down = bin_count < sd_down,
        both = bin_up | bin_down
    )

problematic_bins <- read_counts_f %>%
    group_by(bin_index) %>%
    count(bool_count = both) %>%
    filter(bool_count == TRUE) %>%
    arrange(desc(n))

problematic_bins <- problematic_bins %>%
    filter(n > 1)

problematic_bins %>%
    mutate(chr = ann$Chr[bin_index]) %>%
    ungroup() %>%
    dplyr::count(chr) %>%
    View()


plot(problematic_bins$n)

blacklisted_200kb_varbin <- ann %>%
    filter(ann$geneID %in% unique(problematic_bins$bin_index))

hg38_scaffold_200k_gr <- hg38_scaffold_200k_gr[1:length(hg38_scaffold_200k_gr) %!in% blacklisted_200kb_varbin$geneID]
# Fig. 3 dot plot, as in Functional_Analysis/Fig03_part2_dotplot_results_ViSEAGO.ipynb (cells 2-23):
# GO terms present in the ViSEAGO tables of all 10 carbon sources (LB excluded); mean -log10 p and mean gene frequency
# per category across conditions, min-max normalised; plotted if normalised -log10 p > 3 and frequency > 0.25.
# Bug fixed here: the notebook read the files with list.files() (alphabetical: Acetate, Fructose, Galactose, ...) but
# labelled them in a different order (Glucose, LB, Glucosamine, ...), so the file excluded as "LB" was Fructose and
# LB was kept under the label "Galactose". Here each file is read by its condition name.
# Usage: Rscript a15c_fig3_dotplot.R   -> output/Fig03_dotplot_*.png/svg and output/a15c_fig3_terms.csv
suppressPackageStartupMessages({ library(dplyr); library(reshape2); library(ggplot2); library(svglite) })
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
repo <- Sys.getenv("FA_REPO", dirname(dirname(here)))
conditions <- c("Glucose", "LB", "Glucosamine", "Glycerol", "Xylose", "Mannose", "Galactose", "Fructose", "Acetate",
                "Pyruvate", "Succinate")

read_all <- function(dir, buggy = FALSE) {
  files <- if (buggy) list.files(dir, pattern = "BP_Results_.*\\.xls", full.names = TRUE)
           else file.path(dir, paste0("BP_Results_", conditions, ".xls"))
  do.call(rbind, lapply(seq_along(files), function(i) {
    d <- read.table(files[i], header = TRUE, sep = "\t", quote = "\"", comment.char = "")
    d$Condition <- conditions[i]; d }))
}

dotplot_terms <- function(final_data) {
  sel <- final_data[, c(1, 2, 6, 11, 16, 21, 24)]
  colnames(sel) <- c("GO", "term", "Essential", "Important", "No effect", "Enchancing", "Condition")
  sel <- sel[sel$Condition != "LB", ]
  conds <- unique(sel$Condition)
  common <- unique(sel$GO[sapply(sel$GO, function(go) all(conds %in% sel$Condition[sel$GO == go]))])
  ac <- sel[sel$GO %in% common, ]
  fr <- final_data[, c(1, 2, 4, 9, 14, 19, 24)]
  colnames(fr) <- colnames(sel)
  fr <- fr[fr$Condition != "LB" & fr$GO %in% common, ]
  fr[3:6] <- lapply(fr[3:6], function(v) as.numeric(sub("%.*", "", v)))
  num <- c("Essential", "Important", "No effect")
  pv <- ac %>% group_by(term) %>% summarise(across(all_of(num), ~ mean(as.numeric(.x), na.rm = TRUE)))
  pv <- pv %>% mutate(across(all_of(num), ~ (.x - min(.x)) / (max(.x) - min(.x)) * 10))
  fq <- fr %>% group_by(term) %>% summarise(across(all_of(num), ~ mean(.x, na.rm = TRUE)))
  fq <- fq %>% mutate(across(all_of(num), ~ (.x - min(.x)) / (max(.x) - min(.x))))
  m <- melt(as.data.frame(pv), id.vars = "term", measure.vars = num)
  m$freq <- melt(as.data.frame(fq), id.vars = "term", measure.vars = num)$value
  list(n_common = length(common), plotted = m %>% filter(value > 3, freq > 0.25))
}

plot_save <- function(mp, name) {
  # axis label "Mean-effect", the name of the category in the text (was "No effect"; 21 Oct 2026)
  mp$variable <- factor(mp$variable, levels = c("Essential", "Important", "No effect"),
                        labels = c("Essential", "Important", "Mean-effect"))
  p <- ggplot(mp, aes(x = variable, y = term, color = value, size = freq)) + geom_point() +
    labs(x = "Categories", y = "GO Term", size = "Fraction of genes in the GO term") + theme_minimal() +
    scale_color_viridis_c(name = "-log10 (pvalue)") + theme(axis.line = element_blank()) +
    theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, colour = "black")) +  # centred labels (21 Oct 2026)
    theme(axis.text.y = element_text(colour = "black")) + ylab("") + theme(axis.ticks = element_blank())
  ggsave(file.path(here, "output", paste0(name, ".png")), p, width = 8, height = 8)
  ggsave(file.path(here, "output", paste0(name, ".svg")), p, width = 8, height = 8)
}

runs <- list(
  repo_buggy      = list(dir = file.path(dirname(here), "original_submission", "Viseago_Res"), buggy = TRUE),
  repo_fixed      = list(dir = file.path(dirname(here), "original_submission", "Viseago_Res"), buggy = FALSE),
  corrected       = list(dir = file.path(here, "output", "viseago", "corrected"), buggy = FALSE),
  corrected_no287 = list(dir = file.path(here, "output", "viseago", "no287"), buggy = FALSE))
all_terms <- list()
for (r in names(runs)) {
  res <- dotplot_terms(read_all(runs[[r]]$dir, runs[[r]]$buggy))
  cat(sprintf("%-16s GO terms in all 10 carbon sources: %3d | dots plotted: %3d | terms plotted: %3d\n", r,
              res$n_common, nrow(res$plotted), length(unique(res$plotted$term))))
  all_terms[[r]] <- transform(res$plotted, run = r)
  if (r %in% c("corrected", "repo_buggy")) plot_save(res$plotted, paste0("Fig03_dotplot_", r))
}
write.csv(do.call(rbind, all_terms), file.path(here, "output", "a15c_fig3_terms.csv"), row.names = FALSE)

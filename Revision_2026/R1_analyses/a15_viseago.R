# GO enrichment per condition (Fig. 3, Supplementary Tables S2-S12), as in
# Functional_Analysis/Fig03_part1_run_ViSEAGO.ipynb (same parameters: BP, nodeSize = 20, topGO elim + Fisher,
# cutOff = 0.01, merge cutoff 0.01; background = all E. coli K-12 genes in org.EcK12.eg.db), run in three settings:
#   original  : repo inputs (misaligned Table S1 fitness) + the 287 essential genes without fitness data
#   corrected : regenerated_inputs (correct replicate averages, glucosamine HCl) + the 287 genes   -> new S2-S12
#   no287     : regenerated_inputs without the 287 genes  (R1-M1 sensitivity analysis)
# Usage: Rscript a15_viseago.R [original|corrected|no287 ...]   (default: all three)
# Output: output/viseago/<setting>/BP_Results_<condition>.xls
suppressPackageStartupMessages({
  library(org.EcK12.eg.db)
  library(ViSEAGO)
  library(dplyr)
})

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
repo <- Sys.getenv("FA_REPO", dirname(dirname(here)))
settings <- commandArgs(trailingOnly = TRUE)
if (length(settings) == 0) settings <- c("original", "corrected", "no287")

infiles <- list(
  original  = file.path(dirname(here), "original_submission", "Infiles", "genes_categories"),  # as first submitted
  corrected = file.path(dirname(here), "regenerated_inputs", "genes_categories"),
  no287     = file.path(dirname(here), "regenerated_inputs", "genes_categories")
)
nofit287 <- read.csv(file.path(repo, "Functional_Analysis", "Infiles", "essential_287gdeNofit_otherConds.csv"))

background <- keys(org.EcK12.eg.db, keytype = "ENTREZID")
entrezids <- AnnotationDbi::select(org.EcK12.eg.db, keys = background, columns = c("ENTREZID", "SYMBOL"))
Bioconductor <- ViSEAGO::Bioconductor2GO()
myGENE2GO <- ViSEAGO::annotate("org.EcK12.eg.db", Bioconductor)

gene_sets <- function(file_path, add287) {
  d <- read.csv(file_path)
  es <- d %>% filter(Categ == "Essen") %>% select(gene)
  if (add287) es <- rbind(es, nofit287)
  sets <- list(
    es  = es,
    imp = d %>% filter(Categ == "Important") %>% select(gene),
    ne  = d %>% filter(Categ == "No_effect" | Categ == "No_effect2") %>% select(gene),
    enh = d %>% filter(Categ == "Enhancing") %>% select(gene)
  )
  lapply(sets, function(s) merge(s, entrezids, by.x = "gene", by.y = "SYMBOL", all.x = TRUE)$ENTREZID)
}

conditions <- c("Glucose", "LB", "Glucosamine", "Glycerol", "Xylose", "Mannose", "Galactose", "Fructose", "Acetate",
                "Pyruvate", "Succinate")

for (setting in settings) {
  out <- file.path(here, "output", "viseago", setting)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  for (c in conditions) {
    s <- gene_sets(file.path(infiles[[setting]], paste0("genesfitness_", c, ".csv")), add287 = setting != "no287")
    td <- lapply(s, function(g) ViSEAGO::create_topGOdata(geneSel = g, allGenes = background, gene2GO = myGENE2GO,
                                                          ont = "BP", nodeSize = 20))
    BP_essential <- td$es; BP_important <- td$imp; BP_no_effect <- td$ne; BP_enhance <- td$enh
    elim_BP_essential <- topGO::runTest(BP_essential, algorithm = "elim", statistic = "fisher", cutOff = 0.01)
    elim_BP_important <- topGO::runTest(BP_important, algorithm = "elim", statistic = "fisher", cutOff = 0.01)
    elim_BP_no_effect <- topGO::runTest(BP_no_effect, algorithm = "elim", statistic = "fisher", cutOff = 0.01)
    elim_BP_enhance <- topGO::runTest(BP_enhance, algorithm = "elim", statistic = "fisher", cutOff = 0.01)
    BP_sResults <- ViSEAGO::merge_enrich_terms(
      cutoff = 0.01,
      Input = list(Essential = c("BP_essential", "elim_BP_essential"),
                   Important = c("BP_important", "elim_BP_important"),
                   No_effect = c("BP_no_effect", "elim_BP_no_effect"),
                   Fitness_enhancing = c("BP_enhance", "elim_BP_enhance")))
    ViSEAGO::show_table(BP_sResults, file.path(out, paste0("BP_Results_", c, ".xls")))
    cat(setting, c, "done\n")
  }
}

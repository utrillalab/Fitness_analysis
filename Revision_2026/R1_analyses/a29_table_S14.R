# Supplementary Table S14: GO enrichment of the core mean-effect genes (mean-effect in all 46 conditions, 1,946 genes),
# as in Functional_Analysis/Fig05_p1_ViSEAGO_listag_noeff.ipynb (same input Infiles/core_no_effect.xlsx, same
# parameters: BP, nodeSize = 20, topGO elim + Fisher, cutOff = 0.01, merge cutoff 0.01; background = all E. coli K-12
# genes in org.EcK12.eg.db), rerun with the current GO annotation used for the corrected Tables S2-S12 (a15).
# The core set does not depend on the Table S1 averaging error (recomputed from the 46 conditions: same 1,946 genes).
# Output: output/viseago/core_mean_effect/BP_Results_core_mean_effect.xls (read by a29_table_S14.py)
suppressPackageStartupMessages({
  library(org.EcK12.eg.db)
  library(ViSEAGO)
  library(readxl)
})

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
repo <- Sys.getenv("FA_REPO", dirname(dirname(here)))

background <- keys(org.EcK12.eg.db, keytype = "ENTREZID")
entrezids <- AnnotationDbi::select(org.EcK12.eg.db, keys = background, columns = c("ENTREZID", "SYMBOL"))
Bioconductor <- ViSEAGO::Bioconductor2GO()
myGENE2GO <- ViSEAGO::annotate("org.EcK12.eg.db", Bioconductor)

core <- read_excel(file.path(repo, "Functional_Analysis", "Infiles", "core_no_effect.xlsx"))
g <- merge(data.frame(gene = core$gene_name), entrezids, by.x = "gene", by.y = "SYMBOL", all.x = TRUE)
cat("core genes:", nrow(core), "; mapped to ENTREZID:", sum(!is.na(g$ENTREZID)), "\n")

BP_core <- ViSEAGO::create_topGOdata(geneSel = g$ENTREZID, allGenes = background, gene2GO = myGENE2GO,
                                     ont = "BP", nodeSize = 20)
elim_BP_core <- topGO::runTest(BP_core, algorithm = "elim", statistic = "fisher", cutOff = 0.01)
BP_sResults <- ViSEAGO::merge_enrich_terms(cutoff = 0.01, Input = list(Core_mean_effect = c("BP_core", "elim_BP_core")))
out <- file.path(here, "output", "viseago", "core_mean_effect")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
ViSEAGO::show_table(BP_sResults, file.path(out, "BP_Results_core_mean_effect.xls"))
cat("done\n")

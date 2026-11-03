# Supplementary Table S15: GO enrichment of the 108 recurrent fitness-enhancing candidates (fitness > 0.5 in >= 10 of 46),
# with the parameters of Fig05_p1_ViSEAGO_listag_noeff.ipynb (copied from the enhancing analysis; input Infiles/108_enhancing_10Conditions.xlsx; same
# parameters: BP, nodeSize = 20, topGO elim + Fisher, cutOff = 0.01, merge cutoff 0.01; background = all E. coli K-12
# genes in org.EcK12.eg.db), rerun with the current GO annotation used for the corrected Tables S2-S12 (a15).
# The 108 candidates do not depend on the Table S1 averaging error (they come from the 46 conditions).
# Output: output/viseago/enhancing_108/BP_Results_enhancing_108.xls (compared with the original Table S15 in a30_table_S15_GO.py)
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

core <- read_excel(file.path(repo, "Functional_Analysis", "Infiles", "108_enhancing_10Conditions.xlsx"), skip = 1)
g <- merge(data.frame(gene = core[["gene name"]]), entrezids, by.x = "gene", by.y = "SYMBOL", all.x = TRUE)
cat("enhancing candidates:", nrow(core), "; mapped to ENTREZID:", sum(!is.na(g$ENTREZID)), "\n")

BP_core <- ViSEAGO::create_topGOdata(geneSel = g$ENTREZID, allGenes = background, gene2GO = myGENE2GO,
                                     ont = "BP", nodeSize = 20)
elim_BP_core <- topGO::runTest(BP_core, algorithm = "elim", statistic = "fisher", cutOff = 0.01)
BP_sResults <- ViSEAGO::merge_enrich_terms(cutoff = 0.01, Input = list(Fitness_enhancing = c("BP_core", "elim_BP_core")))
out <- file.path(here, "output", "viseago", "enhancing_108")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
ViSEAGO::show_table(BP_sResults, file.path(out, "BP_Results_enhancing_108.xls"))
cat("done\n")

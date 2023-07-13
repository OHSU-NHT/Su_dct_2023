
library(Seurat)
library(here)

# 1. Create score-defining gene lists ---------------------------------------------

# DCT1 or DCT2 score - Top 10 DEGs in DCT1 or DCT2
DCT1_marker_gene_list <- list(c("Erbb4", "Egf", "Trpm7", "Fgf13", "Col5a2", "Umod", "Ptgfr", "Stk32b", "Rtl4", "Abca13"))
DCT2_marker_gene_list <- list(c("Slc8a1", "Arl15", "Calb1", "Slc2a9", "Phactr1", "Gls", "S100g", "Kl", "Klk1", "Egfem1"))

# DCT score - DCT1 and DCT2 score-defining genes plus Slc12a3 (NCC)
DCT_marker_gene_list <- list(c("Slc12a3", "Erbb4", "Egf", "Trpm7", "Fgf13", "Col5a2", "Umod", "Ptgfr", "Stk32b", "Rtl4", "Abca13", "Slc8a1", "Arl15", "Calb1", "Slc2a9", "Phactr1", "Gls", "S100g", "Kl", "Klk1", "Egfem1"))

# Ca, Mg, or ENaC score - genes that are known to be relavant to Ca transport, Mg transport, or ENaC activity
Ca_cassette <- list(c("Slc8a1", "Vdr", "Trpv5", "Calb1", "S100g", "Ryr2"))
Mg_cassette <- list(c("Slc41a3", "Cnnm2", "Fxyd2", "Prox1", "Umod", "Egf", "Trpm6", "Trpm7")) 
ENaC_score <- list(c("Nr3c2", "Klk1", "Scnn1a", "Scnn1b", "Scnn1g"))

# 2. AddModuleScore in control samples --------------------------------------------
SO <- readRDS(here("CONTROL_DCT.rds"))

SO <- AddModuleScore(object = SO, features = DCT_marker_gene_list, name = "DCT_score")
SO <- AddModuleScore(object = SO, features = DCT1_marker_gene_list, name = "DCT1_score")
SO <- AddModuleScore(object = SO, features = DCT2_marker_gene_list, name = "DCT2_score")
SO <- AddModuleScore(object = SO, features = Ca_cassette, name = "Ca_score")
SO <- AddModuleScore(object = SO, features = Mg_cassette, name = "Mg_score")
SO <- AddModuleScore(object = SO, features = ENaC_score, name = "ENaC_score")

saveRDS(SO, here("CONTROL_DCT.rds"))

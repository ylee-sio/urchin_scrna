library(tidyverse)
library(Seurat)
library(gtools)

# unmodified Strongylocentrotus purpuratus single cell RNA sequencing dataset from the Wessel Lab at Brown University
urchin_0 = readRDS("data_sources/primary/scrna/GSE149221_SpInteg")

# urchin_1 is a renaming of clusters in the original urchin_0 SeuratObject
# it is not a reclustering of urchin_0
urchin_1 = RenameIdents(urchin_0,
                        '0' = "aboral-ectoderm",
                        '1' = "oral_ectoderm_eb",
                        '2' = "ciliated_cells",
                        '3' = "neural_8_to_eg",
                        '4' = "oral_ectoderm_hb",
                        '5' = "aboral_ectoderm/neural",
                        '6' = "endoderm_hb",
                        '7' = "ciliated_cells_m",
                        '8' = "endoderm_endo-msd",
                        '9' = "neural_AnEctoE",
                        '10' = "transient_8_to_eb",
                        '11' = "nsm_pigment_cells",
                        '12' = "oral_ectoderm_64",
                        '13' = "oral_ectoderm_m",
                        '14' = "endoderm_eb",
                        '15' = "vault_protein_containing",
                        '16' = "skeleton_64",
                        '17' = "neural_AnEctoL",
                        '18' = "neural_neuro_prog",
                        '19' = "skeleton_hb",
                        '20' = "transient",
                        '21' = "germline"
                        )

# urchin_2 is a renaming of clusters and a reorganization of cluster names in order to facilitate 
# organization of clusters by germ layer when visualizing the data
urchin_2 = RenameIdents(urchin_0,
                        '6' = "endoderm_hb",
                        '8' = "endoderm_endo-msd",
                        '14' = "endoderm_eb",
                        '11' = "nsm_pigment_cells",
                        '16' = "skeleton_64",
                        '19' = "skeleton_hb",
                        '20' = "transient",
                        '21' = "germline",
                        '0' = "aboral-ectoderm",
                        '1' = "oral_ectoderm_eb",
                        '2' = "ciliated_cells",
                        '3' = "neural_8_to_eg",
                        '4' = "oral_ectoderm_hb",
                        '5' = "aboral_ectoderm/neural",
                        '7' = "ciliated_cells_m",
                        '9' = "neural_AnEctoE",
                        '12' = "oral_ectoderm_64",
                        '13' = "oral_ectoderm_m",
                        '17' = "neural_AnEctoL",
                        '18' = "neural_neuro_prog",
                        '10' = "transient_8_to_eb",
                        '15' = "vault_protein_containing")

# urchin_3 is a renaming of clusters and grouping of the clusters which are from similar cell types
urchin_3 = RenameIdents(urchin_0,
                        '6' = "endoderm",
                        '8' = "endoderm",
                        '14' = "endoderm",
                        '11' = "nsm_pigment_cells",
                        '16' = "skeletal",
                        '19' = "skeletal",
                        '20' = "transient",
                        '21' = "germline",
                        '0' = "ectoderm",
                        '1' = "ectoderm",
                        '2' = "ciliated_cells",
                        '3' = "neural",
                        '4' = "ectoderm",
                        '5' = "neural",
                        '7' = "ciliated_cells",
                        '9' = "neural",
                        '12' = "ectoderm",
                        '13' = "ectoderm",
                        '17' = "neural",
                        '18' = "neural",
                        '10' = "transient",
                        '15' = "vault_protein_containing")


# UNIVERSAL

# obtains list of abbreviated labels used to delineate sample sets, in this case, this also means developmental stages; UNIVERSAL
dev_stage_names_abbreviated_list = unique(urchin_0@meta.data$orig.ident)

# creates list of desired names for use to replace those in dev_stage_names_abbreviated_list; UNIVERSAL
dev_stage_names_list = c(
                    "8 Cell",
                    "64 Cell",
                    "Morula",
                    "Early Blastula",
                    "Hatched Blastula",
                    "Mid Blastula",
                    "Early Gastrula",
                    "Late Gastrula")

# creates list of all "named" loci that exists in the unmodified and modified SeuratObjects; UNIVERSAL
named_locus_list = rownames(urchin_0)[which(rownames(urchin_0) %>% str_sub(1,3) != "LOC")] %>% sort()

# create tibbles for all curated protein families found
abc = read_csv("data_sources/tertiary/protein_groups/abc.csv")
aqp = read_csv("data_sources/tertiary/protein_groups/aqp.csv")
cyp = read_csv("data_sources/tertiary/protein_groups/cyp.csv")
gpcr = read_csv("data_sources/tertiary/protein_groups/gpcr.csv")
nhr = read_csv("data_sources/tertiary/protein_groups/nhr.csv")
receptors_synth = read_csv("data_sources/tertiary/protein_groups/receptor_synth.csv")
slc = read_csv("data_sources/tertiary/protein_groups/slc.csv")
slc2 = read_csv("data_sources/tertiary/protein_groups/slc2.csv")
tlr = read_csv("data_sources/tertiary/protein_groups/tlr.csv")

# URCHIN_2 PROCESSING AND PARSING **********************************************************

# replaces orig.idents in urchin_2 enumerated in dev_stage_names_abbreviated_list with dev_stage_names_list
urchin_2@meta.data$orig.ident = intramutate(urchin_2@meta.data$orig.ident, dev_stage_names_abbreviated_list, dev_stage_names_list)

# separates the urchin_2 by sample sets, i.e developmental stages
urchin_2_cell8 = subset(urchin_2, orig.ident  == "8 Cell")
urchin_2_cell64 = subset(urchin_2, orig.ident == "64 Cell")
urchin_2_morula = subset(urchin_2, orig.ident == "Morula")
urchin_2_eb = subset(urchin_2, orig.ident == "Early Blastula")
urchin_2_hb = subset(urchin_2, orig.ident == "Hatched Blastula")
urchin_2_mb = subset(urchin_2, orig.ident == "Mid Blastula")
urchin_2_eg = subset(urchin_2, orig.ident == "Early Gastrula")
urchin_2_lg = subset(urchin_2, orig.ident == "Late Gastrula")

# creates single list containing each of the sample sets above
urchin_2_devolist = c(
	urchin_2_cell8, 
	urchin_2_cell64, 
	urchin_2_morula, 
	urchin_2_eb, 
	urchin_2_hb, 
	urchin_2_mb, 
	urchin_2_eg, 
	urchin_2_lg
	)

# URCHIN_3 PROCESSING AND PARSING **********************************************************

# replaces orig.idents in urchin_3 enumerated in dev_stage_names_abbreviated_list with dev_stage_names_list

urchin_3@meta.data$orig.ident = intramutate(urchin_3@meta.data$orig.ident, dev_stage_names_abbreviated_list, dev_stage_names_list)

# separates the urchin_3 by sample sets, i.e developmental stages
urchin_3_cell8 = subset(urchin_3, orig.ident  == "8 Cell")
urchin_3_cell64 = subset(urchin_3, orig.ident == "64 Cell")
urchin_3_morula = subset(urchin_3, orig.ident == "Morula")
urchin_3_eb = subset(urchin_3, orig.ident == "Early Blastula")
urchin_3_hb = subset(urchin_3, orig.ident == "Hatched Blastula")
urchin_3_mb = subset(urchin_3, orig.ident == "Mid Blastula")
urchin_3_eg = subset(urchin_3, orig.ident == "Early Gastrula")
urchin_3_lg = subset(urchin_3, orig.ident == "Late Gastrula")

# creates single list containing each of the sample sets above
urchin_3_devolist = c(
	urchin_3_cell8, 
	urchin_3_cell64, 
	urchin_3_morula, 
	urchin_3_eb, 
	urchin_3_hb, 
	urchin_3_mb, 
	urchin_3_eg, 
	urchin_3_lg
	)

# wildtype SeuratObjects
sp48_72 = readRDS("data_sources/primary/scrna/GSE155427_Sp48and72")
sp48_72 = RenameIdents(object = sp48_72, c(
                        '0' = "Ectoderm/Uncharacterized",
                        '1' = "Apical ectoderm",
                        '2' = "Differentiated pigment cells",
                        '3' = "Ciliary band neurons (1)",
                        '4' = "Apical plate/uncharacterized",
                        '5' = "Oral ectoderm/Mouth",
                        '6' = "Lateral ectoderm (right)",
                        '7' = "Ciliary band neurons (2)",
                        '8' = "Aboral ectoderm (1)",
                        '9' = "Aboral ectoderm (2)",
                        '10' = "Apical plate, proneural",
                        '11' = "Skeleton",
                        '12' = "Mesodermal Cells",
                        '13' = "Mitotic Pigment Cells",
                        '14' = "Mid-gut",
                        '15' = "Serotoninergic neurons (apical plate) (1)",
                        '16' = "Ciliary band neurons (3)",
                        '17' = "Serotoninergic neurons (apical plate) (2)")
                        )
levels(sp48_72@meta.data$seurat_clusters) = c(
	"Ectoderm/Uncharacterized",
	"Apical ectoderm",
	"Differentiated pigment cells",
	"Ciliary band neurons (1)",
	"Apical plate/uncharacterized",
	"Oral ectoderm/Mouth",
	"Lateral ectoderm (right)",
	"Ciliary band neurons (2)",
	"Aboral ectoderm (1)",
	"Aboral ectoderm (2)",
	"Apical plate, proneural",
	"Skeleton",
	"Mesodermal Cells",
	"Mitotic Pigment Cells",
	"Mid-gut",
	"Serotoninergic neurons (apical plate) (1)",
	"Ciliary band neurons (3)",
	"Serotoninergic neurons (apical plate) (2)")

sp48_wt = subset(sp48_72, orig.ident == "Sp48")
sp72_wt = subset(sp48_72, orig.ident == "Sp72")


# morpholino SeuratObjects
sp48mo = readRDS("data_sources/primary/scrna/GSE155427_Sp48MO")
sp48mo = RenameIdents(object = sp48mo,
                        c('0' = "Aboral ectoderm (1)",
                        '1' = "Ciliary band neurons (1)",
                        '2' = "Ciliary band neurons (2)",
                        '3' = "Ectoderm-Blastopore",
                        '4' = "Aboral ectoderm (2)",
                        '5' = "Endoderm/Mid-gut",
                        '6' = "Differentiated pigment cells",
                        '7' = "Apical ectoderm",
                        '8' = "Aboral ectoderm (3)",
                        '9' = "Ciliary band neurons (3)",
                        '10' = "Skeleton (1)",
                        '11' = "Endoderm/Foregut",
                        '12' = "Skeleton (2)",
                        '13' = "Coelomic Pouches",
                        '14' = "Serotoninergic neurons")
                        )
levels(sp48mo@meta.data$seurat_clusters) = c(
	"Aboral ectoderm (1)",
	"Ciliary band neurons (1)",
	"Ciliary band neurons (2)",
	"Ectoderm-Blastopore",
	"Aboral ectoderm (2)",
	"Endoderm/Mid-gut",
	"Differentiated pigment cells",
	"Apical ectoderm",
	"Aboral ectoderm (3)",
	"Ciliary band neurons (3)",
	"Skeleton (1)",
	"Endoderm/Foregut",
	"Skeleton (2)",
	"Coelomic Pouches",
	"Serotoninergic neurons")

sp48_gcm_mo = subset(sp48mo, orig.ident == "Sp48gcmMO")
sp48_con_mo = subset(sp48mo, orig.ident == "Sp48controlMO")




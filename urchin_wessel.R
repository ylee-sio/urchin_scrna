urchin_0 = readRDS("data_sources/GSE149221_SpInteg")
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
                        '15' = "wtf",
                        '16' = "skeleton_64",
                        '17' = "neural_AnEctoL",
                        '18' = "neural_neuro_prog",
                        '19' = "skeleton_hb",
                        '20' = "transient",
                        '21' = "germline"
                        )

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
                        '15' = "vault_protein_containing"
)

dev_stage_names_abbreviated_list = unique(urchin_2@meta.data$orig.ident)
dev_stage_names_list = c(
                    "8 Cell",
                    "64 Cell",
                    "Morula",
                    "Early Blastula",
                    "Hatched Blastula",
                    "Mid Blastula",
                    "Early Gastrula",
                    "Late Gastrula")
urchin_2@meta.data$orig.ident = intramutate(urchin_2@meta.data$orig.ident, dev_stage_names_abbreviated_list, dev_stage_names_list)
named_locus_list = rownames(urchin_2)[which(rownames(urchin_2) %>% str_sub(1,3) != "LOC")]

cell8 = subset(urchin_2, orig.ident  == "8 Cell")
cell64 = subset(urchin_2, orig.ident == "64 Cell")
morula = subset(urchin_2, orig.ident == "Morula")
eb = subset(urchin_2, orig.ident == "Early Blastula")
hb = subset(urchin_2, orig.ident == "Hatched Blastula")
mb = subset(urchin_2, orig.ident == "Mid Blastula")
eg = subset(urchin_2, orig.ident == "Early Gastrula")
lg = subset(urchin_2, orig.ident == "Late Gastrula")

devolist = c(cell8, cell64, morula, eb, hb, mb, eg, lg)
print("Loaded the following: urchin_0, urchin_1, cell8, cell64, morula, eb, hb, mb, eg, lg, dev_stage_names, devolist")

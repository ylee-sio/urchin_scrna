urchin_2 = readRDS("data_sources/GSE155427_Sp48and72")
urchin_3 = readRDS("data_sources/GSE155427_Sp48MO")

urchin_2 = RenameIdents(object = urchin_2, c(
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
levels(urchin_2@meta.data$seurat_clusters) = c(
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

u2_48 = subset(urchin_2, orig.ident == "Sp48")
u2_72 = subset(urchin_2, orig.ident == "Sp72")

urchin_3 = RenameIdents(object = urchin_3,
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
levels(urchin_3@meta.data$seurat_clusters) = c(
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

u3_gcmMO = subset(urchin_3, orig.ident == "Sp48gcmMO")
u3_conMO = subset(urchin_3, orig.ident == "Sp48controlMO")

print("Loaded the following: urchin_2, u2_48, u2_72, urchin_3, u3_gcmMO, u3_conMO")
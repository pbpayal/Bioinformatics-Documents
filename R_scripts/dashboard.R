getwd()
setwd("/Users/pbanerjee/Desktop/Dashboard")

f_name = "BioinformaticsInitia_DATA_LABELS_2019-07-17_0911.csv"
data_df = read.csv2(file=f_name, header = TRUE, sep = ",")
head(data_df)
# attach(data_df)
# detach(data_df)
# data_mat = data.matrix(data_df)
# head(data_mat)

data_df_mod = data_df[,c(4:8,11:12,15:20,24,26:35,49:56,59:65)]

# install.packages("reshape")
library(reshape)
# detach("package:dplyr", unload=TRUE) # Since both reshape and dplyr have rename functions
# and to run this coomad just use the rehape package function.

data_df_mod_rename = rename(data_df_mod,c(Record.ID="Record_ID",Survey.Timestamp="Date" ,Project.Title="Project_Title",
                                          Project.Description="Project_Description", First.Name="First_Name", Last.Name="Last_Name",
                                          Principal.Investigator.Name.="PI_Name",
                                          Institution..choice.Childrens.National.Health.System.="CNHS",
                                          Institution..choice.George.Washington.University.="GWU",
                                          Center.Department..choice.Center.for.Cancer...Immunology.Research.="Cancer_&_Immunology",
                                          Center.Department..choice.Center.for.Genetic.Medicine.Research.="Genmed",
                                          Center.Department..choice.Center.for.Neuroscience.Research.="Neuroscience",
                                          Center.Department..choice.Center.for.Translational.Science.="CTS",
                                          Center.Department..choice.Clinical.and.Translational.Science.Institute.at.Childrens.National.="CTSI",
                                          Center.Department..choice.Sheikh.Zayed.Institute.for.Pediatric.Surgical.Innovation.="Sheikh_Zayed",
                                          IDDRC.Funded.Research="IDDRC", Bioinformatics.Services..choice.Data.Analysis.="Data_Analysis",
                                          Bioinformatics.Services..choice.Study.Design.="Study_Design",
                                          Experimental.Design..choice.ChIP.Sequencing.="ChIP",
                                          Experimental.Design..choice.Exome.Sequencing.="Exome",
                                          Experimental.Design..choice.Metagenomics.Microbiome.="Metagenomics",
                                          Experimental.Design..choice.Microarray.="Microarray",
                                          Experimental.Design..choice.RNA.Sequencing.="RNA_Sequencing",
                                          Experimental.Design..choice.Single.Cell.Sequencing.="Single_Cell_Sequencing",
                                          Experimental.Design..choice.T.cell.Receptor.Sequencing.="T-Cell_Receptor_Sequencing",
                                          Experimental.Design..choice.Whole.Genome.Sequencing.="WGS",
                                          Bioinformatics.Analysis..choice.Copy.Number.Variation.="CNV_detection",
                                          Bioinformatics.Analysis..choice.De.Novo.Assembly.="DeNovo",
                                          Bioinformatics.Analysis..choice.Differential.Gene.Expression.="Differential_Gene_expression",
                                          Bioinformatics.Analysis..choice.Fusion.Analysis.="Fusion_Analysis",
                                          Bioinformatics.Analysis..choice.Metagenome.Microbiome.Analysis.="Microbiome",
                                          Bioinformatics.Analysis..choice.Pathway.Analysis.="Pathway_Analysis",
                                          Bioinformatics.Analysis..choice.Single.Nucleotide.Variation.Detection.="SNV_SNP_Detection",
                                          Bioinformatics.Analysis..choice.Structural.Variant.Detection.="Structural_Variant_Detection",
                                          Organism..choice.Fruit.Fly.="Fruit_Fly",
                                          Organism..choice.Human.="Human",
                                          Organism..choice.Mouse.="Mouse",
                                          Organism..choice.Pig.="Pig",
                                          Organism..choice.Rat.="Rat",
                                          Organism..choice.Yeast.="Yeast",
                                          Organism..choice.Zebrafish.="Zebrafish"))

write.table(data_df_mod_rename, "/Users/pbanerjee/Desktop/CBU_Projects.txt", sep="\t", col.names=T)


library(plyr)
data_df_mod_rename$CNHS <- revalue(data_df_mod_rename$CNHS, c("Checked"="CNHS"))
data_df_mod_rename$GWU <- revalue(data_df_mod_rename$GWU , c("Checked"="GWU"))
data_df_mod_rename$`Cancer_&_Immunology` <- revalue(data_df_mod_rename$`Cancer_&_Immunology` , c("Checked"="Cancer_&_Immunology"))
data_df_mod_rename$Genmed <- revalue(data_df_mod_rename$Genmed, c("Checked"="Genmed"))
data_df_mod_rename$Neuroscience <- revalue(data_df_mod_rename$Neuroscience, c("Checked"="Neuroscience"))
data_df_mod_rename$CTS <- revalue(data_df_mod_rename$CTS, c("Checked"="CTS"))
data_df_mod_rename$Sheikh_Zayed <- revalue(data_df_mod_rename$Sheikh_Zayed, c("Checked"="Sheikh_Zayed"))
data_df_mod_rename$CTSI <- revalue(data_df_mod_rename$CTSI, c("Checked"="CTSI"))
data_df_mod_rename$Study_Design <- revalue(data_df_mod_rename$Study_Design, c("Checked"="Study_Design"))
data_df_mod_rename$Data_Analysis <- revalue(data_df_mod_rename$Data_Analysis, c("Checked"="Data_Analysis"))
data_df_mod_rename$RNA_Sequencing <- revalue(data_df_mod_rename$RNA_Sequencing, c("Checked"="RNA_Sequencing"))
data_df_mod_rename$Exome <- revalue(data_df_mod_rename$Exome, c("Checked"="Exome"))
data_df_mod_rename$WGS <- revalue(data_df_mod_rename$WGS, c("Checked"="WGS"))
data_df_mod_rename$Microarray <- revalue(data_df_mod_rename$Microarray, c("Checked"="Microarray"))
data_df_mod_rename$Metagenomics <- revalue(data_df_mod_rename$Metagenomics, c("Checked"="Metagenomics"))
data_df_mod_rename$ChIP <- revalue(data_df_mod_rename$ChIP, c("Checked"="CHiP"))
data_df_mod_rename$Single_Cell_Sequencing <- revalue(data_df_mod_rename$Single_Cell_Sequencing, c("Checked"="Single_Cell_Sequencing"))
data_df_mod_rename$`T-Cell_Receptor_Sequencing` <- revalue(data_df_mod_rename$`T-Cell_Receptor_Sequencing`, c("Checked"="T-Cell_Receptor_Sequencing"))
data_df_mod_rename$DeNovo <- revalue(data_df_mod_rename$DeNovo, c("Checked"="DeNovo"))
data_df_mod_rename$SNV_SNP_Detection <- revalue(data_df_mod_rename$SNV_SNP_Detection, c("Checked"="SNV_SNP_Detection"))
data_df_mod_rename$CNV_detection <- revalue(data_df_mod_rename$CNV_detection, c("Checked"="CNV_detection"))
data_df_mod_rename$Structural_Variant_Detection <- revalue(data_df_mod_rename$Structural_Variant_Detection, c("Checked"="Structural_Variant_Detection"))
data_df_mod_rename$Fusion_Analysis <- revalue(data_df_mod_rename$Fusion_Analysis, c("Checked"="Fusion_Analysis"))
data_df_mod_rename$Differential_Gene_expression <- revalue(data_df_mod_rename$Differential_Gene_expression, c("Checked"="Differential_Gene_expression"))
data_df_mod_rename$Pathway_Analysis <- revalue(data_df_mod_rename$Pathway_Analysis, c("Checked"="Pathway_Analysis"))
data_df_mod_rename$Microbiome <- revalue(data_df_mod_rename$Microbiome, c("Checked"="Microbiome"))
data_df_mod_rename$Human <- revalue(data_df_mod_rename$Human, c("Checked"="Human"))
data_df_mod_rename$Mouse <- revalue(data_df_mod_rename$Mouse, c("Checked"="Mouse"))
data_df_mod_rename$Rat <- revalue(data_df_mod_rename$Rat, c("Checked"="Rat"))
data_df_mod_rename$Fruit_Fly <- revalue(data_df_mod_rename$Fruit_Fly, c("Checked"="Fruit_Fly"))
data_df_mod_rename$Zebrafish <- revalue(data_df_mod_rename$Zebrafish, c("Checked"="Zebrafish"))
data_df_mod_rename$Yeast <- revalue(data_df_mod_rename$Yeast, c("Checked"="Yeast"))
data_df_mod_rename$Pig <- revalue(data_df_mod_rename$Pig, c("Checked"="Pig"))



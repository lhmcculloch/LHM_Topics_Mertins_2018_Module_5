#LHM
#10/26/18
#ATiB Module 5 Code

#I worked with the proteomics data. The objectives I achieved here were:

#1. Cluster the data into three clusters according to the subset of PAM50 genes present in the data.
#I did this using FANNY rather than standard kmeans to avoid imputation/dealing with missing values at this step
#(and because the supplement suggested that the authors used that method with the proteomics data, for similar reasons)

#2. Perform GSEA on the full dataset (after some value removal and imputation) based on the original PAM50 classifications
#for the samples. Our group ultimately decided to combine the LumA and LumB samples in the comparison (giving us three
#subtypes), although I originally did the comparisons using four types as well.

#Necessary packages:
library(cluster)
library(ggplot2)
library(RColorBrewer)
library(bnstruct)
library(biomart)

#Reading in the data:
proteome <- read.table("proteome.txt",header = TRUE, sep="\t",stringsAsFactors=FALSE, quote="") #Used file sent to class
proteome_subtypes <- read.table("subtypes.txt",header = TRUE, sep="\t",stringsAsFactors=FALSE, quote="") #Used file sent to class
pam50_genes <- read.table("PAM50_genes.txt",header = TRUE, sep=",",stringsAsFactors=FALSE, quote="") #List of PAM50 genes

#Narrowing down the proteome data to just the genes of interest:
proteome_pam50 <- data.frame() #Initialize a data frame to pull my data into

for (i in 1:length(pam50_genes)) {
  y <- paste(colnames(pam50_genes[i]), "-", sep = "")
  z <- paste("^", y, sep="")
  x <- subset(proteome, grepl(z, gene))
  proteome_pam50 <- rbind(proteome_pam50, x) 
}

#Getting rid of unwanted samples:

#Get rid of the normal samples
proteome_pam50$X263d3f.I.CPTAC <- NULL
proteome_pam50$c4155b.C.CPTAC <- NULL
proteome_pam50$blcdb9.I.CPTAC <- NULL

#Identifying and getting rid of the duplicate tumor samples:
tumor_names <- colnames(proteome_pam50[2:81])
tumor_names_start <- substr(tumor_names, 1, 7)
duplicated(unlist(strsplit(tumor_names_start, " ")))

#The duplicated values are AO.A12D (1st and 10th on the list), C8.A131 (2nd and 68th on the list), and AO.A12B
#(3rd and 74th on the list). I wasn't sure which replicate was used in the paper's analysis, so I opted to remove
#the first three columns.
proteome_pam50$AO.A12D.01TCGA <- NULL
proteome_pam50$C8.A131.01TCGA <- NULL
proteome_pam50$AO.A12B.01TCGA <- NULL

#Next, going to get rid of the extra isoforms for MELK, MAPT, and MLPH. I'm going to keep the lowest numbered isoform
#in each case. Here, that means that I want to keep MELK-NP_055606 (isoform 1), MAPT-NP_058519 (isoform 1), and
#MLPH-NP_077006 (isoform 1)
#I want to get rid of MELK-NP_001243614 (isoform 2), MAPT-NP_001116539, MAPT-NP_001190181, MAPT-NP_058518, and
#MLPH-NP_001035932
#The rows I want to get rid of are rows 7, 31, 32, 33, 37
proteome_pam50 <-proteome_pam50[-c(7, 31, 32, 33, 37), ] 

#Version with gene names as row names and sample names as column names (gets rid of "gene" column):
prot_pam50_alt <- data.frame(proteome_pam50[2:78], row.names = (proteome_pam50$gene))

#This left 41 genes, while the paper had 35 genes. My next goal was to get rid of any genes that were more than 50%
#NA, per the paper's instructions:
rowSums(is.na(prot_pam50_alt))

#Based on this, if I want to cut genes with more than 50% NA (so more than 38 NA out of 77), I'll want to get rid of
#EXO1-NP_569082, CCNE1-NP_001229, and MYC-NP_002458. That would get me down to 38 genes. These are rows 12, 14, 19:
proteome_pam38 <-prot_pam50_alt[-c(12, 14, 19), ] 

#This still left 38 genes. In an attempt to get rid of three more genes to match the paper's total of 35 genes, I next
#decided to cut out the three remaining genes with the highest NA total. (This ultimately ended up being the equivalent
#of requiring that any included gene have no more than 25% NA.)

#Cutting three more to get to 35 genes would involve getting rid of MELK-NP_055606 (35 NA), FOXC1-NP_001444 (38 NA),
#and BCL2-NP_000624 (20 NA). These are rows 7, 22, 31. Could say that this is getting rid of anything where more than
#a quarter of the samples are NA.

proteome_pam35 <-proteome_pam38[-c(7, 22, 31), ]


#The data was now ready for clustering. I used fanny from the cluster package to do clustering without having to get
#rid of all of the NAs. (Clustering was done with a different package than the one cited in the paper due to
#difficulties.) Using k=3 for clustering.

proteome_pam35_t <- t(proteome_pam35)
prot_pam35_t_fanny_k3 <- fanny(proteome_pam35_t, 3, memb.exp = 1.1)
summary(prot_pam35_t_fanny_k3)
plot(prot_pam35_t_fanny_k3)

#Organizing the clustering info and renumbering the original clusters to match my new numbering system for an easier
#comparison when plotting:
proteome_subtypes_77 <-proteome_subtypes[-c(1, 2, 3, 81, 82, 83), ]

#Add a new column for my new subtype classification that I want to try lining up with my other data. Set Basal = 1, Her2 = 2,
#and LumA and LumB = 3, as that seems most likely to be how the pam35 fanny clusters assigned stuff
proteome_subtypes_77$new_cluster <- 0

#New naming:
for (i in 1:length(rownames(proteome_subtypes_77))) {
  if (proteome_subtypes_77$subtypecode[i] == 0) {
    proteome_subtypes_77$new_cluster[i] <- 1
  }
  if (proteome_subtypes_77$subtypecode[i] == 1) {
    proteome_subtypes_77$new_cluster[i] <- 2
  }
  if (proteome_subtypes_77$subtypecode[i] == 2) {
    proteome_subtypes_77$new_cluster[i] <- 3
  }
  if (proteome_subtypes_77$subtypecode[i] == 3) {
    proteome_subtypes_77$new_cluster[i] <- 3
  }
}

#Plotting the clustering data plus heatmap:
proteome_subtypes_77$Sample <- gsub("-",".",proteome_subtypes_77$Sample) #Changed the names to work better with the heatmap below
l <- proteome_subtypes_77$Sample[order(proteome_subtypes_77$Fanny)]

proteome_subtypes_77$Sample <- factor(proteome_subtypes_77$Sample, levels=l)

prot77_plot_data <- ggplot(proteome_subtypes_77, mapping = aes(Sample, y=2, fill=factor(Fanny)))+
  geom_tile()+
  geom_tile(aes(Sample, y=1, fill=factor(new_cluster))) +
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=rainbow(3)[c(1, 2, 3)])+
  theme(axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank()
  )

prot77_plot_data #This gave the desired clustering rectangles for the samples (sorted by new cluster on top, with
#PAM50 clusters below, and with luminal A and B merged in the PAM50 clusters). Pulled this out and manually placed
#it over my heatmap for viewing.

#Adding a heatmap with samples listed in the same order as in the clustering data:
prot_pam35t_dup <- data.frame(proteome_pam35_t)
prot_pam35t_dup$names <- rownames(prot_pam35t_dup)
prot_pam35t_dup <- prot_pam35t_dup[order(factor(prot_pam35t_dup$names, levels = l)),]
prot_pam35t_dup$names <- NULL

#Tried the heatmap with two slightly different color schemes:

#One palette to try:
my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 299)
heatmap(t(prot_pam35t_dup), Colv=NA, col=my_palette)

#A different palette to try:
RdBuPalette <- brewer.pal(11,"RdBu")
heatmap(t(prot_pam35t_dup), Colv=NA, col=RdBuPalette)



#Preparing the data for GSEA:

#Start with the original proteome spreadsheet and cut it down to just the samples we want:
proteome_77_full <- data.frame(proteome)

#Null the columns I don't want:
proteome_77_full$X263d3f.I.CPTAC <- NULL
proteome_77_full$c4155b.C.CPTAC <- NULL
proteome_77_full$blcdb9.I.CPTAC <- NULL
proteome_77_full$AO.A12D.01TCGA <- NULL
proteome_77_full$C8.A131.01TCGA <- NULL
proteome_77_full$AO.A12B.01TCGA <- NULL

#Moving the gene names to be the row names so the rest of the data frame is just the values
proteome_77_alt <- data.frame(proteome_77_full[2:78], row.names = (proteome_77_full$gene))

#I have 77 samples. I want to eliminate any proteins that have more than 50% NA values. In this case, any with 39+ NA values.

prot_77_filtered <- subset(proteome_77_alt, rowSums(is.na(proteome_77_alt)) <39)

#Next, want to get rid of the NAs. Trying k Nearest Neighbor (knn) imputation, as mentioned in the paper. Using the
#knn.impute method from the bnstruct package without modifications from the default settings.

prot77matrix = as.matrix(prot_77_filtered)
proteome_77_imputed <- knn.impute(prot77matrix)

#Renaming the proteins by their NP numbers (or similar equivalent number):
prot77_imp_df <- data.frame(proteome_77_imputed)
protnames <- rownames(prot77_imp_df)
for (i in 1:length(protnames)) {
  x <- sub("^[^-]*", "", protnames[i])
  x <- substring(x, 2)
  while (grepl("-", x)) { #Added this to deal with the proteins that have multiple hyphens in them.
    x <- sub("^[^-]*", "", x)
    x <- substring(x, 2)
  }
  protnames[i] <- x
}

#Realized I needed gene names instead of NP numbers for use with the GSEA Java app (which is what I decided to use
#for GSEA here). Renamed everything with gene names.
prot77_imp_df_2 <- prot77_imp_df #Cloned this for ease of access when testing things
rownames(prot77_imp_df_2) <- protnames

prot77v3 <- prot77_imp_df_2
prot77v3$names <- rownames(prot77v3)

mart = useDataset("hsapiens_gene_ensembl",useMart("ensembl"))
hh <- listAttributes(mart)
mm <- listFilters(mart)

refseqs <- prot77v3$names
gene_id <- getBM(filters="refseq_peptide",
                 attributes=c("external_gene_name", "refseq_peptide"), values=refseqs,
                 mart=mart)

newNamesProtGeneID <- merge(prot77v3, gene_id,
                            by.x = "names", by.y = "refseq_peptide",
                            all.x = T)

prot77v3forexport <- data.frame(newNamesProtGeneID)
prot77v3forexport$DESCRIPTION <- "na" #Required for the GSEA Java app

#Exporting the data:
write.table(prot77v3forexport,"prot77export.txt",sep="\t", row.names = FALSE)

#Moved the columns around in Excel to match the Broad's description of what the file should look like on their website.
#Imported it into the GSEA app without problems. Also made the .cls file per their instructions in a text editor.
#The rest of the analysis for GSEA was done in the Java app, and consequently is not shown here.
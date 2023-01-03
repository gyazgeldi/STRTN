# TIME-COURSE GENE EXPRESSION ANALYSIS
targetPackages <- c("stringr", "dplyr","ggplot2", "cowplot", "ggbeeswarm", "forcats", "Seurat")
for(package in targetPackages) library(package, character.only = T) 

QC.file <- list.files(pattern = "-QC.txt", full.names = TRUE)
QC <- read.delim(QC.file,header = T,check.names = F)
# Remove outlier
QC <- QC [-c(35),] 

barcode_stages <- read.delim("Example-BarcodesStages.txt", sep="\t", header=F)
barcode_stages$V1 <- gsub('.{11}$', '', barcode_stages$V1)
colnames(barcode_stages) <- c("Barcode", "Stage")
barcode_stages <- barcode_stages [-c(35),]

Mapped_reads <- QC[,c(1,5)]
Spikein_reads <- QC[,c(1,7)]
Spikein_5end_rate <- QC[,c(1,9)]
Mapped_rate <- QC[,c(1,6)]
Coding_5end_rate <- QC[,c(1,12)]

# Mapped reads
dataset <- merge(Mapped_reads, barcode_stages, by = "Barcode")
dataset <- dataset %>% mutate(Developmental_stages=fct_relevel(Stage, 'oocyte', "2cell", '4cell', '8cell', 'morula', 'blastocyst'))
mapped.reads.plot <- ggplot(dataset, aes(x=Developmental_stages, y=Mapped_reads)) + geom_beeswarm() + ggtitle("Mapped reads")

# Spike in reads
dataset <- merge(Spikein_reads , barcode_stages, by = "Barcode")
dataset <- dataset %>% mutate(Developmental_stages=fct_relevel(Stage, 'oocyte', "2cell", '4cell', '8cell', 'morula', 'blastocyst'))
spikein.plot <- ggplot(dataset, aes(x=Developmental_stages, y=Spikein_reads)) + geom_beeswarm() + ggtitle("Spikein reads")

# Spike in 5end rate
dataset <- merge(Spikein_5end_rate, barcode_stages, by = "Barcode")
dataset <- dataset %>% mutate(Developmental_stages=fct_relevel(Stage, 'oocyte', "2cell", '4cell', '8cell', 'morula', 'blastocyst'))
spike5.plot <- ggplot(dataset, aes(x=Developmental_stages, y=`Spikein-5end_rate`)) + geom_beeswarm() + ggtitle("Spikein 5'-end rate")

# Mapped rate
dataset <- merge(Mapped_rate, barcode_stages, by = "Barcode")
dataset <- dataset %>% mutate(Developmental_stages=fct_relevel(Stage, 'oocyte', "2cell", '4cell', '8cell', 'morula', 'blastocyst'))
mapped.rate.plot <- ggplot(dataset, aes(x=Developmental_stages, y=Mapped_rate)) + geom_beeswarm() + ggtitle("Mapped rate")

# Mapped/Spike in
mapped <- (as.numeric(Mapped_reads[,2])-as.numeric(Spikein_reads[,2]))/as.numeric(Spikein_reads[,2])
dataset <- dataset %>% mutate(Developmental_stages=fct_relevel(Stage, 'oocyte', "2cell", '4cell', '8cell', 'morula', 'blastocyst'))
mapped.plot <- ggplot(dataset, aes(x=Developmental_stages, y=mapped)) + geom_beeswarm() + ggtitle("Mapped / Spikein")

# Coding 5end rate
dataset <- merge(Coding_5end_rate, barcode_stages, by = "Barcode")
dataset <- dataset %>% mutate(Developmental_stages=fct_relevel(Stage, 'oocyte', "2cell", '4cell', '8cell', 'morula', 'blastocyst'))
coding5.plot <- ggplot(dataset, aes(x=Developmental_stages, y=`Coding-5end_rate`)) + geom_beeswarm() + ggtitle("Coding 5'-end rate")

# Plotting
pdf(paste0(str_split(QC.file,".txt")[[1]][1],"-BeeswarmPlots.pdf"), width=19, height=9)
plot_grid(mapped.reads.plot, spikein.plot, spike5.plot, mapped.rate.plot, mapped.plot, coding5.plot, nrow = 2, align = "hv")
dev.off()

# scRNA SEQUENCE ANALYSIS USING SEURAT PACKAGE
cts.file <- list.files(pattern = "_byGene-counts.txt$", full.names = TRUE)
cts <- read.delim(cts.file, skip = 1, header=TRUE, sep="\t")
row.names(cts) <- cts[,1]
# Remove features columns
cts <- cts[, -c(1,2,3,4,5,6)]
# Remove outliers
cts <- cts [, -c(1, 47, 48, 11, 27, 35)]
strt.seurat.obj <- CreateSeuratObject(counts = cts, project = "STRT", min.cells = 3, min.features = 200)

# Adding sample developmental stage information
barcode_stages <- read.delim("Example-BarcodesStages.txt", sep="\t", header=F)
# Remove outliers
barcode_stages <- barcode_stages [-c(1, 42, 43, 3, 20, 29),]
strt.seurat.obj@meta.data[["orig.ident"]] <- barcode_stages[,2]

# Spike-in normalization
strt.seurat.obj <- NormalizeData(strt.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
reads.p1 <- cts + 1
reads.p1.spikes <- colSums(reads.p1[1:97,])
strt.seurat.obj[["RNA"]]@data <- as(as.matrix(log(reads.p1/rep(reads.p1.spikes, each=nrow(reads.p1)))), "dgCMatrix")

strt.seurat.obj <- FindVariableFeatures(strt.seurat.obj, selection.method = "vst", nfeatures = 2000)
all.genes <- row.names(strt.seurat.obj)
strt.seurat.obj <- ScaleData(strt.seurat.obj, features = all.genes)

# PCA plotting
strt.seurat.obj <- RunPCA(strt.seurat.obj, features = VariableFeatures(object = strt.seurat.obj), npcs = 33)
strt.seurat.obj <- JackStraw(strt.seurat.obj, num.replicate = 100)
strt.seurat.obj <- ScoreJackStraw(strt.seurat.obj, dims = 1:20)
JackStrawPlot(strt.seurat.obj, dims = 1:20)
ElbowPlot(strt.seurat.obj)
DimPlot(strt.seurat.obj, reduction = "pca", group.by = "orig.ident", label = T)

# UMAP plotting
strt.seurat.obj <- RunUMAP(strt.seurat.obj, dims = 1:4, verbose = F)
DimPlot(strt.seurat.obj,label.size = 5,repel = T, label = T, group.by = "orig.ident")

# Violin plotting for interested genes in pre-implantation development as an example
features <- c("Actb", "Actg1", "Actg2", "Hprt", "Tuba1b", "Tuba1c", "Nanog", "Pou5f1", "Gdf9")
strt.seurat.obj$orig.ident <- factor(strt.seurat.obj$orig.ident, levels = c('oocyte', "2cell", '4cell', '8cell', 'morula', 'blastocyst'))
VlnPlot(object = strt.seurat.obj, features, group.by = "orig.ident", ncol=3)
# Read the count matrix and manipulate dataframe
cts.file <- list.files(pattern = "_byGene-counts.txt$", full.names = TRUE)
cts <- read.delim(cts.file, skip = 1, header=TRUE, sep="\t")
row.names(cts) <- cts[,1]

# Remove columns except gene name and expression values
cts <- cts[, -c(1,2,3,4,5,6)]
scale_factor_df <- cts[1:97,]

# Calculate the scale factor for each sample
scale_factor_df <- as.data.frame(1/colSums(scale_factor_df))
scale_factor_df$samples <- rownames(scale_factor_df)
rownames(scale_factor_df) <- NULL
scale_factor_df <- scale_factor_df[,c(2,1)]
colnames(scale_factor_df) <- c("samples", "scale_factor")
write.table(scale_factor_df, file = "scale_factor_df.txt", sep = "\t", row.names = FALSE)

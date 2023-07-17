# First, we will import our data and store it an an object or data frame(df)
df <- read.table('data.txt', sep ='\t', header = TRUE, row.names=1)

#write the first 10 rows of the dataset:
head(df)
# Check dimension of data (how many rows and columns)
dim(df)

# Check type of data
str(df)

# Extract the sample names (column names)
Samples <- colnames(df)
print(Samples)

Genes  <- rownames(df)
print(Genes)

# remove NA values
df1 <-  na.exclude(df)
# To check the dimension of data
dim(df1)

#Descriptive Statistics and Visualization

# Descriptive statistics, i.e., compute the summary statistics of data
summ <- summary(df1)
# print summary statistics
print(summ)
# Write summary statistics of data  into a file and export it
write.table(summ, file="data.txt", col.names=TRUE, sep="\t")


#Box plot

# Draw boxplot for all samples (for FPKM values)
boxplot(df1)
# Add the main title, axis titles, and color 
boxplot(df1, main="Boxplot for FPKM data", xlab="", ylab="Gene expression (FPKM)", col="red", las=2, cex.axis = 0.65)

log_df <- log(df1+1)
#Check the dimension of data
dim(log_df)

# let’s check the summary statistics of data after log transformation
summ_log <- summary(log_df)
# Print summary statistics of data 
print(summ_log)
# Write summary statistics to a file and export it
write.table(summ_log, file="stat_sum-log-data.txt", col.names=TRUE, sep="\t")


# Draw boxplot for all samples (for FPKM values)
boxplot(log_df, main="Boxplot for log-transformed Data", xlab="", ylab="Gene expression (log[FPKM])", col="red", las=2, cex.axis = 0.7)

Group1 <- log_df[,1:5]

Group1 <- log_df[grep('^Non.malignant', names(log_df))]

# Extract specific samples from data with their name
Non_ML1 <- as.data.frame(log_df[grep('^Non.malignant', names(log_df))])
Claudin1 <- as.data.frame(log_df[grep('^Claudin', names(log_df))])


#Compute the mean for each gene within each group
Non_ML<- rowMeans(Non_ML1)
Claudin <- rowMeans(Claudin1)

#Bind groups together into 2 columns
group <- cbind(Non_ML, Claudin)

#Provide the column names to both groups
colnames(group) <- c("Non-malignant", "Claudin-low")

#Draw boxplot for groups of samples
boxplot(group, main="Boxplot for groups", xlab="Groups", ylab="Gene expression(log[FPKM])", col=c('cyan', 'pink'))

# Extract samples of same group
sample1 <- df1[,1]
sample2 <- df1[,2]
# Combine them column-wise
s1_2 <- cbind(sample1,sample2)
# Draw Scatterplot
plot(s1_2)

# Extract samples of different groups
sample3 <- df1[,1]
sample4 <- df1[,8]
# Combine them column-wise
s3_4 <- cbind(sample3,sample4)
# Draw Scatterplot
plot(s3_4)

# Extract samples of same group
sample1 <- log_df[,1]
sample2 <- log_df[,2]
# Combine them column-wise
s1_2 <- cbind(sample1,sample2)
# Draw Scatter Plot
plot(s1_2)

# Extract samples of different groups
sample3 <- log_df[,1]
sample4 <- log_df[,8]
# Combine them column-wise
s3_4 <- cbind(sample3,sample4)
# Draw Scatter Plot
plot(s3_4)


# draw scatter plot for whole data (with FPKM values)
plot(df1)

# Scatter Plot for log transformed data
plot(log_df)

# Take the 2nd sample
sample1 <- log_df[,2]
#Make a bar plot
barplot(sample1, main="Barplot of Sample", xlab="Genes", ylab="Gene expression (log[FPKM])")

# Extract 2nd sample with first 20 genes
sample_20genes <- log_df[1:20,2]
# Draw barplot
barplot(sample_20genes,  main="Barplot of Sample for 20 genes", xlab="Genes", ylab="Gene expression (log[FPKM])")

# Add colors to the barplot
barplot(sample_20genes,  main="Barplot of Sample for 20 genes", xlab="Genes", ylab="Gene expression (log[FPKM])", col = rainbow(length(sample_20genes)))

# Extract 135th gene from the data
gene_135 <- log_df[135,]
# Convert this data frame into a matrix
gene_135_m <- as.matrix(gene_135)
# Extract gene ids list from the data
Genes  <- rownames(log_df)
# Extract the gene id for 135th gene from list
gene_id <- Genes[135]
# Draw barplot for the selected gene
barplot(gene_135_m,  main= gene_id, xlab="", ylab="Gene expression (log[FPKM])", las = 2, cex.names = 0.7, col = "red")




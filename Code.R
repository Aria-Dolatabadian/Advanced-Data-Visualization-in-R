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


#Convert dataframe (with log values) into a matrix
﻿log_df1 <- as.matrix(log_df) 
#Draw histogram, breaks represents the bins 
﻿hist(log_df1,  main= "Histogram for all samples" , xlab= "samples" , breaks= 80 , col=rainbow( 11 ))


#Extract the sample 
﻿Non_malignant <- as.matrix(log_df[ 1 ])

#Draw histogram 
﻿hist(Non_malignant, col= "blue" )

#Extract sample 
﻿Claudin_low <- as.matrix(log_df[6])

#Plot histogram for claudin- low samples 
﻿hist(Claudin_low, col="blue")

#Extract desired  samples, here, we are extracting sample 1 and sample 8 
﻿sample1 <- log_df1[, 1 ]
sample2 <- log_df1[, 8 ]

#Define colors for histogram 
﻿col1 <- rgb( 0 , 0 , 1 , 0.5 )
col2 <- rgb( 1 , 0 , 0 , 0.5 )


#Draw histograms 
﻿hist(sample1, col=col1, main = "Histograms of two samples", xlab= "", breaks= 200 )
hist(sample2, col=col2, main= "", xlab= "", breaks= 200 , add=T)

#Add density for each sample 
﻿d1 <- density(sample1)
d2 <- density(sample2)

#Plot Density plot 
﻿plot(d2, col="red")
polygon(d2, col=col2, border=col2)
lines(d1, col="blue" )
polygon(d1, col=col1, border=col1)

# Transform dataframe into a matrix
df_mat <- as.matrix(df1)
#Create heatmap
heatmap(df_mat, cexCol= 0.8 )

#Convert into a matrix
df_log_mat <- as.matrix(log_df)

#Plot heatmap
heatmap(df_log_mat, cexCol= 0.8)




#Load libraries 
﻿library(ggplot2) 
﻿library(reshape2)


df <- read.table( "data2.txt", header=TRUE,  sep="\t" , row.names= 1 )


#Remove "NA" values from the data 
﻿df1 <- na.exclude(df) 
#check the dimensions of data 
﻿dim(df1)

#Compute summary statistics of data 
﻿summ = summary(df1)

#print summary statistics 
﻿print(summ)


#Write summary statistics  into a file 
﻿write.table(summ, file= "stat_sum.txt", col.names=TRUE, sep= "\t" )

#Change the dataframe into long data frame format (which will be used by ggplot) 
﻿df_f<- melt(df1)


#Draw Boxplot using ggplot function 
﻿b <- ggplot(df_f, aes(factor(variable), value)) + geom_boxplot(aes(fill = variable))

#Add axis title, legend title and adjust font size 
﻿b + theme(legend.key.size = unit(0.9, "cm"), legend.text = element_text(color= "Black", size= 9), axis.text.x = element_text(color= "Black", size= 9, angle= 90), axis.text.y = element_text(size= 10, angle= 45)) + xlab("Samples") + ylab("Gene Expression(FPKM)") + labs(fill = "Samples")  

# Transform FPKM values into log values 
﻿log_df <- log(df1+1)

#Descriptive statistics  
﻿log_summary <- summary(log_df)


#write into a file 
﻿write.table(log_summary, file= "stat_log_sum.txt" , col.names=TRUE, sep= "\t")

# Reshape or transform data into a long data frame format 
﻿df_log <- melt(log_df)


#Draw ﻿boxplot 
p <- ggplot(df_log, aes(factor(variable), value))  + geom_boxplot(aes(fill = variable))

#Add axis title, legend title and adjust font size 
﻿p + theme(legend.key.size = unit(0.9 , "cm"), legend.text = element_text(color= "Black", size= 9), axis.text.x = element_text( color= "Black", size= 9, angle= 90), axis.text.y = element_text(size= 10, angle= 45)) + xlab("Samples") + ylab("Gene Expression(Log[FPKM])" ) + labs(fill = "Samples")  

#Extract samples belong to each group with name 
﻿Non_ML1<- as.data.frame(log_df[grep('^Non', names(log_df))])
Claudin1 <- as.data.frame(log_df[grep('^Claudin', names(log_df))])

#Compute the mean of samples for each gene within each group 
﻿Non_ML<- rowMeans(Non_ML1)
Claudin <- rowMeans(Claudin1)

#Bind groups together column-wise 
﻿group <- cbind(Non_ML, Claudin)

#Provide the column names to both groups of samples 
﻿colnames(group) <- c("Non-malignant", "Claudin-low")

#Create a long formatted dataframe 
﻿df_group <- melt(group)
head(df_group)


#Draw boxplot  
﻿pp <- ggplot(df_group, aes(factor(Var2), value)) + geom_boxplot(aes(fill = Var2))

#Add axis title, legend title and adjust font size 
﻿pp + theme(legend.key.size = unit(0.9, "cm"), legend.text = element_text(color= "Black", size= 9), axis.text.x= element_text(color= "Black" , size= 9, angle=90), axis.text.y = element_text(size= 10 , angle= 45)) + xlab("Groups") + ylab("Gene Expression (Log[FPKM])") + labs(fill = "Group") 






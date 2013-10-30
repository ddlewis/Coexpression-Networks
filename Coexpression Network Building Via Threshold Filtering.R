#Coexpression Network Mapping by Correlation Threshold

#This program is designed to filter an existing correlation table
#and give you a list of gene pairs with high (or low) correlation values
#the correlation table was originally taken from the ATTED II datbase
#Ted filtered the table so that the columns are genes in mallorie's yeast 1 hybrid network (~225)
#and the rows are every gene in the Arabidopsis genome (~20,000)


#imports coexpession correlation coefficients between 
#all arabidopsis genes and genes identified by Y1H network
corr.coex <- read.table("SourceTargetList.Corrs.tsv", stringsAsFactors = F)
#convert NA to 0
corr.coex[is.na(corr.coex)] = 0
#need to make a new data frame that you can manipulate
corr.coex.list <- corr.coex
#Add another column to the list to make IDing genes easier
corr.coex.list$rowid <- rownames(corr.coex)
#need to make a set of data that has PCC values in a column
#THIS COMMAND TAKES FOREVER-- I ran it overnight and it was done in the morning
#the reshape command will format the table so that each interaction is a pair
#each row ID has already been listed in the column rowid
#each column ID will be listed under the argument timevar = "colid"
#the new column holding the values formerly arranged in a matrix is named by v.names = "pcc"
#for a helpful overview of the reshape command, go to http://www.ats.ucla.edu/stat/r/faq/reshape.htm
corr.coex.list <- reshape(corr.coex.list, varying = y1hnames, v.names = "pcc", timevar = "colid", times = y1hnames, idvar = "rowid", direction = "long", sep = "")
#use the subset command to take the columns that you want (using select)
#and the values that you want PCC >= threshold
seven.corr.coex.list <- subset(corr.coex.list, pcc >= .7, select= c(rowid, colid, pcc))
eight.corr.coex.list <- subset(corr.coex.list, pcc >= .8, select = c(rowid, colid, pcc))
write.csv(eight.corr.coex.list, "0.8 correlation list.csv")
write.csv(seven.corr.coex.list, "0.7 correlation list.csv")
#to subset the data to find genes that are negatively correlated, use pcc <= threshold
negfive.corr.coex.list <- subset(corr.coex.list, pcc <= -.5, select = c(rowid, colid, pcc))
negthree.corr.coex.list <- subset(corr.coex.list, pcc <= -.3, select = c(rowid, colid, pcc))
write.csv(negfive.corr.coex.list, "-.5 correlation list.csv")
write.csv(negthree.corr.coex.list, "-.3 correlation list.csv")

#Translating table, for making sense of the ATG numbers

#translating correlation list to gene function
#sep designates the data separating character, "\t" designates it as a tab
#this is a huge table, you may want to delete it after translating your lists
translate <- read.table("TAIR10_GeneInfo.tsv", header = T, sep = "\t", stringsAsFactors = F)
#get rid of all rows except gene and symbols
translate <- translate[, c("gene", "symbols")]
#make a new column that can be manipulated
gene.symbols <- translate$symbols
#get rid of periods in gene names to facilitate matching between lists
names(gene.symbols) <- sub("\\..*$", "", translate$gene)
#load list --- MAKE SURE TO GET RID OF STRINGS AS FACTORS!!!
eight.corr.coex.list <- read.csv("0.8 correlation list.csv", row.names = 1, stringsAsFactors = F)
#Add two new columns translating the gene names-- find a way to only display one names?
#couldn't think of good column names
eight.corr.coex.list$from <- gene.symbols[eight.corr.coex.list$rowid]
eight.corr.coex.list$to <- gene.symbols[eight.corr.coex.list$colid]
write.csv(eight.corr.coex.list, "translated .8 corr table.csv", row.names = F)

seven.corr.coex.list <- read.csv("0.7 correlation list.csv", stringsAsFactors = F)
seven.corr.coex.list$from <- gene.symbols[seven.corr.coex.list$rowid]
seven.corr.coex.list$to <- gene.symbols[seven.corr.coex.list$colid]
write.csv(seven.corr.coex.list, "translated .7 corr table.csv", row.names = F)

negfive.corr.coex.list <- read.csv("-.5 correlation list.csv", row.names = 1, stringsAsFactors = F)
negfive.corr.coex.list$from <- gene.symbols[negfive.corr.coex.list$rowid]
negfive.corr.coex.list$to <- gene.symbols[negfive.corr.coex.list$colid]
write.csv(negfive.corr.coex.list, "translated -.5 corr table.csv", row.names = F)

negthree.corr.coex.list <- read.csv("-.3 correlation list.csv", row.names = 1, stringsAsFactors = F)
negthree.corr.coex.list$from <- gene.symbols[negthree.corr.coex.list$rowid]
negthree.corr.coex.list$to <- gene.symbols[negthree.corr.coex.list$colid]
write.csv(negthree.corr.coex.list, "translated -.3 corr table.csv", row.names = F)

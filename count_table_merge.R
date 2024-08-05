# PIMANNIE Work
# 30 July 2024


library(tidyverse)
library(corrplot)
library(RColorBrewer)
library(heatmaply)

# PART_1: Create combined CT from individual CTs

# Read in the CTs and arrange by 'sgRNA'
ct_1 <- read_tsv("04_Working/01_49117_1/49117_1.count.txt")
head(ct_1)
colnames(ct_1) <- c('sgRNA', 'Gene', 'ct1_counts')
head(ct_1)
ct_1 <- ct_1 %>%
  arrange(sgRNA)
head(ct_1)

# Read in CT_2 and merge with the first one
ct_2 <- read_tsv("04_Working/02_49117_2/49117_2.count.txt")
head(ct_2)
colnames(ct_2) <- c('sgRNA', 'Gene', 'ct2_counts')
head(ct_2)
ct_2 <- ct_2 %>%
  arrange(sgRNA)
head(ct_2)
ct_all <- full_join(ct_1, ct_2)
head(ct_all)

# Repeat for CT_3
ct_3 <- read_tsv("04_Working/03_49117_3/49117_3.count.txt")
head(ct_3)
colnames(ct_3) <- c('sgRNA', 'Gene', 'ct3_counts')
head(ct_3)
ct_3 <- ct_3 %>%
  arrange(sgRNA)
head(ct_3)
ct_all <- full_join(ct_all, ct_3)
head(ct_all)

# Repeat for CT_4
ct_4 <- read_tsv("04_Working/04_49117_4/49117_4.count.txt")
head(ct_4)
colnames(ct_4) <- c('sgRNA', 'Gene', 'ct4_counts')
head(ct_4)
ct_4 <- ct_4 %>%
  arrange(sgRNA)
head(ct_4)
ct_all <- full_join(ct_all, ct_4)
head(ct_all)

# Repeat for CT_5
ct_5 <- read_tsv("04_Working/05_49117_5/49117_5.count.txt")
head(ct_5)
colnames(ct_5) <- c('sgRNA', 'Gene', 'ct5_counts')
head(ct_5)
ct_5 <- ct_5 %>%
  arrange(sgRNA)
head(ct_5)
ct_all <- full_join(ct_all, ct_5)
head(ct_all)

# Repeat for CT_6
ct_6 <- read_tsv("04_Working/06_49117_6/49117_6.count.txt")
head(ct_6)
colnames(ct_6) <- c('sgRNA', 'Gene', 'ct6_counts')
head(ct_6)
ct_6 <- ct_6 %>%
  arrange(sgRNA)
head(ct_6)
ct_all <- full_join(ct_all, ct_6)
head(ct_all)

# Repeat for CT_7
ct_7 <- read_tsv("04_Working/07_49117_7/49117_7.count.txt")
head(ct_7)
colnames(ct_7) <- c('sgRNA', 'Gene', 'ct7_counts')
head(ct_7)
ct_7 <- ct_7 %>%
  arrange(sgRNA)
head(ct_7)
ct_all <- full_join(ct_all, ct_7)
head(ct_all)

# Repeat for CT_8
ct_8 <- read_tsv("04_Working/08_49117_8/49117_8.count.txt")
head(ct_8)
colnames(ct_8) <- c('sgRNA', 'Gene', 'ct8_counts')
head(ct_8)
ct_8 <- ct_8 %>%
  arrange(sgRNA)
head(ct_8)
ct_all <- full_join(ct_all, ct_8)
head(ct_all)

# Repeat for CT_9
ct_9 <- read_tsv("04_Working/09_49117_9/49117_9.count.txt")
head(ct_9)
colnames(ct_9) <- c('sgRNA', 'Gene', 'ct9_counts')
head(ct_9)
ct_9 <- ct_9 %>%
  arrange(sgRNA)
head(ct_9)
ct_all <- full_join(ct_all, ct_9)
head(ct_all)

# Repeat for CT_10
ct_10 <- read_tsv("04_Working/10_49117_10/49117_10.count.txt")
head(ct_10)
colnames(ct_10) <- c('sgRNA', 'Gene', 'ct10_counts')
head(ct_10)
ct_10 <- ct_10 %>%
  arrange(sgRNA)
head(ct_10)
ct_all <- full_join(ct_all, ct_10)
head(ct_all)



# Write full CT out to file - load this for later work
write_tsv(ct_all, "01_Scripts/ct_all.txt")



# PART_2:  Create Correlation Plot

# Load complete CT from last section
ct_all <- read_tsv("01_Scripts/ct_all.txt")
head(ct_all)

# Need to remove 'sgRNA' and 'Gene' columns first
ct_all_2 <- ct_all[,3:12]
head(ct_all_2)

# Create correlation matrix
M = cor(ct_all_2)
head(round(M, 2))

# Create plot - try various versions. See notes at end of script for web pages on this.
# Are ways to highlight statistically significant pairs - see notes
# NB: Need to look at this in light of PCR sample info - makes more sense.  Need hierarchical clustering to see relationships.


# With hierarchical clustering and setting nos. clusters using 'addrect'
corrplot(M, method="color", outline=T, addgrid.col = "darkgray", order="hclust", addrect=6, rect.col="black", rect.lwd=5, addCoef.col = "black",
         cl.pos="b", tl.col="indianred4", number.cex=0.75)

# Another version of the above - better?
corrplot(M, method="color", outline=T, addgrid.col="darkgrey", order="hclust", addrect=6, rect.col="black", rect.lwd = 5,
         cl.pos="b", tl.col="indianred4", addCoef.col = "white", number.digits = 2, number.cex = 0.75, col = colorRampPalette(c("darkred", "lightgrey", "midnightblue"))(100))


# Re-Do Plot With Renamed_Samples
# ***THINK THIS IS THE BEST ONE***
# Renaming samples according to first- and second-round PCRs, and also by vector - allows easier interpretation of results

# Make a copy of 'ct_all_2'
ct_all_3 <- ct_all_2
head(ct_all_3)

# Rename columns
colnames(ct_all_3) <- c('BFP_1-1', 'BFP_2-1', 'BFP_2-2', 'BFP_3-1', 'BFP_3-2',
                        'Cherry_1-1', 'Cherry_2-1', 'Cherry_2-2', 'Cherry_3-1', 'Cherry_3-2')
head(ct_all_3)

# Create correlation matrix
M2 <- cor(ct_all_3)
head(M2)

# Create correlation plot
corrplot(M2, title = "Correlation Plot: VBC Library Preps", mar=c(0,0,2,0), method="color", outline=T, addgrid.col="darkgrey", order="hclust",
         addrect=6, rect.col="red", rect.lwd = 5, cl.pos="b", tl.col="indianred4", addCoef.col = "white", number.digits = 2, number.cex = 0.75,
         col = colorRampPalette(c("darkred", "lightgrey", "midnightblue"))(100))

# Shows that main determinant of clustering is (a) vector used (BFP vs Cherry) and (b) how many first-round PCRs were used
# Nos second-round PCRs doesn't seem to have as much - if any - of an impact on clustering




# PART_3:  Individual Density Plots
# Try this for ct_5 and ct_10 initially - trying out several conditions

# Using CT_5 - doesn't really show anything useful
plot5 <- ggplot(ct_5, aes(x=ct5_counts)) + theme_classic() +
  geom_density() +
  ggtitle("Distribution of raw read counts for CT5 (BFP_3-2)") +
  xlab('count distribution')
plot5

# CT5 - try logging the counts
# Gives warning message - 93182 rows removed (presumably because of '0' counts)
plot5 <- ggplot(ct_5, aes(x=log10(ct5_counts))) + theme_classic() +
  geom_density() +
  ggtitle("Distribution of log10 read counts for CT5 (BFP_3-2)") +
  xlab('log10 count distribution')
plot5

# CT5 - retry above, but add pseudocount of 1, to avoid problem with '0' values
# It works, but now plot is not hugely informative
plot5 <- ggplot(ct_5, aes(x=log10(ct5_counts +1))) + theme_classic() +
  geom_density() +
  ggtitle("Distribution of log10 (read counts + 1) for CT5 (BFP_3-2)") +
  xlab('log10 (count + 1) distribution')
plot5


# Try above with CT10
head(ct_10)

# Using CT_10 - doesn't really show anything useful
plot10 <- ggplot(ct_10, aes(x=ct10_counts)) + theme_classic() +
  geom_density() +
  ggtitle("Distribution of raw read counts for CT10 (Cherry_3-2)") +
  xlab('count distribution')
plot10

# CT10 - try logging the counts
# Gives warning message - 93182 rows removed (presumably because of '0' counts)
plot10 <- ggplot(ct_10, aes(x=log10(ct10_counts))) + theme_classic() +
  geom_density() +
  ggtitle("Distribution of log10 read counts for CT10 (Cherry_3-2)") +
  xlab('log10 count distribution')
plot10

# CT10 - retry above, but add pseudocount of 1, to avoid problem with '0' values
# It works, but now plot is not hugely informative
plot10 <- ggplot(ct_10, aes(x=log10(ct10_counts +1))) + theme_classic() +
  geom_density() +
  ggtitle("Distribution of log10 (read counts + 1) for CT10 (Cherry_3-2)") +
  xlab('log10 (count + 1) distribution')
plot10



# PART_4:  Combined Density plots
# Need to remodel the combined CT using 'pivot_longer' etc

head(ct_all)
# Make copy of 'ct_all'
ct_all_4 <- ct_all
head(ct_all_4)
# Rename columns
colnames(ct_all_4) <- c('sgRNA', 'Gene', 'BFP_1-1', 'BFP_2-1', 'BFP_2-2', 'BFP_3-1', 'BFP_3-2',
                        'Cherry_1-1', 'Cherry_2-1', 'Cherry_2-2', 'Cherry_3-1', 'Cherry_3-2')
head(ct_all_4)


# Write 'ct_all_4' out to file - load this for later work
write_tsv(ct_all_4, "01_Scripts/ct_all_4.txt")


# Remodel the CT
ct_all_4 <- ct_all_4 %>%
  pivot_longer(cols = 3:12, names_to = 'PCR', values_to = 'Counts')
head(ct_all_4, n = 20)


# Create combined plot
p <- ggplot(ct_all_4, aes(x=log10(Counts))) + theme_classic() +
  geom_density() +
  geom_vline(aes(xintercept=4), color="blue", linetype="dashed") +
  labs(title = "VBC library preps: log10 distribution of read counts") +
  theme(strip.background = element_rect(fill = "grey"),
        plot.title = element_text(size = 20)) +
  facet_wrap(vars(PCR), scales='fixed', ncol = 5)
p




###########################










# Notes

# 1. Web pages on how to create correlation plots
# https://rpubs.com/melike/corrplot
# http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram
# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html

# 2. Interactive heatmap info here: https://www.datanovia.com/en/blog/how-to-create-an-interactive-correlation-matrix-heatmap-in-r/
heatmaply_cor(M, main="", margins=c(80,80,80,80), k_col = 2, k_row = 2)

# 3. These are good, but don't really let you see which ones cluster, as don't have 'r' info
corrplot(M, method="circle", outline=T, addgrid.col = "darkgray")
# As above, but with hierarchical clustering and use of 'addrect' to select nos cluster groups
corrplot(M, method="circle", order='hclust', addrect=6, rect.lwd=5, outline=T, addgrid.col = "darkgray")




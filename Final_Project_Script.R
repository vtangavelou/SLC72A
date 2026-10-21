
#Data Preparation

t_eso_sub <-t(eso_sub)
t_eso_sub <- row_to_names(t_eso_sub, row_number = 1)
#Transposing the data chart to make plotting easier.GeeksForGeeks
#Row to name was utilized by library(janitor) to move barcode to name

t_eso_sub <- as.data.frame(t_eso_sub) #Convert new matrix to df
t_eso_sub[] <- lapply(t_eso_sub, as.numeric) #Convert char to num variables

#ENSG00000003989.18 - *Main Selected Genes*
#ENSG00000002079.14 - V All Genes for Heat Map V
#ENSG00000005893.16 *Gene for Scatter*
#ENSG00000001629.10
#ENSG00000001617.12
#ENSG00000005001.10
#ENSG00000117461.15
#ENSG00000117472.10
#ENSG00000117533.15
#ENSG00000117543.21 


summary(t_eso_sub$ENSG00000003989.18) #Summary analysis for selected gene
sd(t_eso_sub$ENSG00000003989.18) #SD for selected gene
summary(t_eso_sub$ENSG00000005893.16) #Summary analysis for 2nd gene
sd(t_eso_sub$ENSG00000005893.16) #SD for 2nd gene

#Histogram
Histogram <- ggplot(t_eso_sub, aes(x = ENSG00000003989.18)) +
  geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
  scale_x_continuous(labels = scales::comma)
  labs(title = "Count Data of Gene ENSG00000003989.18", 
       x = "Frequency", 
       y = "RNA-Seq Read Count")

#Utilizing ggplot to generate a histogram based on the count
#of the main gene ENSG00000003989.18. The gene shows a heavy
#right-skew in the data.

#Scatterplot
Scatterplot <- ggplot(t_eso_sub, aes(x = ENSG00000003989.18, y = ENSG00000005893.16)) +
  geom_point(alpha = 0.4, color = "blue") +
  scale_x_continuous(labels = scales::comma) +
  labs(title = "Scatter plot of ENSG00000003989.18 and ENSG00000005893.16 RNA-Seq Read Count",
       x = "ENSG00000003989.18",
       y = "ENSG00000005893.16") +
  theme_minimal() 

#Utilizing ggplot to generate a scatterplot between the two selected genes.
#The Scatterplot compares the two counts based on barcode and a relatively
#normal distribution

t_eso_sub_m <- rownames_to_column(t_eso_sub, var = "barcode")
#Converts barcode to column to prepare for combining data frames based on barcode

meta_data$barcode <- gsub("-", ".", meta_data$barcode)
#Converts the "-" in the meta_data df to "." so the sample barcodes will match

t_eso_sub_meta <- list(meta_data, t_eso_sub_m) %>% 
  reduce(full_join, by = "barcode")
#Combines and creates a new data frame based on the barcodes

#Violin Plot
violin_plot <- ggplot(t_eso_sub_meta, aes(x = tissue_type, y = ENSG00000005893.16)) +
  geom_violin(fill = "red", alpha = 0.7) +
  labs(title = "Violin Plot of ENSG00000003989.18 by Tumor Type ", 
       x = "Tumor Type", 
       y = "RNA-Seq Read Count")
#ggplot used to generate a violin plot. Visualizes RNA-seq count of the main gene
#based on whether the sample tissue is a tumor or normal.


#Heatmap Df

library(ComplexHeatmap) #Packages to generate a heatplot
library(circlize)

genes_to_plot <- c("ENSG00000003989.18",
                   "ENSG00000002079.14",
                   "ENSG00000005893.16",
                   "ENSG00000001629.10",
                   "ENSG00000001617.12",
                   "ENSG00000005001.10",
                   "ENSG00000117461.15",
                   "ENSG00000117472.10",
                   "ENSG00000117533.15",
                   "ENSG00000117543.21") #Makes list of genes to plot

counts_subset <- t_eso_sub_meta[, genes_to_plot ]
#Generates new subset based on genes_to_plot list
counts_subset <- t(as.matrix(counts_subset))
#Transposes the new subset as a matrix

column_ha <- HeatmapAnnotation(
  Tumor_Type = t_eso_sub_meta$tissue_type,
  col = list(
    Tumor_Type = c(Tumor = "red", Normal = "blue")
  )
) #Annotation of Heat Map Based on Tumor Type of Tumor or Normal

Heatmap(
  counts_subset,
  name = "Counts",
  top_annotation = column_ha,
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  column_title = "Recorded Specimens",
  row_title = "Genes",
  heatmap_legend_param = list(title = "Count"))

#Generates Heat Map of the selected genes in the subset based on the count
#and identified by whether the tissue sample is normal or tumor.



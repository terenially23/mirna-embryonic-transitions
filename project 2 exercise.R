# Read the file into R
deseq2_results <- read.delim("Galaxy38-[DESeq2_result_file_on_data_34,_data_32,_and_others].tabular", header = FALSE, sep = "\t")

# Check the data to ensure it loaded correctly
head(deseq2_results)

# Renaming multiple columns
colnames(deseq2_results) <- c("GeneID", "BaseMean", "log2FC", "StdErr", "WaldStats", "PValue", "PAdj")

# Remove rows with any missing values
deseq2_results <- na.omit(deseq2_results)

# Get the total number of microRNAs
total_miRNAs <- nrow(deseq2_results)
print(total_miRNAs)

# Filter for significantly differentially expressed microRNAs (PAdj < 0.05)
significant_miRNAs <- deseq2_results[deseq2_results$PAdj < 0.05, ]

# Count the number of significant microRNAs
num_significant_miRNAs <- nrow(significant_miRNAs)
print(paste("Total significant microRNAs:", num_significant_miRNAs))

# Identify up-regulated microRNAs (log2FC > 0)
upregulated_miRNAs <- significant_miRNAs[significant_miRNAs$log2FC > 0, ]
num_upregulated_miRNAs <- nrow(upregulated_miRNAs)

# Identify down-regulated microRNAs (log2FC < 0)
downregulated_miRNAs <- significant_miRNAs[significant_miRNAs$log2FC < 0, ]
num_downregulated_miRNAs <- nrow(downregulated_miRNAs)

# Print the results
print(paste("Number of up-regulated microRNAs:", num_upregulated_miRNAs))
print(paste("Number of down-regulated microRNAs:", num_downregulated_miRNAs))

# Show GeneIDs of up-regulated and down-regulated microRNAs
print(paste("GeneIDs of up-regulated microRNAs:", upregulated_miRNAs$GeneID))
print(paste("GeneIDs of down-regulated microRNAs:", downregulated_miRNAs$GeneID))


#install.packages("ggplot2")
# Assuming your data frame is named 'miRNA_data'
# Load necessary libraries
library(ggplot2)

# Assuming significant miRNAs are those with Padj < 0.05 and |log2FC| > 1
# Create a new column that categorizes miRNAs into upregulated, downregulated, or not significant
deseq2_results$Expression <- ifelse(deseq2_results$PAdj < 0.05 & deseq2_results$log2FC > 1, "Significantly Upregulated",
                                ifelse(deseq2_results$PAdj < 0.05 & deseq2_results$log2FC < -1, "Significantly Downregulated", "Not Significant"))

# Volcano Plot with Color-Coding and Labels
ggplot(deseq2_results, aes(x = log2FC, y = -log10(PAdj), color = Expression, label = ifelse(Expression != "Not Significant", GeneID, ""))) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significantly Upregulated" = "red", "Significantly Downregulated" = "blue")) +
  labs(,
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  ) +
  theme_minimal() +
  theme(legend.position = "top",   
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Fold-change threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Significance threshold line
  geom_text(aes(label = ifelse(Expression != "Not Significant", GeneID, "")), 
            vjust = 1.5, hjust = 0.5, size = 3, show.legend = FALSE)  # Labels without legend

ggsave("volcano_plot_high_res.png", width = 10, height = 7, dpi = 300)



##############################

# install.packages("readxl")
# install.packages("dplyr")
# install.packages("openxlsx")

# Load necessary libraries
library(readxl)
library(dplyr)
library(openxlsx)  # For saving to Excel if needed

# List all files in the folder that have .xlsx extension
file_list <- list.files(pattern = "*.xlsx", full.names = TRUE)

# Create an empty list to store filtered data from each file, with each element named by the file
filtered_data_list <- list()

# Create a data frame to store the counts before and after filtering
summary_counts <- data.frame(
  miRNA = character(),
  genes_before = integer(),
  genes_after = integer(),
  stringsAsFactors = FALSE
)

# Loop through each file
for (file in file_list) {
  # Read the data
  data <- read_excel(file)
  
  # Clean up column names: remove whitespace and special characters
  colnames(data) <- gsub("\\s+", "_", trimws(colnames(data)))  # Replace any whitespace/newlines with underscore
  
  # Count genes before filtering
  genes_before <- nrow(data)
  
  # Step 1: Convert all columns to character and replace "N/A" with NA
  data <- data %>%
    mutate(across(everything(), as.character)) %>%         # Convert all columns to character
    mutate(across(everything(), ~ na_if(., "N/A")))        # Replace "N/A" with NA
  
  # Step 2: Convert specific columns back to numeric if needed
  data <- data %>%
    mutate(Aggregate_PCT = as.numeric(Aggregate_PCT),
           Cumulative_weighted_context_score = as.numeric(Cumulative_weighted_context_score),
           `3P-seq_tags_+_5` = as.numeric(`3P-seq_tags_+_5`))  # Add any other numeric columns here
  
  # Apply the filtering criteria (Moderate Confidence with Balanced Coverage)
  filtered_data <- data %>%
    filter((!is.na(Aggregate_PCT) & Aggregate_PCT >= 0.5) |  
             (is.na(Aggregate_PCT) & Cumulative_weighted_context_score <= -0.3))  
  
  # Count genes after filtering
  genes_after <- nrow(filtered_data)
  
  # Extract only the gene list column (assuming the column is named `Target_gene`)
  gene_list <- filtered_data %>% select(Target_gene) %>% pull()  # Convert to vector format
  
  # Save the gene list to a new text file
  gene_list_file_name <- paste0("gene_list_", tools::file_path_sans_ext(basename(file)), ".txt")
  writeLines(gene_list, con = gene_list_file_name)
  
  # Add a column to identify the miRNA (based on the file name)
  filtered_data$miRNA <- basename(file)
  
  # Store the filtered data in the list with the file name (or miRNA name) as the key
  filtered_data_list[[basename(file)]] <- filtered_data
  
  # Save the filtered data to a new Excel file
  # output_file_name <- paste0("filtered_", basename(file))
  # write.xlsx(filtered_data, file = output_file_name)
  
  # Add counts to the summary data frame
  summary_counts <- rbind(summary_counts, data.frame(
    miRNA = basename(file),
    genes_before = genes_before,
    genes_after = genes_after
  ))
}

########################

install.packages("readr")
install.packages("stringr")
install.packages("gridExtra")
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(readr)
library(gridExtra)
library(stringr)

# List all files in the folder that contain "flyenrichr" and have a .txt extension
file_list <- list.files(pattern = "flyenrichr-mf*\\.txt$", full.names = TRUE)

# Create an empty list to store individual plots
plot_list <- list()

# Generate labels (A, B, C, ...) for each plot
labels <- LETTERS[1:length(file_list)]

# Loop through each file to create a plot for each gene
for (i in seq_along(file_list)) {  # Ensure 'i' is properly defined here
  file <- file_list[i]
  
  # Read the data
  data <- read_tsv(file, show_col_types = FALSE)  # Suppress column spec messages
  
  # Clean up column names
  colnames(data) <- gsub("\\s+", "_", trimws(colnames(data)))
  
  # Convert Combined_Score to numeric if it isn't already
  data <- data %>%
    mutate(Combined_Score = as.numeric(Combined_Score))

  # Remove "-flyenrichr.txt" from the file name for the title
  title_text <- sub("-flyenrichr-mf\\.txt$", "", basename(file))
  
  # Filter out NA combined scores and select top 10 terms
  top_data <- data %>%
    filter(!is.na(Combined_Score)) %>%
    arrange(desc(Combined_Score)) %>%
    slice(1:5)  %>% # Select only the top 10 rows
    mutate(Term = str_wrap(Term, width = 30))  # Wrap text to 15 characters per line
  
  
  # Create a bar plot for this gene with an alphabet label in the top-left corner
  p <- ggplot(top_data, aes(x = reorder(Term, Combined_Score), y = Combined_Score, fill = Combined_Score)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # Flip for horizontal bars
    labs(
      title = paste(title_text),  # Use file name as title
      x = "Functional Enrichment Term",
      y = "Combined Score"
    ) +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +  # Color gradient
    theme_minimal(base_size = 10) +
    theme(legend.position = "none", axis.text.y = element_text(size = 8), plot.title = element_text(hjust = 0.5))   # Remove legend
    # annotate("text", x = 0.5, y = 0.5, label = labels[i], hjust = 1.5, vjust = 1.5, size = 6, fontface = "bold")  # Add label
  
  # Add the plot to the list
  plot_list[[basename(file)]] <- p
}

# Save the plot to a PNG file
png("enrichment_plots_mf.png", width = 3200, height = 4600, res = 300)  # Adjust resolution as needed
grid.arrange(grobs = plot_list, ncol = 2, nrow = 6)
dev.off()  # Close the PNG device

# Group plots into sets of 2 and save each group to a separate PNG
# plots_per_file <- 2  # Number of plots per file
# num_files <- ceiling(length(plot_list) / plots_per_file)  # Calculate the number of files needed


####################

# install.packages("cowplot")
library(cowplot)
# install.packages("magick")
library(magick)

# Load pathway and function plots
pathway_plot <- ggdraw() + draw_image("enrichment_plots_bp.png")
function_plot <- ggdraw() + draw_image("enrichment_plots_mf.png")

# Combine into a single plot with labels
combined_plot <- plot_grid(
  pathway_plot, function_plot,
  ncol = 1,  # Arrange vertically
  labels = c("A", "B"),  # Add subfigure labels
  label_size = 12
)

# Save combined plot
ggsave("combined_enrichment_plot.png", combined_plot, width = 10, height = 15, dpi = 300, bg="white")

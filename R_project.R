# Install and load necessary packages
install.packages(c("rentrez", "tm", "stringr"))
library(rentrez)
library(tm)
library(stringr)

# Load your gene list from the TSV file on Desktop
gene_file <- "<yourpathhere>/gene_symbols.tsv"
genes_df <- read.delim(gene_file, header = TRUE, stringsAsFactors = FALSE)

# Extract gene symbols from the data frame
genes <- genes_df$Symbol  # Adjust column name if needed

# Define search terms for conditions
condition1 <- "circadian rhythm"
condition2 <- "autism"

# Set up a data frame to store results
results <- data.frame(Gene = character(), PMID = character(), Title = character(), Abstract = character(), stringsAsFactors = FALSE)

# Define path for intermediate results
intermediate_file <- "/<yourpathhere>/PubMed_Gene_Circadian_Autism_Intermediate.csv"

# Load intermediate results if they already exist
if (file.exists(intermediate_file)) {
  results <- read.csv(intermediate_file, stringsAsFactors = FALSE)
  processed_genes <- unique(results$Gene)  # Get genes already processed
  genes <- setdiff(genes, processed_genes)  # Exclude processed genes
} else {
  # If no intermediate file exists, initialize it
  write.csv(results, intermediate_file, row.names = FALSE)
}

# Batch processing for genes to reduce API calls
batch_size <- 5
gene_batches <- split(genes, ceiling(seq_along(genes) / batch_size))

for (gene_batch in gene_batches) {

  # Construct search query for a batch of genes
  gene_query <- paste(gene_batch, collapse = " OR ")
  query <- paste(gene_query, condition1, condition2, sep = " AND ")

  # Search PubMed with combined query
  search_results <- entrez_search(db="pubmed", term=query, retmax=10)  # Limiting to 10 results per query

  # Pause to avoid server overload
  Sys.sleep(1)

  # Process results if any are found
  if (search_results$count > 0) {

    # Fetch abstracts for articles
    fetched_abstracts <- entrez_fetch(db="pubmed", id=search_results$ids, rettype="abstract", retmode="text")
    abstracts_split <- strsplit(fetched_abstracts, "\\n\\n")[[1]]

    for (i in seq_along(search_results$ids)) {
      pmid <- search_results$ids[i]
      title <- str_extract(abstracts_split[i], "^.*\\.")  # Extract title up to the first period
      abstract <- sub("^.*?\\.", "", abstracts_split[i])  # Remove title to get the abstract text

      for (gene in gene_batch) {
        # Append each gene and result to the results data frame
        results <- rbind(results, data.frame(Gene = gene, PMID = pmid, Title = title, Abstract = abstract, stringsAsFactors = FALSE))
      }
    }

    # Write updated results to the intermediate file to save progress
    write.csv(results, intermediate_file, row.names = FALSE, append = TRUE)
  }
}

# Print the resulting data frame with PubMed IDs and abstracts
print(results)

# OPTIONAL: Save final results to a separate CSV
write.csv(results, "/<yourpathhere>/PubMed_Gene_Circadian_Autism_Final_Results.csv", row.names = FALSE)

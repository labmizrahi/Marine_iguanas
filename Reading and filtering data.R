---
  # --- Load Data ---
  # Note: Ensure data files are located in a 'data/' subdirectory relative to this script
  otu_raw  <- read.csv('otu_table location')
  meta_raw <- read.csv("metadata location")[1:112, ]
  taxa_raw <- read.csv("taxa_table_location")
  
  # --- Clean Metadata ---
  # Create consistent sample IDs
  rownames(meta_raw) <- paste0("sa", meta_raw$Serial..)
  meta_df <- sample_data(meta_raw)
  
  # --- Clean OTU Table ---
  # Standardize column names to match metadata sample IDs
  colnames(otu_raw) <- gsub(pattern = "X", replacement = "sa", x = colnames(otu_raw))
  otu_tab <- otu_table(otu_raw, taxa_are_rows = TRUE)
  
  # --- Clean Taxonomy ---
  rownames(taxa_raw) <- taxa_raw$Feature.ID
  tax_tab <- tax_table(as.matrix(taxa_raw))
  
  # --- Construct Phyloseq Object ---
  ps_iguana <- phyloseq(otu_tab, meta_df, tax_tab)
  
  # Basic validation (optional logging)
  cat("Initial Dataset Summary:\n")
  print(ps_iguana)
  # Define filter parameters
  min_count <- 2
  min_sample_frac <- 0.05
  
  # Create filter function
  flt <- phyloseq::genefilter_sample(
    ps_iguana, 
    filterfun_sample(function(x) x >= min_count), 
    A = min_sample_frac * nsamples(ps_iguana)
  )
  
  # Prune taxa
  ps_iguana_filt <- prune_taxa(flt, ps_iguana)
  
  # Update the main object to the filtered version
  ps_iguana <- ps_iguana_filt
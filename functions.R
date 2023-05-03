### Spatial -----------------------------------------------------------------------

## Read in spatial data
read_spatial <- function(spatial_dir, project_name, coord_dir){
  gene_exp <- Read10X(spatial_dir)
  slide_seq <- CreateSeuratObject(gene_exp, project = project_name, assay = "Spatial")

  barcode_xy <- read.table(coord_file, row.names = 1)
  colnames(barcode_xy) <- c("x", "y")

  slide_seq[['image']] <- new(
      Class = 'SlideSeq',
      assay = "Spatial",
      coordinates = barcode_xy
  )
}

## Returns all beads around the vicinity of a spatial object ----------------------
# seurat_spatial: spatial seurat object
# vicinity: in microns
# center_coords: dataframe of center xy coordinates to find vicinity around- a subset of spatial@images$image@coordinates
find_surrounding_beads <- function(seurat_spatial, center_coords_df,  vicinity = 20){
  full_coords_df <- seurat_spatial@images$image@coordinates
  full_coords_df$barcode <- rownames(full_coords_df)
  center_coords_df$barcode <- rownames(center_coords_df)
  
  surr_beads <- c()
  
  center_coords_df <- center_coords_df %>%
    mutate(left = x - vicinity,
         right = x + vicinity,
         top = y + vicinity, 
         down = y - vicinity)
  
  for(i in c(1:nrow(center_coords_df))){
    left = center_coords_df[i,"left"]
    right = center_coords_df[i,"right"]
    top = center_coords_df[i,"top"]
    down = center_coords_df[i,"down"]
    for(j in c(1:nrow(full_coords_df))){
      if(full_coords_df[j,"x"] > left & full_coords_df[j,"x"] < right &
         full_coords_df[j,"y"] < top & full_coords_df[j,"y"] > down){
        surr_beads <- append(surr_beads, full_coords_df[j,"barcode"])}
    }
  }
  
  unique_surr_beads <- unique(surr_beads)
  return(unique_surr_beads)
}


### Single cell RNA seq ----------------------------------------------------------------

## Loads (raw) 10x data into a seurat object -------------------------
load_data <- function(data_dir, name){
  raw_data <- Read10X(data.dir = data_dir)
  
  # Create Seurat object with gene and hashtags
  seurat <- CreateSeuratObject(counts = raw_data[["Gene Expression"]], project = name)
  seurat[["HTO"]] <- CreateAssayObject(counts = raw_data[["Antibody Capture"]])
  
  # Calculate mitochondrial and blood contamination percentage
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")
  seurat[["percent.hbb"]] <- PercentageFeatureSet(seurat, pattern = "^Hbb-")
  
  # Remove cells with 0 barcoding and very lax filtering (in order to cut down on the object size)
  seurat <- subset(seurat, subset = nFeature_HTO > 0 & nFeature_RNA > 50 & nCount_RNA < 50000)
  
  return(seurat)
}

## Demultiplexes seurat object -----------------------------------------------------
demultiplex <- function(seurat){
  # Normalize RNA count and HTO data, use centered log-ratio (CLR) transformation for HTO
  seurat <- NormalizeData(seurat, assay = "HTO", normalization.method = "CLR")
  
  # Demultiplex cells based on HTO enrichment
  seurat <- HTODemux(seurat, assay = "HTO", positive.quantile = 0.99)
  
  # Rename levels
  seurat$hash.ID <- factor(seurat$hash.ID, levels = c('HTO-A', 'HTO-B', 'HTO-C', 'HTO-D', 'HTO-E', 'HTO-F', 'Negative', 'Doublet'))
  
  table1 <- table(seurat$HTO_classification.global)
  # Show the signals for each hashtag in each 'sample'
  plot2 <- RidgePlot(seurat, assay = "HTO", features = rownames(seurat[["HTO"]]), ncol = 3)
  # Print global classification result
  plot3 <- table(seurat$hash.ID)
  
  print(table1)
  print(plot2)
  print(plot3)
  return(seurat)
 }


## Density plot showing the no. of cells split by a metadata col -------------------
plot_cell_count_by_metadata <- function(seurat, metadata_col){
  # Extract unique values from the metadata column
  md_factors <- seurat@meta.data[[metadata_col]]
  md_factors <- unique(md_factors)
  
  # Create dataframe with unique values as row names
  countDf <- data.frame(matrix(ncol = 2, nrow = 0))
  # Count number of cells associated with each metadata value
  for (md in md_factors){
    count <- c(md, dim(subset(seurat@meta.data, seurat@meta.data[[metadata_col]] == md))[1])
    countDf <- rbind(countDf, count)
  }
  
  # Rename counts dataframe
  colnames(countDf) <- c("metadata", "count")
  countDf[,2] <- as.numeric(countDf[,2])
  
  # Plot barplot showing the number of cells associated with each metadata column, along with the mean number of cells in the entire object
  plot <- ggplot(countDf, aes(x = metadata, y = count)) +
    theme_minimal() +
    #coord_cartesian(ylim=c(10000, 28000)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = count)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    geom_hline(yintercept = mean(countDf$count), color = "red")
  
  return(plot)
}

## Density plot showing the distribution of numerical columns found in metadata ----
plot_distribution_by_metadata <- function(seurat,  metadata_col, lowerSD = 1, upperSD = 3, sk6Metrics = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.hbb")){
  # Create a dataframe holding the mean and standard deviation of the given distribution
  sdTable <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(sdTable) <- c("mean", "SD", paste0("lowerSD", lowerSD), paste0("upperSD", upperSD))
  
  for(metric in sk6Metrics){
    meanVal <- mean(as.numeric(seurat@meta.data[[metric]]))
    sdVal <- sd(as.numeric(seurat@meta.data[[metric]]))
    sdTable[metric,] <- c(meanVal, sdVal, meanVal - lowerSD*sdVal, meanVal + upperSD*sdVal)
  }
  
  # Remove any metrics with NA values- they can't be plotted with an intercept line.
  sdTable <- na.omit(sdTable)
  
  # Get remaining metrics
  remaining_metrics <- rownames(sdTable)
  
  # Create list of plots to return
  plot_list <- list()
  
  plot_list <- append(plot_list, list(sdTable))
  
  # Create plots for remaining metrics and store in list
  for(metric in remaining_metrics){
    plot <- ggplot(seurat@meta.data, aes(x=seurat@meta.data[,metric], fill = seurat@meta.data[,metadata_col])) + 
  	  geom_density(alpha = 0.2, show.legend = FALSE) + 
  	  scale_x_log10() + 
  	  theme_classic() +
  	  ylab("Density") +
      xlab(paste0(metric, " (log scaled)")) +
      facet_wrap(as.formula(paste0('.~', metadata_col))) +
      geom_vline(xintercept = as.numeric(sdTable[metric,c(3,4)]), color = "blue")
    print(plot)
  }
}

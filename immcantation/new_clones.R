suppressPackageStartupMessages(library(shazam))

# Get the sample ID from command line argument
sampleID <- commandArgs(trailingOnly = TRUE)[1]

# Define the loci to process
loci <- c("H", "IGK", "IGL")

# Iterate over each locus for the given sample
for (locus in loci) {
  # Construct the input file path
  inputFileName <- paste0("./presto/changeo/", sampleID, "/igblast_output_", locus, ".fmt7.tsv")
  
  # Check if the input file exists
  if (!file.exists(inputFileName)) {
    next  # Skip this locus if file does not exist
  }
  
  # Read the input file
  changeoTable <- read.table(inputFileName, sep="\t", header=TRUE)
  
  # Initialize threshold
  threshold <- 0.16
  
  # Try to calculate the threshold using Hamming distance
  tryCatch({
    dist_ham <- distToNearest(changeoTable, sequenceColumn="junction", vCallColumn="v_call", jCallColumn="j_call", model="ham", normalize="len", nproc=11)
    output <- findThreshold(dist_ham$dist_nearest, method="density")
    calculated_threshold <- output@threshold
    
    if (!is.na(calculated_threshold)) {
      threshold <- calculated_threshold
    }
  }, error=function(e) {
    # Use default threshold if Hamming distance calculation fails
    # Already set to 0.16
  })
  
  print(paste("Threshold for sample", sampleID, "locus", locus, ":", threshold))
  
  # Run DefineClones.py with the calculated/default threshold
  system(paste0("DefineClones.py -d ", inputFileName, " --mode allele --act set --nproc 11 --model ham --norm len --dist ", threshold), intern = FALSE, ignore.stdout = FALSE)
  print(paste("Sample", sampleID, "locus", locus, "finished"))
}
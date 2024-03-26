suppressPackageStartupMessages(library(shazam))

# Read sample names file
samplesFile <- commandArgs(trailingOnly = TRUE)[1]
samples <- read.table(samplesFile, header=FALSE)$V1
samples <- trimws(samples)
# Iterate over each sample
for(sampleID in samples) {
  inputFileName <- paste0("./presto/", sampleID, "/changeo/", sampleID, "_merged-changeo-rename.tsv")

  # Check if the Change-O file exists
  if (!file.exists(inputFileName)) {
    next  # Skip this sampleID if file does not exist
  }

  # Read the input Change-O file
  changeoTable <- read.table(inputFileName, sep="\t", header=TRUE)

  # Initialize threshold
  threshold <- 0.16

  # Try to calculate the threshold using Hamming distance
  tryCatch({
    dist_ham <- distToNearest(changeoTable, sequenceColumn="junction", vCallColumn="v_call", jCallColumn="j_call", model="ham", normalize="len", nproc=5)
    output <- findThreshold(dist_ham$dist_nearest, method="density")
    calculated_threshold <- output@threshold
    if (!is.na(calculated_threshold)) {
      threshold <- calculated_threshold
    }
  }, error=function(e) {
    # Use default threshold if Hamming distance calculation fails
    # Already set to 0.16
  })

  print(paste("Threshold for sample", sampleID, ":", threshold))

  # Run DefineClones.py with the calculated/default threshold
  system(paste0("DefineClones.py -d ", inputFileName, " --mode allele --act set --nproc 11 --model ham --norm len --dist ", threshold), intern = FALSE, ignore.stdout = FALSE)
}

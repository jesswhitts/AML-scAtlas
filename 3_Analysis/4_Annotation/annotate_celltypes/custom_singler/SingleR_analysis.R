## Custom SingleR Annotations

# Load Packages
library(scater)
library(SingleR)
library(SingleCellExperiment)
library(dplyr)
library(parallel)
cores <- as.numeric(detectCores())

dir.create('outs')
dir.create('plots')

# Read in reference dataset
vgalen <- readRDS('/data/stemcell/jwhittle/ref/van-galen-celltypes/vgalen_sce.rds')
vgalen <- logNormCounts(vgalen)

# Seeing the available labels in this dataset.
table(vgalen$CellType)

# Read Test Dataset
print("Reading SCE...")
sce <- readRDS('../outs/sce_hvg.rds')


print("Running SingleR...")

# Run on Individual Cells - n of 50
pred.grun <- SingleR(test=sce, ref=vgalen, labels=vgalen$CellType, 
                      de.method="wilcox", de.n=50)

# Save Results
save(pred.grun, file = "outs/pred.grun.RData")

singler_scores <- pred.grun$scores %>% as_tibble() %>% dplyr::mutate(assigned_score = NA)

for ( i in seq_len(nrow(singler_scores)) ) {
  pred.grun$assigned_score[i] <- singler_scores[[pred.grun$labels[i]]][i]}

write.csv(pred.grun, paste0("outs/SingleR_custom_predictions_full.csv"))

png('plots/score_heatmap.png')
plotScoreHeatmap(pred.grun)
dev.off()

png("plots/delta_distribution.png")
plotDeltaDistribution(pred.grun)


print(sessionInfo())

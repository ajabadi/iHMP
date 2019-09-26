library(MultiAssayExperiment)
library(data.table)
## input
metab_files <- list.files(path = 'data/merged_Metabolomics/', full.names = TRUE)
wgs_files <- list.files(path = 'data/Merged_WGS/', full.names = TRUE)
metadata_file <- 'data/metadata.tsv'
prot_files <- list.files('data/Proteomics/', full.names = TRUE)

metadata <- fread(metadata_file)
## ----------- metabolomics ----------- 
## read data.table as data.frame with compound rownames
metab_list <- lapply(metab_files, function(x) as.matrix(data.frame(fread(x), row.names = 1)))
## name the list with file name
names(metab_list) <- gsub(".tsv", "", basename(metab_files))
View(metab_list[[1]][1:10,1:10])

## are all features unique?
sum(rownames(metab_list[[1]]) %in% rownames(metab_list[[2]]))
sum(rownames(metab_list[[1]]) %in% rownames(metab_list[[3]]))
## yes

## samples?
length(colnames(metab_list[[1]]))
sum(colnames(metab_list[[1]]) %in% colnames(metab_list[[2]]))
sum(colnames(metab_list[[1]]) %in% colnames(metab_list[[3]]))
## all the same

## so it seems they are same samples (phew!)
## let's MAE them

metab_mae <- MultiAssayExperiment(experiments = ExperimentList(metab_list))

## are all features unique?
sum(rownames(metab_list[[1]]) %in% rownames(metab_list[[2]]))
sum(rownames(metab_list[[1]]) %in% rownames(metab_list[[3]]))
## yes

## samples?
length(colnames(metab_list[[1]]))
sum(colnames(metab_list[[1]]) %in% colnames(metab_list[[2]]))
sum(colnames(metab_list[[1]]) %in% colnames(metab_list[[3]]))
## all the same

## so it seems they are same samples (phew!)
## let's merge them
metab_all <- Reduce(f = rbind, x = metab_list)
saveRDS(metab_all, file = 'output/merged_metab_all-samples.rds')


# metab_all <- readRDS('output/merged_metab_all-samples.rds')

## ----------- WGS ----------- 
## read data.table as data.frame with compound rownames
wgs_list <- lapply(wgs_files, function(x) data.frame(fread(x), row.names = 1))
View(wgs_list[[2]][1:10,1:10])

## are all features unique?
sum(rownames(wgs_list[[1]]) %in% rownames(wgs_list[[2]]))
sum(rownames(wgs_list[[1]]) %in% rownames(wgs_list[[3]]))
## samples?
length(colnames(wgs_list[[1]]))
sum(colnames(wgs_list[[1]]) %in% colnames(wgs_list[[2]]))
sum(colnames(wgs_list[[1]]) %in% colnames(wgs_list[[3]]))

## let's merge them
wgs_all <- Reduce(f = rbind, x = wgs_list)
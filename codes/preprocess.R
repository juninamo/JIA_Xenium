source("scripts/JIA_Xenium/utils.R")

### n_jobs <- 4
## n_perms <- 50
n_jobs = commandArgs(trailingOnly=TRUE)[1] %>% as.integer()
print(paste0("N cores is ", n_jobs))
n_perms = commandArgs(trailingOnly=TRUE)[2] %>% as.integer()
print(paste0("N permutation is ", n_perms))

celltype="ALL"
resolution <- switch(celltype,
                     "ALL" = "0.20",
                     "T" = "1.00",
                     "M" = "0.20",
                     "B" = "0.40",
                     "Tissue" = "0.60"
)

file_list = list.files("/pl/active/yomogidalab/JIA_Xenium/Xenium_Rawdata/", pattern = "output.*", full.names = TRUE)
sb_values <- stringr::str_extract(file_list, "SB[0-9]+")
print(sb_values)

for(sample_id in 1:length(file_list)){
  
  data_path = file_list[sample_id]
  data <- ReadXenium(data_path, outs = c("matrix", "microns"), type = c("centroids", "segmentations"))
  
  segmentations.data <- list(
    centroids = CreateCentroids(data$centroids),
    segmentation = CreateSegmentation(data$segmentations))
  coords <- CreateFOV(
    coords = segmentations.data, 
    type = c("segmentation", "centroids"), 
    molecules = data$microns, 
    assay = "Spatial")
  xenium.obj <- CreateSeuratObject(
    counts = data$matrix[["Gene Expression"]], 
    assay = "Spatial")
  xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
  xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
  xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
  xenium.obj[["fov"]] <- coords
  rm(data); gc(); gc()
  rm(coords); gc(); gc()
  rm(segmentations.data); gc(); gc()
  
  xenium.obj <- subset(xenium.obj, subset = nCount_Spatial > 0)
  
  print(dim(xenium.obj@assays$Spatial$counts))
  xenium.obj@assays$Spatial$counts[1:5,1:5]
  head(xenium.obj@meta.data)
  summary(xenium.obj@meta.data$nCount_Spatial)
  summary(xenium.obj@meta.data$nFeature_Spatial)
  
  # Harmony with Seurat SCTransform: https://github.com/immunogenomics/harmony/issues/41
  xenium.obj <- subset(xenium.obj, subset = nFeature_Spatial > 200 & nFeature_Spatial < 5000 & nCount_Spatial < 25000)
  xenium.obj <- SCTransform(xenium.obj, assay = "Spatial",
                            conserve.memory = TRUE, vst.flavor="v2")
  
  assign(paste0("xenium.obj_", sample_id), xenium.obj)
  
}
rm(xenium.obj); gc(); gc()

ste <- merge(xenium.obj_1, y = list(xenium.obj_2, xenium.obj_3, xenium.obj_4, xenium.obj_5, xenium.obj_6, xenium.obj_7, xenium.obj_8, xenium.obj_9), 
             add.cell.ids = sb_values)
names(ste@images) = sb_values
ste@meta.data$slide = dplyr::case_when(
  stringr::str_split(rownames(ste@meta.data),pattern="_",simplify=TRUE)[,1] == "SB4"  ~ "slide1",
  stringr::str_split(rownames(ste@meta.data),pattern="_",simplify=TRUE)[,1] == "SB14" ~ "slide1",
  stringr::str_split(rownames(ste@meta.data),pattern="_",simplify=TRUE)[,1] == "SB10" ~ "slide2",
  stringr::str_split(rownames(ste@meta.data),pattern="_",simplify=TRUE)[,1] == "SB19" ~ "slide2",
  stringr::str_split(rownames(ste@meta.data),pattern="_",simplify=TRUE)[,1] == "SB25" ~ "slide3",
  stringr::str_split(rownames(ste@meta.data),pattern="_",simplify=TRUE)[,1] == "SB21" ~ "slide3",
  stringr::str_split(rownames(ste@meta.data),pattern="_",simplify=TRUE)[,1] == "SB29" ~ "slide4",
  stringr::str_split(rownames(ste@meta.data),pattern="_",simplify=TRUE)[,1] == "SB30" ~ "slide4",
  stringr::str_split(rownames(ste@meta.data),pattern="_",simplify=TRUE)[,1] == "SB62" ~ "slide4"
)
table(ste@meta.data$slide)
ste@meta.data$sample_id = stringr::str_split(rownames(ste@meta.data),pattern="_",simplify=TRUE)[,1]
table(ste@meta.data$sample_id)
head(ste@meta.data)

var_genes = SelectIntegrationFeatures(list(xenium.obj_1,xenium.obj_2,xenium.obj_3,xenium.obj_4,xenium.obj_5,xenium.obj_6,xenium.obj_7,xenium.obj_8,xenium.obj_9), nfeatures = 2000, verbose = TRUE)
rm(xenium.obj_1); gc(); gc()
rm(xenium.obj_2); gc(); gc()
rm(xenium.obj_3); gc(); gc()
rm(xenium.obj_4); gc(); gc()
rm(xenium.obj_5); gc(); gc()
rm(xenium.obj_6); gc(); gc()
rm(xenium.obj_7); gc(); gc()
rm(xenium.obj_8); gc(); gc()
rm(xenium.obj_9); gc(); gc()


ste[["Spatial"]] <- JoinLayers(ste[["Spatial"]])

ste <- SCTransform(ste, assay = "Spatial",
                   conserve.memory = TRUE, vst.flavor="v2")
VariableFeatures(ste[["SCT"]]) <- var_genes
set.seed(1234)
ste <- RunPCA(ste, assay = "SCT", verbose = FALSE)

for (name in names(ste@images)) {
  ste@images[[name]]@key <- paste0(name, "_")  # Update the Key to match the name (SB4_, SB10_, etc.)
}

meta = read.csv("/path/to/jinamo/JIA_Xenium/meta.csv")
meta = dplyr::left_join(ste@meta.data, meta, by = "sample_id")
ste <- AddMetaData(ste, meta)

print("Harmony...")
set.seed(1234)
ste <- ste %>% 
  RunHarmony(c("slide","sample_id"), 
             assay.use = "SCT",
             plot_convergence = TRUE)


ste <- FindNeighbors(ste, 
                     reduction = "harmony", 
                     k.param = 30,
                     dims = 1:30)

ste <- RunUMAP(ste, 
               reduction = "harmony", 
               n.neighbors = 30L,
               min.dist = 0.3,
               dims = 1:30)

print("Clustering...")
# Louvain
snn_pcs <- BuildSNNSeurat(ste[["harmony"]]@cell.embeddings[,1:30], 
                          nn.eps = 0)

resolution_list <- c(0.2, 0.4, 0.6, 0.8, 1.0)
ids_cos <- Reduce(cbind, parallel::mclapply(resolution_list, function(res_use) {
  Seurat:::RunModularityClustering(SNN = snn_pcs, 
                                   modularity = 1, 
                                   resolution = res_use, 
                                   algorithm = 3, 
                                   n.start = 10, 
                                   n.iter = 10, random.seed = 0, print.output = FALSE, 
                                   temp.file.location = NULL, edge.file.name = NULL)    
}, mc.cores = min(16, length(resolution_list))))
ids_cos %<>% data.frame()
colnames(ids_cos) <- sprintf("res_%.2f", resolution_list)

ids_cos$res_0.20 <- as.character(ids_cos$res_0.20)
ids_cos$res_0.40 <- as.character(ids_cos$res_0.40)
ids_cos$res_0.60 <- as.character(ids_cos$res_0.60)
ids_cos$res_0.80 <- as.character(ids_cos$res_0.80)
ids_cos$res_1.00 <- as.character(ids_cos$res_1.00)

rownames(ids_cos) = rownames(ste@meta.data)
ids_cos <- ids_cos %>%
  dplyr::mutate(across(everything(), ~ factor(.x, levels = sort(unique(.x)))))
ste <- AddMetaData(ste, ids_cos)
head(ste@meta.data)

ste@meta.data$seurat_clusters = ste@meta.data[,paste0("res_",resolution)]
Idents(ste) = ste@meta.data[,paste0("res_",resolution)]

ste <- CellCycleScoring(ste, 
                        s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, 
                        set.ident = TRUE)

Idents(ste) = ste@meta.data[,paste0("res_",resolution)]

## harmonoized data [hPC1-20]
lisi_res_harmony <- rbind(
  lisi::compute_lisi(ste[["harmony"]]@cell.embeddings, 
                     ste@meta.data,
                     c("slide", "seurat_clusters", "sample_id"))
) %>% 
  tidyr::gather(key, val, slide, seurat_clusters, sample_id)
print(head(lisi_res_harmony))
saveRDS(lisi_res_harmony, file=paste0("/path/to/jinamo/JIA_Xenium/output/LISI_harmony.rds"))

## pre-harmonoized data [PC1-20]
lisi_res_preharmony <- rbind(
  lisi::compute_lisi(ste[["pca"]]@cell.embeddings, 
                     ste@meta.data,
                     c("slide", "seurat_clusters", "sample_id"))
) %>% 
  tidyr::gather(key, val, slide, seurat_clusters, sample_id)
print(head(lisi_res_preharmony))
saveRDS(lisi_res_preharmony, file=paste0("/path/to/jinamo/JIA_Xenium/output/LISI_PCA.rds"))


if(length(colnames(ste@meta.data)[grepl("score",colnames(ste@meta.data))] > 0)){
  for(score_name in colnames(ste@meta.data)[grepl("score",colnames(ste@meta.data))]){
    ste[[score_name]] <- NULL
  }
}

feature_lt <- list(
  c("CD4","PDCD1","ICOS","CXCL13"),# Tph
  c("CD4","FOXP3","CTLA4","CCR4"), # Treg
  c("CD4","CXCR3","CYB561"), # Th1
  c("CD4","CCR4","GALT","GATA3"), # Th2
  c("CD4","RORC","CCR6"), # Th17
  c("CD4","CXCR5","FGFR2","SGPP2"), # Tfh
  c("CD8B","GZMK"), 
  c("CD8B","GZMB"),
  
  c("CR2", "FCER2", "NR4A1"), # CD21, CD23: activaion
  c("TCL1A","IGHD"), # Naive
  c("CD27", "TNFRSF13B","S100A10","S100A4"), # memory
  c("IGHG1","IGHA1","SDC1","CXCR4"), # CD138: plasma cells
  c("RPN2","XBP1","PRDX4", "MZB1", "CD38"), # plasmablasts
  c("ITGAX", "TBX21", "ZEB2", "EMP3","FGR", "FCRL5"), # "CD11c", "TBET": ABC
  
  c("PDPN", "THY1", "ICAM1","VCAM1","CD82","CXCL9","CCL19","CD34","TNFSF13B","PDPN", "THY1", "ICAM1","VCAM1","PRG4","FN1"),
  
  c("CD3E","CD4","CD8A"),
  c("CD19", "MS4A1", "IGHG1","IGHA1","SDC1", "XBP1","MZB1"),
  c("NCAM1","FCER1G","KLRF1", "NKG7","GNLY","PRF1"),
  c("TYROBP","CD14","CD163","CLEC1B", "ACTB","FTL", "FCGR3B")
)

score_names = c("Tph.score",
                "Treg.score",
                "Th1.score",
                "Th2.score",
                "Th17.score",
                "Tfh.score",
                "GZMK.score",
                "GZMB.score",
                
                "Activation.score",
                "Naive.score",
                "Memory.score",
                "Plasma.score",
                "Plasmablast.score",
                "ABC.score",
                
                "Fibroblast.score",
                
                "T.score",
                "B.plasma.score",
                "NK.score",
                "Myeloid.score"
)

ste <- AddModuleScore(
  object = ste,
  features = feature_lt,
  name = score_names
)

colnames(ste@meta.data) = gsub("\\.score\\d+", ".score", colnames(ste@meta.data))

sce = SingleCellExperiment(assays = list(counts=ste@assays$Spatial$counts), colData=ste@meta.data)
sce <- scDblFinder(sce, samples="sample_id")

library(scuttle)
sce <- logNormCounts(sce)

library(scran)
dec <- modelGeneVar(sce)
hvgs <- getTopHVGs(dec, n=2000)

set.seed(1003)
scores <- computeDoubletDensity(sce, subset.row=hvgs)
classes = sce$scDblFinder.class
names(scores) = names(classes) = colnames(ste@assays$Spatial$counts)
saveRDS(scores, paste0("/path/to/jinamo/JIA_Xenium/output/scDblFinder_Scores.rds"))
saveRDS(classes, paste0("/path/to/jinamo/JIA_Xenium/output/scDblFinder_Class.rds"))

rm(list = c("sce"))
gc(reset = TRUE)
gc(reset = TRUE)

scores = readRDS(paste0("/path/to/jinamo/JIA_Xenium/output/scDblFinder_Scores.rds"))
classes = readRDS(paste0("/path/to/jinamo/JIA_Xenium/output/scDblFinder_Class.rds"))

ste <- AddMetaData(ste, data.frame(cell = names(classes),
                                   doublet_score = scores,
                                   doublet_class = classes) %>% 
                     magrittr::set_colnames(c("cell","doublet_score", "doublet_class")))

table(ste@meta.data$doublet_class)
table(ste@meta.data$doublet_class, ste@meta.data$sample_id)

ste@meta.data$log_doublet_score = log(ste@meta.data$doublet_score + 1)

# Auto-correlation by genes
ste = spatial_autocorr(
  ste, 
  neighbors.k = 30, 
  connectivity_key = "nn",
  genes = NULL,
  mode = "moran",
  transformation = TRUE, 
  #n_perms = n_perms,
  n_perms = 50,
  two_tailed = FALSE,
  corr_method = "BH",
  assay = "SCT",
  attr = "data",
  seed = 1938493,
  copy = FALSE
)

saveRDS(ste, file = paste0("/path/to/jinamo/JIA_Xenium/output/Xenium_merged.rds"))

ste@meta.data$celltype = dplyr::case_when(
  ste@meta.data[,paste0("res_",resolution)] %in% c("1","2","4","6","9","10","12","13") ~ "Tissue", #13: Muskcle (MYH2,PYGM)
  ste@meta.data[,paste0("res_",resolution)] %in% c("0","7","11","14","15","16") ~ "M",
  ste@meta.data[,paste0("res_",resolution)] %in% c("3") ~ "T",
  ste@meta.data[,paste0("res_",resolution)] %in% c("5","8") ~ "B"
)
table(ste@meta.data$celltype)

saveRDS(rownames(ste@meta.data)[ste@meta.data$celltype=="Tissue"],
        "/path/to/jinamo/JIA_Xenium/output/Tissue_cellids.rds")
saveRDS(rownames(ste@meta.data)[ste@meta.data$celltype=="M"],
        "/path/to/jinamo/JIA_Xenium/output/M_cellids.rds")
saveRDS(rownames(ste@meta.data)[ste@meta.data$celltype=="T"],
        "/path/to/jinamo/JIA_Xenium/output/T_cellids.rds")
saveRDS(rownames(ste@meta.data)[ste@meta.data$celltype=="B"],
        "/path/to/jinamo/JIA_Xenium/output/B_cellids.rds")


celltype="T"

ste = readRDS(paste0("/path/to/jinamo/JIA_Xenium/output/Xenium_merged.rds"))

resolution <- switch(celltype,
                     "ALL" = "0.20",
                     "T" = "1.00",
                     "M" = "0.20",
                     "B" = "0.40",
                     "Tissue" = "0.60"
)

cell_ids = readRDS(paste0("/path/to/jinamo/JIA_Xenium/output/",celltype,"_cellids.rds"))
ste = subset(ste, cells=cell_ids)
table(ste@meta.data$seurat_clusters)

ste <- SCTransform(ste, assay = "Spatial",
                   conserve.memory = TRUE, vst.flavor="v2")

VariableFeatures(ste[["SCT"]]) <- rownames(ste[["SCT"]]$data)
set.seed(1234)
ste <- RunPCA(ste, assay = "SCT", verbose = FALSE)

print("Harmony...")
set.seed(1234)
ste <- ste %>% 
  RunHarmony(c("slide","sample_id"), 
             assay.use = "SCT",
             plot_convergence = TRUE)


ste <- FindNeighbors(ste, 
                     reduction = "harmony", 
                     k.param = 30,
                     dims = 1:30)

ste <- RunUMAP(ste, 
               reduction = "harmony", 
               n.neighbors = 30L,
               min.dist = 0.3,
               dims = 1:30)

print("Clustering...")
# Louvain
snn_pcs <- BuildSNNSeurat(ste[["harmony"]]@cell.embeddings[,1:30], 
                          nn.eps = 0)

resolution_list <- c(0.2, 0.4, 0.6, 0.8, 1.0)
ids_cos <- Reduce(cbind, parallel::mclapply(resolution_list, function(res_use) {
  Seurat:::RunModularityClustering(SNN = snn_pcs, 
                                   modularity = 1, 
                                   resolution = res_use, 
                                   algorithm = 3, 
                                   n.start = 10, 
                                   n.iter = 10, random.seed = 0, print.output = FALSE, 
                                   temp.file.location = NULL, edge.file.name = NULL)    
}, mc.cores = min(16, length(resolution_list))))
ids_cos %<>% data.frame()
colnames(ids_cos) <- sprintf("res_%.2f", resolution_list)

ids_cos$res_0.20 <- as.character(ids_cos$res_0.20)
ids_cos$res_0.40 <- as.character(ids_cos$res_0.40)
ids_cos$res_0.60 <- as.character(ids_cos$res_0.60)
ids_cos$res_0.80 <- as.character(ids_cos$res_0.80)
ids_cos$res_1.00 <- as.character(ids_cos$res_1.00)

rownames(ids_cos) = rownames(ste@meta.data)
ids_cos <- ids_cos %>%
  dplyr::mutate(across(everything(), ~ factor(.x, levels = sort(unique(.x)))))
ste <- AddMetaData(ste, ids_cos %>%
                     magrittr::set_colnames(paste0(colnames(.),"_",celltype)))
head(ste@meta.data)

ste@meta.data$seurat_clusters = ste@meta.data[,paste0("res_",resolution,"_",celltype)]
Idents(ste) = factor(ste@meta.data$seurat_clusters, levels = as.character(sort(as.integer(unique(ste@meta.data$seurat_clusters)))-1))
table(ste@meta.data$seurat_clusters)

## harmonoized data [hPC1-20]
lisi_res_harmony <- rbind(
  lisi::compute_lisi(ste[["harmony"]]@cell.embeddings, 
                     ste@meta.data,
                     c("slide", "seurat_clusters", "sample_id"))
) %>% 
  tidyr::gather(key, val, slide, seurat_clusters, sample_id)
print(head(lisi_res_harmony))
saveRDS(lisi_res_harmony, file=paste0("/path/to/jinamo/JIA_Xenium/output/LISI_harmony_",celltype,".rds"))

## pre-harmonoized data [PC1-20]
lisi_res_preharmony <- rbind(
  lisi::compute_lisi(ste[["pca"]]@cell.embeddings, 
                     ste@meta.data,
                     c("slide", "seurat_clusters", "sample_id"))
) %>% 
  tidyr::gather(key, val, slide, seurat_clusters, sample_id)
print(head(lisi_res_preharmony))
saveRDS(lisi_res_preharmony, file=paste0("/path/to/jinamo/JIA_Xenium/output/LISI_PCA_",celltype,".rds"))

table(ste@meta.data$doublet_class, ste@meta.data$seurat_clusters)
table(ste@meta.data$doublet_class, ste@meta.data$sample_id)


saveRDS(ste, file = paste0("/path/to/jinamo/JIA_Xenium/output/Xenium_",celltype,".rds"))


# AMP RA-synovium paper (Nature 2023)
syno.rna = readRDS("/path/to/jinamo/scRNA_ATACseq/data/RA_syno/syno.rna_Fan.rds")
ra_celltype <- switch(celltype,
                      "T" = c("T cell","NK"),
                      "M" = "Myeloid",
                      "B" = "B cell",
                      "Tissue" = c("Fibroblast","Endothelial"))

head(syno.rna@meta.data)
table(syno.rna@meta.data$ct_subtype)
table(syno.rna@meta.data$ct)
syno.rna = subset(syno.rna, ct %in% ra_celltype)
table(syno.rna@meta.data$ct)

ra = syno.rna@assays$RNA$data
jia = ste@assays$SCT$data

shared_mrna = intersect(rownames(jia),
                        rownames(ra))
length(shared_mrna)
print(shared_mrna)

ra = ra[shared_mrna,]
jia = jia[shared_mrna,]

colnames(ra) = paste0("RA:",colnames(ra))
colnames(jia) = paste0("JIA:",colnames(jia))

print("mRNAs in RA scRNA-seq;")
dim(ra)

print("mRNAs in JIA Xenium;")
dim(jia)

assay_list_indirect = list(
  RA = ra,
  JIA = jia
)

lapply(assay_list_indirect, dim)
lapply(assay_list_indirect, class)

stab_indirect = stabMap(assay_list_indirect,
                        reference_list = c("RA"),
                        maxFeatures = Inf,
                        ncomponentsReference = 30,
                        ncomponentsSubset = 30,
                        plot = FALSE)

set.seed(1234)
stab_indirect_umap <- uwot::umap(as.data.frame(stab_indirect), 
                                 n_neighbors = 30, 
                                 metric = "cosine", 
                                 min_dist = .1, 
                                 n_threads = 4,
                                 verbose = TRUE)
stab_indirect_umap = stab_indirect_umap %>% 
  magrittr::set_colnames(c("UMAP1","UMAP2")) %>%
  as.data.frame()


train_data <- stab_indirect[grepl("^RA:",rownames(stab_indirect)),]  
train_labels <- syno.rna@meta.data$ct_subtype
test_data <- stab_indirect[!grepl("^RA:",rownames(stab_indirect)),]


knn_pred <- class::knn(train = train_data, test = test_data, cl = train_labels, k = 30, prob = TRUE)
knn_prob = attributes(knn_pred)$prob
test_data = data.frame(test_data,
                       ct_subtype = knn_pred,
                       knn_prob = knn_prob)
train_data = data.frame(train_data,
                        ct_subtype = train_labels,
                        knn_prob = 1)
stab_indirect = rbind(train_data,test_data)

stab_indirect = cbind(stab_indirect,stab_indirect_umap)
head(stab_indirect)
tail(stab_indirect)

# Louvain
snn_pcs <- BuildSNNSeurat(stab_indirect[,1:30], nn.eps = 0)

resolution_list <- c(0.2, 0.4, 0.6, 0.8, 1.0)
ids_cos <- base::Reduce(cbind, parallel::mclapply(resolution_list, function(res_use) {
  Seurat:::RunModularityClustering(SNN = snn_pcs, modularity = 1, 
                                   resolution = res_use, algorithm = 3, n.start = 10, 
                                   n.iter = 10, random.seed = 0, print.output = FALSE, 
                                   temp.file.location = NULL, edge.file.name = NULL)    
}, mc.cores = min(16, length(resolution_list))))
ids_cos %<>% data.frame()
colnames(ids_cos) <- sprintf("res_%.2f", resolution_list)

stab_indirect$res_0.20_stab <- as.character(ids_cos$res_0.20)
stab_indirect$res_0.40_stab <- as.character(ids_cos$res_0.40)
stab_indirect$res_0.60_stab <- as.character(ids_cos$res_0.60)
stab_indirect$res_0.80_stab <- as.character(ids_cos$res_0.80)
stab_indirect$res_1.00_stab <- as.character(ids_cos$res_1.00)

saveRDS(stab_indirect, paste0("/path/to/jinamo/JIA_Xenium/output/StabMap_Obj_RAtoJIA_",celltype,".rds"))



celltype="B"

ste = readRDS(paste0("/path/to/jinamo/JIA_Xenium/output/Xenium_merged.rds"))

resolution <- switch(celltype,
                     "ALL" = "0.20",
                     "T" = "1.00",
                     "M" = "0.20",
                     "B" = "0.40",
                     "Tissue" = "0.60"
)

cell_ids = readRDS(paste0("/path/to/jinamo/JIA_Xenium/output/",celltype,"_cellids.rds"))
ste = subset(ste, cells=cell_ids)
ste
table(ste@meta.data$seurat_clusters)

ste <- SCTransform(ste, assay = "Spatial",
                   conserve.memory = TRUE, vst.flavor="v2")

VariableFeatures(ste[["SCT"]]) <- rownames(ste[["SCT"]]$data)
set.seed(1234)
ste <- RunPCA(ste, assay = "SCT", verbose = FALSE)

print("Harmony...")
set.seed(1234)
ste <- ste %>% 
  RunHarmony(c("slide","sample_id"), 
             assay.use = "SCT",
             plot_convergence = TRUE)


ste <- FindNeighbors(ste, 
                     reduction = "harmony", 
                     k.param = 30,
                     dims = 1:30)

ste <- RunUMAP(ste, 
               reduction = "harmony", 
               n.neighbors = 30L,
               min.dist = 0.3,
               dims = 1:30)

print("Clustering...")
# Louvain
snn_pcs <- BuildSNNSeurat(ste[["harmony"]]@cell.embeddings[,1:30], 
                          nn.eps = 0)

resolution_list <- c(0.2, 0.4, 0.6, 0.8, 1.0)
ids_cos <- Reduce(cbind, parallel::mclapply(resolution_list, function(res_use) {
  Seurat:::RunModularityClustering(SNN = snn_pcs, 
                                   modularity = 1, 
                                   resolution = res_use, 
                                   algorithm = 3, 
                                   n.start = 10, 
                                   n.iter = 10, random.seed = 0, print.output = FALSE, 
                                   temp.file.location = NULL, edge.file.name = NULL)    
}, mc.cores = min(16, length(resolution_list))))
ids_cos %<>% data.frame()
colnames(ids_cos) <- sprintf("res_%.2f", resolution_list)

ids_cos$res_0.20 <- as.character(ids_cos$res_0.20)
ids_cos$res_0.40 <- as.character(ids_cos$res_0.40)
ids_cos$res_0.60 <- as.character(ids_cos$res_0.60)
ids_cos$res_0.80 <- as.character(ids_cos$res_0.80)
ids_cos$res_1.00 <- as.character(ids_cos$res_1.00)

rownames(ids_cos) = rownames(ste@meta.data)
ids_cos <- ids_cos %>%
  dplyr::mutate(across(everything(), ~ factor(.x, levels = sort(unique(.x)))))
ste <- AddMetaData(ste, ids_cos %>%
                     magrittr::set_colnames(paste0(colnames(.),"_",celltype)))
head(ste@meta.data)

ste@meta.data$seurat_clusters = ste@meta.data[,paste0("res_",resolution,"_",celltype)]
Idents(ste) = factor(ste@meta.data$seurat_clusters, levels = as.character(sort(as.integer(unique(ste@meta.data$seurat_clusters)))-1))
table(ste@meta.data$seurat_clusters)

## harmonoized data [hPC1-20]
lisi_res_harmony <- rbind(
  lisi::compute_lisi(ste[["harmony"]]@cell.embeddings, 
                     ste@meta.data,
                     c("slide", "seurat_clusters", "sample_id"))
) %>% 
  tidyr::gather(key, val, slide, seurat_clusters, sample_id)
print(head(lisi_res_harmony))
saveRDS(lisi_res_harmony, file=paste0("/path/to/jinamo/JIA_Xenium/output/LISI_harmony_",celltype,".rds"))

## pre-harmonoized data [PC1-20]
lisi_res_preharmony <- rbind(
  lisi::compute_lisi(ste[["pca"]]@cell.embeddings, 
                     ste@meta.data,
                     c("slide", "seurat_clusters", "sample_id"))
) %>% 
  tidyr::gather(key, val, slide, seurat_clusters, sample_id)
print(head(lisi_res_preharmony))
saveRDS(lisi_res_preharmony, file=paste0("/path/to/jinamo/JIA_Xenium/output/LISI_PCA_",celltype,".rds"))

table(ste@meta.data$doublet_class, ste@meta.data$seurat_clusters)
table(ste@meta.data$doublet_class, ste@meta.data$sample_id)


saveRDS(ste, file = paste0("/path/to/jinamo/JIA_Xenium/output/Xenium_",celltype,".rds"))

# AMP RA-synovium paper (Nature 2023)
syno.rna = readRDS("/path/to/jinamo/scRNA_ATACseq/data/RA_syno/syno.rna_Fan.rds")
ra_celltype <- switch(celltype,
                      "T" = c("T cell","NK"),
                      "M" = "Myeloid",
                      "B" = "B cell",
                      "Tissue" = c("Fibroblast","Endothelial"))

head(syno.rna@meta.data)
table(syno.rna@meta.data$ct_subtype)
table(syno.rna@meta.data$ct)
syno.rna = subset(syno.rna, ct %in% ra_celltype)
table(syno.rna@meta.data$ct)

ra = syno.rna@assays$RNA$data
jia = ste@assays$SCT$data

shared_mrna = intersect(rownames(jia),
                        rownames(ra))
length(shared_mrna)
print(shared_mrna)

ra = ra[shared_mrna,]
jia = jia[shared_mrna,]

colnames(ra) = paste0("RA:",colnames(ra))
colnames(jia) = paste0("JIA:",colnames(jia))

print("mRNAs in RA scRNA-seq;")
dim(ra)

print("mRNAs in JIA Xenium;")
dim(jia)

assay_list_indirect = list(
  RA = ra,
  JIA = jia
)

lapply(assay_list_indirect, dim)
lapply(assay_list_indirect, class)

stab_indirect = stabMap(assay_list_indirect,
                        reference_list = c("RA"),
                        maxFeatures = Inf,
                        ncomponentsReference = 30,
                        ncomponentsSubset = 30,
                        plot = FALSE)

set.seed(1234)
stab_indirect_umap <- uwot::umap(as.data.frame(stab_indirect), 
                                 n_neighbors = 30, 
                                 metric = "cosine", 
                                 min_dist = .1, 
                                 n_threads = 4,
                                 verbose = TRUE)
stab_indirect_umap = stab_indirect_umap %>% 
  magrittr::set_colnames(c("UMAP1","UMAP2")) %>%
  as.data.frame()


train_data <- stab_indirect[grepl("^RA:",rownames(stab_indirect)),]  
train_labels <- syno.rna@meta.data$ct_subtype
test_data <- stab_indirect[!grepl("^RA:",rownames(stab_indirect)),]


knn_pred <- class::knn(train = train_data, test = test_data, cl = train_labels, k = 30, prob = TRUE)
knn_prob = attributes(knn_pred)$prob
test_data = data.frame(test_data,
                       ct_subtype = knn_pred,
                       knn_prob = knn_prob)
train_data = data.frame(train_data,
                        ct_subtype = train_labels,
                        knn_prob = 1)
stab_indirect = rbind(train_data,test_data)

stab_indirect = cbind(stab_indirect,stab_indirect_umap)
head(stab_indirect)
tail(stab_indirect)

# Louvain
snn_pcs <- BuildSNNSeurat(stab_indirect[,1:30], nn.eps = 0)

resolution_list <- c(0.2, 0.4, 0.6, 0.8, 1.0)
ids_cos <- base::Reduce(cbind, parallel::mclapply(resolution_list, function(res_use) {
  Seurat:::RunModularityClustering(SNN = snn_pcs, modularity = 1, 
                                   resolution = res_use, algorithm = 3, n.start = 10, 
                                   n.iter = 10, random.seed = 0, print.output = FALSE, 
                                   temp.file.location = NULL, edge.file.name = NULL)    
}, mc.cores = min(16, length(resolution_list))))
ids_cos %<>% data.frame()
colnames(ids_cos) <- sprintf("res_%.2f", resolution_list)

stab_indirect$res_0.20_stab <- as.character(ids_cos$res_0.20)
stab_indirect$res_0.40_stab <- as.character(ids_cos$res_0.40)
stab_indirect$res_0.60_stab <- as.character(ids_cos$res_0.60)
stab_indirect$res_0.80_stab <- as.character(ids_cos$res_0.80)
stab_indirect$res_1.00_stab <- as.character(ids_cos$res_1.00)

saveRDS(stab_indirect, paste0("/path/to/jinamo/JIA_Xenium/output/StabMap_Obj_RAtoJIA_",celltype,".rds"))



celltype="M"

ste = readRDS(paste0("/path/to/jinamo/JIA_Xenium/output/Xenium_merged.rds"))

resolution <- switch(celltype,
                     "ALL" = "0.20",
                     "T" = "1.00",
                     "M" = "0.20",
                     "B" = "0.40",
                     "Tissue" = "0.60"
)

cell_ids = readRDS(paste0("/path/to/jinamo/JIA_Xenium/output/",celltype,"_cellids.rds"))
ste = subset(ste, cells=cell_ids)
table(ste@meta.data$seurat_clusters)

ste <- SCTransform(ste, assay = "Spatial",
                   conserve.memory = TRUE, vst.flavor="v2")

VariableFeatures(ste[["SCT"]]) <- rownames(ste[["SCT"]]$data)
set.seed(1234)
ste <- RunPCA(ste, assay = "SCT", verbose = FALSE)

print("Harmony...")
set.seed(1234)
ste <- ste %>% 
  RunHarmony(c("slide","sample_id"), 
             assay.use = "SCT",
             plot_convergence = TRUE)


ste <- FindNeighbors(ste, 
                     reduction = "harmony", 
                     k.param = 30,
                     dims = 1:30)

ste <- RunUMAP(ste, 
               reduction = "harmony", 
               n.neighbors = 30L,
               min.dist = 0.3,
               dims = 1:30)

print("Clustering...")
# Louvain
snn_pcs <- BuildSNNSeurat(ste[["harmony"]]@cell.embeddings[,1:30], 
                          nn.eps = 0)

resolution_list <- c(0.2, 0.4, 0.6, 0.8, 1.0)
ids_cos <- Reduce(cbind, parallel::mclapply(resolution_list, function(res_use) {
  Seurat:::RunModularityClustering(SNN = snn_pcs, 
                                   modularity = 1, 
                                   resolution = res_use, 
                                   algorithm = 3, 
                                   n.start = 10, 
                                   n.iter = 10, random.seed = 0, print.output = FALSE, 
                                   temp.file.location = NULL, edge.file.name = NULL)    
}, mc.cores = min(16, length(resolution_list))))
ids_cos %<>% data.frame()
colnames(ids_cos) <- sprintf("res_%.2f", resolution_list)

ids_cos$res_0.20 <- as.character(ids_cos$res_0.20)
ids_cos$res_0.40 <- as.character(ids_cos$res_0.40)
ids_cos$res_0.60 <- as.character(ids_cos$res_0.60)
ids_cos$res_0.80 <- as.character(ids_cos$res_0.80)
ids_cos$res_1.00 <- as.character(ids_cos$res_1.00)

rownames(ids_cos) = rownames(ste@meta.data)
ids_cos <- ids_cos %>%
  dplyr::mutate(across(everything(), ~ factor(.x, levels = sort(unique(.x)))))
ste <- AddMetaData(ste, ids_cos %>%
                     magrittr::set_colnames(paste0(colnames(.),"_",celltype)))
head(ste@meta.data)

ste@meta.data$seurat_clusters = ste@meta.data[,paste0("res_",resolution,"_",celltype)]
Idents(ste) = factor(ste@meta.data$seurat_clusters, levels = as.character(sort(as.integer(unique(ste@meta.data$seurat_clusters)))-1))
table(ste@meta.data$seurat_clusters)

## harmonoized data [hPC1-20]
lisi_res_harmony <- rbind(
  lisi::compute_lisi(ste[["harmony"]]@cell.embeddings, 
                     ste@meta.data,
                     c("slide", "seurat_clusters", "sample_id"))
) %>% 
  tidyr::gather(key, val, slide, seurat_clusters, sample_id)
print(head(lisi_res_harmony))
saveRDS(lisi_res_harmony, file=paste0("/path/to/jinamo/JIA_Xenium/output/LISI_harmony_",celltype,".rds"))

## pre-harmonoized data [PC1-20]
lisi_res_preharmony <- rbind(
  lisi::compute_lisi(ste[["pca"]]@cell.embeddings, 
                     ste@meta.data,
                     c("slide", "seurat_clusters", "sample_id"))
) %>% 
  tidyr::gather(key, val, slide, seurat_clusters, sample_id)
print(head(lisi_res_preharmony))
saveRDS(lisi_res_preharmony, file=paste0("/path/to/jinamo/JIA_Xenium/output/LISI_PCA_",celltype,".rds"))

table(ste@meta.data$doublet_class, ste@meta.data$seurat_clusters)
table(ste@meta.data$doublet_class, ste@meta.data$sample_id)


saveRDS(ste, file = paste0("/path/to/jinamo/JIA_Xenium/output/Xenium_",celltype,".rds"))

# AMP RA-synovium paper (Nature 2023)
syno.rna = readRDS("/path/to/jinamo/scRNA_ATACseq/data/RA_syno/syno.rna_Fan.rds")
ra_celltype <- switch(celltype,
                      "T" = c("T cell","NK"),
                      "M" = "Myeloid",
                      "B" = "B cell",
                      "Tissue" = c("Fibroblast","Endothelial"))

head(syno.rna@meta.data)
table(syno.rna@meta.data$ct_subtype)
table(syno.rna@meta.data$ct)
syno.rna = subset(syno.rna, ct %in% ra_celltype)
table(syno.rna@meta.data$ct)

ra = syno.rna@assays$RNA$data
jia = ste@assays$SCT$data

shared_mrna = intersect(rownames(jia),
                        rownames(ra))
length(shared_mrna)
print(shared_mrna)

ra = ra[shared_mrna,]
jia = jia[shared_mrna,]

colnames(ra) = paste0("RA:",colnames(ra))
colnames(jia) = paste0("JIA:",colnames(jia))

print("mRNAs in RA scRNA-seq;")
dim(ra)

print("mRNAs in JIA Xenium;")
dim(jia)

assay_list_indirect = list(
  RA = ra,
  JIA = jia
)

lapply(assay_list_indirect, dim)
lapply(assay_list_indirect, class)

stab_indirect = stabMap(assay_list_indirect,
                        reference_list = c("RA"),
                        maxFeatures = Inf,
                        ncomponentsReference = 30,
                        ncomponentsSubset = 30,
                        plot = FALSE)

set.seed(1234)
stab_indirect_umap <- uwot::umap(as.data.frame(stab_indirect), 
                                 n_neighbors = 30, 
                                 metric = "cosine", 
                                 min_dist = .1, 
                                 n_threads = 4,
                                 verbose = TRUE)
stab_indirect_umap = stab_indirect_umap %>% 
  magrittr::set_colnames(c("UMAP1","UMAP2")) %>%
  as.data.frame()


train_data <- stab_indirect[grepl("^RA:",rownames(stab_indirect)),]  
train_labels <- syno.rna@meta.data$ct_subtype
test_data <- stab_indirect[!grepl("^RA:",rownames(stab_indirect)),]


knn_pred <- class::knn(train = train_data, test = test_data, cl = train_labels, k = 30, prob = TRUE)
knn_prob = attributes(knn_pred)$prob
test_data = data.frame(test_data,
                       ct_subtype = knn_pred,
                       knn_prob = knn_prob)
train_data = data.frame(train_data,
                        ct_subtype = train_labels,
                        knn_prob = 1)
stab_indirect = rbind(train_data,test_data)

stab_indirect = cbind(stab_indirect,stab_indirect_umap)
head(stab_indirect)
tail(stab_indirect)

# Louvain
snn_pcs <- BuildSNNSeurat(stab_indirect[,1:30], nn.eps = 0)

resolution_list <- c(0.2, 0.4, 0.6, 0.8, 1.0)
ids_cos <- base::Reduce(cbind, parallel::mclapply(resolution_list, function(res_use) {
  Seurat:::RunModularityClustering(SNN = snn_pcs, modularity = 1, 
                                   resolution = res_use, algorithm = 3, n.start = 10, 
                                   n.iter = 10, random.seed = 0, print.output = FALSE, 
                                   temp.file.location = NULL, edge.file.name = NULL)    
}, mc.cores = min(16, length(resolution_list))))
ids_cos %<>% data.frame()
colnames(ids_cos) <- sprintf("res_%.2f", resolution_list)

stab_indirect$res_0.20_stab <- as.character(ids_cos$res_0.20)
stab_indirect$res_0.40_stab <- as.character(ids_cos$res_0.40)
stab_indirect$res_0.60_stab <- as.character(ids_cos$res_0.60)
stab_indirect$res_0.80_stab <- as.character(ids_cos$res_0.80)
stab_indirect$res_1.00_stab <- as.character(ids_cos$res_1.00)

saveRDS(stab_indirect, paste0("/path/to/jinamo/JIA_Xenium/output/StabMap_Obj_RAtoJIA_",celltype,".rds"))



celltype="Tissue"

ste = readRDS(paste0("/path/to/jinamo/JIA_Xenium/output/Xenium_merged.rds"))

resolution <- switch(celltype,
                     "ALL" = "0.20",
                     "T" = "1.00",
                     "M" = "0.20",
                     "B" = "0.40",
                     "Tissue" = "0.60"
)

cell_ids = readRDS(paste0("/path/to/jinamo/JIA_Xenium/output/",celltype,"_cellids.rds"))
ste = subset(ste, cells=cell_ids)
table(ste@meta.data$seurat_clusters)

ste <- SCTransform(ste, assay = "Spatial",
                   conserve.memory = TRUE, vst.flavor="v2")

VariableFeatures(ste[["SCT"]]) <- rownames(ste[["SCT"]]$data)
set.seed(1234)
ste <- RunPCA(ste, assay = "SCT", verbose = FALSE)

print("Harmony...")
set.seed(1234)
ste <- ste %>% 
  RunHarmony(c("slide","sample_id"), 
             assay.use = "SCT",
             plot_convergence = TRUE)


ste <- FindNeighbors(ste, 
                     reduction = "harmony", 
                     k.param = 30,
                     dims = 1:30)

ste <- RunUMAP(ste, 
               reduction = "harmony", 
               n.neighbors = 30L,
               min.dist = 0.3,
               dims = 1:30)

print("Clustering...")
# Louvain
snn_pcs <- BuildSNNSeurat(ste[["harmony"]]@cell.embeddings[,1:30], 
                          nn.eps = 0)

resolution_list <- c(0.2, 0.4, 0.6, 0.8, 1.0)
ids_cos <- Reduce(cbind, parallel::mclapply(resolution_list, function(res_use) {
  Seurat:::RunModularityClustering(SNN = snn_pcs, 
                                   modularity = 1, 
                                   resolution = res_use, 
                                   algorithm = 3, 
                                   n.start = 10, 
                                   n.iter = 10, random.seed = 0, print.output = FALSE, 
                                   temp.file.location = NULL, edge.file.name = NULL)    
}, mc.cores = min(16, length(resolution_list))))
ids_cos %<>% data.frame()
colnames(ids_cos) <- sprintf("res_%.2f", resolution_list)

ids_cos$res_0.20 <- as.character(ids_cos$res_0.20)
ids_cos$res_0.40 <- as.character(ids_cos$res_0.40)
ids_cos$res_0.60 <- as.character(ids_cos$res_0.60)
ids_cos$res_0.80 <- as.character(ids_cos$res_0.80)
ids_cos$res_1.00 <- as.character(ids_cos$res_1.00)

rownames(ids_cos) = rownames(ste@meta.data)
ids_cos <- ids_cos %>%
  dplyr::mutate(across(everything(), ~ factor(.x, levels = sort(unique(.x)))))
ste <- AddMetaData(ste, ids_cos %>%
                     magrittr::set_colnames(paste0(colnames(.),"_",celltype)))
head(ste@meta.data)

ste@meta.data$seurat_clusters = ste@meta.data[,paste0("res_",resolution,"_",celltype)]
Idents(ste) = factor(ste@meta.data$seurat_clusters, levels = as.character(sort(as.integer(unique(ste@meta.data$seurat_clusters)))-1))
table(ste@meta.data$seurat_clusters)

## harmonoized data [hPC1-20]
lisi_res_harmony <- rbind(
  lisi::compute_lisi(ste[["harmony"]]@cell.embeddings, 
                     ste@meta.data,
                     c("slide", "seurat_clusters", "sample_id"))
) %>% 
  tidyr::gather(key, val, slide, seurat_clusters, sample_id)
print(head(lisi_res_harmony))
saveRDS(lisi_res_harmony, file=paste0("/path/to/jinamo/JIA_Xenium/output/LISI_harmony_",celltype,".rds"))

## pre-harmonoized data [PC1-20]
lisi_res_preharmony <- rbind(
  lisi::compute_lisi(ste[["pca"]]@cell.embeddings, 
                     ste@meta.data,
                     c("slide", "seurat_clusters", "sample_id"))
) %>% 
  tidyr::gather(key, val, slide, seurat_clusters, sample_id)
print(head(lisi_res_preharmony))
saveRDS(lisi_res_preharmony, file=paste0("/path/to/jinamo/JIA_Xenium/output/LISI_PCA_",celltype,".rds"))

table(ste@meta.data$doublet_class, ste@meta.data$seurat_clusters)
table(ste@meta.data$doublet_class, ste@meta.data$sample_id)


saveRDS(ste, file = paste0("/path/to/jinamo/JIA_Xenium/output/Xenium_",celltype,".rds"))

# AMP RA-synovium paper (Nature 2023)
syno.rna = readRDS("/path/to/jinamo/scRNA_ATACseq/data/RA_syno/syno.rna_Fan.rds")
ra_celltype <- switch(celltype,
                      "T" = c("T cell","NK"),
                      "M" = "Myeloid",
                      "B" = "B cell",
                      "Tissue" = c("Fibroblast","Endothelial"))


head(syno.rna@meta.data)
table(syno.rna@meta.data$ct_subtype)
table(syno.rna@meta.data$ct)
syno.rna = subset(syno.rna, ct %in% ra_celltype)
table(syno.rna@meta.data$ct)

ra = syno.rna@assays$RNA$data
jia = ste@assays$SCT$data

shared_mrna = intersect(rownames(jia),
                        rownames(ra))
length(shared_mrna)
print(shared_mrna)

ra = ra[shared_mrna,]
jia = jia[shared_mrna,]

colnames(ra) = paste0("RA:",colnames(ra))
colnames(jia) = paste0("JIA:",colnames(jia))

print("mRNAs in RA scRNA-seq;")
dim(ra)

print("mRNAs in JIA Xenium;")
dim(jia)

assay_list_indirect = list(
  RA = ra,
  JIA = jia
)

lapply(assay_list_indirect, dim)
lapply(assay_list_indirect, class)

stab_indirect = stabMap(assay_list_indirect,
                        reference_list = c("RA"),
                        maxFeatures = Inf,
                        ncomponentsReference = 30,
                        ncomponentsSubset = 30,
                        plot = FALSE)

set.seed(1234)
stab_indirect_umap <- uwot::umap(as.data.frame(stab_indirect), 
                                 n_neighbors = 30, 
                                 metric = "cosine", 
                                 min_dist = .1, 
                                 n_threads = 4,
                                 verbose = TRUE)
stab_indirect_umap = stab_indirect_umap %>% 
  magrittr::set_colnames(c("UMAP1","UMAP2")) %>%
  as.data.frame()


train_data <- stab_indirect[grepl("^RA:",rownames(stab_indirect)),]  
train_labels <- syno.rna@meta.data$ct_subtype
test_data <- stab_indirect[!grepl("^RA:",rownames(stab_indirect)),]


knn_pred <- class::knn(train = train_data, test = test_data, cl = train_labels, k = 30, prob = TRUE)
knn_prob = attributes(knn_pred)$prob
test_data = data.frame(test_data,
                       ct_subtype = knn_pred,
                       knn_prob = knn_prob)
train_data = data.frame(train_data,
                        ct_subtype = train_labels,
                        knn_prob = 1)
stab_indirect = rbind(train_data,test_data)

stab_indirect = cbind(stab_indirect,stab_indirect_umap)
head(stab_indirect)
tail(stab_indirect)

# Louvain
snn_pcs <- BuildSNNSeurat(stab_indirect[,1:30], nn.eps = 0)

resolution_list <- c(0.2, 0.4, 0.6, 0.8, 1.0)
ids_cos <- base::Reduce(cbind, parallel::mclapply(resolution_list, function(res_use) {
  Seurat:::RunModularityClustering(SNN = snn_pcs, modularity = 1, 
                                   resolution = res_use, algorithm = 3, n.start = 10, 
                                   n.iter = 10, random.seed = 0, print.output = FALSE, 
                                   temp.file.location = NULL, edge.file.name = NULL)    
}, mc.cores = min(16, length(resolution_list))))
ids_cos %<>% data.frame()
colnames(ids_cos) <- sprintf("res_%.2f", resolution_list)

stab_indirect$res_0.20_stab <- as.character(ids_cos$res_0.20)
stab_indirect$res_0.40_stab <- as.character(ids_cos$res_0.40)
stab_indirect$res_0.60_stab <- as.character(ids_cos$res_0.60)
stab_indirect$res_0.80_stab <- as.character(ids_cos$res_0.80)
stab_indirect$res_1.00_stab <- as.character(ids_cos$res_1.00)

saveRDS(stab_indirect, paste0("/path/to/jinamo/JIA_Xenium/output/StabMap_Obj_RAtoJIA_",celltype,".rds"))



ste = readRDS(paste0("/path/to/jinamo/JIA_Xenium/output/Xenium_merged.rds"))
head(ste@meta.data)

stab_indirect = data.frame()
for(celltype in c("T","M","B","Tissue")){
  stab_indirect = rbind(stab_indirect,
                        readRDS(paste0("/path/to/jinamo/JIA_Xenium/output/StabMap_Obj_RAtoJIA_",celltype,".rds")))
  
}

stab_indirect %<>%
  dplyr::mutate(dataset = dplyr::case_when(
    grepl("^RA",rownames(stab_indirect)) ~ "RA-scRNAseq",
    grepl("^JIA",rownames(stab_indirect)) ~ "JIA-Xenium"),
    cell = gsub("RA:","",rownames(stab_indirect)) %>% gsub("JIA:","",.))
table(stab_indirect$dataset)
table(stab_indirect$ct_subtype)

ste <- AddMetaData(ste, 
                   stab_indirect[stab_indirect$dataset=="JIA-Xenium",c("cell",paste0("ct_subtype"))] %>% magrittr::set_colnames(c("cell","ct_subtype")) %>% magrittr::set_rownames(.$cell))
table(ste@meta.data$ct_subtype)
head(ste@meta.data)

meta_all = data.frame()
for(celltype in c("T","M","B","Tissue")){
  resolution <- switch(celltype,
                       "ALL" = "0.20",
                       "T" = "1.00",
                       "M" = "0.20",
                       "B" = "0.40",
                       "Tissue" = "0.60"
  )
  meta_all = rbind(meta_all,
                   readRDS(paste0("/path/to/jinamo/JIA_Xenium/output/Xenium_",celltype,".rds"))@meta.data[,c("cell",paste0("res_",resolution,"_",celltype))] %>% 
                     magrittr::set_colnames(c("cell","res_subtype")) %>%
                     dplyr::mutate(res_subtype = paste0(celltype,"-",res_subtype)))
}
ste <- AddMetaData(ste, 
                   meta_all %>% magrittr::set_rownames(.$cell))
head(ste@meta.data)

celltype="ALL"
resolution <- switch(celltype,
                     "ALL" = "0.20",
                     "T" = "1.00",
                     "M" = "0.20",
                     "B" = "0.40",
                     "Tissue" = "0.60"
)

cluster_col = "res_subtype" 

Idents(ste) = factor(ste@meta.data[,cluster_col])
DefaultAssay(ste) = "SCT"

cluster_df = read.csv("/path/to/jinamo/JIA_Xenium/cluster_annotation.csv")
cluster_df_ = dplyr::left_join(ste@meta.data,
                               cluster_df %>% dplyr::rename(res_subtype = cluster,
                                                            res_subtype_anno = annotation),
                               by="res_subtype")
ste = AddMetaData(ste,cluster_df_)

ste@meta.data$celltype = dplyr::case_when(
  ste@meta.data[,paste0("res_",resolution)] %in% c("1","2","4","6","9","10","12","13") ~ "Tissue", #13: Muskcle (MYH2,PYGM)
  ste@meta.data[,paste0("res_",resolution)] %in% c("0","7","11","14","15","16") ~ "M",
  ste@meta.data[,paste0("res_",resolution)] %in% c("3") ~ "T",
  ste@meta.data[,paste0("res_",resolution)] %in% c("5","8") ~ "B"
)
table(ste@meta.data$celltype)

colnames(ste@meta.data)

ste@meta.data$res_subtype = paste0(ste@meta.data$res_subtype,":",ste@meta.data$res_subtype_anno)

Idents(ste) = factor(ste@meta.data[,cluster_col])
DefaultAssay(ste) = "SCT"


ste <- nhood_enrichment(
  ste,
  cluster_key = "res_subtype", 
  neighbors.k = 30, 
  connectivity_key = "nn", 
  transformation = TRUE,
  n_perms = n_perms, seed = 1234, n_jobs = n_jobs
)

# Local Co-occur score 
ste@meta.data[,grep("^cooccur_local_",colnames(ste@meta.data))] = NULL
print("Local Co-occur score for fine-scaled cell types:")
cluster_pair_df = expand.grid(
  cluster1 = unique(ste@meta.data[,"res_subtype"]),
  cluster2 = unique(ste@meta.data[,"res_subtype"]),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(cluster1 != cluster2) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    sorted = list(sort(c(cluster1, cluster2)))
  ) %>%
  dplyr::mutate(
    cluster1_sorted = sorted[[1]],
    cluster2_sorted = sorted[[2]]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-sorted) %>%
  dplyr::distinct(cluster1_sorted, cluster2_sorted, .keep_all = TRUE) %>%
  as.data.frame() %>%
  dplyr::select(cluster1 , cluster2) %>%
  dplyr::arrange(cluster1, cluster2) 

result_list <- mclapply(
  X        = seq_len(nrow(cluster_pair_df)),
  FUN      = function(j) {
    cluster_x <- cluster_pair_df$cluster1[j]
    cluster_y <- cluster_pair_df$cluster2[j]
    
    df_tmp <- cooccur_local(
      ste,
      cluster_x        = cluster_x,
      cluster_y        = cluster_y,
      connectivity_key = "nn",
      cluster_key      = "res_subtype",
      sample_key       = "sample_id",
      neighbors.k      = 30, 
      radius           = 10
    )
    
    return(df_tmp)
  },
  mc.cores = n_jobs
)
good_idx <- sapply(result_list, function(x) !inherits(x, "try-error"))
result_lists <- result_list[good_idx]
cooccur_local_df <- do.call(cbind, result_lists)
colnames(cooccur_local_df) = colnames(cooccur_local_df) %>%
  gsub(":", "--", .) %>%
  gsub(" ", "XXX", .) %>%
  gsub("\\+", "YYY", .) %>%
  gsub("ï", "i", ., fixed = TRUE) %>%
  gsub("TPH_TFH", "TPHTFH", .)
if(!is.null(cooccur_local_df)){
  ste = AddMetaData(ste,cooccur_local_df)
}

# Local Co-occur score for broad cell types
print("Local Co-occur score for broad cell types:")
cluster_pair_df = expand.grid(
  cluster1 = unique(ste@meta.data[,"celltype"]),
  cluster2 = unique(ste@meta.data[,"celltype"]),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(cluster1 != cluster2) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    sorted = list(sort(c(cluster1, cluster2)))
  ) %>%
  dplyr::mutate(
    cluster1_sorted = sorted[[1]],
    cluster2_sorted = sorted[[2]]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-sorted) %>%
  dplyr::distinct(cluster1_sorted, cluster2_sorted, .keep_all = TRUE) %>%
  as.data.frame() %>%
  dplyr::select(cluster1 , cluster2) %>%
  dplyr::arrange(cluster1, cluster2) 

result_list <- mclapply(
  X        = seq_len(nrow(cluster_pair_df)),
  FUN      = function(j) {
    cluster_x <- cluster_pair_df$cluster1[j]
    cluster_y <- cluster_pair_df$cluster2[j]
    
    df_tmp <- cooccur_local(
      ste,
      cluster_x        = cluster_x,
      cluster_y        = cluster_y,
      connectivity_key = "nn",
      cluster_key      = "celltype",
      sample_key       = "sample_id",
      neighbors.k      = 30, 
      radius           = 10
    )
    
    return(df_tmp)
  },
  mc.cores = n_jobs
)
good_idx <- sapply(result_list, function(x) !inherits(x, "try-error"))
result_lists <- result_list[good_idx]
cooccur_local_df <- do.call(cbind, result_lists)
if(!is.null(cooccur_local_df)){
  ste = AddMetaData(ste,cooccur_local_df)
}

saveRDS(ste, file = paste0("/path/to/jinamo/JIA_Xenium/output/Xenium_merged_res_subtype.rds"))



cur_levels <- levels(Idents(ste))
df <- data.frame(level = cur_levels, stringsAsFactors = FALSE)
df$prefix <- sub("^([^-]+)-.*", "\\1", df$level)
df$num <- as.integer(sub("^[^-]+-(\\d+):.*", "\\1", df$level))
prefix_order <- c("T", "B", "M", "Tissue")
df$prefix_factor <- factor(df$prefix, levels = prefix_order)
df_sorted <- df[order(df$prefix_factor, df$num), ]
new_levels <- df_sorted$level
Idents(ste) <- factor(Idents(ste), levels = new_levels)

ste@meta.data$new_cluster = factor(stringr::str_split(ste@meta.data$res_subtype,pattern=":",simplify=TRUE)[,2], 
                                   levels = unique(stringr::str_split(levels(Idents(ste)),pattern=":",simplify=TRUE)[,2])) 
Idents(ste) = ste@meta.data$new_cluster

ord = levels(Idents(ste))
ord = c(ord[!grepl("venule|arteriolar|lymphatic|mixed",ord)],
        ord[grepl("venule|arteriolar|lymphatic|mixed",ord)])
Idents(ste) = factor(Idents(ste), levels = ord)
ord = levels(Idents(ste))
ord = c(ord[!grepl("mixed",ord)],
        ord[grepl("mixed",ord)])
Idents(ste) = factor(Idents(ste), levels = ord)
levels(Idents(ste))

cluster_col = "new_cluster"

# Define Niches
ste <- BuildNicheAssay_(object = ste, 
                        group.by = cluster_col,
                        neighbors.k = 30, 
                        connectivity_key = "nn",
                        niches.k = 30)

ste <- nhood_enrichment(
  ste,
  cluster_key = cluster_col, 
  neighbors.k = 30, 
  connectivity_key = "nn", 
  transformation = TRUE,
  n_perms = n_perms, seed = 1234, n_jobs = n_jobs
)

# Local Co-occur score 
ste@meta.data[,grep("^cooccur_local_",colnames(ste@meta.data))] = NULL
print("Local Co-occur score for fine-scaled cell types:")
cluster_pair_df = expand.grid(
  cluster1 = unique(ste@meta.data[,cluster_col]),
  cluster2 = unique(ste@meta.data[,cluster_col]),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(cluster1 != cluster2) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    sorted = list(sort(c(cluster1, cluster2)))
  ) %>%
  dplyr::mutate(
    cluster1_sorted = sorted[[1]],
    cluster2_sorted = sorted[[2]]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-sorted) %>%
  dplyr::distinct(cluster1_sorted, cluster2_sorted, .keep_all = TRUE) %>%
  as.data.frame() %>%
  dplyr::select(cluster1 , cluster2) %>%
  dplyr::arrange(cluster1, cluster2) 

result_list <- mclapply(
  X        = seq_len(nrow(cluster_pair_df)),
  FUN      = function(j) {
    cluster_x <- cluster_pair_df$cluster1[j]
    cluster_y <- cluster_pair_df$cluster2[j]
    
    df_tmp <- cooccur_local(
      ste,
      cluster_x        = cluster_x,
      cluster_y        = cluster_y,
      connectivity_key = "nn",
      cluster_key      = cluster_col,
      sample_key       = "sample_id",
      neighbors.k      = 30, 
      radius           = 10
    )
    
    return(df_tmp)
  },
  mc.cores = n_jobs
)
good_idx <- sapply(result_list, function(x) !inherits(x, "try-error"))
result_lists <- result_list[good_idx]
cooccur_local_df <- do.call(cbind, result_lists)
colnames(cooccur_local_df) = colnames(cooccur_local_df) %>%
  gsub(":", "--", .) %>%
  gsub(" ", "XXX", .) %>%
  gsub("\\+", "YYY", .) %>%
  gsub("ï", "i", ., fixed = TRUE) %>%
  gsub("TPH_TFH", "TPHTFH", .)
if(!is.null(cooccur_local_df)){
  ste = AddMetaData(ste,cooccur_local_df)
}

# Local Co-occur score for broad cell types
print("Local Co-occur score for broad cell types:")
cluster_pair_df = expand.grid(
  cluster1 = unique(ste@meta.data[,"celltype"]),
  cluster2 = unique(ste@meta.data[,"celltype"]),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(cluster1 != cluster2) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    sorted = list(sort(c(cluster1, cluster2)))
  ) %>%
  dplyr::mutate(
    cluster1_sorted = sorted[[1]],
    cluster2_sorted = sorted[[2]]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-sorted) %>%
  dplyr::distinct(cluster1_sorted, cluster2_sorted, .keep_all = TRUE) %>%
  as.data.frame() %>%
  dplyr::select(cluster1 , cluster2) %>%
  dplyr::arrange(cluster1, cluster2) 

result_list <- mclapply(
  X        = seq_len(nrow(cluster_pair_df)),
  FUN      = function(j) {
    cluster_x <- cluster_pair_df$cluster1[j]
    cluster_y <- cluster_pair_df$cluster2[j]
    
    df_tmp <- cooccur_local(
      ste,
      cluster_x        = cluster_x,
      cluster_y        = cluster_y,
      connectivity_key = "nn",
      cluster_key      = "celltype",
      sample_key       = "sample_id",
      neighbors.k      = 30, 
      radius           = 10
    )
    
    return(df_tmp)
  },
  mc.cores = n_jobs
)
good_idx <- sapply(result_list, function(x) !inherits(x, "try-error"))
result_lists <- result_list[good_idx]
cooccur_local_df <- do.call(cbind, result_lists)
if(!is.null(cooccur_local_df)){
  ste = AddMetaData(ste,cooccur_local_df)
}


radius_values <- seq(0,500,5) # μm
co_occ_results <- mclapply(
  X = radius_values,
  FUN = function(r) {
    calc_co_occurrence_for_radius(
      seurat_obj          = ste,
      radius              = r,
      sample_key          = "sample_id", 
      cluster_key         = cluster_col, 
      k                   = 30
    )
  },
  mc.cores = n_jobs
)
names(co_occ_results) <- paste0("radius=", radius_values)
ste@misc[[paste0("cooccur_results")]] <- co_occ_results

saveRDS(ste, file = paste0("/path/to/jinamo/JIA_Xenium/output/Xenium_merged_final.rds"))
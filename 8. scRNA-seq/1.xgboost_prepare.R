library(stringr)
load(file = "train_vali_full_obj.Rdata")
load("scRNA_no_umap.Rdata")
vali_scRNA_obj[["percent.mt"]] <- PercentageFeatureSet(vali_scRNA_obj, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(vali_scRNA_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vali_scRNA_obj <- subset(vali_scRNA_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
vali_scRNA_obj <- NormalizeData(vali_scRNA_obj)
vali_scRNA_obj <- FindVariableFeatures(vali_scRNA_obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(vali_scRNA_obj)
vali_scRNA_obj <- ScaleData(vali_scRNA_obj, features = all.genes)
vali_scRNA_obj <- RunPCA(vali_scRNA_obj, features = VariableFeatures(object = vali_scRNA_obj))
# ElbowPlot(vali_scRNA_obj,ndims = 30)

vali_scRNA_obj <- FindNeighbors(object = vali_scRNA_obj,dims = 1:30)
vali_scRNA_obj <- FindClusters(object = vali_scRNA_obj, resolution = 0.5)
vali_scRNA_obj <- RunUMAP(object = vali_scRNA_obj,dims = 1:30)

vali_umap <- as.data.frame(Embeddings(object = vali_scRNA_obj, reduction = "umap"))
vali_matrix <- as.data.frame(GetAssayData(object = vali_scRNA_obj, slot = "data"))
# 
# 
train_scRNA_obj[["percent.mt"]] <- PercentageFeatureSet(train_scRNA_obj, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(train_scRNA_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
train_scRNA_obj <- subset(train_scRNA_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
train_scRNA_obj <- NormalizeData(train_scRNA_obj)
train_scRNA_obj <- FindVariableFeatures(train_scRNA_obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(train_scRNA_obj)
train_scRNA_obj <- ScaleData(train_scRNA_obj, features = all.genes)
train_scRNA_obj <- RunPCA(train_scRNA_obj, features = VariableFeatures(object = train_scRNA_obj))
ElbowPlot(train_scRNA_obj,ndims = 30)

train_scRNA_obj <- FindNeighbors(object = train_scRNA_obj,dims = 1:30)
train_scRNA_obj <- FindClusters(object = train_scRNA_obj, resolution = 0.5)
train_scRNA_obj <- RunUMAP(object = train_scRNA_obj,dims = 1:30)
train_umap <- as.data.frame(Embeddings(object = train_scRNA_obj, reduction = "umap"))
train_matrix <- as.data.frame(GetAssayData(object = train_scRNA_obj, slot = "data"))

AML_191 <- read.table("AML_LS_RESULTS.txt",sep = "\t",header = F)[,1]

train_matrix_LS <- train_matrix[AML_LS,]%>%t()%>%as.data.frame()
train_matrix_LS$sample_type <- train_scRNA_obj@meta.data[["sample_type"]]

vali_matrix_LS <- vali_matrix[AML_LS,]%>%t()%>%as.data.frame()
vali_matrix_LS$sample_type <- vali_scRNA_obj@meta.data$sample_type

vali_scRNA_obj@meta.data$sample_type <- ifelse(grepl("./GRCh38.RNA/NL[12]", vali_scRNA_obj@meta.data$orig.ident), "normal", 
                      ifelse(grepl("./GRCh38.RNA/PT[1-8]A", vali_scRNA_obj@meta.data$orig.ident), "AML", NA))
vali_scRNA_obj <- subset(vali_scRNA_obj, subset = sample_type%in%c("AML","normal"))

save(train_umap,vali_umap,file = "umap.Rdata")
save(train_matrix_LS,vali_matrix_LS,file = "train_vali_191_matrix.Rdata")
save(train_scRNA_obj,vali_scRNA_obj,file = "train_vali_full_obj.Rdata")

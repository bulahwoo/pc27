library(ggplot2)
library(ggrepel)
library(tidyr)
library(dplyr)
library(Seurat)
library(rmarkdown)
library(DropletUtils)
library(patchwork)
library(RColorBrewer)
library(gt)
library(tximport)
setwd("/media/nguyen/Data1/github/bulahwoo/pc27/data/06_starsoloFiltered/")

`%ni%` <- Negate(`%in%`)

ggp_title <- expression(italic(P.)~italic(californicus)~"testes (PC_27)")

ggp_theme_bw_square_01 <-
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.title = element_text(color="black"),
        axis.text = element_text(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color="black", linewidth = 0.5), aspect.ratio = 1)

ggp_theme_bw_square_02 <-
  ggp_theme_bw_square_01 +
  theme(legend.title = element_text(size=10),
        legend.background=element_blank(),
        legend.justification=c(1, 0.85))

# load dge
dge_raw <- ReadSTARsolo("/media/nguyen/Data1/mao/scseq/select_polyt/starsolo_original/Solo.out/Gene/raw/")
dge <- ReadSTARsolo("/media/nguyen/Data1/mao/scseq/select_polyt/starsolo_original/Solo.out/Gene/filtered/")

# barcodeRanks()
 # This is a filtered matrix processed by `STARsolo` using `emptyDrops()`.

# emptyDrops()
 # This is a filtered matrix processed by `STARsolo` using `emptyDrops()`.

# a function for statistics
scRNAstat <- function(SEURATOBJ) {
  output <- c(
    dim(SEURATOBJ),
    mean(SEURATOBJ@meta.data$nCount_RNA),
    median(SEURATOBJ@meta.data$nCount_RNA),
    mean(SEURATOBJ@meta.data$nFeature_RNA),
    median(SEURATOBJ@meta.data$nFeature_RNA)
  )
  return(output)
}

# pre emptyDrops dge
so_pre_emptydrops <-
  CreateSeuratObject(counts=dge_raw) %>%
  PercentageFeatureSet(pattern = "^agat|^rrn", col.name = "percent.mt")
stat_pre_emptydrops <- scRNAstat(so_pre_emptydrops)
plt_vln_pre_emptydrops <- VlnPlot(so_pre_emptydrops, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ggsave("./files/plt_vln_pre_emptydrops.png", plt_vln_pre_emptydrops, height=12, width=18, units="cm")

# post emptyDrops dge
so_post_emptydrops <-
  CreateSeuratObject(counts=dge) %>%
  PercentageFeatureSet(pattern = "^agat|^rrn", col.name = "percent.mt")
stat_post_emptydrops <-scRNAstat(so_post_emptydrops)
plt_vln_post_emptydrops <- VlnPlot(so_post_emptydrops, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ggsave("./files/plt_vln_post_emptydrops.png", plt_vln_post_emptydrops, height=12, width=18, units="cm")

# Seurat filter dge 01
so_filter01 <-
  CreateSeuratObject(counts=dge, min.cells = 3) %>%
  PercentageFeatureSet(pattern = "^agat|^rrn", col.name = "percent.mt")
stat_filter01 <- scRNAstat(so_filter01)
plt_vln_filter01 <- VlnPlot(so_filter01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ggsave("./files/plt_vln_filter01.png", plt_vln_filter01, height=12, width=18, units="cm")

# Seurat filter dge 02
so_filter02 <-
  CreateSeuratObject(counts=dge, min.cells = 3, min.features = 200) %>%
  PercentageFeatureSet(pattern = "^agat|^rrn", col.name = "percent.mt")
stat_filter02 <- scRNAstat(so_filter02)
plt_vln_filter02 <- VlnPlot(so_filter02, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ggsave("./files/plt_vln_filter02.png", plt_vln_filter02, height=12, width=18, units="cm")

# Seurat filter dge 03
so_filter03 <-
  CreateSeuratObject(counts=dge, min.cells = 3, min.features = 200) %>%
  PercentageFeatureSet(pattern = "^agat|^rrn", col.name = "percent.mt") %>%
  subset(., nFeature_RNA > 200 & nFeature_RNA < 2500)
stat_filter03 <- scRNAstat(so_filter03)
plt_vln_filter03 <- VlnPlot(so_filter03, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ggsave("./files/plt_vln_filter03.png", plt_vln_filter03, height=12, width=18, units="cm")

# Seurat filter dge 04
so_filter04 <-
  CreateSeuratObject(counts=dge, min.cells = 3, min.features = 200) %>%
  PercentageFeatureSet(pattern = "^agat|^rrn", col.name = "percent.mt") %>%
  subset(., percent.mt < 5)
stat_filter04 <- scRNAstat(so_filter04)
plt_vln_filter04 <- VlnPlot(so_filter04, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ggsave("./files/plt_vln_filter04.png", plt_vln_filter04, height=12, width=18, units="cm")

# Seurat filter dge
so_filter <-
  CreateSeuratObject(counts=dge, min.cells = 3, min.features = 200) %>%
  PercentageFeatureSet(pattern = "^agat|^rrn", col.name = "percent.mt") %>%
  subset(., nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
stat_filter <- scRNAstat(so_filter)
plt_vln_filter <- VlnPlot(so_filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ggsave("./files/plt_vln_filter.png", plt_vln_filter, height=12, width=18, units="cm")

# Concatenate stat output
summary_stat <-
  rbind.data.frame(stat_pre_emptydrops,
                   stat_post_emptydrops,
                   stat_filter01,
                   stat_filter02,
                   stat_filter03,
                   stat_filter04,
                   stat_filter, make.row.names = F)

colnames(summary_stat) <- c("gene", "cell", "mean_umi", "median_umi", "mean_gene", "median_gene")

write.table(summary_stat, "./files/summary_stat.tsv", quote=F, sep="\t", row.names = F)

# STARsolo summary
summary_star <- read.table("/media/nguyen/Data1/mao/scseq/select_polyt/starsolo_original/Solo.out/Gene/Summary.csv", sep="\t")
write.table(summary_star, "./files/summary_star.tsv", quote=F, sep="\t", row.names=F, col.names=c("Measurement", "Value"))

# UMAP
so_pc27 <-
  CreateSeuratObject(counts=dge, min.cells = 3, min.features = 200) %>%
  PercentageFeatureSet(pattern = "^agat|^rrn", col.name = "percent.mt") %>% subset(., nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) %>%
  SCTransform(vars.to.regress="percent.mt") %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  RunUMAP(dims = 1:30) %>%
  FindClusters()
df_umap <- so_pc27@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(color=so_pc27@meta.data$seurat_clusters)
num_cluster <- length(unique(df_umap$color))
my_color <- sample(c(brewer.pal(name="Set2", n=8), brewer.pal(name="Dark2", n=8)), num_cluster, replace=F)

plt_umap <-
  ggplot(df_umap) +
  geom_point(aes(x=umap_1, y=umap_2, color=color), size=0.8) +
  geom_text_repel(data=df_umap %>% group_by(color) %>% summarise(q1=quantile(umap_1, 0.5), q2=quantile(umap_2, 0.5)),
                  aes(x=q1, y=q2, label=paste(LETTERS[1:num_cluster], summary(so_pc27@active.ident), sep="_")), size=6) +
  labs(title=ggp_title,
       x="UMAP_1",
       y="UMAP_2") +
  scale_color_manual(values = my_color, name="clusters", labels=LETTERS[1:num_cluster]) +
  guides(color = guide_legend(override.aes = list(size = 5))) + ggp_theme_bw_square_02

ggsave("./files/umap.png", plt_umap, height=15, width=15, units="cm")

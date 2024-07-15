Single-cell RNA sequencing of *Phymatopus californicus* samples
================
Bulah Wu
July 15, 2024

## Genome assembly and annotation

The *P. californicus* genome was sequenced in [Petr Nguyen’s
lab](https://www.entu.cas.cz/en/departments/department-of-ecology-and-conservation-biology/laboratory-of-evolutionary-genetics/)
using Oxford Nanopore technology and assembled with
[Flye](https://github.com/mikolmogorov/Flye). Annotation of the genome
was performed with [BRAKER3](https://github.com/Gaius-Augustus/BRAKER).
The RNA-seq data utilized include the publicly accessible
[SRR1021622](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1021622) and
additional data provided by [Michal
Zurovec](https://www.entu.cas.cz/en/staff/profile/269-michalzurovec/).
Mitochondrial genome was identified by aligning the genome assembly
against the mt genomes of *Bombyx mori* and *Yponomeuta evonymella*
using [minimap2](https://github.com/lh3/minimap2). The annotation of mt
genome was carried out with
[mitos](https://usegalaxy.org.au/root?tool_id=toolshed.g2.bx.psu.edu/repos/iuc/mitos2/mitos2/2.1.3+galaxy0)
on the Galaxy platform.

## Library preparation

The scRNA-seq library was prepared with the RNAdia 2.0 kit (Dolomite
Bio), following the manufacturer’s instructions. Approximately 8,000
cells were captured, and the size distribution of the library was
examined using the Bioanalyzer system.

<div class="figure" style="text-align: center">

<img src="data/bioanalyzer_pc27.png" alt="Size distribution assessed by Bioanalyzer" width="60%" />
<p class="caption">
Size distribution assessed by Bioanalyzer
</p>

</div>

## Quality control checks on raw reads

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) was
used for quality check. Output can be found
[here](/media/nguyen/Data1/github/bulahwoo/scRNAseq_Phymatopus_californicus/data/fastqc_pc27_read1_pre.html)
and
[here](/media/nguyen/Data1/github/bulahwoo/scRNAseq_Phymatopus_californicus/data/fastqc_pc27_read2_pre.html).

## Drop-seq protocol

The digital expression matrix using the Drop-seq protocol developed by
the McCarroll Lab (<https://github.com/broadinstitute/Drop-seq>). Here
are the programs used:

    #!/bin/bash
    #PBS -N pc27mito
    #PBS -l select=1:ncpus=16:mem=140gb:scratch_ssd=800gb
    #PBS -l walltime=24:00:00
    #PBS -m ae

    #export TMPDIR=$SCRATCHDIR
    DATADIR=/auto/plzen1/home/bulah/mao/scseq/dropseq/pc27mito
    echo "||| $HOSTNAME $SCRATCHDIR ||| (pc27mito) ||| ssh $HOSTNAME \"scp -rp bulah@$HOSTNAME:$SCRATCHDIR/* \$HOSTNAME:\$SCRATCHDIR\" |||" >> $DATADIR/jobid

    module add mambaforge picard star samtools
    mamba activate /auto/plzen1/home/bulah/.conda/envs/openjdk2201

    export ID=pc27mito
    export ASSEMBLY=/auto/plzen1/home/bulah/mao/braker3/phymatopus_californicus/genome_mao/repeatmodeler/phycal_genome_masked2.fa
    export ANNOTATION=/auto/plzen1/home/bulah/mao/braker3/phymatopus_californicus/genome_mao/braker3/phycal_mito.gtf
    export READ1=/auto/plzen1/home/bulah/mao/scseq/X201SC24042690-Z01-F001/X201SC24042690-Z01-F001/01.RawData/PC_27/PC_27_1.fq.gz
    export READ2=/auto/plzen1/home/bulah/mao/scseq/X201SC24042690-Z01-F001/X201SC24042690-Z01-F001/01.RawData/PC_27/PC_27_2.fq.gz

    mkdir -p $SCRATCHDIR/dropseq/genome $SCRATCHDIR/dropseq/project_${ID}/{tmp,raw,output_${ID}}
    export DROPSEQDIR=$SCRATCHDIR/dropseq
     export GENOMEDIR=$DROPSEQDIR/genome
     export WORKDIR=$DROPSEQDIR/project_${ID}
      export TMPDIR=$WORKDIR/tmp
      export RAWDIR=$WORKDIR/raw
      export OUTPUTDIR=$WORKDIR/output_${ID}
    export PATH=$PATH:$DROPSEQDIR/dropseq-3.0.1

    cd $DROPSEQDIR

    if [ ! -d "$DROPSEQDIR"/dropseq-3.0.1 ]; then
      wget https://github.com/broadinstitute/Drop-seq/releases/download/v3.0.1/dropseq-3.0.1.zip
      unzip dropseq-3.0.1.zip
      rm -rf dropseq-3.0.1.zip
    fi

    if [ ! -f "$GENOMEDIR"/genome.fa ]; then
      cp -rp $ASSEMBLY $GENOMEDIR/genome.fa
      STAR --runMode genomeGenerate --genomeDir $GENOMEDIR --genomeFastaFiles $GENOMEDIR/genome.fa --genomeSAindexNbases 13 --runThreadN $PBS_NCPUS
    fi

    if [ ! -f "$GENOMEDIR"/genome.gtf ]; then
      cp -rp $ANNOTATION $GENOMEDIR/genome.gtf
    fi

    if [ ! -f "$GENOMEDIR"/genome.dict ]; then
      picard CreateSequenceDictionary -R $GENOMEDIR/genome.fa --TMP_DIR $TMPDIR
    fi

    cp -rp $READ1 $RAWDIR/read1.fq.gz
    cp -rp $READ2 $RAWDIR/read2.fq.gz

    picard FastqToSam -F1 $RAWDIR/read1.fq.gz -F2 $RAWDIR/read2.fq.gz -O $OUTPUTDIR/unaligned.bam -SM ${ID} --TMP_DIR $TMPDIR

    TagBamWithReadSequenceExtended -SUMMARY $OUTPUTDIR/unaligned_tagged_Cellular.bam_summary.txt -BASE_RANGE 3-14 -BASE_QUALITY 10 -BARCODED_READ 1 -DISCARD_READ false -TAG_NAME XC \
    -NUM_BASES_BELOW_QUALITY 1 -INPUT $OUTPUTDIR/unaligned.bam -OUTPUT $OUTPUTDIR/unaligned_tagged_Cell.bam -TMP_DIR $TMPDIR

    TagBamWithReadSequenceExtended -SUMMARY $OUTPUTDIR/unaligned_tagged_Molecular.bam_summary.txt -BASE_RANGE 15-28 -BASE_QUALITY 10 -BARCODED_READ 1 -DISCARD_READ true -TAG_NAME XM \
    -NUM_BASES_BELOW_QUALITY 1 -INPUT $OUTPUTDIR/unaligned_tagged_Cell.bam -OUTPUT $OUTPUTDIR/unaligned_tagged_CellMolecular.bam -TMP_DIR $TMPDIR

    FilterBam -TAG_REJECT XQ -INPUT $OUTPUTDIR/unaligned_tagged_CellMolecular.bam -OUTPUT $OUTPUTDIR/unaligned_tagged_filtered.bam -SUMMARY $OUTPUTDIR/unaligned_tagged_filtered.bam_summary.txt -TMP_DIR $TMPDIR

    TrimStartingSequence -OUTPUT_SUMMARY $OUTPUTDIR/adapter_trimming_report.txt -SEQUENCE AAGCAGTGGTATCAACGCAGAGTGAATGGG -MISMATCHES 0 -NUM_BASES 5 -INPUT $OUTPUTDIR/unaligned_tagged_filtered.bam -OUTPUT $OUTPUTDIR/unaligned_tagged_trimmed_smart.bam -TMP_DIR $TMPDIR

    PolyATrimmer -OUTPUT $OUTPUTDIR/unaligned_mc_tagged_polyA_filtered.bam -OUTPUT_SUMMARY $OUTPUTDIR/polyA_trimming_report.txt -MISMATCHES 0 -NUM_BASES 6 -NEW true -INPUT $OUTPUTDIR/unaligned_tagged_trimmed_smart.bam -TMP_DIR $TMPDIR

    picard SamToFastq --INPUT $OUTPUTDIR/unaligned_mc_tagged_polyA_filtered.bam --FASTQ $OUTPUTDIR/unaligned_mc_tagged_polyA_filtered.fastq --TMP_DIR $TMPDIR

    STAR --genomeDir $GENOMEDIR --outFileNamePrefix $OUTPUTDIR/star --readFilesIn $OUTPUTDIR/unaligned_mc_tagged_polyA_filtered.fastq --runThreadN $PBS_NCPUS

    picard SortSam --INPUT $OUTPUTDIR/starAligned.out.sam --OUTPUT $OUTPUTDIR/aligned.sorted.bam --SORT_ORDER queryname --TMP_DIR $TMPDIR

    picard MergeBamAlignment --REFERENCE_SEQUENCE $GENOMEDIR/genome.fa --UNMAPPED_BAM $OUTPUTDIR/unaligned_mc_tagged_polyA_filtered.bam --ALIGNED_BAM $OUTPUTDIR/aligned.sorted.bam \
    --INCLUDE_SECONDARY_ALIGNMENTS false --PAIRED_RUN false --CLIP_ADAPTERS false --OUTPUT $OUTPUTDIR/merged.bam --TMP_DIR $TMPDIR

    TagReadWithGeneFunction -INPUT $OUTPUTDIR/merged.bam -OUTPUT $OUTPUTDIR/function_tagged.bam -ANNOTATIONS_FILE $GENOMEDIR/genome.gtf -TMP_DIR $TMPDIR

    DetectBeadSubstitutionErrors -INPUT $OUTPUTDIR/function_tagged.bam -OUTPUT $OUTPUTDIR/substitution_repaired.bam -MIN_UMIS_PER_CELL 20 -OUTPUT_REPORT $OUTPUTDIR/substitution_error_report.txt -NUM_THREADS $PBS_NCPUS -TMP_DIR $TMPDIR

    DetectBeadSynthesisErrors -INPUT $OUTPUTDIR/substitution_repaired.bam -MIN_UMIS_PER_CELL 20 -OUTPUT_STATS $OUTPUTDIR/synthesis_error_stats.txt -SUMMARY $OUTPUTDIR/synthesis_error_summary.txt \
    -PRIMER_SEQUENCE AAGCAGTGGTATCAACGCAGAGTAC -REPORT $OUTPUTDIR/synthesis_error_report.txt -CREATE_INDEX true -OUTPUT $OUTPUTDIR/final.bam -NUM_THREADS $PBS_NCPUS -TMP_DIR $TMPDIR

    DigitalExpression -INPUT $OUTPUTDIR/final.bam -OUTPUT $OUTPUTDIR/dge_c10k.txt.gz -SUMMARY $OUTPUTDIR/dge_c10k.summary.txt -NUM_CORE_BARCODES 10000 -TMP_DIR $TMPDIR

    DigitalExpression -INPUT $OUTPUTDIR/final.bam -OUTPUT $OUTPUTDIR/dge_t20.txt.gz -SUMMARY $OUTPUTDIR/dge_t20.summary.txt -MIN_NUM_TRANSCRIPTS_PER_CELL 20 -TMP_DIR $TMPDIR

    BamTagHistogram -INPUT $OUTPUTDIR/final.bam -OUTPUT $OUTPUTDIR/cell_readcounts.txt.gz -TAG XC -TMP_DIR $TMPDIR

    GatherMolecularBarcodeDistributionByGene -INPUT $OUTPUTDIR/final.bam -OUTPUT $OUTPUTDIR/gmbd.txt.gz -NUM_CORE_BARCODES 10000 -TMP_DIR $TMPDIR

## Knee-plot analysis

The y axis indicates the “cumulative fraction of uniquely mapped reads.”

``` r
pc27=read.table("/media/nguyen/Data1/mao/scseq/dropseq/pc27mito/cell_readcounts.txt.gz", header=F, stringsAsFactors=F)
csum_pc27=cumsum(pc27$V1)
df_pc27 <- cbind.data.frame(xvalue=1:length(csum_pc27), yvalue=csum_pc27/max(csum_pc27))
ggplot(df_pc27, aes(xvalue, yvalue)) +
  geom_point(size=0.1, color="cornflowerblue") + scale_x_continuous(limits = c(0,50000))+
  #geom_hline(aes(yintercept=df_pc27 %>% filter(xvalue==10000) %>% pull(yvalue)), color="brown", linetype=2)+
  geom_hline(aes(yintercept=df_pc27 %>% filter(xvalue==8000) %>% pull(yvalue)), color="brown", linetype=2)+
  #geom_vline(aes(xintercept=10000), color="orange", linetype=3)+
  geom_vline(aes(xintercept=8000), color="orange", linetype=3)+
  labs(title=expression(italic(P.)~italic(californicus)~"testis (PC_27)"), x="Cell barcodes sorted by number of reads [descending]", y="Cumulative fraction of reads") +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.title = element_text(color="black"),
        axis.text = element_text(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 1, color="black"), aspect.ratio = 1)
```

## DropletUtils analysis

[Nikos Konstantinides](https://konstantinides-lab.com) suggested using
[DropletUtils](https://doi.org/doi:10.18129/B9.bioc.DropletUtils) to
identify the knee and inflection points. The y axis indicates the “total
UMI count for each barcode.”

``` r
for_row_names <- read.table("/media/nguyen/Data1/mao/scseq/dropseq/pc27mito/dge_c10k.txt.gz", header=T, stringsAsFactors=F)$GENE
m_pc27 <- read.table("/media/nguyen/Data1/mao/scseq/dropseq/pc27mito/dge_c10k.txt.gz", header=T, stringsAsFactors=F, row.names = for_row_names)[,-1]
br.out <- barcodeRanks(m_pc27)
o <- order(br.out$rank)
# metadata(br.out)$knee: 273
# metadata(br.out)$inflection : 114
which(br.out$total==273)[1]
min_rank <- br.out$rank[3905] # 7321
which(br.out$total==114)[1]
max_rank <- br.out$rank[7873] # 9401
ggplot()+
  geom_point(aes(x=br.out$rank, y=br.out$total+1), size=0.5, alpha=0.5)+
  geom_line(aes(x=br.out$rank[o],y=br.out$fitted[o]), color="magenta")+
  geom_hline(aes(yintercept=metadata(br.out)$knee), color="dodgerblue", linetype=2)+
  geom_hline(aes(yintercept=metadata(br.out)$inflection), color="brown", linetype=2)+
  geom_vline(aes(xintercept=min_rank), color="orange", linetype=3)+
  geom_vline(aes(xintercept=max_rank), color="orange", linetype=3)+
  scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10')+
  labs(title=expression(italic(P.)~italic(californicus)~"testis (PC_27)"),
       x="Cell barcodes sorted by number of counts [descending]",
       y="Total UMI count for each barcode") +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.title = element_text(color="black"),
        axis.text = element_text(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 1, color="black"), aspect.ratio = 1)
```

## UMAP plot

We used [Seurat](https://satijalab.org/seurat) to generate the UMAP
plot.

``` r
mtx_pc27 <- read.table("/media/nguyen/Data1/mao/scseq/dropseq/pc27mito/dge_c10k.txt.gz", header = TRUE, row.names = 1, colClasses =c("character", rep("numeric", 10000)))
so_pc27 <- CreateSeuratObject(counts = mtx_pc27, min.cells = 3, min.features = 200, project = "pc27") %>%
           PercentageFeatureSet(pattern = "^agat|^rrn", col.name = "percent.mt") %>%
           SCTransform(vars.to.regress = "percent.mt") %>%
           RunPCA() %>%
           FindNeighbors(dims = 1:30) %>%
           RunUMAP(dims = 1:30) %>%
           FindClusters()
df_umap <- so_pc27@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(color=so_pc27@meta.data$seurat_clusters)
my_color <- c(brewer.pal(name="Set2", n=8),brewer.pal(name="Dark2", n=8))[c(1,3,2,4,5:8,12,14,15)]
ggplot(df_umap) +
  geom_point(aes(x=umap_1, y=umap_2, color=color), size=0.8) +
  geom_text_repel(data=df_umap %>% group_by(color) %>% summarise(q1=quantile(umap_1, 0.5), q2=quantile(umap_2, 0.5)),
                  aes(x=q1, y=q2, label = LETTERS[1:11]), size=8) +
  labs(title=expression(italic(P.)~italic(californicus)~"testis (PC_27)"),
       x="UMAP_1",
       y="UMAP_2") +
  #scale_color_brewer(palette = "Set2", name="clusters", labels=LETTERS[1:11]) +
  scale_color_manual(values = my_color, name="clusters", labels=LETTERS[1:11]) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.title = element_text(color="black"),
        axis.text = element_text(color="black"),
        legend.title = element_text(size=10),
        legend.background=element_blank(),
        legend.justification=c(1, 0.85),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 1, color="black"), aspect.ratio = 1)
```

## Session info

    ## R version 4.3.3 (2024-02-29)
    ## Platform: x86_64-conda-linux-gnu (64-bit)
    ## Running under: Ubuntu 22.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /home/nguyen/miniforge-pypy3/envs/seurat510/lib/libopenblasp-r0.3.27.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Europe/Prague
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] RColorBrewer_1.1-3          patchwork_1.2.0            
    ##  [3] DropletUtils_1.22.0         SingleCellExperiment_1.24.0
    ##  [5] SummarizedExperiment_1.32.0 Biobase_2.62.0             
    ##  [7] GenomicRanges_1.54.1        GenomeInfoDb_1.38.8        
    ##  [9] IRanges_2.36.0              S4Vectors_0.40.2           
    ## [11] BiocGenerics_0.48.1         MatrixGenerics_1.14.0      
    ## [13] matrixStats_1.3.0           rmarkdown_2.27             
    ## [15] Seurat_5.1.0                SeuratObject_5.0.2         
    ## [17] sp_2.1-4                    dplyr_1.1.4                
    ## [19] tidyr_1.3.1                 ggrepel_0.9.5              
    ## [21] ggplot2_3.5.1              
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] jsonlite_1.8.8            magrittr_2.0.3           
    ##   [3] spatstat.utils_3.0-5      zlibbioc_1.48.2          
    ##   [5] vctrs_0.6.5               ROCR_1.0-11              
    ##   [7] DelayedMatrixStats_1.24.0 spatstat.explore_3.2-7   
    ##   [9] RCurl_1.98-1.14           S4Arrays_1.2.1           
    ##  [11] htmltools_0.5.8.1         Rhdf5lib_1.24.2          
    ##  [13] rhdf5_2.46.1              SparseArray_1.2.4        
    ##  [15] sctransform_0.4.1         parallelly_1.37.1        
    ##  [17] KernSmooth_2.23-24        htmlwidgets_1.6.4        
    ##  [19] ica_1.0-3                 plyr_1.8.9               
    ##  [21] plotly_4.10.4             zoo_1.8-12               
    ##  [23] igraph_2.0.3              mime_0.12                
    ##  [25] lifecycle_1.0.4           pkgconfig_2.0.3          
    ##  [27] Matrix_1.6-5              R6_2.5.1                 
    ##  [29] fastmap_1.2.0             GenomeInfoDbData_1.2.11  
    ##  [31] fitdistrplus_1.1-11       future_1.33.2            
    ##  [33] shiny_1.8.1.1             digest_0.6.36            
    ##  [35] colorspace_2.1-0          tensor_1.5               
    ##  [37] dqrng_0.4.1               RSpectra_0.16-1          
    ##  [39] irlba_2.3.5.1             beachmat_2.18.1          
    ##  [41] progressr_0.14.0          fansi_1.0.6              
    ##  [43] spatstat.sparse_3.1-0     httr_1.4.7               
    ##  [45] polyclip_1.10-6           abind_1.4-5              
    ##  [47] compiler_4.3.3            withr_3.0.0              
    ##  [49] BiocParallel_1.36.0       fastDummies_1.7.3        
    ##  [51] highr_0.11                R.utils_2.12.3           
    ##  [53] HDF5Array_1.30.1          MASS_7.3-60              
    ##  [55] DelayedArray_0.28.0       tools_4.3.3              
    ##  [57] lmtest_0.9-40             httpuv_1.6.15            
    ##  [59] future.apply_1.11.2       goftest_1.2-3            
    ##  [61] R.oo_1.26.0               glue_1.7.0               
    ##  [63] rhdf5filters_1.14.1       nlme_3.1-165             
    ##  [65] promises_1.3.0            grid_4.3.3               
    ##  [67] Rtsne_0.17                cluster_2.1.6            
    ##  [69] reshape2_1.4.4            generics_0.1.3           
    ##  [71] gtable_0.3.5              spatstat.data_3.1-2      
    ##  [73] R.methodsS3_1.8.2         data.table_1.15.4        
    ##  [75] utf8_1.2.4                XVector_0.42.0           
    ##  [77] spatstat.geom_3.2-9       RcppAnnoy_0.0.22         
    ##  [79] RANN_2.6.1                pillar_1.9.0             
    ##  [81] stringr_1.5.1             limma_3.58.1             
    ##  [83] spam_2.10-0               RcppHNSW_0.6.0           
    ##  [85] later_1.3.2               splines_4.3.3            
    ##  [87] lattice_0.22-6            survival_3.7-0           
    ##  [89] deldir_2.0-4              tidyselect_1.2.1         
    ##  [91] locfit_1.5-9.10           scuttle_1.12.0           
    ##  [93] miniUI_0.1.1.1            pbapply_1.7-2            
    ##  [95] knitr_1.48                gridExtra_2.3            
    ##  [97] edgeR_4.0.16              scattermore_1.2          
    ##  [99] xfun_0.45                 statmod_1.5.0            
    ## [101] stringi_1.8.4             lazyeval_0.2.2           
    ## [103] yaml_2.3.9                evaluate_0.24.0          
    ## [105] codetools_0.2-20          tibble_3.2.1             
    ## [107] cli_3.6.3                 uwot_0.2.2               
    ## [109] xtable_1.8-4              reticulate_1.38.0        
    ## [111] munsell_0.5.1             Rcpp_1.0.12              
    ## [113] globals_0.16.3            spatstat.random_3.2-3    
    ## [115] png_0.1-8                 parallel_4.3.3           
    ## [117] dotCall64_1.1-1           sparseMatrixStats_1.14.0 
    ## [119] bitops_1.0-7              listenv_0.9.1            
    ## [121] viridisLite_0.4.2         scales_1.3.0             
    ## [123] ggridges_0.5.6            crayon_1.5.3             
    ## [125] leiden_0.4.3.1            purrr_1.0.2              
    ## [127] rlang_1.1.4               cowplot_1.1.3

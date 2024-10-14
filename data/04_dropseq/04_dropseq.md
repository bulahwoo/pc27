Single-cell RNA sequencing of *Phymatopus californicus* sample PC27 with
Drop-seq core computational protocol
================
Andrea Elizabeth Acurio Armas, Bulah Wu, Petr Nguyen  
October 14, 2024

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
[mitos2](https://usegalaxy.org.au/root?tool_id=toolshed.g2.bx.psu.edu/repos/iuc/mitos2/mitos2/2.1.9+galaxy0)
on the Galaxy platform.

## Library preparation

The scRNA-seq library was prepared with the RNAdia 2.0 kit (Dolomite
Bio), following the manufacturer’s instructions. Approximately 8,000
cells were captured, and the size distribution of the library was
examined using the Bioanalyzer system.

<div class="figure" style="text-align: center">

<img src="../shared/misc/bioanalyzer_pc27.png" alt="Size distribution assessed by Bioanalyzer" width="60%" />
<p class="caption">
Size distribution assessed by Bioanalyzer
</p>

</div>

## Quality control checks on raw reads

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) was
used for quality check. Output can be found
[here](../shared/fastqc/pc27_read1/fastqc.md) and
[here](../shared/fastqc/pc27_read2/fastqc.md).

## Drop-seq protocol

The digital expression matrix using the Drop-seq core computational
protocol developed by the McCarroll Lab
(<https://github.com/broadinstitute/Drop-seq>). Briefly, cell barcodes
(CBs) and unique molecular identifiers (UMIs) were identified from raw
reads. Reads with CBs/UMIs of low quality bases were removed, and 5’
primer and 3’ polyA sequences were trimmed. The reads were then aligned
to the reference genome using STAR aligner. Substitution and indel
errors in barcodes were repaired before the digital expression matrix
was created.

    #!/bin/bash
    #PBS -N dropseq_original
    #PBS -l select=1:ncpus=16:mem=140gb:scratch_ssd=500gb
    #PBS -l walltime=48:00:00
    #PBS -m ae

    #RANDOM=1
    #export SEEDID=1
    #export PERCENT=$(echo "scale=2; ${SEEDID}/100" | bc -l)
    export ID=dropseq_original
    DATADIR=/storage/plzen1/home/bulah/mao/scseq/select_polyt/dropseq_original
    echo "||| $HOSTNAME $SCRATCHDIR ||| (${ID}) ||| ssh $HOSTNAME \"scp -rp bulah@$HOSTNAME:$SCRATCHDIR/* \$HOSTNAME:\$SCRATCHDIR\" |||" >> $DATADIR/jobid

    module add mambaforge picard star samtools
    mamba activate /storage/plzen1/home/bulah/.conda/envs/openjdk2201
    ssh kirke59.meta.zcu.cz "scp -rp bulah@kirke59.meta.zcu.cz:/scratch.ssd/bulah/job_4671652.pbs-m1.metacentrum.cz/* $HOSTNAME:$SCRATCHDIR"

    #export ASSEMBLY=/storage/plzen1/home/bulah/mao/braker3/phymatopus_californicus/genome_mao/repeatmodeler/phycal_genome_masked2.fa
    #export ANNOTATION=/storage/plzen1/home/bulah/mao/braker3/phymatopus_californicus/genome_mao/braker3/edit02_phycal_mito.gtf
    #export READ1=/storage/plzen1/home/bulah/mao/scseq/X201SC24042690-Z01-F001/X201SC24042690-Z01-F001/01.RawData/PC_27/PC_27_1.fq.gz
    #export READ2=/storage/plzen1/home/bulah/mao/scseq/X201SC24042690-Z01-F001/X201SC24042690-Z01-F001/01.RawData/PC_27/PC_27_2.fq.gz

    mkdir -p $SCRATCHDIR/{genome,raw,tmp}
    export GENOMEDIR=$SCRATCHDIR/genome
    export TMPDIR=$SCRATCHDIR/tmp
    export RAWDIR=$SCRATCHDIR/raw
    export PATH=$PATH:$SCRATCHDIR/dropseq-3.0.1

    cd $SCRATCHDIR

    if [ ! -d "$SCRATCHDIR"/dropseq-3.0.1 ]; then
      wget https://github.com/broadinstitute/Drop-seq/releases/download/v3.0.1/dropseq-3.0.1.zip
      unzip dropseq-3.0.1.zip
      rm -rf dropseq-3.0.1.zip
    fi

    if [ ! -f "$GENOMEDIR"/genome.fa ]; then
      cp -rp $ASSEMBLY $GENOMEDIR/genome.fa
      cp -rp $ANNOTATION $GENOMEDIR/genome.gtf
      STAR --runMode genomeGenerate --genomeDir $GENOMEDIR --genomeFastaFiles $GENOMEDIR/genome.fa --sjdbGTFfile $GENOMEDIR/genome.gtf --genomeSAindexNbases 13 --runThreadN $PBS_NCPUS
    fi

    if [ ! -f "$GENOMEDIR"/genome.dict ]; then
      picard CreateSequenceDictionary -R $GENOMEDIR/genome.fa --TMP_DIR $TMPDIR
    fi

    #cp -rp $READ1 $RAWDIR/original_read1.fq.gz
    #cp -rp $READ2 $RAWDIR/original_read2.fq.gz

    #picard FastqToSam -F1 $RAWDIR/original_read1.fq.gz -F2 $RAWDIR/original_read2.fq.gz -O $RAWDIR/unaligned.bam -SM pc27 --TMP_DIR $TMPDIR
    #rm -rf $RAWDIR/original_read1.fq.gz $RAWDIR/original_read2.fq.gz

    for ((i=1; i<=1; i++)); do

      mkdir -p $SCRATCHDIR/$ID
      export OUTPUTDIR=$SCRATCHDIR/$ID

    #  SEED=$RANDOM
    #  samtools view -h --subsample-seed $SEED --subsample $PERCENT $RAWDIR/unaligned.bam > $OUTPUTDIR/unaligned.bam

    #  md5sum $OUTPUTDIR/unaligned.bam > $OUTPUTDIR/subread_md5.txt

      ln -s $RAWDIR/unaligned.bam $OUTPUTDIR/unaligned.bam

      TagBamWithReadSequenceExtended -SUMMARY $OUTPUTDIR/unaligned_tagged_Cellular.bam_summary.txt -BASE_RANGE 3-14 -BASE_QUALITY 10 -BARCODED_READ 1 -DISCARD_READ false -TAG_NAME XC \
    -NUM_BASES_BELOW_QUALITY 1 -INPUT $OUTPUTDIR/unaligned.bam -OUTPUT $OUTPUTDIR/unaligned_tagged_Cell.bam -TMP_DIR $TMPDIR

      rm -rf $OUTPUTDIR/unaligned.bam

      TagBamWithReadSequenceExtended -SUMMARY $OUTPUTDIR/unaligned_tagged_Molecular.bam_summary.txt -BASE_RANGE 15-28 -BASE_QUALITY 10 -BARCODED_READ 1 -DISCARD_READ true -TAG_NAME XM \
    -NUM_BASES_BELOW_QUALITY 1 -INPUT $OUTPUTDIR/unaligned_tagged_Cell.bam -OUTPUT $OUTPUTDIR/unaligned_tagged_CellMolecular.bam -TMP_DIR $TMPDIR

      rm -rf $OUTPUTDIR/unaligned_tagged_Cell.bam

      FilterBam -TAG_REJECT XQ -INPUT $OUTPUTDIR/unaligned_tagged_CellMolecular.bam -OUTPUT $OUTPUTDIR/unaligned_tagged_filtered.bam -SUMMARY $OUTPUTDIR/unaligned_tagged_filtered.bam_summary.txt -TMP_DIR $TMPDIR

      rm -rf $OUTPUTDIR/unaligned_tagged_CellMolecular.bam

      TrimStartingSequence -OUTPUT_SUMMARY $OUTPUTDIR/adapter_trimming_report.txt -SEQUENCE AAGCAGTGGTATCAACGCAGAGTGAATGGG -MISMATCHES 0 -NUM_BASES 5 -INPUT $OUTPUTDIR/unaligned_tagged_filtered.bam -OUTPUT $OUTPUTDIR/unaligned_tagged_trimmed_smart.bam -TMP_DIR $TMPDIR

      rm -rf $OUTPUTDIR/unaligned_tagged_filtered.bam

      PolyATrimmer -OUTPUT $OUTPUTDIR/unaligned_mc_tagged_polyA_filtered.bam -OUTPUT_SUMMARY $OUTPUTDIR/polyA_trimming_report.txt -MISMATCHES 0 -NUM_BASES 6 -NEW true -INPUT $OUTPUTDIR/unaligned_tagged_trimmed_smart.bam -TMP_DIR $TMPDIR

      rm -rf $OUTPUTDIR/unaligned_tagged_trimmed_smart.bam

      picard SamToFastq --INPUT $OUTPUTDIR/unaligned_mc_tagged_polyA_filtered.bam --FASTQ $OUTPUTDIR/unaligned_mc_tagged_polyA_filtered.fastq --TMP_DIR $TMPDIR

      STAR --genomeDir $GENOMEDIR --outFileNamePrefix $OUTPUTDIR/star --readFilesIn $OUTPUTDIR/unaligned_mc_tagged_polyA_filtered.fastq --runThreadN $PBS_NCPUS

      rm -rf $OUTPUTDIR/unaligned_mc_tagged_polyA_filtered.fastq

      picard SortSam --INPUT $OUTPUTDIR/starAligned.out.sam --OUTPUT $OUTPUTDIR/aligned.sorted.bam --SORT_ORDER queryname --TMP_DIR $TMPDIR

      rm -rf $OUTPUTDIR/starAligned.out.sam

      picard MergeBamAlignment --REFERENCE_SEQUENCE $GENOMEDIR/genome.fa --UNMAPPED_BAM $OUTPUTDIR/unaligned_mc_tagged_polyA_filtered.bam --ALIGNED_BAM $OUTPUTDIR/aligned.sorted.bam \
    --INCLUDE_SECONDARY_ALIGNMENTS false --PAIRED_RUN false --CLIP_ADAPTERS false --OUTPUT $OUTPUTDIR/merged.bam --TMP_DIR $TMPDIR

      rm -rf $OUTPUTDIR/unaligned_mc_tagged_polyA_filtered.bam $OUTPUTDIR/aligned.sorted.bam

      TagReadWithGeneFunction -INPUT $OUTPUTDIR/merged.bam -OUTPUT $OUTPUTDIR/function_tagged.bam -ANNOTATIONS_FILE $GENOMEDIR/genome.gtf -TMP_DIR $TMPDIR

      rm -rf $OUTPUTDIR/merged.bam

      DetectBeadSubstitutionErrors -INPUT $OUTPUTDIR/function_tagged.bam -OUTPUT $OUTPUTDIR/substitution_repaired.bam -MIN_UMIS_PER_CELL 20 -OUTPUT_REPORT $OUTPUTDIR/substitution_error_report.txt -NUM_THREADS $PBS_NCPUS -TMP_DIR $TMPDIR

      rm -rf $OUTPUTDIR/function_tagged.bam

      DetectBeadSynthesisErrors -INPUT $OUTPUTDIR/substitution_repaired.bam -MIN_UMIS_PER_CELL 20 -OUTPUT_STATS $OUTPUTDIR/synthesis_error_stats.txt -SUMMARY $OUTPUTDIR/synthesis_error_summary.txt \
    -PRIMER_SEQUENCE AAGCAGTGGTATCAACGCAGAGTAC -REPORT $OUTPUTDIR/synthesis_error_report.txt -CREATE_INDEX true -OUTPUT $OUTPUTDIR/final.bam -NUM_THREADS $PBS_NCPUS -TMP_DIR $TMPDIR

      rm -rf $OUTPUTDIR/substitution_repaired.bam

      DigitalExpression -INPUT $OUTPUTDIR/final.bam -OUTPUT $OUTPUTDIR/dge_t20.txt.gz -SUMMARY $OUTPUTDIR/dge_t20.summary.txt -MIN_NUM_TRANSCRIPTS_PER_CELL 20 -TMP_DIR $TMPDIR

      BamTagHistogram -INPUT $OUTPUTDIR/final.bam -OUTPUT $OUTPUTDIR/cell_readcounts.txt.gz -TAG XC -TMP_DIR $TMPDIR

      GatherMolecularBarcodeDistributionByGene -INPUT $OUTPUTDIR/final.bam -OUTPUT $OUTPUTDIR/mb_dist_by_gene_t20.txt.gz -MIN_NUM_TRANSCRIPTS_PER_CELL 20 -TMP_DIR $TMPDIR

    done

## Knee-plot analysis

In the dot plot, the x-axis represents the cell barcodes (organized by
the number of reads, arranged from highest to lowest), and the y-axis
shows the cumulative fraction of uniquely mapped reads. The transition
from beads sampling cellular RNA to beads sampling ambient RNA is marked
by the inflection point ([Macosko et al.,
2015](http://dx.doi.org/10.1016/j.cell.2015.05.002)).

``` r
pc27=read.table("/media/nguyen/Data1/mao/scseq/select_polyt/dropseq_original/cell_readcounts.txt.gz", header=F, stringsAsFactors=F)
csum_pc27=cumsum(pc27$V1)
df_pc27 <- cbind.data.frame(xvalue=1:length(csum_pc27), yvalue=csum_pc27/max(csum_pc27))
ggplot(df_pc27, aes(xvalue, yvalue)) +
  geom_point(size=0.1, color="cornflowerblue") + scale_x_continuous(limits = c(0,100000))+
  #geom_hline(aes(yintercept=df_pc27 %>% filter(xvalue==8000) %>% pull(yvalue)), color="brown", linetype=2)+
  #geom_vline(aes(xintercept=8000), color="orange", linetype=3)+
  #annotate("text", x=15000, y=0.25, label="(8000, 0.2998)")+ # 0.2998019
  labs(title=expression(italic(P.)~italic(californicus)~"PC27"), x="Cell barcodes sorted by number of reads [descending]", y="Cumulative fraction of reads") +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.title = element_text(color="black"),
        axis.text = element_text(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 1, color="black"), aspect.ratio = 1)
```

![](04_dropseq_files/figure-gfm/knee-plot-1.png)<!-- -->

The plot indicates the top 8000 cells contribute nearly 30% of the total
uniquely mapped reads.

## DropletUtils analysis

We used
[DropletUtils](https://doi.org/doi:10.18129/B9.bioc.DropletUtils) to
identify the knee and inflection points. Here the x-axis indicates the
cell barcodes (organized by the number of reads, arranged from highest
to lowest) and the y-axis the total UMI count for each barcode.

``` r
dge <- read.table("/media/nguyen/Data1/mao/scseq/select_polyt/dropseq_original/dge_t20.txt.gz", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
br.out <- barcodeRanks(dge)
o <- order(br.out$rank)
metadata(br.out) # knee and inflection are values of br.out$total (y-axis)
```

    ## $knee
    ## [1] 1471
    ## 
    ## $inflection
    ## [1] 759

``` r
knee <- metadata(br.out)$knee %>% round()
inflection <- metadata(br.out)$inflection %>% round()
# which br.out$rank value (x-axis) corresponds to knee value (y-axis)
knee_rank <- br.out %>% as_tibble() %>% filter(total==knee) %>% arrange(rank) %>% dplyr::slice(1) %>% pull(rank)
# which br.out$rank value (x-axis) corresponds to inflection value (y-axis)
inflection_rank <- br.out %>% as_tibble() %>% filter(total==inflection) %>% arrange(rank) %>% dplyr::slice(1) %>% pull(rank)
ggplot()+
  geom_point(aes(x=br.out$rank, y=br.out$total+1), color="grey50", size=0.5, alpha=0.5)+
  geom_line(aes(x=br.out$rank[o],y=br.out$fitted[o]), color="magenta")+
  geom_hline(aes(yintercept=knee), color="dodgerblue", linetype=2)+
  geom_hline(aes(yintercept=inflection), color="brown", linetype=2)+
  geom_vline(aes(xintercept=knee_rank), color="orange", linetype=3)+
  geom_vline(aes(xintercept=inflection_rank), color="orange", linetype=3)+
  annotate("text", x=11, y=1100, label=paste0("(", knee_rank, ", ", knee, ")"))+
  annotate("text", x=1100, y=110, label=paste0("(", inflection_rank, ", ", inflection, ")"))+
  #scale_x_continuous(trans='log10', limits=c(1, 10^(floor(log10(inflection_rank)) + 2)), breaks=c(10^(1:(floor(log10(inflection_rank)) + 2))), labels = scales::number)+
  #scale_y_continuous(trans='log10', limits=c(1, 10^(floor(log10(knee)) + 2)), breaks=c(10^(1:(floor(log10(knee)) + 2))), labels = scales::number)+
  scale_x_continuous(limits=c(1, 50000), breaks=c(1,10,100,1000,10000, 40000), trans='log10', labels = scales::number)+
  scale_y_continuous(limits=c(1, 5000), breaks=c(1, 10, 50, 100, 500, 1000, 4000), trans='log10', labels = scales::number)+
  labs(title=ggp_title,
       x="Cell barcodes sorted by number of counts [descending]",
       y="Total UMI count for each barcode") +
  ggp_theme_bw_square_01
```

![](04_dropseq_files/figure-gfm/barcodeRanks-1.png)<!-- -->

We used DropletUtils::emptyDrops() to identify the empty droplets.

``` r
set.seed(100)
e.out <- emptyDrops(dge, niters=50000)
is.cell <- e.out$FDR <= 0.01
sum(is.cell, na.rm=TRUE)
```

    ## [1] 1587

In total 5919 cells are retained.

``` r
table(Limited=e.out$Limited, Significant=is.cell)
```

    ##        Significant
    ## Limited FALSE  TRUE
    ##   FALSE 34135   639
    ##   TRUE      0   948

The zero in the above table indicates no entry with false positives
frequency above the threshold 0.001 (Significant FALSE) can be achieved
by increasing the number of permutations (Limited TRUE). So the
parameter *niters=50000* is good.

``` r
dge <- dge[,which(e.out$FDR<=0.01)]
ggplot(e.out %>% as.data.frame() %>% filter(!is.na(LogProb)) %>% mutate(FDR_fill=ifelse(FDR>0.01, "cornflowerblue", "orange")))+
  geom_point(aes(x=Total, y=-LogProb, color=FDR_fill), size=1, alpha=0.5)+
  scale_color_manual(values=c("cornflowerblue", "orange"), labels=c('FDR > 0.01', 'FDR <= 0.01'))+
  labs(title=ggp_title,
       x="Total UMI count",
       y="-Log Probability") +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.title = element_text(color="black"),
        axis.text = element_text(color="black"),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 1, color="black"), aspect.ratio = 1)
```

![](04_dropseq_files/figure-gfm/emptyDrops3-1.png)<!-- -->

The retained droplets (orange color) should have large UMI counts and/or
large negative log probabilities.

## Analysis using Seurat

[Seurat v5](https://satijalab.org/seurat) was used to analyze the data.

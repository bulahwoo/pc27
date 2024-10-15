pc27 analysis using STARsolo raw matrix
================
Andrea Elizabeth Acurio Armas, Bulah Wu, Petr Nguyen  
October 15, 2024

### FastQC

Output can be found [here (read
1)](../shared/fastqc/pc27_read1/fastqc.md) and [here (read
2)](../shared/fastqc/pc27_read2/fastqc.md).

 

### barcodeRanks()

<div class="figure" style="text-align: center">

<img src="./files/barcoderanks.png" alt="barcodeRanks() output" width="60%" />
<p class="caption">
barcodeRanks() output
</p>

</div>

 

### emptyDrops()

|       | FALSE | TRUE |
|:------|------:|-----:|
| FALSE | 36941 |  318 |
| TRUE  |     0 |  635 |

953 cells are identified.

 

### vlnplot()

- The raw matrix
  <div class="figure" style="text-align: center">

  <img src="./files/plt_vln_pre_emptydrops.png" alt="pre-emptydrops" width="60%" />
  <p class="caption">
  pre-emptydrops
  </p>

  </div>

|  | Gene | Cell | Mean UMI/Cell | Median UMI/Cell | Mean Gene/Cell | Median Gene/Cell |
|:---|---:|---:|---:|---:|---:|---:|
| Raw | 11320 | 1277682 | 10.53559 | 1 | 7.912049 | 1 |

 

- After emptyDrops()
  <div class="figure" style="text-align: center">

  <img src="./files/plt_vln_post_emptydrops.png" alt="post-emptydrops" width="60%" />
  <p class="caption">
  post-emptydrops
  </p>

  </div>

|  | Gene | Cell | Mean UMI/Cell | Median UMI/Cell | Mean Gene/Cell | Median Gene/Cell |
|:---|---:|---:|---:|---:|---:|---:|
| emptyDrops | 11320 | 953 | 1287.112 | 920 | 659.4848 | 595 |

 

- STEP 1: filter genes detected in \< 3 cells
  <div class="figure" style="text-align: center">

  <img src="./files/plt_vln_filter01.png" alt="seurat filter step 1" width="60%" />
  <p class="caption">
  seurat filter step 1
  </p>

  </div>

|       | Gene | Cell | Mean UMI/Cell | Median UMI/Cell | Mean Gene/Cell | Median Gene/Cell |
|:------|-----:|-----:|--------------:|----------------:|---------------:|-----------------:|
| Step1 | 6920 |  953 |      1285.474 |             919 |       657.8835 |              594 |

 

- STEP 2: filter cells that contain \< 200 genes detected
  <div class="figure" style="text-align: center">

  <img src="./files/plt_vln_filter02.png" alt="seurat filter step 2" width="60%" />
  <p class="caption">
  seurat filter step 2
  </p>

  </div>

|       | Gene | Cell | Mean UMI/Cell | Median UMI/Cell | Mean Gene/Cell | Median Gene/Cell |
|:------|-----:|-----:|--------------:|----------------:|---------------:|-----------------:|
| Step2 | 6875 |  839 |      1436.646 |            1123 |       727.5542 |              688 |

 

- STEP 3.1: following STEP 2, filter cells that contain ≤ 200 genes or ≥
  2500 genes detected
  <div class="figure" style="text-align: center">

  <img src="./files/plt_vln_filter03.png" alt="seurat filter step 3.1" width="60%" />
  <p class="caption">
  seurat filter step 3.1
  </p>

  </div>

|  | Gene | Cell | Mean UMI/Cell | Median UMI/Cell | Mean Gene/Cell | Median Gene/Cell |
|:---|---:|---:|---:|---:|---:|---:|
| Step3.1 | 6875 | 836 | 1440.894 | 1128.5 | 729.4533 | 690 |

 

- STEP 3.2: following STEP 2, filter cells that contain ≥ 5%
  mitochondrial counts
  <div class="figure" style="text-align: center">

  <img src="./files/plt_vln_filter04.png" alt="seurat filter step 3.2" width="60%" />
  <p class="caption">
  seurat filter step 3.2
  </p>

  </div>

|  | Gene | Cell | Mean UMI/Cell | Median UMI/Cell | Mean Gene/Cell | Median Gene/Cell |
|:---|---:|---:|---:|---:|---:|---:|
| Step3.2 | 6875 | 613 | 1588.168 | 1385 | 789.0261 | 799 |

 

- STEP 3: following STEP 2, filter cells that contain ≤ 200 genes or ≥
  2500 genes detected, and filter cells that contain ≥ 5% mitochondrial
  counts
  <div class="figure" style="text-align: center">

  <img src="./files/plt_vln_filter.png" alt="seurat filter step 3" width="60%" />
  <p class="caption">
  seurat filter step 3
  </p>

  </div>

|       | Gene | Cell | Mean UMI/Cell | Median UMI/Cell | Mean Gene/Cell | Median Gene/Cell |
|:------|-----:|-----:|--------------:|----------------:|---------------:|-----------------:|
| Step3 | 6875 |  612 |      1590.408 |          1387.5 |       789.9902 |            799.5 |

 

- Summary

|  | Gene | Cell | Mean UMI/Cell | Median UMI/Cell | Mean Gene/Cell | Median Gene/Cell |
|:---|---:|---:|---:|---:|---:|---:|
| Raw | 11320 | 1277682 | 10.53559 | 1.0 | 7.912049 | 1.0 |
| emptyDrops | 11320 | 953 | 1287.11228 | 920.0 | 659.484785 | 595.0 |
| Step1 | 6920 | 953 | 1285.47429 | 919.0 | 657.883526 | 594.0 |
| Step2 | 6875 | 839 | 1436.64601 | 1123.0 | 727.554231 | 688.0 |
| Step3.1 | 6875 | 836 | 1440.89354 | 1128.5 | 729.453349 | 690.0 |
| Step3.2 | 6875 | 613 | 1588.16803 | 1385.0 | 789.026101 | 799.0 |
| Step3 | 6875 | 612 | 1590.40850 | 1387.5 | 789.990196 | 799.5 |

 

### UMAP

<div class="figure" style="text-align: center">

<img src="./files/umap.png" alt="umap" width="60%" />
<p class="caption">
umap
</p>

</div>

 

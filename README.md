# dREG.HD

Refining TRE regions from dREG by imputing DHS.

## Abstract:


## Installation Instructions:

dREG.HD is only available for the Linux and Mac OSX. The source code can be downloaded from this repository (https://github.com/Danko-Lab/dREG.HD.git) . 

### Required software and packages
    
1. R (http://www.r-project.org/)
    
2. dREG package (https://github.com/Danko-Lab/dREG).
    
3. bedtools (https://github.com/arq5x/bedtools2/)
    
4. bedGraphToBigWig command in UCSC Blat application binaries (http://hgdownload.cse.ucsc.edu/admin/exe/)
    
5. Extra R Package: snowfall, data.table.
    
### Install dREG.HD

Please install the required R package before you install dREG.HD package. After the  installation of `dREG`, `snowfall` and `data.table` package, please install the #dREG.HD# as following steps.

```
git clone https://github.com/Danko-Lab/dREG.HD.git

cd dREG.HD

R CMD INSTALL dREG.HD

```

##Usage instructions##

dREG.HD uses 4 files as input, and outputs one file. Input files include the PRO-seq read distributions on the plus and minus strand (which are separate files), a bed file estimated by dREG model and its SVR model pre-trained by dREG package.

>PRO-seq files are required to be in the bigWig format standard created by the UCSC (more information can be found here: http://genome.ucsc.edu/goldenPath/help/bigWig.html).

>The SVR model is included in dREG package (under dREG_model/asvm.RData). Users are advised to use that when possible.

To use dREG.HD, type: 

```
bash run_dREG-HD.bsh dREG_bed plus_strand.bw minus_strand.bw HD_model.rdata chromInfo [nthreads] [GPU]

dREG_bed        -- the dREG peaks in bed format.
plus_strand.bw	-- PRO-seq data (plus  strand) formatted as a bigWig file.
minus_strand.bw	-- PRO-seq data (minus strand) formatted as a bigWig file.
HD_model.rdata	-- The path to the RData file containing the pre-trained dREG-HD SVM model.
chromInfo       -- the chromInfo file required for generating bigwig file of imputed DNase-I signal.
[nthreads]	    -- [optional, default=1] The number of threads to use when evaluating dREG-HD sites.
[GPU]	        -- [optional, gpu or blank, default is blank] indicating whether GPU package can be used.

```

dREG.HD needs the package `Rgtsvm` to run on GPU. This SVM package on GPU platform can be downloaded from (https://github.com/Danko-Lab/Rgtsvm.git)

dREG.HD is an R package, and that provides some additional flexibility for users familiar with R. Please check the example code in the package. 









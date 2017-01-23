# dREG.HD

Refining TRE regions from dREG by imputing DHS.

## Abstract:

The profile of nascent transcription from precision run-on and sequencing (PRO-seq) contains precise information on the location of active transcriptional regulatory elements (TREs). The dREG package uses support vector regression (SVR) to identify active TREs solely from PRO-seq data. Here, we introduce its successor dREG.HD, which improves resolution by imputing smoothed DNase-I hypersensitivity signal and calling putative DNase-I hypersensitive sites (DHS). Briefly, we used an epsilon-support vector regression (SVR) with a Gaussian kernel to map the distribution of PRO-seq reads to smoothed DNase-I signal intensities. Training was conducted on randomly chosen positions within dREG peaks extended by 200bp on either side. Selection of feature vectors was optimized based on Pearson correlation coefficients between the imputed and experimental DNase-I score over the validation set.  PRO-seq data was normalized by sequencing depth and further scaled such that the maximum value of any prediction dataset is within 90 percentile of the training examples.  We chose a step size to be 60bp and extending 30 steps on each direction. The final model was trained using matched DNase-I and PRO-seq data in K562 cells.  

Next we identified peaks in the imputed DNase-I hypersensitivity profile by fitting the imputed DNase-I signal using a cubic spline and identifying local maxima.  We optimized two free parameters that control the (1) smoothness of spline curve fitting, and (2) threshold on the imputed DNase-I signal intensity.  Parameters were optimized to achieve an appropriate trade-off between FDR and sensitivity on the testing K562 dataset. Parameters were tuned using a grid optimization over free parameters. Testing the optimized dREG.HD (including both DNase-I imputation and peak calling) on GM12878, a GRO-seq dataset completely held out from model training and parameter optimization, revealed 82% sensitivity for DNase-I peaks within dREG sites at a 10% false discovery rate (FDR). dREG.HD will output peaks called under both the  relaxed condition (FDR=16%) and stringent condition (FDR=10%).


## Installation Instructions:

dREG.HD is only available for the Linux and Mac OSX. The source code can be downloaded from this repository (https://github.com/Danko-Lab/dREG.HD.git) . 

### Required software and packages
    
1. R (http://www.r-project.org/)
    
2. dREG package (https://github.com/Danko-Lab/dREG).
    
3. bedtools (https://github.com/arq5x/bedtools2/)
    
4. bedGraphToBigWig command in UCSC Blat application binaries (http://hgdownload.cse.ucsc.edu/admin/exe/)
    
5. Extra R Package: snowfall, data.table.
    
### Install dREG.HD

Please install the required R package before you install dREG.HD package. After the  installation of `dREG`, `snowfall` and `data.table` package, please install the **dREG.HD** as following steps.

```
git clone https://github.com/Danko-Lab/dREG.HD.git

cd dREG.HD

R CMD INSTALL dREG.HD

```

##Usage instructions##

dREG.HD uses 4 files as input, and outputs one file. Input files include the PRO-seq read distributions on the plus and minus strand (which are separate files), a bed file estimated by dREG model and its SVR model pre-trained by dREG package.

>PRO-seq files are required to be in the bigWig format standard created by the UCSC (more information can be found here: http://genome.ucsc.edu/goldenPath/help/bigWig.html).

>The SVR model is included in dREG package (under inst/extdata/dREG_HD.model.rdata). Users are advised to use that when possible.


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

An example is 

bash ./run_dREG-HD.bsh ./dREG.HD/inst/extdata/k562.chr21.predictions.bed ./dREG.HD/inst/extdata/K562.chr21.plus.bw ./dREG.HD/inst/extdata/K562.chr21.minus.bw ./dREG.HD/inst/extdata/dREG_HD.model.rdata ./dREG.HD/inst/extdata/chromInfo.hg19 14 GPU


dREG.HD needs the package `Rgtsvm` to run on GPU. This SVM package on GPU platform can be downloaded from (https://github.com/Danko-Lab/Rgtsvm.git)

dREG.HD is an R package, and that provides some additional flexibility for users familiar with R. Please check the example code in the package. 


>dREG.HD will output three files under the same directory of the input files, with the prefix being "dREG_bed". The details of these files are listed below.

```

${prefix}_imputedDnase.bw	      -- the imputed DNase-I signal.
${prefix}_dREG_HD_relaxed.bed	  -- dREG.HD peaks called under relaxed condition (FDR=16%)
${prefix}_dREG_HD_stringent.bed   -- dREG.HD peaks called under stringent condition (FDR=10%)

```






require(dREG.HD)

## Process command arguments
args <- commandArgs(trailingOnly=TRUE)

## Load the model.  Do this before loading ps_plus_path, just in case those are saved in the model file.
model_path <- args[4]
#load(dreg_model) ## Should have (by default) gdm and asvm.

## Read PRO-seq data.
bed_path <- args[1]
ps_plus_path  <- args[2]
ps_minus_path <- args[3]
ncores <- as.integer(args[5])
if (is.na(ncores)) ncores <- 1;

use_rgtsvm <- FALSE;
use_gpu <- toupper(as.character(args[6]))
if (!is.na(use_gpu) && use_gpu=="GPU") use_rgtsvm <- TRUE;


if (use_rgtsvm)
{
   if(!requireNamespace("Rgtsvm"))
   stop("Rgtsvm has not been installed fotr GPU computing.");
}

load(model_path);

## Now running dREG-HD
t <- system.time( dREG_HD(bed_path= bed_path, bigwig_plus = ps_plus_path, bigwig_minus = ps_minus_path, model= model, ncores = ncores, use_rgtsvm= use_rgtsvm))

cat("Running time [User]:", t[1], "[System]:", t[2], "[Elapsed]:", t[3], "\n");



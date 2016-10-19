\name{dREG_HD}
\alias{dREG_HD}
\title{Run dREG-HD}
\description{Run dREG-HD}
\usage{
dREG_HD(bed_path, 
      bigwig_plus, 
      bigwig_minus, 
      chromInfo, 
      model,
      ncores=1, 
      use_rgtsvm = FALSE )
}
\arguments{
  \item{bed_path}{path to the dREG peaks in bed format}
  \item{bigwig_plus}{the path the the plus strand of bigwig file}
  \item{bigwig_minus}{the path the the minus strand of bigwig file}
  \item{chromInfo}{the file specifies chromosome sizes}
  \item{model}{the dREG-HD model trained by dREG-HD}
  \item{ncores}{Number of CPU cores in parallel computing}
  \item{use_rgtsvm}{use GPU or not}
}

\examples{

bed_path <-system.file("extdata","k562.chr21.predictions.bed", package="dREG.HD")
bigwig_plus<-system.file("extdata","K562.chr21.plus.bw", package="dREG.HD")
bigwig_minus<-system.file("extdata","K562.chr21.minus.bw", package="dREG.HD")
chromInfo <-system.file("extdata","chromInfo.hg19", package="dREG.HD")
model_path <-system.file("extdata","dREG_HD.model.rdata", package="dREG.HD")
load(model_path)

dREG_HD(bed_path= bed_path, bigwig_plus = bigwig_plus, bigwig_minus =bigwig_minus,chromInfo= chromInfo, model=model)



}
\keyword{ function }

#utility functions

#get the chromosome length
get.chrom.length<-function(chr, chrom.info.table){
	return(chrom.info.table[(chrom.info.table[,1]==chr),2])
}

#extract to construct bedgraph and then to bigwig
tobigwig<-function(filename, temp.bg, chromInfo){
	#get bedgraph

	bedgraph.sorted=tempfile()

	options(scipen =99) # not to use scientific notation when writing out

	#write bedgraph formatted dataframes to tempfile
  	#write.table(bedgraph,file= bedgraph.file,quote=F,sep="\t",col.names=F,row.names=F)
  	command=paste("LC_COLLATE=C sort -k1,1 -k2,2n", temp.bg, "| uniq >", bedgraph.sorted,sep=" ")
  	#browser()
  	 try(system(command))

  	command=paste("bedGraphToBigWig", bedgraph.sorted, chromInfo ,filename, sep=" ")
  	#cat(command,"\n")
    try(system(command))
    #unlink(bedgraph.file)
    unlink(bedgraph.sorted)
}

#betools merge
bedTools.merge<-function(bed)
{
  #create temp files
  a.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out

  #write bed formatted dataframes to tempfile
  write.table(bed,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)

  # create the command string and call the command using system()
  command=paste("LC_COLLATE=C sort -k1,1 -k2,2n",a.file,"| mergeBed -i stdin >",out,sep=" ")
  #cat(command,"\n")
  try(system(command))

  res=read.table(out,header=F)
  unlink(a.file);unlink(out)
  return(res)
}

#function that calls bedtools and operate on two bed dataframes
bedTools.2in<-function(functionstring="bedIntersect",bed1,bed2,opt.string="")
{
  #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out

  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)

  # create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  # cat(command,"\n")
  try(system(command))

  res=read.table(out,header=F)
  unlink(a.file);unlink(b.file);unlink(out)
  return(res)
}


#dREG HD peak calling functions
#scan for continuous peak that passes minimun peak length defined
#returns a matrix in the form of matrix(start, end)
scan_for_peak<-function(potential.peak.positions){
	length.tot<-length(potential.peak.positions)
	if (length.tot==0) return(cbind(numeric(0),numeric(0)))
	peak.starts<-c(1)
	peak.ends<-c()
	i<-1
	while (i< length.tot){
		if(potential.peak.positions[i+1]!=potential.peak.positions[i]+1){
		#break point found
		peak.starts<-c(peak.starts,i+1)
		peak.ends<-c(peak.ends,i)
		}
		i<-i+1
	}
	peak.ends<-c(peak.ends ,length.tot)
	return (cbind(potential.peak.positions [peak.starts[!is.na(peak.ends)]], (potential.peak.positions [peak.ends[!is.na(peak.ends)]]+1)))
}



D1.criterion.filter<-function(peaks, spl.model, start){
	peak.pass.or.not<-c() #true false vector
	for (i in c(1:nrow(peaks))){
		position.pass.or.not<-c()
		for (j in c((peaks[i,1]-start+2):(peaks[i,2]-start+1))){
			D1.val<-predict(spl.model,x=c(j, (j-1)),deriv=1)$y
			if(prod(D1.val)<=0) position.pass.or.not <-c(position.pass.or.not, TRUE)
			else position.pass.or.not <-c(position.pass.or.not, FALSE)
		}
		peak.pass.or.not<-c(peak.pass.or.not, sum(position.pass.or.not)>0) #or ==1? whichever is faster
	}
	return (matrix(peaks [peak.pass.or.not,],byrow=FALSE, ncol=2))
}


split_peak<-function(predicted_data, knots.ratio, background){
	#print(attr(predicted_data,"start"))
	full.vec<-seq(attr(predicted_data,"start"),(attr(predicted_data,"end")-1))
	positions.ONEbased<-full.vec-full.vec[1]+1
	#chose degree of freedom proportional to peak length
	df<-as.integer(length(full.vec)/knots.ratio)+3


	#print(length(positions.ONEbased))
	#print(length(predicted_data))

	dnase.spl<-smooth.spline(positions.ONEbased,predicted_data,df=df)
	D0.val<-predict(dnase.spl,x=c(1:length(full.vec)))$y
	D2.val<-predict(dnase.spl,x=c(1:length(full.vec)),deriv=2)$y
	#get the positions below 2nd derivative cutoff
	potential.peak.positions.D2<-full.vec[D2.val<0]
	#get the positions above the Dnase value background level
	 potential.peak.positions.D0 <- full.vec[predicted_data> background]


	potential.peak.positions<-intersect (potential.peak.positions.D2, potential.peak.positions.D0)
	#scan for continuous peak that passes minimun peak length defined
	peaks<-scan_for_peak(potential.peak.positions)

	#further filtering the peak using "change of sign criterion" of 1st order derivative
	if(!is.null(peaks)){
		if(nrow(peaks)!=0){
    			peaks.D1.filtered <- D1.criterion.filter(peaks, dnase.spl, start=attr(predicted_data,"start"))
			peak.dataframe<-cbind.data.frame(rep(attr(predicted_data,"chrom"), nrow(peaks.D1.filtered)), peaks.D1.filtered)
			return (peak.dataframe)
		}
	}

}


#for regions of outlier (too high) PRO-seq reads, where rarely found in training example, we do scaling by a factor to pull everything back to the distribution of traning examples. This aims at fixing the chunk region in prediction.

fix_distribution<-function(dat){
	#training.max reference determined by 150K+ training examples over true +ve dREG regions (both dnase and grocap +ve) 90% percentile of max(training example input vectors)

	training.max<-4.341949
	dat.max<-max(dat)
	if(dat.max <=training.max) return(dat)
	else return(dat/(dat.max/training.max))
}





#dREG_HD_pred functions

#return the input data
dREG_HD_get_dat <-function( bed_line,  zoom, bigwig_plus, bigwig_minus,total){
	chrom<-bed_line[[1]]
	start<-as.integer(bed_line[[2]])
	end<-as.integer(bed_line[[3]])-1
	positions<-c(start:end)

	#scaling using total read depth
	dat_unscaled <- .Call("get_genomic_data_R", as.character(rep(chrom,length(positions))), as.integer(positions), as.character(bigwig_plus), as.character(bigwig_minus), zoom, FALSE, PACKAGE= "dREG")
	dat_unscaled <-unlist (dat_unscaled)/(total/1E6)
	dat<- cbind(rep(),t(matrix(dat_unscaled, ncol=NROW(positions))))
	dat<-fix_distribution(dat)

	stopifnot(nrow(dat)==length(positions))

	cbind.data.frame(chrom, positions, positions+1, dat)

}


#@param bed is a data.frame from the block of bed_file
run_dREG_HD_pred<-function(gdm, bed,bigwig_plus,bigwig_minus, model,total, temp.bg, ncores, use_rgtsvm) {

	stopifnot(NROW(gdm@window_sizes) == NROW(gdm@half_nWindows))
	zoom<- list(as.integer(gdm@window_sizes), as.integer(gdm@half_nWindows))


	if (nrow(bed)<ncores) blocks=nrow(bed)
	else blocks= ncores

	line.cutoff<-as.integer(seq(from=0, to=nrow(bed),length.out= blocks+1))


	if(use_rgtsvm){

		cpu.fun<-function(idx, line.cutoff, dREG_bed_ext, zoom, bigwig_plus, bigwig_minus,total){
			requireNamespace("dREG")

			do.call(rbind.data.frame,apply(dREG_bed_ext[c((line.cutoff[idx]+1):line.cutoff[idx+1]),],MARGIN=1,FUN= dREG_HD_get_dat,zoom= zoom, bigwig_plus= bigwig_plus, bigwig_minus= bigwig_minus,total= total))
		}

		sfInit(parallel = TRUE, cpus = blocks, type = "SOCK" )
		sfExport("blocks","line.cutoff","bed","zoom","bigwig_plus","bigwig_minus","total");
			#sfExport("blocks","line.cutoff","bed","zoom","bigwig_plus","bigwig_minus","total","dREG_HD_get_dat","fix_distribution");

		dat<-do.call(rbind.data.frame,sfLapply(x=1:blocks,fun= cpu.fun, line.cutoff= line.cutoff, dREG_bed_ext= bed,zoom= zoom,bigwig_plus = bigwig_plus, bigwig_minus = bigwig_minus, total = total))
		sfStop()
		#sfRemoveAll()


		pos<-dat[,1:3]
		ret <- Rgtsvm::predict.gtsvm(model,dat[,4:ncol(dat)])
		rm(dat)
		gc(verbose=TRUE, reset=TRUE)
		stopifnot(nrow(pos)==sum(bed$V3-bed$V2))

		options(scipen =99)
		write.table(cbind.data.frame(pos,ret),file=temp.bg,quote=F,sep="\t",col.names=F,row.names=F,append = TRUE)
	}
	else{

		cpu.fun<-function(idx, line.cutoff, dREG_bed_ext, zoom, bigwig_plus, bigwig_minus,total,model)	{
			requireNamespace("dREG")
			requireNamespace("e1071")
			dat <-do.call(rbind.data.frame,apply(dREG_bed_ext[c((line.cutoff[idx]+1):line.cutoff[idx+1]),],MARGIN=1,FUN= dREG_HD_get_dat,zoom= zoom, bigwig_plus= bigwig_plus, bigwig_minus= bigwig_minus,total= total))
			pos<-dat[,1:3]
			ret <- predict(model,dat[,4:ncol(dat)])
			rm(dat)
			gc(verbose=TRUE, reset=TRUE)
			cbind.data.frame(pos,ret)

		}

		sfInit(parallel = TRUE, cpus = blocks, type = "SOCK" )
		sfExport("blocks","line.cutoff","bed","zoom","bigwig_plus","bigwig_minus","total")
		dat<-do.call(rbind.data.frame,sfLapply(x=1:blocks,fun= cpu.fun, line.cutoff= line.cutoff, dREG_bed_ext= bed,zoom= zoom,bigwig_plus = bigwig_plus, bigwig_minus = bigwig_minus, total = total,model=model))
		sfStop()
		sfRemoveAll()

		stopifnot(nrow(dat)==sum(bed$V3-bed$V2))
		options(scipen =99)

		write.table(dat,file=temp.bg,quote=F,sep="\t",col.names=F,row.names=F,append = TRUE)
	}

	gc(verbose=TRUE, reset=TRUE)
	return(NULL)
}


split.bed.evenly<-function(bed){
	GPU<-1.8
	tot.num.examples <-as.integer(GPU*1024*1024*1024/8/123)
	line.cutoff<-c(0)
	current.row<-1
	while(current.row<=nrow(bed)){
			current.num.examples <-bed[current.row,3]-bed[current.row,2]
		while(current.num.examples <= tot.num.examples && current.row<=nrow(bed)){
			current.row<-current.row+1
			current.num.examples<-current.num.examples+(bed[current.row,3]-bed[current.row,2])

		}
		line.cutoff<-c(line.cutoff,current.row-1)
	}

	return(line.cutoff)

}

get.chromosome.info <- function(bw.plus, bw.minus)
{
    chrom <- rbind( cbind( bw.plus$chroms, bw.plus$chromSizes), cbind( bw.minus$chroms, bw.minus$chromSizes) );
    chr.size <- unlist( lapply( unique(chrom[,1]), function(chr){max( as.numeric( chrom[which(chrom[,1]==chr),2])) } ) );
	return(data.frame( V1=unique(chrom[,1]), V2=chr.size ));
}

#main function runs prediction for blocks of dREG_bed in parallel
#dREG_HD<-function(bed_path, bigwig_plus, bigwig_minus, #chromInfo, model, ncores=1, use_rgtsvm=FALSE){
dREG_HD<-function(bed_path, bigwig_plus, bigwig_minus, model, ncores=1, use_rgtsvm=FALSE){

	#Step1: imputing Dnase-I signal in parallel mode
	cat("running dREG-HD on ", bed_path,"\n");

	bw.plus <- load.bigWig(bigwig_plus);
	bw.minus <-  load.bigWig(bigwig_minus);
	total<-abs(bw.plus$mean*bw.plus$basesCovered)+abs(bw.minus$mean*bw.minus$basesCovered);
	unload.bigWig(bw.plus);
	unload.bigWig(bw.minus);

	ext=200 #impute on extended dREG site
	dREG_bed<-read.table(bed_path);

	#chrom.info.table<-read.table(chromInfo);
	chrom.info.table <- get.chromosome.info( bw.plus, bw.minus );
	chromInfo <- tempfile("chrom.info.");
	write.table( chrom.info.table, file=chromInfo, quote=F, row.names=F, col.names=F, sep="\t");

	dREG_bed_ext<-cbind.data.frame(dREG_bed$V1, apply(cbind(0,(dREG_bed$V2-ext)), MARGIN=1, FUN=max), apply(cbind(sapply(as.character(dREG_bed$V1), FUN= get.chrom.length, chrom.info.table),(dREG_bed$V3+ext)), MARGIN=1, FUN=min));
	dREG_bed_ext<-bedTools.merge(dREG_bed_ext);

	if(use_rgtsvm) class(model)<-"gtsvm"
	gdm<-genomic_data_model(60,30);

	line.cutoff<-split.bed.evenly(dREG_bed_ext);
	blocks<-length(line.cutoff)-1;

	cat("number of total blocks=",blocks,"\n")

	cpu.fun<-function(idx, line.cutoff, dREG_bed_ext,bigwig_plus, bigwig_minus, gdm, model,total, temp.bg, ncores, use_rgtsvm){
		print(idx);
		dREG_bed_ext_part<-dREG_bed_ext[c((line.cutoff[idx]+1):line.cutoff[idx+1]),]
		run_dREG_HD_pred(gdm=gdm, bed= dREG_bed_ext_part,bigwig_plus= bigwig_plus, bigwig_minus= bigwig_minus, model= model,total= total, temp.bg= temp.bg, ncores= ncores, use_rgtsvm= use_rgtsvm)
	  	rm(list=ls());
		gc(verbose=TRUE, reset=TRUE);
		return(NULL);
	  }

	temp.bg=tempfile();
	lapply(c(1:blocks),FUN= cpu.fun,line.cutoff= line.cutoff, dREG_bed_ext= dREG_bed_ext,bigwig_plus= bigwig_plus, bigwig_minus= bigwig_minus, gdm=gdm, model=model,total=total, temp.bg= temp.bg,ncores= ncores, use_rgtsvm= use_rgtsvm)

	#Step2: generate bigwig file from "returned_pred_data" file
	bw.filename<-paste(bed_path,"_imputedDnase.bw",sep="")
	tobigwig(filename=bw.filename, temp.bg = temp.bg, chromInfo= chromInfo)
	unlink(temp.bg);
	rm(model);
	gc(verbose=TRUE, reset=TRUE);

	#Step three generate dREG HD peaks
	imputed_dnase.bw <- load.bigWig(bw.filename);
	returned_pred_data <- bed.step.bpQuery.bigWig(bw= imputed_dnase.bw,bed= dREG_bed_ext,step=1);
	unload.bigWig(imputed_dnase.bw);

	line.cutoff=as.integer(seq(from=0, to=length(returned_pred_data),length.out= ncores+1 ))
	returned_pred_data_list <- list();
	for( i in 1:(length(line.cutoff)-1) )
    	returned_pred_data_list[[i]] <- returned_pred_data[ c((line.cutoff[i]+1):line.cutoff[i+1])];

	rm(returned_pred_data);
	gc(verbose=TRUE, reset=TRUE);

	cpu.fun<-function( pred_data, knots.ratio, background){
	    rbindlist(lapply( pred_data, FUN= split_peak, knots.ratio= knots.ratio, background= background));
	}

	#relaxed mode
	print("calling peaks under relaxed condition");

	# mclapply cause much more memory (ncore times ) are allocated potentially.
	#dREG_HD_bed <- rbindlist( mclapply( returned_pred_data_list, FUN= cpu.fun, knots.ratio= 397.4,background=0.02, mc.cores = ncores ) );

	sfInit(parallel = TRUE, cpus = ncores, type = "SOCK" );
	dREG_HD_bed <- rbindlist( sfLapply( x = returned_pred_data_list,fun= cpu.fun, knots.ratio= 397.4,background=0.02 ) );
	sfStop();

	dREG_HD.filename <- paste(bed_path,"_dREG_HD_relaxed.bed",sep="");

	dREG_HD_bed.intersected<-bedTools.2in ("bedtools intersect -u", dREG_HD_bed, dREG_bed);
	dREG_HD_bed.intersected.merged<-bedTools.merge(dREG_HD_bed.intersected);
	write.table(dREG_HD_bed.intersected.merged,file=dREG_HD.filename,sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE);

	rm(list=c("dREG_HD_bed","dREG_HD_bed.intersected","dREG_HD_bed.intersected.merged"));
	gc(verbose=TRUE, reset=TRUE);

	print("calling peaks under stringent condition");

	#stringent mode
	#dREG_HD_bed <-rbindlist( mclapply( X=1:ncores, FUN= cpu.fun, knots.ratio= 1350, background= 0.02723683,mc.cores = ncores));

	sfInit(parallel = TRUE, cpus = ncores, type = "SOCK" );
	dREG_HD_bed <- rbindlist( sfLapply(x = returned_pred_data_list,fun= cpu.fun, knots.ratio= 1350, background= 0.02723683 ) );
	sfStop();

	dREG_HD.filename <- paste(bed_path,"_dREG_HD_stringent.bed",sep="");

	dREG_HD_bed.intersected<-bedTools.2in ("bedtools intersect -u", dREG_HD_bed, dREG_bed);
	dREG_HD_bed.intersected.merged<-bedTools.merge(dREG_HD_bed.intersected);
	write.table(dREG_HD_bed.intersected.merged,file=dREG_HD.filename,sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE);
}


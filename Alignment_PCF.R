#!/usr/bin Rscript


library("optparse")

option_list = list(
  make_option(c("-F", "--R1_L1"), type="character", default=NULL,
              help="R1 read fasta for lane", metavar="character"),
  make_option(c("-f", "--R2_L1"), type="character", default=NULL,
              help="R2 read fasta for lane", metavar="character"),
  make_option(c("-t", "--tools_dir"), type="character", default="~/Scratch",
              help="Tools DIR location", metavar="character"),
  make_option(c("-c", "--cores"), type="integer", default=6,
              help="Number of cores", metavar="integer"),
  make_option(c("-o", "--output_dir"), type="character", default=".",
              help="Output DIR location", metavar="character"),
  make_option(c("-l", "--R2_L2"), type="character", default="",
              help="R2 read fasta for lane 2", metavar="character"),
  make_option(c("-q", "--seq"), type="character", default="tg",
              help="Sequencing method. Options [tg/wgs]", metavar="character"),
  make_option(c("-Q", "--tg_ver"), type="character", default="v2",
              help="Targeted Panel Version. Only if targeted sequencing. Options [v2/v3]", metavar="character"),
  make_option(c("-L", "--R1_L2"), type="character", default="",
              help="R1 read fasta for lane 2. Chromosomes", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default="hg19",
              help="Genome build to use for reference. Options[hg19/hg38]", metavar="character"),
  make_option(c("-T", "--tmp_dir"), type="character", default="~/Scratch/tmp",
              help="Tmp directory location", metavar="character"),
  make_option(c("-v", "--verbose"), type="character", default="TRUE",
              help="Extra verbose", metavar="character")
);




thisFile <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


print(opt)

if (is.null(opt$R1_L1)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}


script.dir=dirname(thisFile())
if(opt$seq=="tg"){

	if(opt$tg_ver=="v2"){

		input_loc=paste0(script.dir,"/panel/v2")

	}else if(opt$tg_ver=="v3"){
			
		input_loc=paste0(script.dir,"/panel/v3")

	}else{	
		stop("Wrong input for argument -Q. Available options are [v2/v3]")
	}
	prim_tg_il=paste0(input_loc,"/prim_tg.interval_list")
	cap_tg_il=paste0(input_loc,"/cap_tg.interval_list")
	off_tg_bed=paste0(input_loc,"/off_tg.bed")
	prim_tg_bed=paste0(input_loc,"/prim_tg.bed")
	
}else if(opt$seq!="wgs"&opt$seq!="tg"){

  stop("Wrong input for argument -q. Available options are [wgs/tg]")

}

if (opt$genome=="hg19"){

	input_loc=paste0(script.dir,"/hg19")
	ref_genome=paste0(input_loc,"/reference/hs37d5.fa")

} else if (opt$genome=="hg38"){

	input_loc=paste0(script.dir,"/hg38")
	ref_genome=paste0(input_loc,"/reference/ucsc.hg38.fa")
}

dbsnp=paste0(input_loc,"/database/00-common_all.vcf.gz")



sample_name=ULPwgs::intersect_sample_name(ULPwgs::get_sample_name(opt$R1_L1),ULPwgs::get_sample_name(opt$R2_L1))
if (opt$R1_L2!=""){
	sample_name_L2=ULPwgs::intersect_sample_name(ULPwgs::get_sample_name(opt$R1_L2),ULPwgs::get_sample_name(opt$R2_L2))
	sample_name_merged=ULPwgs::intersect_sample_name(sample_name,sample_name_L2)
}else{
	sample_name_merged=sample_name
}
sep="/"

if(opt$output_dir==""){
    sep=""
  }


if (!any(grepl(paste0(opt$output_dir,sep,sample_name_merged,"_DONE"),list.dirs()))){


if (!any(grepl(paste0(opt$output_dir,sep,sample_name_merged,"_trimmed"),list.dirs()))){
	tictoc::tic("FastQC Pre-Trimming:")
	if(opt$R1_L2!=""){
		ULPwgs::fastqc(bin_path=paste0(opt$tools_dir,"/",formals(ULPwgs::fastqc)$bin_path),file_R1=opt$R1_L1,
		file_R2=opt$R2_L1,n_cores=opt$cores,output_dir=paste0(opt$output_dir,sep,sample_name_merged,"_FastQC_reports/L1"),verbose=eval(parse(text= opt$verbose)))
		ULPwgs::fastqc(bin_path=paste0(opt$tools_dir,"/",formals(ULPwgs::fastqc)$bin_path),file_R1=opt$R1_L2,
               	file_R2=opt$R2_L2,n_cores=opt$cores,output_dir=paste0(opt$output_dir,sep,sample_name_merged,"_FastQC_reports/L2"),verbose=eval(parse(text= opt$verbose)))
	}else{
	        ULPwgs::fastqc(bin_path=paste0(opt$tools_dir,"/",formals(ULPwgs::fastqc)$bin_path),file_R1=opt$R1_L1,
               	file_R2=opt$R2_L1,n_cores=opt$cores,output_dir=opt$output_dir,verbose=eval(parse(text = opt$verbose)))
}
	tictoc::toc()
}else{
	print("Skipping FastQC Pre-Trimming")
}


if (!any(grepl(paste0(opt$output_dir,sep,sample_name_merged,"-trimmed_FastQC_reports"),list.dirs()))){
	tictoc::tic("Trimming:")
	if(opt$R1_L2!=""){
		ULPwgs::trimming(bin_path=paste0(opt$tools_dir,"/",formals(ULPwgs::trimming)$bin_path),
		file_R1=opt$R1_L1,file_R2=opt$R2_L1,n_cores=opt$cores,output_dir=paste0(opt$output_dir,sep,sample_name_merged,"_trimmed/L1"),
		verbose=eval(parse(text=opt$verbose)))
		ULPwgs::trimming(bin_path=paste0(opt$tools_dir,"/",formals(ULPwgs::trimming)$bin_path),
		file_R1=opt$R1_L2,file_R2=opt$R2_L2,n_cores=opt$cores,output_dir=paste0(opt$output_dir,sep,sample_name_merged,"_trimmed/L2"),
		verbose=eval(parse(text=opt$verbose)))
	}else{
		ULPwgs::trimming(bin_path=paste0(opt$tools_dir,"/",formals(ULPwgs::trimming)$bin_path),
               	file_R1=opt$R1_L1,file_R2=opt$R2_L1,n_cores=opt$cores,output_dir=opt$output_dir,verbose=eval(parse(text=opt$verbose)))
	}
	tictoc::toc()
}else{
	print("Skipping Trimming")
}
if (!any(grepl(paste0(opt$output_dir,sep,sample_name_merged,"-trimmed_BAM"),list.dirs()))){
	tictoc::tic("FastQC Post-Trimming:")
	if(opt$R1_L2!=""){
		ULPwgs::fastqc(bin_path=paste0(opt$tools_dir,"/",formals(ULPwgs::fastqc)$bin_path),
		file_R1=paste0(opt$output_dir,sep,sample_name_merged,"_trimmed","/L1/",sample_name,"_trimmed/",sample_name,"-trimmed-pair1.fastq.gz"),
		file_R2=paste0(opt$output_dir,sep,sample_name_merged,"_trimmed","/L1/",sample_name,"_trimmed/",sample_name,"-trimmed-pair2.fastq.gz"),	
		n_cores=opt$cores,output_dir=paste0(opt$output_dir,sep,sample_name_merged,"-trimmed_FastQC_reports/L1"),verbose=eval(parse(text = opt$verbose)))
		ULPwgs::fastqc(bin_path=paste0(opt$tools_dir,"/",formals(ULPwgs::fastqc)$bin_path),
               	file_R1=paste0(opt$output_dir,sep,sample_name_merged,"_trimmed","/L2/",sample_name_L2,"_trimmed/",sample_name_L2,"-trimmed-pair1.fastq.gz"),
               	file_R2=paste0(opt$output_dir,sep,sample_name_merged,"_trimmed","/L2/",sample_name_L2,"_trimmed/",sample_name_L2,"-trimmed-pair2.fastq.gz"),
               	n_cores=opt$cores,output_dir=paste0(opt$output_dir,sep,sample_name_merged,"-trimmed_FastQC_reports/L2"),verbose=eval(parse(text = opt$verbose)))
	}else{
		ULPwgs::fastqc(bin_path=paste0(opt$tools_dir,"/",formals(ULPwgs::fastqc)$bin_path),
               	file_R1=paste0(opt$output_dir,sep,sample_name,"_trimmed","/",sample_name,"-trimmed-pair1.fastq.gz"),
               	file_R2=paste0(opt$output_dir,sep,sample_name,"_trimmed","/",sample_name,"-trimmed-pair2.fastq.gz"),
               	n_cores=opt$cores,output_dir=opt$output_dir,verbose=eval(parse(text = opt$verbose)))
}
	tictoc::toc()
}else{
	print("Skipping FastQC Post-Trimming")
}

if (!any(grepl(paste0(opt$output_dir,sep,sample_name_merged,"-trimmed_SORTED.BAM"),list.dirs()))){
	tictoc::tic("Alignment:")
	if(opt$R1_L2!=""){

		#LANE 1
		ULPwgs::alignment(bin_path=paste0(opt$tools_dir,"/",formals(ULPwgs::alignment)$bin_path),
		bin_path2=paste0(opt$tools_dir,"/",formals(ULPwgs::alignment)$bin_path2),ref_genome=ref_genome,
		file_R1=paste0(opt$output_dir,sep,sample_name_merged,"_trimmed","/L1/",sample_name,"_trimmed/",sample_name,"-trimmed-pair1.fastq.gz"),
		file_R2=paste0(opt$output_dir,sep,sample_name_merged,"_trimmed","/L1/",sample_name,"_trimmed/",sample_name,"-trimmed-pair2.fastq.gz"),
		n_cores=opt$cores,output_dir=paste0(opt$output_dir,sep,sample_name_merged,"-trimmed_BAM/L1"),
		verbose=eval(parse(text=opt$verbose)))
		
		#LANE 2
		ULPwgs::alignment(bin_path=paste0(opt$tools_dir,"/",formals(ULPwgs::alignment)$bin_path),
               	bin_path2=paste0(opt$tools_dir,"/",formals(ULPwgs::alignment)$bin_path2),ref_genome=ref_genome,
               	file_R1=paste0(opt$output_dir,sep,sample_name_merged,"_trimmed","/L2/",sample_name_L2,"_trimmed/",sample_name_L2,"-trimmed-pair1.fastq.gz"),
               	file_R2=paste0(opt$output_dir,sep,sample_name_merged,"_trimmed","/L2/",sample_name_L2,"_trimmed/",sample_name_L2,"-trimmed-pair2.fastq.gz"),
               	n_cores=opt$cores,output_dir=paste0(opt$output_dir,sep,sample_name_merged,"-trimmed_BAM/L2"),
		verbose=eval(parse(text = opt$verbose)))
		
		#MERGE LANES
		
		ULPwgs::merge_bams(bin_path=paste0(opt$tools_dir,"/",formals(ULPwgs::concatenate_bams)$bin_path),
		bams=c(paste0(opt$output_dir,sep,sample_name_merged,
		"-trimmed_BAM","/L1/",sample_name,"-trimmed_BAM/",
		sample_name,"-trimmed.bam"),paste0(opt$output_dir,sep,sample_name_merged,"-trimmed_BAM","/L2/",sample_name_L2,"-trimmed_BAM/",
		sample_name_L2,"-trimmed.bam")),output_name=paste0(opt$output_dir,sep,sample_name_merged,"-trimmed_BAM/",sample_name_merged,"-trimmed"),threads=opt$cores)
	}else{
		#NO LANE
		ULPwgs::alignment(bin_path=paste0(opt$tools_dir,"/",formals(ULPwgs::alignment)$bin_path),
               	bin_path2=paste0(opt$tools_dir,"/",formals(ULPwgs::alignment)$bin_path2),ref_genome=ref_genome,
               	file_R1=paste0(opt$output_dir,sep,sample_name,"_trimmed","/",sample_name,"-trimmed-pair1.fastq.gz"),
               	file_R2=paste0(opt$output_dir,sep,sample_name,"_trimmed","/",sample_name,"-trimmed-pair2.fastq.gz"),
               	n_cores=opt$cores,output_dir=opt$output_dir,verbose=eval(parse(text = opt$verbose)))

}
	tictoc::toc()
}else{
	print("Skipping Alignment")
}

if (!any(grepl(paste0(opt$output_dir,sep,sample_name_merged,"-trimmed_SORTED.RMDUP.SORTED.BAM"),list.dirs()))){
	tictoc::tic("Sorting and Indexing Pre-Duplicate Removal:")
	
	ULPwgs::sort_and_index(bin_path=paste0(opt$tools_dir,"/",formals(ULPwgs::sort_and_index)$bin_path),
	file=paste0(opt$output_dir,sep,sample_name_merged,"-trimmed_BAM","/",sample_name_merged,"-trimmed.bam"),
	output_dir=opt$output_dir,verbose=eval(parse(text = opt$verbose)),threads=opt$cores,coord_sort=TRUE)
	tictoc::toc()
}else{
	print("Skipping Sorting and Indexing Pre-Duplicate Removal")
}


### Switched to a multi-threadead version of MarkDuplicates using GATK. Legacy code.

##if (!any(grepl(paste0(opt$output_dir,sep,sample_name_merged,"-trimmed_SORTED.RMDUP.SORTED.BAM"),list.dirs()))){
##	tictoc::tic("Remove Duplicates:")
##	ULPwgs::remove_duplicates(bin_path=paste0(opt$tools_dir,"/",formals(ULPwgs::remove_duplicates)$bin_path),
##	file=paste0(opt$output_dir,sep,sample_name_merged,"-trimmed_SORTED.BAM","/",sample_name_merged,"-trimmed.SORTED.bam"),hnd=5000,ram=8,
##	tmp_dir=opt$tmp_dir,output_dir=opt$output_dir)
##	tictoc::toc()
## }else{
## 	print("Skipping Remove Duplicates")
## }


## if (!any(grepl(paste0(opt$output_dir,sep,sample_name_merged,"-trimmed_alignQC_report"),list.dirs()))){
##       	tictoc::tic("Sorting and Indexing Post-Duplicate Removal:")
##       	ULPwgs::sort_and_index(bin_path=paste0(opt$tools_dir,"/",formals(ULPwgs::sort_and_index)$bin_path),
##      	file=paste0(opt$output_dir,sep,sample_name_merged,"-trimmed_RMDUP.SORTED.BAM","/",sample_name_merged,"-trimmed.RMDUP.SORTED.bam"),
##		output_dir=opt$output_dir,verbose=eval(parse(text = opt$verbose)),threads=opt$cores,coord_sort=TRUE)
##       	tictoc::toc()
## }else{
##       	print("Skipping Sorting and Indexing Post-Duplicate Removal:")
##}



if (!any(grepl(paste0(opt$output_dir,sep,sample_name_merged,"-trimmed_RECAL.SORTED.RMDUP.SORTED.BAM"),list.dirs()))){
	tictoc::tic("Remove Duplicates:")
 	ULPwgs::remove_duplicates_gatk(bin_path=paste0(opt$tools_dir,"/",formals(ULPwgs::remove_duplicates_gatk)$bin_path),
	file=paste0(opt$output_dir,sep,sample_name_merged,"-trimmed_SORTED.BAM","/",sample_name_merged,"-trimmed.SORTED.bam"),
	tmp_dir=opt$tmp_dir,threads=opt$cores,output_dir=opt$output_dir,verbose=eval(parse(text = opt$verbose)))
     	tictoc::toc()
}else{
     print("Skipping Remove Duplicates")
}


if (!any(grepl(paste0(opt$output_dir,sep,sample_name_merged,"-trimmed_alignQC_reports"),list.dirs()))){
       	tictoc::tic("Recalibrating Base Quality:")
       	ULPwgs::recalibrate_bq(bin_path=paste0(opt$tools_dir,"/",formals(ULPwgs::recalibrate_bq)$bin_path),bin_path2=paste0(opt$tools_dir,"/",	
	formals(ULPwgs::recalibrate_bq)$bin_path2),bin_path3=paste0(opt$tools_dir,"/",
	formals(ULPwgs::recalibrate_bq)$bin_path3),ref_genome=ref_genome,
	bam=paste0(opt$output_dir,sep,sample_name_merged,
	"-trimmed_SORTED.RMDUP.SORTED.BAM","/",sample_name_merged,"-trimmed.SORTED.RMDUP.SORTED.bam"),
	snpdb=dbsnp,threads=opt$cores,
	output_dir=opt$output_dir,verbose=eval(parse(text = opt$verbose)))
       	tictoc::toc()
}else{
       	print("Skipping Base Recalibration :")
}



if (!any(grepl(paste0(opt$output_dir,sep,sample_name_merged,"_DONE"),list.dirs()))){
	tictoc::tic("Alignment Quality Control:")
	if(opt$seq=="wgs"){
		ULPwgs::qc_metrics(bin_path=paste0(opt$tools_dir,"/",formals(ULPwgs::qc_metrics)$bin_path),
		bin_path2=paste0(opt$tools_dir,"/",formals(ULPwgs::qc_metrics)$bin_path2),
		bin_path3=paste0(opt$tools_dir,"/",formals(ULPwgs::qc_metrics)$bin_path3),
		bam=paste0(opt$output_dir,sep,sample_name_merged,"-trimmed_RECAL.SORTED.RMDUP.SORTED.BAM","/",
		sample_name_merged,"-trimmed.SORTED.RECAL.SORTED.RMDUP.SORTED.bam"),
		ram=8,tmp_dir=opt$tmp_dir,ref_genome=ref_genome,verbose=eval(parse(text = opt$verbose)),output_dir=opt$output_dir)


	}else if (opt$seq=="tg"){
		ULPwgs::qc_metrics(bin_path=paste0(opt$tools_dir,"/",formals(ULPwgs::qc_metrics)$bin_path),
        	bin_path2=paste0(opt$tools_dir,"/",formals(ULPwgs::qc_metrics)$bin_path2),
		bin_path3=paste0(opt$tools_dir,"/",formals(ULPwgs::qc_metrics)$bin_path3),
        	bam=paste0(opt$output_dir,sep,sample_name_merged,"-trimmed_RECAL.SORTED.RMDUP.SORTED.BAM",
		"/",sample_name_merged,"-trimmed.SORTED.RECAL.SORTED.RMDUP.SORTED.bam"),
        	ram=8,bi=cap_tg_il,ti=prim_tg_il,off_tar=off_tg_bed,
        	on_tar=prim_tg_bed,tmp_dir=opt$tmp_dir,ref_genome=ref_genome,verbose=eval(parse(text = opt$verbose)),output_dir=opt$output_dir)	
	}
	tictoc::toc()
}else{
	print("Skipping Alignment Quality Control")
	}
}


dir.create(paste0(opt$output_dir,sep,sample_name_merged,"_DONE"))
print("Cleaning directory...")
system(paste("rm -rf" ,paste0(opt$output_dir,sep,sample_name,"_trimmed")))
##system(paste("rm -rf" ,paste0(opt$output_dir,sep,sample_name,"-trimmed_BAM")))
system(paste("rm -rf" ,paste0(opt$output_dir,sep,sample_name,"-trimmed_SORTED.BAM/*.bam*")))
system(paste("rm -rf" ,paste0(opt$output_dir,sep,sample_name,"-trimmed_SORTED.RMDUP.SORTED.BAM/*.bam*")))

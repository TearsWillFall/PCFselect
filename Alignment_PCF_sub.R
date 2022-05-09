##!/usr/bin Rscript


library("optparse")
library("tidyverse")

option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default=NULL,
              help="FastQ Input Directory", metavar="character"),
  make_option(c("-r", "--rds"), type="logical", default=TRUE,
              help="Get files from RDS", metavar="logical"),
  make_option(c("-t", "--tools_dir"), type="character", default="~/Scratch",
              help="Tools DIR location", metavar="character"),
  make_option(c("-c", "--cores"), type="integer", default=6,
              help="Number of cores", metavar="integer"),
  make_option(c("-o", "--output_dir"), type="character", default=".",
              help="Output DIR location", metavar="character"),
  make_option(c("-q", "--seq"), type="character", default="tg",
              help="Sequencing method. Options [tg/wgs]", metavar="character"),
  make_option(c("-Q", "--tg_ver"), type="character", default="v2",
              help="Targeted Panel Version. Only if targeted sequencing. Options [v2/v3]", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default="hg19",
              help="Genome build to use for reference. Options[hg19/hg38]", metavar="character"),
  make_option(c("-T", "--tmp_dir"), type="character", default="~/Scratch/tmp",
              help="Tmp directory location", metavar="character"),
  make_option(c("-v", "--verbose"), type="character", default="TRUE",
              help="Extra verbose", metavar="character")

);



opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$input_dir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

files=system(paste0("ssh ssh.rd.ucl.ac.uk \' find ",opt$input_dir,"|grep .gz$| egrep \"fq|fasta\"  |sort \'"),intern=TRUE)

files=data.frame(file=files,filename=unlist(lapply(files,FUN=function(x){basename(x)})),id=unlist(lapply(files,FUN=function(x){basename(dirname(x))})))
files= files %>% group_by(id) %>% mutate(pos=1:n()) %>% mutate(read=ifelse(pos%%2==0,2,1),lane=ceiling(pos/2))

lapply(unique(files$id),FUN=function(x){tmp_files=files %>% filter(id==x);
	system(paste0("mkdir ",opt$output_dir,"/",x),wait=TRUE)
	
	lapply(tmp_files$file,FUN=function(y){
		if(opt$rds){
			system(paste0("scp ssh.rd.ucl.ac.uk:",y," ",opt$output_dir,"/",x),wait=TRUE)
		
		}else{
			system(paste0("cp ",y," ",opt$output_dir,"/",x),wait=TRUE)
		}
	})
	R1_L1=tmp_files %>% filter(read==1,lane==1)
	R2_L1=tmp_files %>% filter(read==2,lane==1)
	R1_L2=tmp_files %>% filter(read==1,lane==2)
        R2_L2=tmp_files %>% filter(read==2,lane==2)

	system(paste("qsub -N ",x,"-pe smp ",opt$cores," -wd ",paste0(opt$output_dir,"/",x)," Alignment_PCF.sh ",paste0(opt$output_dir,"/",x,"/",R1_L1$filename),
	paste0(opt$output_dir,"/",x,"/",R2_L1$filename), 
	paste0(opt$output_dir,"/",x), opt$genome,  opt$seq, opt$tg_ver, opt$tmp_dir,  opt$tools_dir, opt$cores, opt$verbose, 
	ifelse(!is.null(R1_L2$filename)&!is.null(R2_L2$filename),paste0(R1_L2$filename, " ",R2_L2$filename),"")))	
	

        print(paste("qsub -N ",x,"-pe smp ",opt$cores," -wd ",paste0(opt$output_dir,"/",x)," Alignment_PCF.sh ",paste0(opt$output_dir,"/",x,"/",R1_L1$filename),
        paste0(opt$output_dir,"/",x,"/",R2_L1$filename),
        paste0(opt$output_dir,"/",x), opt$genome,  opt$seq, opt$tg_ver,
        opt$tmp_dir,  opt$tools_dir, opt$cores, opt$verbose,
        ifelse(!is.null(R1_L2$filename)&!is.null(R2_L2$filename),paste0(paste0(opt$output_dir,"/",x,"/",R1_L2$filename), " ",paste0(opt$output_dir,"/",x,"/",R2_L2$filename)))))

	
})



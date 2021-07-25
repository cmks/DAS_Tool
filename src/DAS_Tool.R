#!/usr/bin/env Rscript

#
# DAS Tool for genome-resolved metagenomics
# by Christian MK Sieber (cmksieber@gmail.com)
#
# Please cite https://doi.org/10.1101/107789
#

doc <- "DAS Tool

Usage:
  DAS_Tool [options] -i <contig2bin> -c <contigs_fasta> -o <outputbasename>
  DAS_Tool -i <contig2bin> -c <contigs_fasta> -o <outputbasename> [--labels=<labels>] [--proteins=<proteins_fasta>] [--threads=<threads>] [--search_engine=<search_engine>] [--score_threshold=<score_threshold>] [--dbDirectory=<dbDirectory> ] [--megabin_penalty=<megabin_penalty>] [--duplicate_penalty=<duplicate_penalty>] [--write_bin_evals] [--create_plots] [--write_bins] [--write_unbinned] [--resume] [--debug]
  DAS_Tool [--version]
  DAS_Tool [--help]

Options:
   -i --bins=<contig2bin>                   Comma separated list of tab separated contigs to bin tables.
   -c --contigs=<contigs>                   Contigs in fasta format.
   -o --outputbasename=<outputbasename>     Basename of output files.
   -l --labels=<labels>                     Comma separated list of binning prediction names.
   --search_engine=<search_engine>          Engine used for single copy gene identification (blast/diamond/usearch) [default: diamond].
   -p --proteins=<proteins>                 Predicted proteins (optional) in prodigal fasta format (>contigID_geneNo).
                                            Gene prediction step will be skipped.
   --write_bin_evals                        Write evaluation of input bin sets.
   --write_bins                             Export bins as fasta files.
   --write_unbinned                         Export unbinned contigs as fasta file (--write_bins needs to be set).
   -t --threads=<threads>                   Number of threads to use [default: 1].
   --score_threshold=<score_threshold>      Score threshold until selection algorithm will keep selecting bins (0..1) [default: 0.5].
   --duplicate_penalty=<duplicate_penalty>  Penalty for duplicate single copy genes per bin (weight b).
                                            Only change if you know what you are doing (0..3) [default: 0.6].
   --megabin_penalty=<megabin_penalty>      Penalty for megabins (weight c). Only change if you know what you are doing (0..3) [default: 0.5].
   --dbDirectory=<dbDirectory>              Directory of single copy gene database [default: db].
   --customDbDir=<customDbDir>              Directory of custom single copy gene database(s). One directory per set. Multiple sets can be comma separated.
   --useCustomDbOnly                        Only use custom single copy gene database(s).
   --customDbFormat                         Format of custom single copy gene database(s) [default: anvio].
   --resume                                 Use existing predicted single copy gene files from a previous run.
   --debug                                  Write debug information to log file.
   -v --version                             Print version number and exit.
   -h --help                                Show this.


Please cite: Sieber et al., 2018, Nature Microbiology (https://doi.org/10.1038/s41564-018-0171-1).
"

if(length(commandArgs(trailingOnly = TRUE)) == 0L) {
   docopt:::help(doc)
   quit()
}

##
## Load packages, define functions
##

suppressMessages(suppressWarnings(library(data.table, warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(magrittr, warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(docopt, warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(rhmmer, warn.conflicts = F, quietly = T)))

cherry_pick <- function(binTab,scgTab,contigTab,output_basename,score_threshold,duplicate_penalty,megabin_penalty,write_unbinned,write_bin_evals){
   
   # thresholds:
   score_threshold <- max(score_threshold,-42)
   internal_ratio_threshold <- 0.0
   internal_score_threshold <- min(0.0,score_threshold)
   
   # set keys
   setkey(binTab,'contig_id')
   setkey(scgTab,'contig_id')
   setkey(contigTab,'contig_id')
   
   # join tables
   binTabScg <- binTab[scgTab]
   binTabContig <- contigTab[binTab]
   
   # score bins
   binTabEval <- score_bins(bin_tab_scg = binTabScg, bin_tab_contig = binTabContig, b = duplicate_penalty,c = megabin_penalty)
   
   if(write_bin_evals){
      write.table(binTabEval[,.(bin=bin_id,
                                bin_set=binner_name,
                                unique_SCGs=uniqueSCG,
                                redundant_SCGs=duplicatedSCG,
                                SCG_set=protein_set,
                                size=binSize,
                                contigs=nContig,
                                N50=contigN50,
                                bin_score=score,
                                SCG_completeness=round(completeness,2)*100,
                                SCG_redundancy=round(contamination,2)*100)],
                  paste(output_basename,'_allBins.eval',sep=''),sep='\t',col.names=T,row.names=F,quote=F)
   }
   
   max_score <- max(binTabEval[,score])
   
   if(max_score < score_threshold){
      write.log(message = paste0('No bins with bin-score >', score_threshold, ' found. Adjust score_threshold to report bins with lower quality.\n Aborting.\n'),filename = logfile)
      return(-1)
   }
   
   append <- F
   sub_bins <- c()
   # write.log(message = paste0('starting bin selection from ', dim(bin_summary)[1],' bins'),filename = logfile,write_to_stdout = T)
   while(max_score > internal_score_threshold){
      
      # identify highest scoring bin
      topBin <- binTabEval[ 1 ]
      
      # identify contigs of highest scoring bin
      topContig2Bin <- binTabContig[bin_id == topBin[,bin_id],.(contig_id,bin_id)]
      
      
      if(max_score >= score_threshold && max(binTabEval[,completeness]) > internal_ratio_threshold){
         
         # write contig2bins of highest scoring bin
         fwrite(topContig2Bin,paste0(output_basename,'_DASTool_contig2bin.tsv'),sep='\t',col.names=F,row.names = F,quote=F,append = append)
         
         # write summary of highest scoring bin
         fwrite(topBin[,.(bin=bin_id,
                          bin_set=binner_name,
                          unique_SCGs=uniqueSCG,
                          redundant_SCGs=duplicatedSCG,
                          SCG_set=protein_set,
                          size=binSize,
                          contigs=nContig,
                          N50=contigN50,
                          bin_score=score,
                          SCG_completeness=round(completeness,2)*100,
                          SCG_redundancy=round(contamination,2)*100)],
                paste0(output_basename,'_DASTool_summary.tsv'),sep='\t',col.names=(!append),row.names = F,quote=F,append = append)
         
         append <- T
      }else{
         # write unbinned contigs if bin-score is below threshold
         if(write_unbinned){
            topContig2Bin[,bin_id:='unbinned'] %>% 
               fwrite(paste0(output_basename,'_DASTool_contig2bin.tsv'),sep='\t',col.names=F,row.names = F,quote=F,append = append)
         }
      }
      
      # remove contigs of highest scoring bin
      affected_bins <- unique(binTabContig[ contig_id %in% topContig2Bin[,contig_id], bin_id ])
      binTabContig[ bin_id %in% affected_bins, bin_id:= paste(gsub('_sub$','',bin_id),'sub',sep='_')]
      binTabContig <- binTabContig[ !contig_id %in% topContig2Bin[,contig_id] ]
      binTabScg[ bin_id %in% affected_bins, bin_id:= paste(gsub('_sub$','',bin_id),'sub',sep='_')]
      binTabScg <- binTabScg[ !contig_id %in% topContig2Bin[,contig_id] ]
      
      # update scores of remaining bins
      binTabEval <- score_bins(bin_tab_scg = binTabScg, bin_tab_contig = binTabContig, b = duplicate_penalty,c = megabin_penalty)
      
      
      if(nrow(binTabEval) == 0){
         max_score <- -Inf
      }else{
         max_score <- max(binTabEval[,score])
      }
   }
   
   # write remaining unbinned contigs
   if(write_unbinned){
      binTabContig[,.(contig_id)] %>% 
         .[,bin_id:='unbinned'] %>% 
         fwrite(paste0(output_basename,'_DASTool_contig2bin.tsv'),sep='\t',col.names=F,row.names = F,quote=F,append = append)
      
   }
   
   # write.log(message = paste0('bin selection complete: ',reported_bins,' bins above score threshold selected.'),filename = logfile,write_to_stdout = T)
   
}


calc_N50 <- function(contig_id,contig_length){
   
   seq_len <- data.table(contig_id,contig_length) %>% 
      .[!duplicated(contig_id)] %>% 
      .[order(contig_length,decreasing = T)] %>% 
      .[,contig_length]
   
   N50 <- seq_len[cumsum(seq_len) > sum(seq_len)/2][1] 
   
   return(N50)
}


calc_bins_size <- function(contig_id,contig_length){
   
   binSize <- data.table(contig_id,contig_length) %>%
      .[!duplicated(contig_id)] %>% 
      .[,contig_length] %>% 
      sum()
   
   return(binSize)
}


calc_duplicates <- function(protein_names){
   tmp <- table(protein_names)
   
   return(max(length(tmp[tmp>1]),0,na.rm = T))
}


calc_bin_score <- function(b,c,protein_set_size,uniqueSCG,duplicatedSCG,sumSCG,additionalSCG){
   
   bin_score <- (uniqueSCG/protein_set_size) - b*(duplicatedSCG / uniqueSCG) - c*(additionalSCG /  protein_set_size)
   
   return(bin_score)
}


score_bins <- function(bin_tab_scg,bin_tab_contig,b=.6,c=.5){
   bin_tab_scg_eval <- bin_tab_scg[,.(uniqueSCG=length(unique(protein_name)),
                                      duplicatedSCG=calc_duplicates(protein_name),
                                      sumSCG=.N),by=c('bin_id','protein_set','binner_name','protein_set_size')] %>% 
      .[,additionalSCG:= (sumSCG - uniqueSCG)] %>% 
      setkey(bin_id,binner_name)
   
   bin_tab_contig_eval <- bin_tab_contig[,.(binSize=calc_bins_size(contig_id,contig_length),
                                            contigN50=calc_N50(contig_id,contig_length),
                                            nContig=.N),by=c('bin_id','binner_name')] %>% 
      setkey(bin_id,binner_name)
   
   bin_tab_eval <- bin_tab_contig_eval[bin_tab_scg_eval] %>% 
      .[,score:=calc_bin_score(b=.6,c=.5,protein_set_size,uniqueSCG,duplicatedSCG,sumSCG,additionalSCG)] %>% 
      .[,completeness:=uniqueSCG/protein_set_size] %>% 
      .[,contamination:=duplicatedSCG/protein_set_size] %>% 
      .[ order(score,contigN50,decreasing = T)] %>% 
      .[.[, .I[which.max(score)], by=bin_id]$V1]
   
   return(bin_tab_eval)
}

exit <- function() {
   invokeRestart("abort") 
} 

write.log <- function(message,append=T,filename,write_to_file=T,type='none'){
   
   if( 'character' %in% class(message)){
      if(write_to_file){
         cat(message,'\n',file=filename,sep=' ',append = append)
      }
      if(type == 'cat'){
         cat(message,'\n')
      }
      if(type == 'warning'){
         cat('Warning: ',message,'\n')
         # warning(paste0(message,'\n'))
      }
      if(type == 'stop'){
         cat('Error: ',message,'\n')
         exit()
      }
   }else if("data.frame" %in% class(message)){
      write.table(message,file=filename,sep='\t',col.names=F,row.names = F,quote=F,append = append)
      cat('\n',file=filename,sep=' ',append = append)
   }
}


##
## Run DAS Tool
##
version <- 'DAS Tool 1.1.4\n'

arguments <- docopt(doc, version = version)
# print(arguments)

# Init log file
logFile <- paste0(arguments$outputbasename,'_DASTool.log')
cat('\n')
write.log(version,filename = logFile,append = F,write_to_file = T,type = 'cat')


write.log(paste(Sys.time()),filename = logFile,append = T,write_to_file = T)

write.log('\nParameters:',filename = logFile,append = T,write_to_file = T)
argsTab <- data.table(opt=names(arguments),args=arguments) %>% 
   .[grepl('^--',opt)]
write.log(argsTab,filename = logFile,append = T,write_to_file = T)


threads <- max(1,as.numeric(arguments$threads))

# Check options
## search engine %in% diamond, usearch, blastp?
if(!tolower(arguments$search_engine) %in% c('diamond', 'usearch', 'blastp')){
   write.log(paste0('Unknown argument for --search_engine: ',arguments$search_engine,'\n',
                    'Defaulting to diamond'),filename = logFile,append = T,write_to_file = T,type = 'warning')
   searchEngine <- 'diamond'
}else{
   searchEngine <- tolower(arguments$search_engine)
}

# Check dependencies
dependencies <- Sys.which(c("prodigal", "diamond", "pullseq", "ruby", "usearch", "blastp")) %>% 
   data.table(dependency=names(.),path=.)

min_dependencies <- c("pullseq", "ruby",searchEngine)
if(is.null(arguments$proteins)){
   min_dependencies <- c(min_dependencies,'prodigal')
}

if(any(dependencies[ path == '', dependency ] %in% min_dependencies)){
   
   write.log(paste0('Cannot find dependencies: ', 
                    paste0(dependencies[ path == '', dependency ],collapse = ', ')),filename = logFile,append = T,write_to_file = T,type = 'stop')
}

write.log('\nDependencies:',filename = logFile,append = T,write_to_file = T)
write.log(dependencies,filename = logFile,append = T,write_to_file = T,type = 'cat')


# Get script directory:
cmdArgs <- commandArgs(trailingOnly = FALSE)
scriptDir <- dirname(normalizePath(sub("--file=", "", cmdArgs[grep("--file=", cmdArgs)])))
scriptDir <- gsub('\\/src$','',scriptDir)
# print(scriptDir)

# Set default db directory if not defined:
dbDirectory <- ifelse(arguments$dbDirectory == 'db',
                       paste0(scriptDir,'/','db/'),
                       arguments$dbDirectory)


##
## Check files and directories
##
binSets <- unlist(strsplit(arguments$bins,','))

# Input files existence
## Contig2bin tables:
for(binSet in binSets){
   if(!file.exists(binSet)){
      write.log(paste0('File does not exist: ',binSet),filename = logFile,append = T,write_to_file = T,type = 'stop')
   }
}
## Contig2bin tables checks:
### - duplicated contig names - does it matter?
### - do contig names match assembly contigs?


## Assembly:
if(!file.exists(arguments$contigs)){
   write.log(paste0('File does not exist: ',arguments$contigs),filename = logFile,append = T,write_to_file = T,type = 'stop')
   # stop(paste0('File does not exist: ',arguments$contigs))
}
if(file.size(arguments$contigs) == 0){
   write.log(paste0('File is empty: ',arguments$contigs),filename = logFile,append = T,write_to_file = T,type = 'stop')
   # stop(paste0('File does not exist: ',arguments$contigs))
}
tmpfile <- file(arguments$contigs)
if(summary( tmpfile )$class == 'gzfile'){
   write.log(paste0('Unsupported file: Assembly file is compressed. Please extract: ',arguments$contigs),filename = logFile,append = T,write_to_file = T,type = 'stop')
}
close(tmpfile)

# Check bin set labels
if( !is.null(arguments$labels) ){
   binSetLabels <- unlist(strsplit(arguments$labels,','))
   
   if(length(binSetLabels) != length(binSets)){
      write.log(paste0('Number of bin sets (', length(binSets),') is different from number of labels (', length(binSetLabels),')'),
                filename = logFile,append = T,write_to_file = T,type = 'warning')
      # warning('Number of bin sets (', length(binSets),') is different from number of labels (', length(binSetLabels),')')
      binSetLabels <- basename(binSets)
   }
   if(any(duplicated(binSetLabels))){
      write.log('Non-unique bin set labels given', filename = logFile,append = T,write_to_file = T,type = 'warning')
      # warning('Duplicated bin set labels given')
      binSetLabels <- paste0(binSetLabels,'.',sprintf("%02d",c(1:length(binSets))))
   }
}else{
   binSetLabels <- paste0('binner.',sprintf("%02d",c(1:length(binSets))))
}

binSetToLabel <- data.table(inputNo=c(1:length(binSets)),binSetLabel=binSetLabels,binSet=binSets)

# Check contig2bin files:
binTab <- data.table()
for(i in 1:nrow(binSetToLabel)){
   # print(binSetToLabel[i,binSet])
   if(file.size(binSetToLabel[i,binSet]) == 0){
      ## Warn if contig2bin table is empty, but keep going:
      write.log(paste0('File is empty: ',binSetToLabel[i,binSet]), filename = logFile,append = T,write_to_file = T,type = 'warning')
      # warning(paste0('File is empty: ',binSetToLabel[i,binSet]))
   }else{
      ## Read contig2bin table and concatenate:
      dt <- fread(binSetToLabel[i,binSet],header=F,col.names=c('contig_id','bin_id')) %>% 
         .[,binner_name:=binSetToLabel[i,binSetLabel]]
      binTab <- rbind(binTab,dt)
   }
}

## Stop if all contig2bin tables are empty:
if(nrow(binTab) == 0){
   write.log('No bins provided. Check contig2bin files.', filename = logFile,append = T,write_to_file = T,type = 'error')
   # stop('No bins provided. Check contig2bin files.')
}

# Check for duplicate bin-IDs
# TODO: Check this!
if(nrow(unique(binTab[,.(binner_name,bin_id)])) > length(unique(binTab[,bin_id]))){
   binTab[,bin_id:=paste0(binner_name,'_',bin_id)]
   write.log(paste0('Non-unique bin-ids given. Renaming bin-ids.'), filename = logFile,append = T,write_to_file = T,type = 'warning')
}


# print(binSetToLabel[,.(binSetLabel,binSet)])

# Existence and permissions of output directory
if(! dir.exists(dirname(arguments$outputbasename))){
   write.log(paste0('Output directory does not exist. Creating: ', dirname(arguments$outputbasename)), 
             filename = logFile,append = T,write_to_file = T,type = 'warning')
   # warning('Output directory does not exist. Creating: ', dirname(arguments$outputbasename))
   system(paste0('mkdir ',dirname(arguments$outputbasename)))
}

# Existence of database directory
if( !file.exists(paste0(dbDirectory,'bac.all.faa')) || 
    !file.exists(paste0(dbDirectory,'arc.all.faa')) || 
    !file.exists(paste0(dbDirectory,'bac.scg.faa')) || 
    !file.exists(paste0(dbDirectory,'arc.scg.faa')) || 
    !file.exists(paste0(dbDirectory,'bac.scg.lookup')) || 
    !file.exists(paste0(dbDirectory,'arc.scg.lookup'))){
   write.log(paste0('Database directory does not exist or is incomplete:', dbDirectory,'\n
                  Attempting to extract ',scriptDir,'/db.zip'), 
             filename = logFile,append = T,write_to_file = T,type = 'warning')
   if(file.exists(paste0(scriptDir,'/db.zip'))){
      system(paste0('unzip -o ',scriptDir,'/db.zip -d ',dbDirectory))
   }else{
      write.log(paste0('File does not exist: ',scriptDir,'/db.zip'), 
                filename = logFile,append = T,write_to_file = T,type = 'error')
      # stop(paste0('File does not exist: ',scriptDir,'/db.zip'))
   }
}

##
## Calculate contig lengths
##
write.log('Analyzing assembly',filename = logFile,append = T,write_to_file = T,type = 'cat')
system(paste0('bash ', scriptDir,'/src/contig_length.sh ',
              arguments$contigs,' ',
              arguments$outputbasename,'.seqlength'))

contigTab <- fread(paste0(arguments$outputbasename,'.seqlength'),header=F,col.names=c('contig_id','contig_length')) %>% 
   .[,contig_id:=gsub(' .*','',contig_id)]

stopifnot(nrow(contigTab)>0)
# Check if bin-set contigs are in assembly:
missingContigs <- binTab[!contig_id %in% contigTab[,contig_id]]
if( nrow(missingContigs) > 0){
   write.log(paste0('Error: Contigs of contig2bin files not found in assembly: ', nrow(missingContigs)), 
             filename = logFile,append = T,write_to_file = T,type = 'cat')
   write.log(head(missingContigs[,.(contig_id)]), 
             filename = logFile,append = T,write_to_file = T,type = 'cat')
   write.log('...', 
             filename = logFile,append = T,write_to_file = T,type = 'cat')
   exit()
}

##
## Predict single copy genes
##

# Run prodigal
if(!is.null(arguments$proteins)){
   proteins <- arguments$proteins
   write.log(paste0('Skipping gene prediction, using protein fasta file: ',proteins),filename = logFile,append = T,write_to_file = T,type = 'cat')
}else{
   write.log('Predicting genes',filename = logFile,append = T,write_to_file = T,type = 'cat')
   
   proteins <- paste0(arguments$outputbasename,'_proteins.faa')

   if(threads == 1){
      system(paste('prodigal -i ',arguments$contigs,
                   ' -a ',proteins,
                   ' -p meta -m -q > /dev/null 2>&1'))
   }else{
      system(paste0('bash ', scriptDir,'/src/prodigal_parallel.sh ',
                    arguments$contigs,' ',
                    proteins, ' ',
                    threads))
   }
}

# Identify single copy genes
scgTab <- data.table()
## Predict bacterial SCGs
write.log('Annotating single copy genes',filename = logFile,append = T,write_to_file = T,type = 'cat')
system(paste0('ruby ', scriptDir,'/src/scg_blank_',searchEngine,'.rb ',
              searchEngine,' ',
              proteins,' ',
              dbDirectory,'/bac.all.faa ',
              dbDirectory,'/bac.scg.faa ',
              dbDirectory,'/bac.scg.lookup ',
              threads,
              ifelse(arguments$debug,paste0(' 2>&1 | tee -a ',logFile),' > /dev/null 2>&1')))
system(paste0('mv ',proteins,'.scg ',proteins,'.bacteria.scg'))

## Predict archaeal SCGs
system(paste0('ruby ', scriptDir,'/src/scg_blank_',searchEngine,'.rb ',
              searchEngine,' ',
              proteins,' ',
              dbDirectory,'/arc.all.faa ',
              dbDirectory,'/arc.scg.faa ',
              dbDirectory,'/arc.scg.lookup ',
              threads,
              ifelse(arguments$debug,paste0(' 2>&1 | tee -a ',logFile),' > /dev/null 2>&1')))
system(paste0('mv ',proteins,'.scg ',proteins,'.archaea.scg'))

## Stop if no single copy genes were predicted:
if(file.size(paste0(proteins,'.bacteria.scg')) == 0 || file.size(paste0(proteins,'.archaea.scg')) == 0){
   write.log('No single copy genes predicted', 
             filename = logFile,append = T,write_to_file = T,type = 'error')
   # stop('No single copy genes predicted')
}else{
   scgTab <- rbind(fread(paste0(proteins,'.bacteria.scg'),header=F,col.names=c('protein_id','protein_name')) %>% 
                            .[,protein_set:='bacteria'] %>% 
                            .[,contig_id:=gsub('_[0-9]+$','',protein_id)] %>% 
                            .[,protein_set_size:=51],
                   fread(paste0(proteins,'.archaea.scg'),header=F,col.names=c('protein_id','protein_name')) %>% 
                            .[,protein_set:='.archaea'] %>% 
                            .[,contig_id:=gsub('_[0-9]+$','',protein_id)] %>% 
                            .[,protein_set_size:=38])
   stopifnot(nrow(scgTab)>0)
}


# Identify custom single copy genes
if(!is.null(arguments$customDbDir)){
   
   customSets <- unlist(strsplit(arguments$customDbDir,','))
   customSets <- unlist(strsplit(customDbDir,','))
   if(any(duplicated(basename(customSets)))){
      write.log('Custom SCG set directory names are not unique', 
                filename = logFile,append = T,write_to_file = T,type = 'error')
   }
   
   for(scgSet in customSets){
      # TODO: check if all files are present, hmms are extracted, hmmpress was applied
      
      # scgSet <- customSets[1]
      setName <- basename(scgSet)
      
      # determine set size
      scgSetSize <- fread(paste0(scgSet,'/genes.txt'),header = T,sep='\t') %>% 
         nrow()
      
      # determine noise cutoff
      scgNoiseCutoff <- fread(paste0(scgSet,'/noise_cutoff_terms.txt'),header = F,sep='\t')$V1
      
      # run HMMER on custom set
      hmmscan_out <- paste0(arguments$outputbasename,'_',setName,'_hmmscan.tblout')
      system(paste0('hmmscan --tblout ', hmmscan_out ,' ',scgNoiseCutoff,' --cpu ',arguments$threads,' "',
                    scgSet,'/genes.hmm ','" "',proteins,'">/dev/null 2>&1'))
      
      # parse HMMER result
      # TODO: account for empty hmmscan results
      hmmscan_tab <- read_tblout(hmmscan_out) %>%
         data.table() %>% 
         setnames(old = c('domain_name','query_name'),new=c('protein_name','protein_id')) %>% 
         .[.[, .I[which.min(sequence_evalue)], by=protein_id]$V1]
      scgTab <- rbind(scgTab,
                      hmmscan_tab[,.(protein_id,
                                     protein_name,
                                     protein_set=setName,
                                     contig_id=gsub('_[0-9]+$','',protein_id),
                                     protein_set_size=scgSetSize)])
      
   }
   
   
}


##
## Run bin selection
##
write.log('Dereplicating, aggregating, and scoring bins',filename = logFile,append = T,write_to_file = T,type = 'cat')
bin_evaluations <- cherry_pick(binTab=binTab,
                                scgTab=scgTab,
                                contigTab=contigTab,
                                output_basename=arguments$outputbasename,
                                score_threshold=as.numeric(arguments$score_threshold),
                                duplicate_penalty=as.numeric(arguments$duplicate_penalty),
                                megabin_penalty=as.numeric(arguments$megabin_penalty),
                                write_unbinned=arguments$write_unbinned,
                                write_bin_evals=arguments$write_bin_evals)


##
## Extract bins
##
if(arguments$write_bins){
   write.log('Writing bins',filename = logFile,append = T,write_to_file = T,type = 'cat')
   binDir <- paste0(dirname(arguments$outputbasename),'/DASTool_bins')
   if(!dir.exists(binDir)){
      system(paste0('mkdir ',binDir))
   }
   system(paste0('bash ', scriptDir,'/src/extract_bins.sh ',arguments$outputbasename,'_DASTool_contig2bin.tsv ',arguments$contigs,' ',binDir))
}


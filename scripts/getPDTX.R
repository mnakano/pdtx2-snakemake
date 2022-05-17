#!/usr/bin/env Rscript
library(PharmacoGx)
library(CoreGx)
library(SummarizedExperiment)
library(data.table)
library(Biobase)
library(reshape2)
library(dplyr)

## temporarily adding for mapping, once we rerun the beadarray normalization file, we can move it earlier in pipeline
if(!require("illuminaHumanv4.db", quietly=TRUE)){
  BiocManager::install("illuminaHumanv4.db", ask=FALSE)
}
library(illuminaHumanv4.db)
# library(biocompute)

# options(stringsAsFactors=FALSE)
# print("Retrieving selection")
# args = commandArgs(trailingOnly=TRUE)

# processedDataGithubPrefix <- "/pfs/downloadPDTXBreastData"
# sensitivityDataPrefix <- "/pfs/normalizeAndComputePDTXSens"
# annotationRepoPrefix <- '/pfs/downAnnotations'
# out_dir <- ''

# processedDataGithubPrefix <- "~/Documents/pfs/downloadBreatPDTXData/"
# sensitivityDataPrefix <- "~/Documents/pfs/normalizeAndComputePDTXSens/"
# annotationRepoPrefix <- '~/Documents/pfs/annotation/'

# rnaseq_select <- args
# print(rnaseq_select)
# rnaseq_results <- list()
# ORCESTRA_ID = tail(rnaseq_select, n=1)

# cnv_select <-  grep('cnv', rnaseq_select)
# mutation_select <-  grep('mutation', rnaseq_select)
# microarray_select <-  grep('microarray', rnaseq_select)
# microarray_v_select <- tail(args,n=3)[2]
# fusion_select <-  grep('fusion', rnaseq_select)

# tools <- grep(pattern = 'Kallisto|Salmon', x = rnaseq_select)
# tools <- rnaseq_select[tools]
# tools <- gsub("-", "_", tools)
# transcriptome <- grep(pattern = 'Gencode|Ensembl', x = rnaseq_select)
# transcriptome <- rnaseq_select[transcriptome]
# tool_path = expand.grid(a = tools,b = transcriptome)
# tool_path = paste0(tool_path$a, "_",tool_path$b)

# print(tool_path)

# version <- head(args)[3]
# drug_version <- head(args)[4]
# print(version)
# print(drug_version)

# standardize <- args[grep("filtered", args)]

#standardize drug concentration range function
standardizeRawDataConcRange <- function(sens.info, sens.raw){
  unq.drugs <- unique(sens.info$drugid)
  
  conc.m <- data.table(melt(sens.raw[,,1], as.is=TRUE))
  conc.m[,drugid := sens.info$drugid[match(Var1, rownames(sens.info))]]
  conc.ranges <- conc.m[,.(l = min(value, na.rm=T), r = max(value, na.rm=T)), c("drugid", "Var1")]
  conc.ranges[,Var1 := NULL]
  conc.ranges <- conc.ranges[,unique(.SD), drugid]	
  # conc.ranges[,N := .N, drugid]
  conc.ranges.disj <- conc.ranges[, {sq <- sort(unique(c(l,r))); 
  l = sq[seq(1,length(sq)-1)];
  r = sq[seq(2,length(sq))];
  .(l=l,r=r)}, drugid]
  ## Function below returns all consecutive ranges of ints between 1 and N
  returnConsInts <- function(N) {
    stopifnot(N>0)
    unlist(sapply(seq(1,N), function(ii) return(sapply(seq(ii, N), function(jj) return(seq(ii,jj))))), recursive=FALSE)
  }
  rangeNoHoles <- function(indicies, lr.tbl){
    if(length(indicies) == 1) return(TRUE)
    sq <- seq(indicies[1], indicies[length(indicies)]-1)
    all(lr.tbl[["l"]][sq+1] <= lr.tbl[["r"]][sq])
  }
  per.drug.range.indicies <- sapply(conc.ranges.disj[,.N,drugid][,N], returnConsInts)
  
  names(per.drug.range.indicies) <- conc.ranges.disj[,unique(drugid)] ## checked this: conc.ranges.disj[,.N,drugid][,drugid] == conc.ranges.disj[,unique(drugid)]
  
  
  # Check if there are any holes in the chosen range combination
  per.drug.range.indicies <- sapply(names(per.drug.range.indicies), function(drug){
    
    lr.tbl <- conc.ranges.disj[drugid == drug]
    per.drug.range.indicies[[drug]][sapply(per.drug.range.indicies[[drug]], rangeNoHoles, lr.tbl = lr.tbl)]
    
  })
  per.drug.range.indicies.2 <- sapply(names(per.drug.range.indicies), function(drug){
    
    lr.tbl <- conc.ranges.disj[drugid == drug]
    res <- t(sapply(per.drug.range.indicies[[drug]], function(x) return(c(lr.tbl[x[1],l], lr.tbl[x[length(x)],r]))))
    colnames(res) <- c("l", "r")
    res <- data.frame(res)
    res <- cbind(drugid = drug, res)
  }, simplify=FALSE)
  per.drug.range.indicies.dt <- rbindlist(per.drug.range.indicies.2)
  
  conc.ranges <- conc.m[,.(l = min(value, na.rm=T), r = max(value, na.rm=T)), c("drugid", "Var1")]
  setkey(conc.m, Var1)
  conc.m <- na.omit(conc.m)
  setkey(conc.m, drugid, Var1, value)
  setkey(conc.ranges, drugid, l, r)
  # tic()
  ## NOTE:: Data.table used for maximum speed. Probably possible to do this more intelligently by 
  ## NOTE:: being aware of which conditions overlap, but its fast enough right now as it is.
  chosen.drug.ranges <- lapply(unq.drugs, function(drug){
    num.points.in.range <- apply(per.drug.range.indicies.dt[drugid==drug, .(l,r)], 1, function(rng){
      conc.m[drugid==drug][conc.ranges[drugid==drug][l<=rng["l"]][r>=rng["r"],Var1], on="Var1"][value >= rng["l"]][value <= rng["r"],.N]
      # conc.m[drugid==drug][, Var1]
    })
    max.ranges <- per.drug.range.indicies.dt[drugid==drug][which(num.points.in.range==max(num.points.in.range))]
    max.ranges[which.max(log10(r) - log10(l)), ]
  })
  # toc()
  names(chosen.drug.ranges) <- sapply(chosen.drug.ranges, `[[`, "drugid")
  removed.experiments <- unlist(lapply(unq.drugs, function(drug){
    rng <- unlist(chosen.drug.ranges[[drug]][,.(l,r)])
    exp.out.range <- conc.ranges[drugid==drug][l>rng["l"] | r<rng["r"],Var1]
    return(exp.out.range)
  }))
  
  sens.raw[removed.experiments,,] <- NA_real_
  conc.ranges.kept <- conc.ranges[!Var1 %in% removed.experiments]
  
  for(drug in unq.drugs){
    rng <- unlist(chosen.drug.ranges[[drug]][,.(l,r)])
    myx <- conc.ranges.kept[drugid==drug,Var1]
    doses <- sens.raw[myx, ,"Dose"]
    which.remove <- (doses < rng["l"] | doses > rng["r"])
    sens.raw[myx, ,"Dose"][which(which.remove,arr.ind=TRUE)] <- NA_real_
    sens.raw[myx, ,"Viability"][which(which.remove,arr.ind=TRUE)] <- NA_real_
    
    ## Annotate sens info with chosen range
    sens.info[sens.info$drugid==drug,"chosen.min.range"] <- rng["l"]
    sens.info[sens.info$drugid==drug,"chosen.max.range"] <- rng["r"]
  }
  sens.info$rm.by.conc.range <- FALSE
  sens.info[removed.experiments,"rm.by.conc.range"] <- TRUE
  
  return(list("sens.info" = sens.info, sens.raw = sens.raw))
}


#filter noisy curves from PSet (modified function to take into account standardized conc range)
filterNoisyCurves2 <- function(pSet, epsilon=25 , positive.cutoff.percent=.80, mean.viablity=200, nthread=1) {
  acceptable <- mclapply(rownames(sensitivityInfo(pSet)), function(xp) {
    #for(xp in rownames(sensitivityInfo(pSet))){
    drug.responses <- as.data.frame(apply(pSet@sensitivity$raw[xp , ,], 2, as.numeric), stringsAsFactors=FALSE)
    if (!all(is.na(drug.responses))){
      
      
      drug.responses <- drug.responses[complete.cases(drug.responses), ]
      doses.no <- nrow(drug.responses)
      drug.responses[,"delta"] <- .computeDelta(drug.responses$Viability)
      
      delta.sum <- sum(drug.responses$delta, na.rm = TRUE)
      
      max.cum.sum <- .computeCumSumDelta(drug.responses$Viability)
      
      if ((table(drug.responses$delta < epsilon)["TRUE"] >= (doses.no * positive.cutoff.percent)) &
          (delta.sum < epsilon) &
          (max.cum.sum < (2 * epsilon)) &
          (mean(drug.responses$Viability) < mean.viablity)) {
        return (xp)
      }
    }
    
  }, mc.cores=nthread)
  acceptable <- unlist(acceptable)
  noisy <- setdiff(rownames(sensitivityInfo(pSet)), acceptable)
  return(list("noisy"=noisy, "ok"=acceptable))
}

.computeDelta <- function(xx ,trunc = TRUE) {
  xx <- as.numeric(xx)
  if(trunc)
  {
    return(c(pmin(100, xx[2:length(xx)]) - pmin(100, xx[1:length(xx)-1]), 0))
  }else{
    return(c(xx[2:length(xx)] - xx[1:length(xx)-1]), 0)
  }
}

#' @importFrom utils combn
.computeCumSumDelta <- function(xx, trunc = TRUE) {
  xx <- as.numeric(xx)
  if(trunc) {
    xx <- pmin(xx, 100)
  }
  tt <- t(combn(1:length(xx), 2 , simplify = TRUE))
  tt <- tt[which(((tt[,2] - tt[,1]) >= 2) == TRUE),]
  if (is.null(nrow(tt))){
    tt <- matrix(tt, ncol = 2)
  }
  cum.sum <- unlist(lapply(1:nrow(tt), function(x){xx[tt[x,2]]-xx[tt[x,1]]}))
  return(max(cum.sum))
}   

matchToIDTable <- function(ids,tbl, column, returnColumn="unique.cellid") {
  sapply(ids, function(x) {
    myx <- grep(paste0("((///)|^)",Hmisc::escapeRegex(x),"((///)|$)"), tbl[,column])
    if(length(myx) > 1){
      stop("Something went wrong in curating ids, we have multiple matches")
    }
    if(length(myx) == 0){return(NA_character_)}
    return(tbl[myx, returnColumn])
  })
}



.converteSetToSE <- function(eSets) {
  
  SEfinal <- lapply(eSets, function(eSet) {
        if (grepl("^rna$", Biobase::annotation(eSet))) {
            rownames(eSet) <- Biobase::fData(eSet)$EnsemblGeneId
        }
        SE <- SummarizedExperiment::SummarizedExperiment(assays = SimpleList(as.list(Biobase::assayData(eSet))),
            rowData = S4Vectors::DataFrame(Biobase::fData(eSet),
                rownames = rownames(Biobase::fData(eSet))), colData = S4Vectors::DataFrame(Biobase::pData(eSet),
                rownames = rownames(Biobase::pData(eSet))), metadata = list(experimentData = eSet@experimentData,
                annotation = Biobase::annotation(eSet), protocolData = Biobase::protocolData(eSet)))
        assayNames(SE) <- assayDataElementNames(eSet)
        mDataType <- Biobase::annotation(eSet)
        molecularProfilesSlot(cSet)[[mDataType]] <- SE
    })
  #setNames(pSet@molecularProfiles, names(eSets))
  return(SEfinal)
}

args <- commandArgs(trailingOnly = TRUE)
processedDataGithubPrefix <- args[1]
sensitivityDataPrefix <- args[2]
annotationRepoPrefix <- args[3]
out_dir <- args[4]
filename <- args[5]

standardize <- args[grep("filtered", args)]

## If you want to run this, update to correct path:
drug.annotation.all <- read.csv(file.path(annotationRepoPrefix,"drugs_with_ids.csv"))


cell <- read.csv(file.path(processedDataGithubPrefix,"cell.csv"), header=TRUE)
drug <- read.csv(file.path(processedDataGithubPrefix,"raw_drug.csv"), header=TRUE)


drug$unique.drugid <- matchToIDTable(drug$DRUG_NAME, drug.annotation.all, "PDTX.drugid", "unique.drugid")

rownames(drug) <- drug$unique.drugid
rownames(cell) <- cell$unique.cellid

tmp <- read.csv(file.path(processedDataGithubPrefix,"cell_annotation_all.csv"), header=TRUE)
curationCell <- read.csv(file.path(processedDataGithubPrefix,"cell_annotation_all.csv"), header=TRUE, row.names=1)
curationDrug <- drug.annotation.all[,c("unique.drugid", "PDTX.drugid")]
curationDrug <- curationDrug[complete.cases(curationDrug),]
rownames(curationDrug) <- curationDrug$unique.drugid


info <- readRDS(file.path(sensitivityDataPrefix,"info.rds"))
profiles <- readRDS(file.path(sensitivityDataPrefix,"profiles.rds"))
raw <- readRDS(file.path(sensitivityDataPrefix,"raw.rds"))

## Map Info Table 

info$drugid <- matchToIDTable(info$drugid, drug.annotation.all, "PDTX.drugid", "unique.drugid")

eSet <- readRDS(file.path(processedDataGithubPrefix,"final_eset.Rda"))
Biobase::annotation(eSet) <- "rna"

## adding mapping to ENSG

# list.files(annotationRepoPrefix)

# load(file.path(annotationRepoPrefix, "Gencode.v33.annotation.RData"))

# symbol1.matches <- match(Biobase::fData(eSet)$symbol, features_gene$gene_name)
# symbol2.matches <- match(Biobase::fData(eSet)$symbol2, features_gene$gene_name)


## checking the mismatches 

# Biobase::fData(eSet)[which(!symbol1.matches == symbol2.matches),]

ids <- Biobase::fData(eSet)$IlluminaID
ensg2 <- unlist(mget(as.character(ids), illuminaHumanv4ENSEMBLREANNOTATED, ifnotfound=NA))
ensg <- unlist(mget(as.character(ids), illuminaHumanv4ENSEMBL, ifnotfound=NA))



Biobase::fData(eSet)$EnsemblGeneId <- ensg2

# eSet <- eSet[!is.na(ensg2),]


# z <- .converteSetToSE(list("rna"=eSet))


final_PSet <- PharmacoSet(name="BreastPDTX",
                    molecularProfiles=list(rna=as(eSet, "SummarizedExperiment")),
                    cell=cell,
                    drug=drug,
                    sensitivityInfo=info,
                    sensitivityRaw=raw,
                    sensitivityProfiles=profiles,
                    curationCell=curationCell,
                    curationDrug=curationDrug,
                    datasetType="sensitivity",
                    verify=TRUE)




if (length(standardize) > 0){
  
  noisy_out <- filterNoisyCurves2(final_PSet)
  print("filter done")
  final_PSet@sensitivity$profiles[noisy_out$noisy, ] <- NA
  
} else {
  print("unfiltered PSet")
  
}

final_PSet@annotation$version <- 2     


saveRDS(final_PSet, file=paste0(out_dir, filename))

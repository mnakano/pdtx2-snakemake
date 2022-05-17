library(PharmacoGx)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input.dir <- args[1]
output.dir <- args[2]

# input.dir <- "/pfs/downloadPDTXPublished/"
# output.dir <- "/pfs/out/"

# input.dir <- "~/Documents/pfs/downloadPDTXPublished/"
# output.dir <- "~/Documents/pfs/normalizeAndComputePDTXSens/"

raw_drug <- data.frame(fread(file.path(input.dir,"RawDataDrugsSingleAgents.txt")))
raw_drug <- raw_drug[, -c(1, 3, 15, 16)]
raw_drug <- raw_drug[order(raw_drug$DRUG_ID), ]

x <- aggregate(raw_drug[, 3:7],
                    by=list(raw_drug$ID, raw_drug$DRUG_ID),
                    FUN=median) # same thing as removing duplicates but this way we make sure the rownames match
                                # exactly with the rownames in profiles and everything else

x$num <- c()
for (i in 1:nrow(x)) {
  x$num[i] <- length(unique(x[i, 3:7]))
}

cellid <- paste(x$Group.1)


info <- data.frame(cellid=cellid, drugid=x$Group.2, nbr.conc.tested=x$num,
                  Dose1.uM=x$D1_CONC, Dose2.uM=x$D2_CONC, Dose3.uM=x$D3_CONC,
                  Dose4.uM=x$D4_CONC, Dose5.uM=x$D5_CONC)
rownames(info) <- paste("drugid", info$drugid, info$cellid, sep="_")

saveRDS(info, file=file.path(output.dir,"info.rds"))

raw_drug <-  data.frame(fread(file.path(input.dir,"RawDataDrugsSingleAgents.txt")))
raw_drug <- raw_drug[, -c(1, 3)]
raw_drug <- raw_drug[order(raw_drug$DRUG_ID), ]
# Normalize raw intensity values to produce viability values
raw_drug[, 8:12] <- (raw_drug[, 8:12] - raw_drug$Blank)/(raw_drug$Control - raw_drug$Blank) * 100

# info <- readRDS("~/Desktop/BreastPDTX/data/results/sensitivity/info.rds")

med_intensity <- aggregate(raw_drug[, 8:12],
                          by=list(raw_drug$ID, raw_drug$DRUG_ID),
                          FUN=median)
med_intensity <- cbind(med_intensity, info[, 4:8])

rownames(med_intensity) <- paste("drugid", med_intensity$Group.2, med_intensity$Group.1, sep="_")
saveRDS(med_intensity, file.path(output.dir,"dose_viability.rds"))

dose_viability <- med_intensity
dose_viability <- dose_viability[, -c(1, 2)]

dose <- as.matrix(dose_viability[, 6:10])
viability <- as.matrix(dose_viability[, 1:5])

dose <- dose[,seq(ncol(dose), 1)]
viability <- viability[,seq(ncol(viability), 1)]


raw <- array(c(dose, viability), dim=c(2550, 5, 2),
            dimnames=list(rownames(dose_viability),
                          sprintf("doses%d", seq(1, 5)),
                          c("Dose", "Viability")))

saveRDS(raw, file=file.path(output.dir,"raw.rds"))

# ic50_recomputed <- c()
# for (i in 1:nrow(med_intensity)) {
#   zz <- computeIC50(med_intensity[i, 8:12], med_intensity[i, 3:7])
  
#   ic50_recomputed <- rbind(ic50_recomputed, zz)
# }

# auc_recomputed <- c()
# for (i in 1:nrow(med_intensity)) {
#   zz <- computeAUC(med_intensity[i, 8:12], med_intensity[i, 3:7])
  
#   auc_recomputed <- rbind(auc_recomputed, zz)
# }


recomputed <- PharmacoGx:::.calculateFromRaw(raw)


parTable <- do.call(rbind,recomputed[[3]])
# print(head(rownames(parTable)))
# print(str(recomputed[[3]]))
profiles <- data.frame("aac_recomputed" = as.numeric(unlist(recomputed[[1]]))/100, 
        "ic50_recomputed" = as.numeric(unlist(recomputed[[2]])), 
        "HS" = as.numeric(unlist(parTable[,1])),
        "E_inf" = as.numeric(unlist(parTable[,2])),
        "EC50" = as.numeric(unlist(parTable[,3]))) 

rownames(profiles) <- rownames(raw)


published <- data.frame(fread(file.path(input.dir,"DrugResponsesAUCSamples.txt")))
rownames(published) <- paste("drugid", published$Drug, published$ID, sep="_")
stopifnot(setequal(rownames(published), rownames(info)))

published <- published[rownames(info), ]

profiles$auc_published <- published$AUC

profiles$ic50_published <- as.numeric(published$iC50)

saveRDS(profiles, file=file.path(output.dir,"profiles.rds"))
# QTL mapping of the a chicken intercross: wild-type x white leghorn, F8

# Run R/4.1.1

# = variables ================================================================ #


wd <- "" # CHANGE THIS TO THE PROJECT PATH

# only keep QTL w/ lod > 3, thresholds are appiled later
thrSex <- 3
thrNonsex <- 3
minprob <- 0.8

# = libraries ================================================================ #

library(qtl2)
library(plyr)

# = functions ================================================================ #

get_args <- function(){
  args <- commandArgs(TRUE)
  
  if (length(args) == 0) print("No arguments given")
  
  if (length(args) > 0){
    for (i in seq_along(args)){
      # add quotes to args or the parse will not work
      VAR <- sub("=", "='", args[i])  # behind =
      VAR <- sub("$", "'", VAR)       # last char
      
      eval(parse(text = VAR), envir = .GlobalEnv)
    }
    print(paste("num =", num))
  }
}
get.markers <- function(dat) {
  markers <- find_marker(pseudomarker_map, chr = as.vector(dat$chr), interval = c(as.vector(dat$ci_lo), as.vector(dat$ci_hi)))

  d <- grep("^c", markers)
  if (length(d) > 0) markers <- markers[-grep("^c", markers)]

  markerpos <- find_markerpos(cross, markers = markers)

  if (length(markers) > 2) {
    start <- markers[1]
    end <- markers[length(markers)]
    # pos <- markers[2]
    pos <- markers[which.min(abs(dat$pos - markerpos$gmap))]
  }
  if (length(markers) == 2) {
    start <- markers[1]
    end <- markers[2]

    x <- c(as.vector(dat$ci_lo), as.vector(dat$ci_hi))
    pos <- markers[which.min(abs(x - as.vector(dat$pos)))]
  }
  if (length(markers) == 1) {
    start <- markers[1]
    end <- markers[1]
    pos <- markers[1]
  }
  l1 <- list()
  l1$start <- start
  l1$end <- end
  l1$pos <- pos
  l1
}
make.table <- function(dat) {
  markers <- get.markers(dat)

  out <- dat[, 3:8]
  out$phenotype <- methbin
  out$marker <- markers$pos
  out$ci_lo_marker <- markers$start
  out$ci_hi_marker <- markers$end
  out$ci <- paste0("chr", out$chr, ":", out$ci_lo, "-", out$ci_hi)
  out$qtlchr <- methbinChrom

  out[, c("phenotype", "marker", "lod", "ci", "sexint", "chr", "qtlchr", "pos", "ci_lo", "ci_hi", "ci_lo_marker", "ci_hi_marker")]
}
calc.probabilites <- function(d, minprob) {
  # Account for the probability of the imputed genotype with the maxmarg-function.
  # This will only keep the QTL that furfill minprob value (80%) 0.8 * 124 = 100 w/o NAs

  if (d[3] == "X") {
    tmp <- maxmarg(prob, pseudomarker_map, chr = d[3], pos = as.numeric(d[4]), minprob = minprob)
  } else {
    tmp <- maxmarg(prob, pseudomarker_map, chr = d[3], pos = as.numeric(d[4]))
  }
  length(which(is.na(as.vector(tmp))))
}

# = code ===================================================================== #

# load datasets
#load(paste0(wd, "rdata/f8_meth_hypo_cross.Rdata")) # 1050176 phenotypes
cross_url <- url("f8_meth_hypo_cross.Rdata")

load(cross_url)
close(cross_url)

# Subdivided between 20 runs
phenos <- colnames(cross$pheno)

get_args() # gets command line arguments from R session starter

#max <- round(length(phenos) / 20)
#range <- split(phenos, ceiling( seq_along(phenos) / max) )
#phenos <- range[[num]]
#print(paste0("First: ", phenos[1], " | Last: ", phenos[length(phenos)]))


# psuedomarkers between my markers into the genetic map
pseudomarker_map <- insert_pseudomarkers(cross$gmap, step = 1)

# genotype probability in all of the markers
prob <- calc_genoprob(cross, pseudomarker_map, error_prob = 0.001, quiet = F)

# allele probability needed if genome scan of additive allele model is to be performed
allele_prob <- genoprob_to_alleleprob(prob)

# Sex choromosome covar
Xcovar <- get_x_covar(cross)

# covar
# pull phenotypes and covariates; ensure that covariates have names attribute
sex <- match(cross$covar$sex, c("male", "female"))
sex <- sub(1, 0, sex)
sex <- sub(2, 1, sex)
sex <- as.numeric(sex)
batch <- as.numeric(cross$covar$batch)
risky <- as.numeric(cross$covar$risky)
covar <- cbind(sex, batch, risky)
rownames(covar) <- rownames(cross$covar)

intcovar <- sex
names(intcovar) <- rownames(cross$covar)

l1 <- list()
for (i in seq_along(phenos)){
  methbin <- phenos[i]
  methbinChrom <- sub("chr", "", strsplit(methbin, split = "_")[[1]][1])

  scan1sex <- scan1(
    genoprobs = prob, 
    pheno = cross$pheno[, methbin], 
    addcovar = covar, 
    intcovar = intcovar, 
    Xcovar = Xcovar
  )
  
  scan1nonsex <- scan1(
    genoprobs = prob, 
    pheno = cross$pheno[, methbin], 
    addcovar = covar, 
    Xcovar = Xcovar
  )

  peaksSex <- find_peaks(
    scan1sex, 
    pseudomarker_map, 
    threshold = thrSex, 
    drop = 1.8
  )

  peaksNonsex <- find_peaks(
    scan1nonsex, 
    pseudomarker_map, 
    threshold = thrNonsex, 
    drop = 1.8
  )

  if (nrow(peaksSex) > 0) peaksSex <- peaksSex[apply(peaksSex, 1, calc.probabilites, minprob) < 20, ]
  if (nrow(peaksNonsex) > 0) peaksNonsex <- peaksNonsex[apply(peaksNonsex, 1, calc.probabilites, minprob) < 20, ]

  sex1 <- ifelse(nrow(peaksSex) > 0, TRUE, FALSE)
  nonsex1 <- ifelse(nrow(peaksNonsex) > 0, TRUE, FALSE)

  if (sex1 | nonsex1) {
    if (sex1 & nonsex1) {
      peaksSex$sexint <- "sex"
      tmp1 <- subset(peaksSex, lod == max(peaksSex$lod))

      peaksNonsex$sexint <- "nonsex"
      tmp2 <- subset(peaksNonsex, lod == max(peaksNonsex$lod))

      tmp <- rbind(tmp1, tmp2)
    } else {
      if (sex1) {
        peaksSex$sexint <- "sex"
        tmp <- subset(peaksSex, lod == max(peaksSex$lod))
      }
      if (nonsex1) {
        peaksNonsex$sexint <- "nonsex"
        tmp <- subset(peaksNonsex, lod == max(peaksNonsex$lod))
      }
    }
    out <- adply(tmp, 1, make.table)

    # Only continue with sex chromosome mQTL
    out <- subset(out, chr == "X" | qtlchr == "Z")

    if (nrow(out) == 0) next

    if (nrow(out) > 1){
      # Identify if sexint or non-sexint is to keep, compare lod scores
      sex1 <- subset(out, sexint %in% "sex")
      nonsex1 <- subset(out, sexint %in% "nonsex")
      out <- if (sex1$lod - 1 > nonsex1$lod) sex1 else nonsex1
    }
    
    # extract the QTL effects
    chrom <- out$chr
    sexint <- out$sexint
    pos <- out$pos

    fit <- fit1(
      genoprobs = pull_genoprobpos(prob, pseudomarker_map, chr = chrom, pos = pos),
      pheno = cross$pheno[, methbin],
      addcovar = covar,
      #intcovar = intcovar,
      intcovar = if (sexint == "sex") intcovar else NULL,
      se = TRUE,

      # add additive and dominante effects
      #contrasts = cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0))
      contrasts = cbind(
        mu = c(1,1,1,1,1), 
        a = c(-1, 0, 1, 0, 0), 
        d = c(0, 1, 0, 0, 0),
        AX = c(0, 0, 0, 1, -1),
        BX = c(0, 0, 0, -1 ,1))
    )
    out <- cbind(out, coef=t(fit$coef), SE=t(fit$SE))
    #se <- data.frame(attributes(coef1)["SE"])
    #se <- se[out$marker, ]
    #out <- cbind(out, t(coef1[out$marker, ]), se)

    # Save output
    if (nrow(out) > 1) l1[[i]] <- out[, 3:ncol(out)]
  }

  print(paste0(Sys.time(), " | scan # ", i, " done..."))
}

d <- Reduce(function(x, y) rbind(x, y), l1)

name <- paste0("methQTL_scan", num, ".Rds")
filename <- file.path(wd, "data", name)
saveRDS(d, file = filename)
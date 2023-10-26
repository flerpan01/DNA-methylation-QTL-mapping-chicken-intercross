# R/4.1.1

# Import the QTL-mapping files generated previous, bin/qtlscan.R
# Filter by thresholds
# Summaries into one file
# Plot all QTL for visual inspection, to remove outlier driven QTL

# = variables ================================================================ #

#setwd() # CHANGE THIS TO THE PROJECT PATH
cross_url <- url("https://www.dropbox.com/scl/fi/o4rnvbd6b1n0kd5ctbm6m/f8_meth_hypo_cross.Rdata?rlkey=aozwijnieg196alubm68uq01k&dl=1")

# = libraries ================================================================ #

library(qtl2)
library(ggplot2)

# = functions ================================================================ #

plotter <- function(d, pdfname, type) {
  col <- c("red", "tomato3", "yellow3", "snow3", "grey95")
  x_lims <- c("RJF/-", "RJF/RJF", "RJF/WL", "WL/WL", "WL/-")
  
  pdf(paste0(pdfname, ".pdf"))
  
  for (i in 1:nrow(d)) {
    x <- as.vector(maxmarg(prob, map, minprob = 0.8, chr = d$chr[i], pos = d$pos[i], return_char = TRUE))
    y <- as.vector(cross$pheno[, d$phenotype[i]])
    NAs <- !is.na(x) & !is.na(y)

    y <- y[NAs]
    x <- x[NAs]
    x <- sub("AA", "RJF/RJF", x)
    x <- sub("AB", "RJF/WL", x)
    x <- sub("BB", "WL/WL", x)
    x <- sub("AY", "RJF/-", x)
    x <- sub("BY", "WL/-", x)
    x <- factor(x)

    fac <- x_lims %in% levels(x)

    if(type=="point") {
        title <- paste0(d$type[i], " mQTL | chr", d$chr[i], "@", d$pos[i], "cM | LOD: ", round(d$lod[i], 2), "\n", "Methbin: ", d$phenotype[i], " | marker:", d$marker[i], "\n", "CI: ", d$ci[i], "cM", "\n", "i=",i)
        p <- ggplot() + 
          geom_jitter(aes(x=x, y=y, fill = x), shape = 21, width = 0.2, size = 3, color = "black") + 
          theme_classic() + 
          theme(legend.position="none") + 
          labs(title = title) + 
          scale_x_discrete(limits = x_lims[fac])

    }
    if(type=="box") {
        title <- paste0(d$type[i], " mQTL | chr", d$chr[i], "@", d$pos[i], "cM | LOD: ", round(d$lod[i], 2), "\n", "Methbin: ", d$phenotype[i], " | marker:", d$marker[i], "\n", "CI: ", d$ci[i], "cM")
        p <- ggplot() + geom_boxplot(aes(x=x, y=y), fill = col[fac]) + theme_classic() + theme(legend.position="none") + labs(x = "Genotype", y = "Normalised DNA methylation levels", title = title) + scale_x_discrete(limits = x_lims[fac])
    }

    print(p)

    show(i)
  }
  dev.off()
}

# = code ===================================================================== #

pmap <- read.csv(file.path("data", "f8_pmap_galgal6.csv"))

files <- list.files(
  file.path("data"),
  pattern = "methQTL_scan",
  full.names = TRUE
)

files <- data.frame(
  files = files,
  stringsAsFactors = FALSE
)

# import Rds
dat <- lapply(files$files, readRDS)
dat <- Reduce(function(x,y) rbind(x,y), dat)

# Filter out type of QTL; cis or trans
foo <- function(x) {
  col1 <- which(names(x) %in% "chr")
  col2 <- which(names(x) %in% "qtlchr")
  chr <- ifelse(x[col1] == "X", "Z", x[col1])
  qtlchr <- x[col2]
  ifelse(chr == qtlchr, "cis", "trans")
}
dat$type <- apply(dat[, c("chr", "qtlchr")], 1, foo)

# check distance of markers to assess if cis or trans
for (i in 1:nrow(dat)) {
  type <- dat$type[i]
  if (type %in% "cis") {
    methPos <- as.numeric(strsplit(dat$phenotype[i], split = "_")[[1]][2])
    m1 <- dat$ci_lo_marker[i]
    m2 <- dat$ci_hi_marker[i]

    pos1 <- subset(pmap, marker %in% m1)
    pos2 <- subset(pmap, marker %in% m2)

    dat$type[i] <- ifelse(pos1$pos <= methPos & methPos <= pos2$pos, "cis", "trans")
  }
  print(i)
}


# Threshold from permutation tests

# trans autosome
# sex sign; 7.6
# sex sugg; 6.31
# nonsex sign; 5.88
# nonsex sugg; 4.77

# trans sex chromsome (loci on sex chromosome)
# sex sign; 7.7
# sex sugg; 5.92
# nonsex sign; 7.68
# nonsex sugg; 5.93

# cis
# Sex sign; 5.73
# Sex sugg; 4.87
# Nonsex sign; 4.29
# Nonsex sugg; 3.58

sexcis <- subset(dat, sexint == "sex" & type == "cis" & lod >= 4.87)
nonsexcis <- subset(dat, sexint == "nonsex" & type == "cis" & lod >= 3.58)
sextrans <- subset(dat, sexint == "sex" & type == "trans" & lod >= 5.92)
nonsextrans <- subset(dat, sexint == "nonsex" & type == "trans" & lod >= 5.93)

dat <- rbind(sexcis, sextrans, nonsexcis, nonsextrans)
dat$lodthr <- "suggestive"

rows <- which(dat$type == "cis" & dat$sexint == "sex" & dat$lod >= 5.73)
dat$lodthr[rows] <- "significant"

rows <- which(dat$type == "cis" & dat$sexint == "nonsex" & dat$lod >= 4.29)
dat$lodthr[rows] <- "significant"

rows <- which(dat$type == "trans" & dat$sexint == "sex" & dat$lod >= 7.7)
dat$lodthr[rows] <- "significant"

rows <- which(dat$type == "trans" & dat$sexint == "nonsex" & dat$lod >= 7.68)
dat$lodthr[rows] <- "significant"

rownames(dat) <- NULL

# plot QTL to spot outliers
# also remove QTL w/ genotypes < 10 observations
# this step is done visual
# save the index of the QTL to be imported in the next script: outlier_check.R

load(cross_url)

# Make map
map <- insert_pseudomarkers(cross$gmap, step = 1)
# genotype probability in all of the markers
prob <- calc_genoprob(cross, map, error_prob = 0.001, cores = 1, quiet = F)

plotter(dat, file.path("tmp", "tmp"), "point")

saveRDS(dat, file = file.path("data", "f8_mqtl_table_chrz_unfiltered.Rds"))
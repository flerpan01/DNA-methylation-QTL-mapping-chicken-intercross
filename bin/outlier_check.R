# R/4.1.1

# After visual inspection for outliers 
# import text file with index (i) as identifier
# remove these from the dataset
# Also, QTL with < 10 observations per genotypes were removed

# = variables ================================================================ #

#setwd() # CHANGE THIS TO THE PROJECT PATH
text_file_with_outliers_index <- "methQTL_toRemove.txt"
cross_url <- url("https://www.dropbox.com/scl/fi/o4rnvbd6b1n0kd5ctbm6m/f8_meth_hypo_cross.Rdata?rlkey=aozwijnieg196alubm68uq01k&dl=1")

# = libraries ================================================================ #

library(qtl2)
library(ggplot2)
library(openxlsx)

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

dat <- readRDS(file.path("data", "f8_mqtl_table_chrz_unfiltered.Rds"))
toRemove <- read.table(file.path("tmp", text_file_with_outliers_index), header=T)
toRemove <- unname(unlist(toRemove))

d <- dat[-toRemove, ]

load(cross_url)

# Make map
map <- insert_pseudomarkers(cross$gmap, step = 1)
# genotype probability in all of the markers
prob <- calc_genoprob(cross, map, error_prob = 0.001, cores = 1, quiet = F)

# plot the solid QTL
plotter(d, paste0("img/mQTL_plots1.0"), "box")
write.csv(d, file = file.path("data", "f8_mqtl_table_chrz.csv"), row.names=F, quote=F)
write.xlsx(d, file = file.path("data", "mQTL_table1.0.xlsx"), row.names=F)
Filtre_GMODatabank <- function (host, gmoDB, gmoDB_Output) {

	CDShost <- read.csv(host, fill = FALSE, header = TRUE, sep = ",", row.names = 1)
	CDSgmoDB <- read.csv(gmoDB, fill = FALSE, header = TRUE, sep = ",", row.names = 1)
	
	Medl4m2 <- median(CDShost[,1], na.rm = FALSE)
	Medl3m1 <- median(CDShost[,6], na.rm = FALSE)
	Medl9m7 <- median(CDShost[,7], na.rm = FALSE)

	NewGmoDatabank <- subset(CDSgmoDB, BrayCurtis_l4_m2Prop > Medl4m2 | CodonUsage > Medl3m1 | BrayCurtis_l9_m7Freq > Medl9m7)
	 if (dim(NewGmoDatabank)[1] == 0) {
                NewGmoDatabank <- CDSgmoDB
        }
	write.table(NewGmoDatabank, file=gmoDB_Output, sep=",", quote=FALSE)
}

Filtre_GMODatabank(snakemake@input[[1]], snakemake@input[[2]], snakemake@output[[1]])

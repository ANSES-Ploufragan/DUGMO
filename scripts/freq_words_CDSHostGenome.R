#### Normalization by the sum of the counts - calculation of frequencies for the CDSs of host genome
#### Selection of the most over-represented or under-represented words by percentage

normalisation <- function(dataRmes, pourcentSurRep, pourcentSousRep, freqGenRef) {
    
	#Selection only a percentage of the total words corresponding to over-represented or under-represented words
	DataRmes <- read.table(dataRmes, sep="|", header=TRUE, fill=FALSE)
	nblignesSurRep <- trunc(as.numeric(pourcentSurRep)/100*nrow(DataRmes)) 
	nblignesSousRep <- trunc(as.numeric(pourcentSousRep)/100*nrow(DataRmes))
	DataSurRep <- tail(DataRmes, nblignesSurRep)
	DataSousRep <- head(DataRmes, nblignesSousRep)
	Data <- rbind(DataSousRep, DataSurRep)
	#Frequencies calcultation
	sumcomptage <- apply(as.matrix(Data[,2]), 2, sum)
	freqmot <- lapply(as.numeric(Data[,2]), function(x) x/(sumcomptage))
	#Output file
	write.table(cbind(as.matrix(Data[,1]), round(as.numeric(freqmot), 9), as.matrix(Data[,5]), as.matrix(Data[,2])), quote = FALSE, file=freqGenRef, row.names=FALSE, col.names=c("mot", "frequence", "score", "count"), sep="\t")
}

normalisation(snakemake@input[[1]], snakemake@params[["pourcentSurRep"]], snakemake@params[["pourcentSousRep"]], snakemake@output[[1]])

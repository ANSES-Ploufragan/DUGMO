ML <- function(cds, gmo, t25, PredRFPos, PredRFProba, PredRLPos, PredRLProba) {

# Import learning data
cds <- read.csv(cds, fill = FALSE, header = TRUE, sep = ",", row.names = 1)
gmo <- read.csv(gmo, fill = FALSE, header = TRUE, dec = ".", sep = ",", row.names = 1)
cds[, 10] <- 0
gmo[, 10] <- 1
colnames(cds) <- c("BrayCurtis_l4m2Prop", "Length", "MeanScore_l4m2Prop", "DensityNuc_l4m2Prop", "GCpourcent_l4m2Prop", "BrayCurtis_l3m1Prop", "BrayCurtis_l9m7Freq", "MeanScore_l9m7Freq", "DensityNuc_l9m7Freq", "OGM")
colnames(gmo) <- c("BrayCurtis_l4m2Prop", "Length", "MeanScore_l4m2Prop", "DensityNuc_l4m2Prop", "GCpourcent_l4m2Prop", "BrayCurtis_l3m1Prop", "BrayCurtis_l9m7Freq", "MeanScore_l9m7Freq", "DensityNuc_l9m7Freq", "OGM")

gmo$OGM <- as.factor(gmo$OGM)
cds$OGM <- as.factor(cds$OGM)
Data      <- rbind.data.frame(cds, gmo)
Data$OGM  <- as.factor(Data$OGM)
print(summary(Data))

#packages
library(caret)
library(datasets)
library(randomForest)

# import prediction data
Data.new <- read.csv(t25, fill = FALSE, header = TRUE, sep = ",", row.names = 1)
Data.new <- Data.new[,1:9] 
colnames(Data.new) <- c("BrayCurtis_l4m2Prop", "Length", "MeanScore_l4m2Prop", "DensityNuc_l4m2Prop", "GCpourcent_l4m2Prop", "BrayCurtis_l3m1Prop", "BrayCurtis_l9m7Freq", "MeanScore_l9m7Freq", "DensityNuc_l9m7Freq")

# Create model with learning data + prediction using Random Forest
set.seed(1234)
mod.rf  <- train(OGM ~., data=Data, method='rf', preProcess = c("center", "scale"), tuneLength = 12)
pred.rf <- predict(mod.rf, Data.new, "prob")
res <- as.data.frame(pred.rf)
res.rf <- data.frame(Data.new, res)
new.res.rf <- res.rf[ which(res.rf$X1 > 0.50), ]
write.table(new.res.rf, file=PredRFPos, sep=",")
write.table(res.rf, file=PredRFProba, sep=",")

# Create model with learning data + prediction using régression logistique
mod.reglog <- train(OGM ~., data=Data, method='glm', family = binomial, preProcess = c("center", "scale"), tuneLength = 6)
pred.reglog <- predict(mod.reglog, Data.new, "prob")
res <- as.data.frame(pred.reglog)
res.reglog <- data.frame(Data.new, res)
new.res.reglog <- res.reglog[ which(res.reglog$X1 > 0.50), ]
write.table(new.res.reglog, file=PredRLPos, sep=",")
write.table(res.reglog, file=PredRLProba, sep=",")

}

ML(snakemake@input[[1]], snakemake@input[[2]], snakemake@input[[3]], snakemake@output[[1]], snakemake@output[[2]], snakemake@output[[3]], snakemake@output[[4]])

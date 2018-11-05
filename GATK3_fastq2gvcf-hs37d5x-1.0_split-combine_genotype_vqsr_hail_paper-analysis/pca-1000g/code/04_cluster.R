#!/usr/bin/env Rscript
library(mclust)

scores = read.table("../02_MGRB_1000G_pca.scores.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

pred_vars = c("PC1", "PC2", "PC3", "PC4")

model_superpop = MclustDA(scores[scores$superPopulation != "MGRB",pred_vars], scores$superPopulation[scores$superPopulation != "MGRB"], modelType = "EDDA")
summary(model_superpop)
plot(model_superpop, "scatterplot")
model_eur = MclustDA(scores[scores$superPopulation == "EUR",pred_vars], scores$population[scores$superPopulation == "EUR"], modelType = "EDDA")
summary(model_eur)
plot(model_eur, "scatterplot")

temp = predict(model_superpop, scores[,pred_vars])
scores$pred.superPop = temp$classification
temp = temp$z
colnames(temp) = paste("pred.superPop.", colnames(temp), sep = "")
scores = cbind(scores, temp)

temp = predict(model_eur, scores[,pred_vars])
scores$pred.eurPop = temp$classification
scores$pred.eurPop[scores$pred.superPop != "EUR"] = NA
temp = temp$z
temp[scores$pred.superPop != "EUR",] = NA
colnames(temp) = paste("pred.eurPop.", colnames(temp), sep = "")
scores = cbind(scores, temp)

scores$pred.NFE = FALSE
scores$pred.NFE[scores$pred.superPop == "EUR" & scores$pred.eurPop != "FIN"] = TRUE
scores$prob.NFE = 1 - scores$pred.eurPop.FIN

table(scores$pred.superPop, scores$superPopulation)
table(scores$pred.NFE, scores$superPopulation)

write.table(scores, "../04_MGRB_1000G_pca.scores_clustered.tsv", col.names = TRUE, row.names = FALSE)

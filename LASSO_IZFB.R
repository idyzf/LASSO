library(openxlsx)
rm(list = ls())
## load clinical data
clinical <- read.xlsx("../clinical.xlsx")
rownames(clinical) <- clinical$ID

## load molecular data
data <- read.xlsx("../ruvcorrected.xlsx")
rownames(data) <- data$Name
data <- data[,-1]
data <- as.data.frame(data)


order <- rownames(data)
clinical <- clinical[order,]
banco <- data.frame(clinical,data)

## y sera a variavle a ser predita
y <- clinical$status

## x serao as variaveis testadas para predizer (genes)
x <- as.matrix(data)

library(glmnet)

#perform k-fold cross-validation to find optimal lambda value
cv_model <- cv.glmnet(x, y, alpha = 1)

#find optimal lambda value that minimizes test MSE
best_lambda <- cv_model$lambda.min

#produce plot of test MSE by lambda value
png("cv_model.png")
plot(cv_model)
dev.off()

#find coefficients of best model

fit <- glmnet(x, y, alpha = 1, maxit = 10000000, lambda = best_lambda)
coef(fit)


tmp_coeffs <- coef(fit)
teste <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)


write.xlsx(as.data.frame(teste), "selected_genes.xlsx", rowNames = TRUE)

#use fitted best model to make predictions
y_predicted <- predict(fit, s = best_lambda, newx = x)

#find SST and SSE
sst <- sum((y - mean(y))^2)
sse <- sum((y_predicted - y)^2)

#find R-Squared
rsq <- 1 - sse/sst
rsq

#######################################################################
patients <- rownames(data)

SCORE <- y_predicted
colnames(SCORE) <- "SCORE"

mediana <- median(SCORE)			
groups <- clinical[patients,]			
data_kaplan <- data.frame(groups, SCORE)
data_kaplan$RISK[data_kaplan$SCORE >= mediana] <- "High_score"
data_kaplan$RISK[data_kaplan$SCORE < mediana] <- "Low_score"

write.xlsx(data_kaplan,"riskscore.xlsx", rownames = TRUE)

################# kaplan meier
library(openxlsx)
library(survminer)
library(survival)

time <- clinical$time
status<- clinical$status

png("survival_score.png", width = 110, height = 140, res = 400, units = "mm") ####
fit_k <- survfit(Surv(time, status) ~ RISK, data = data_kaplan) ####
ggsurvplot(fit_k, data = data_kaplan, censor.shape = "|", censor.size = 3,
           xlim = c(0,78),
           pval = TRUE,
           risk.table = TRUE,
           legend.title = "",
           palette = c("red","darkblue"),
           legend.labs = c("High_score", "Low_score"),
           pval.coord = c(40,0.92),
           risk.table.y.text = FALSE,
           title = "Risk Score" #####
)
dev.off()

stats <- summary(fit_k)$table
write.xlsx(as.data.frame(stats), "survival_stats.xlsx", rowNames=TRUE)

#library(pROC)
#pROC_obj <- roc(y, SCORE,
 #               smoothed = TRUE,
   #             ci=TRUE, ci.alpha=0.95, stratified=FALSE,
    #            plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
      #          print.auc=TRUE, show.thres=TRUE)
#sens.ci <- ci.se(pROC_obj)
#plot(sens.ci, type="shape", col="lightblue")

#plot.roc(pROC_obj,  type="shape", col = "red", title = "pROC package", lty=1,lwd=2,xlim = c(1.8, 0))
library(PRROC)
PRROC_obj <- roc.curve(scores.class0 = SCORE, weights.class0=y,
                       curve=TRUE)
plot(PRROC_obj)

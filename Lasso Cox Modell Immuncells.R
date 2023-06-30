# lasso OS

#install.packages("glmnet")
#install.packages("survival")

library(glmnet)
library(tidyverse)
library(survival)
library(openxlsx)
library(maxstat)
library(WriteXLS)

df <- drop_na(df)

rt  <- df
rt <- rt%>%
        tibble::column_to_rownames("Patient_ID")

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$OS_5y_months,rt$OS_5y_event))

fit <- glmnet(x, y, family = "cox", maxit = 100000)
plot(fit, xvar = "lambda", label = FALSE, lw=2)

#cvfit <- cv.glmnet(x, y, family="cox", maxit = 100000)
cvfit <- cv.glmnet(x, y, family="cox")
best_lambda <- cvfit$lambda.min
best_lambda
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
plot(cvfit)
# rt$best.lambda    code zeile wirklich nötig?
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoCelltype=row.names(coef)[index]
geneCoef=cbind(Celltype=lassoCelltype,Coef=actCoef)
write.table(geneCoef,file="CoeF-OS.txt",sep="\t",quote=F,row.names=F)

riskScore=predict(cvfit, newx = x, s = "lambda.min" ,type="link")
outCol=c("OS_5y_months","OS_5y_event",lassoCelltype)
outTab=cbind(rt[,outCol],riskScore=as.vector(riskScore))
write.table(cbind(id=rownames(outTab),outTab),
            file="lasso-Risk-Cox-OS-Immunstatus.txt",
            sep="\t",
            quote=F,
            row.names=F)


# best cutoff

stat <- maxstat.test(Surv(outTab$OS_5y_months,outTab$OS_5y_event)~outTab$riskScore, data = outTab, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk=as.vector(ifelse(riskScore>cutoff,"high","low"))
write.table(cbind(id=rownames(outTab),outTab,risk),
            file="lasso-Risk-Cox-best-cut-OS-Immunstatus.txt",
            sep="\t",
            quote=F,
            row.names=F)

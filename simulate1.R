library(lme4)
library(dplyr)
library(did)
library(ggplot2)

n <- 10
t <- 10
ob <- n*t
id <- rep(1:n, each=t)
t <- rep(1:t, n)


mm1 <- c(0)
mm2 <- c(0)
mm3 <- c(0)
cs1 <- c(0)
cs2 <- c(0)
cs3 <- c(0)

for (i in 1:10000) {

  no <- as.matrix(rnorm(n,0,1))
  no1 <- rep(no, each=10)
  
treat <- ifelse(id >25,1,0)
eta <- id + no1
d <- ifelse(treat==1 & id < 40, 1, ifelse(treat==1 & id >= 40,2,0))

y <- ifelse(d==1 & t >4, 5 +t+eta+ rnorm(ob,0,1),
            ifelse(d==2 & t > 4, 10+t+eta+rnorm(ob,0,1), t+eta+rnorm(ob,0,1)))
ttt <- ifelse(treat==1,5,0)
ei <- ifelse(treat==1,t-ttt,-6)

data <- data.frame(id, t, treat, d, ttt, ei,y)
event <- lmer(y~treat + factor(t)+factor(id) +(1|ei),data=data, REML =FALSE)

# Getting attribute-cluster random-intercepts

u0 <- ranef(event,condVar=TRUE)
u <- u0$ei
ue <-sqrt(attr(u0[[1]], "postVar")[1,,])

# Calculating shrinkage values
ran_var <- VarCorr(event)$ei[,1]
ran_re <- attr(VarCorr(event),"sc")^2
n <- table(ei)
nj <- as.numeric(n)
shrinkage <- ran_var/(ran_var + ran_re/nj)

# Stack data and calculate ATT

result <- data.frame(u,ue, shrinkage, ran_var, ran_re)
colnames(result)[colnames(result)=="X.Intercept."] <- "u"
result$ei <- as.numeric(rownames(result))
result$u1 <- result$u/result$shrinkage
result$pp <- mean(result$u1[result$ei<0 & result$ei >-6], na.rm=TRUE)
result$ATT <- result$u - result$pp

mm1[i] = result$ATT[6]
mm2[i] = result$ATT[7]
mm3[i] = result$ATT[8]

out <- att_gt(yname="y", tname="t", idname="id", gname="ttt", xformla=NULL,data=data)
cs1[i] = out$att[4]
cs2[i] = out$att[5]
cs3[i] = out$att[6]
}

par(mfrow=c(1,2))
hist(mm1, breaks =100)
hist(cs1, breaks=100)




# Adding attribute clusters in the model to optimise the results
for (i in 1:1000) {
  
treat <- ifelse(id >4,1,0)
d1 <- ifelse(treat==1, runif(ob,0,1),0)
d1 <- ave(d1, id, FUN = mean)
d <- ifelse(treat==1 & d1 >.5, 1, ifelse(treat==1 & d1<=.5 ,2,0))
x <- rnorm(ob,0,1)
y <- ifelse(d==1 & t >4, 5+x +t+ rnorm(ob,0,1),
            ifelse(d==2 & t > 4, 10+x+t+rnorm(ob,0,1), x+t+rnorm(ob,0,1)))
ttt <- ifelse(treat==1,5,0)
ei <- ifelse(treat==1,t-ttt,-6)

data <- data.frame(id, t, treat, d, ttt, ei,y)
event <- lmer(y~treat +x+ factor(t) +(1|ei) + (1|ei:d),data=data, REML =FALSE)

# Getting attribute-cluster random-intercepts

u0 <- ranef(event,condVar=TRUE)
u <- u0$ei
colnames(u)[colnames(u)=="(Intercept)"] <- "u"
u$ei <- as.numeric(rownames(u))
ue <-sqrt(attr(u0[[2]], "postVar")[1,,])
u <- data.frame(u,ue)

u_ac <- u0$`ei:d`
u_ac$id <- rownames(u_ac)
u_ac$id <- strsplit(u_ac$id, ":")
id1 <- do.call(rbind, u_ac$id)
u_ac$ei <- as.numeric(id1[,1])
u_ac$d <- as.numeric(id1[,2])
ue_ac <- sqrt(attr(u0[[1]], "postVar")[1,,])
colnames(u_ac)[colnames(u_ac)=="(Intercept)"] <- "u_ac"
u_ac <- data.frame(u_ac, ue_ac)

u_comb <- merge(u, u_ac,by="ei", all = FALSE)

# Calculating shrinkage values
ran_var <- VarCorr(event)$ei[,1]
ran_var_ac <- VarCorr(event)$`ei:d`[,1]
ran_re <- attr(VarCorr(event),"sc")^2

data$n <- 1
data$nj <- ave(data$n, data$ei, FUN = sum)
data$njj <- ave(data$n, data$ei, data$d, FUN = sum)
data$tag <- !duplicated(data[,c("ei","d")])
data_sub <- subset(data, tag=="TRUE")
data_sub <- data_sub[,c("d", "ei", "nj", "njj")]
data_sub <- data.frame(data_sub, ran_var, ran_var_ac, ran_re)

# Stack data and calculate ATT

result <- merge(u_comb, data_sub, by=c("ei", "d"), all = FALSE)
result$shrinkage_1 <- result$ran_var / (result$ran_var + result$ran_re/result$nj)
result$shrinkage_2 <- result$ran_var_ac / (result$ran_var_ac + result$ran_re / result$njj)

result$u_1 <- ifelse(result$shrinkage_1==0,0,result$u / result$shrinkage_1)
result$u_2 <- result$u_ac / result$shrinkage_2

result$u_f <- result$u_1 + result$u_2

result$pp1 <- mean(result$u_f[result$ei<0 & result$d==1], na.rm=TRUE)
result$pp2 <- mean(result$u_f[result$ei<0 & result$d==2], na.rm=TRUE)
result$pp <- ifelse(result$d==1, result$pp1,result$pp2)

result$ATT <- result$u_f - result$pp
result$ATE <- ave(result$ATT, result$ei, FUN = mean)

mm1[i] = result$ATE[result$ei==0]
mm2[i] = result$ATE[result$ei==1]
mm3[i] = result$ATE[result$ei==2]

out <- att_gt(yname="y", tname="t", idname="id", gname="ttt", xformla=~x,data=data, control_group="notyettreated")
cs1[i] = out$att[4]
cs2[i] = out$att[5]
cs3[i] = out$att[6]
}


mm_b1 <- seq(min(mm1), max(mm1), length=100)
d_mm_b1 <- dnorm(mm_b1, mean = mean(mm1), sd=sd(mm1))
cs_b1 <- seq(min(cs1), max(cs1), length=100)
d_cs_b1 <- dnorm(cs_b1, mean=mean(cs1), sd=sd(cs1))

plot(mm_b1,d_mm_b1)
lines(mm_b1, d_mm_b1, col="red") +
lines(cs_b1, d_cs_b1, col="blue")

mm_b1 <- seq(min(data_simu$mm1, na.rm = TRUE), max(data_simu$mm1, na.rm=TRUE), length=100)
d_mm_b1 <- dnorm(mm_b1, mean = mean(data_simu$mm1, na.rm=TRUE), sd=sd(data_simu$mm1, na.rm=TRUE))
cs_b1 <- seq(min(data_simu$cs1, na.rm = TRUE), max(data_simu$cs1, na.rm = TRUE), length=100)
d_cs_b1 <- dnorm(cs_b1, mean=mean(data_simu$cs1, na.rm=TRUE), sd=sd(data_simu$cs1, na.rm = TRUE))

ggplot() +
  geom_line(aes(x=mm_b1, y=d_mm_b1, col="mmdid")) + 
  geom_line(aes(x=cs_b1, y=d_cs_b1, col="csdid"))




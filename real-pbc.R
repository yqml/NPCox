rm(list = ls())

library(NPCox)
setwd('C:/Users/qysdu/Desktop/npcox')
source('npcox funcs.R')
source('yqfuc.R')
print(Sys.time())

data(pbc)

dta  = na.omit(pbc[,c('time', 'status', 'age', "albumin", "bili", "protime", "edema")])
dta[,'status'] = sign(dta[,'status'])
dta[,c("bili", "protime", "albumin")] = log(dta[,c("bili", "protime", "albumin")])
colnames(dta) = c('time', 'status', 'age', "log(albumin)", "log(bili)", "log(protime)", "edema")


## Semi PH model:
print(Sys.time())
# res_semi  = spcox(cva_cons = dta[,3:4], cva_time = dta[,5:7],
#               delta = dta$status, obstime = dta$time, SE = F, bandwidth = 783)	
res_semi  = spcox(cva_cons = as.matrix(dta[,3:7]), cva_time = NULL,
              delta = dta$status, obstime = dta$time, SE = F, bandwidth = 783)
print(Sys.time())
semi = cbind(floor(res_semi$constant_coef*1000)/1000, 
             floor(res_semi$constant_coef_SEE*1000)/1000)
colnames(semi) = c('esti','see')
rownames(semi) = colnames(dta)[c(3:7)]
semi

## Comparative analysis 
library(survival)
coxres = coxph(Surv(time, status) ~ age + `log(albumin)` + `log(bili)`+ `log(protime)` + edema, dta )
coxres

if(0){
  cva_cons = as.matrix(dta[,3:7]); cva_time = NULL;
  delta = dta$status; obstime = dta$time; SE = F; bandwidth = 783; resamp = 100
}

## Setting 1: bandwidth = 500, resamp = 100
res1  = npcox(cva = dta[,3:7], delta = dta$status, obstime = dta$time, SE = F, bandwidth = 500)	
print(Sys.time())
## Setting 2: bandwidth = FALSE (Auto selected: 783), resamp = 100
res2  = npcox(cva = dta[,3:7], delta = dta$status, obstime = dta$time, SE = T, bandwidth = FALSE)
print(Sys.time())
## Setting 3: bandwidth = 500, resamp = 200
res3  = npcox(cva = dta[,3:7], delta = dta$status, obstime = dta$time, SE = T, bandwidth = 500, 
              resamp = 200)
print(Sys.time())
# write.table(data.frame(res1$temporal_coef,res1$temporal_coef_SEE,
#                        res1$data),  file = paste("C:/Users/qysdu/Desktop/npcox/nppbc1.txt", sep = "" ))
# write.table(data.frame(res2$temporal_coef,res2$temporal_coef_SEE,
#                        res2$data),  file = paste("C:/Users/qysdu/Desktop/npcox/nppbc2.txt", sep = "" ))
# write.table(data.frame(res3$temporal_coef,res3$temporal_coef_SEE,
#                        res3$data),  file = paste("C:/Users/qysdu/Desktop/npcox/nppbc3.txt", sep = "" ))


# print(Sys.time())

## 画图
real1 = read.table(file = paste("C:/Users/qysdu/Desktop/npcox/nppbc1.txt", sep = "" ))
real2 = read.table(file = paste("C:/Users/qysdu/Desktop/npcox/nppbc2.txt", sep = "" ))
real3 = read.table(file = paste("C:/Users/qysdu/Desktop/npcox/nppbc3.txt", sep = "" ))

namm  = c('age',"log(albumin)", "log(bili)", "log(protime)", "edema") 
colnames(real1) = c( namm, paste(namm, '_SEE', sep = ''), 
                     paste('dat_',namm, sep = ""), 'obstime', 'delta')
colnames(real2) = colnames(real1)
colnames(real3) = colnames(real1)

trf = function(real1, band){
  r1  = real1
  lis = list()# "temporal_coef", "temporal_coef_SEE", "bandwidth", "inconverged", "Time points (not converged)"
  lis$temporal_coef = as.matrix(r1[,c(1:5)])
  lis$temporal_coef_SEE = as.matrix(r1[,c(6:10)])
  lis$data = as.matrix(r1[,c(11:17)])
  lis$bandwidth = band
  
  return(lis)
}

r1 = trf(real1, 500)
r2 = trf(real2, 783)
r3 = trf(real3, 500)

par(mfcol = c(5,3), mar = c(2.5,4.5,2,2))
npplot(r1, xrange = c(500,3000))
npplot(r2, xrange = c(783,3000))
npplot(r3, xrange = c(500,3000))


# write.table(res_rec, file = paste("C:/Users/qysdu/Desktop/npcox/real1.txt", sep = "" ))


# cva = as.matrix(dta[,3:7]); delta = dta$status; 
# obstime = dta$time; SE = F; resamp = 100; bandwidth = 500; 
# dta1 = dta[which(pbc$time < 3800),]
# par(mfrow = c(1,2))
# hist(dta$time); hist(dta1$time)

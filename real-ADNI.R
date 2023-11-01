# rm(list = ls())
library(NPCox)
setwd('C:/Users/qysdu/Desktop/npcox')
source('npcox funcs.R')
source('yqfuc.R')

# data extraction
f  = function(){
  library("mvtnorm")
  library(Hmisc)
  library(ADNIMERGE)
  
  df  = adnimerge
  # write.csv(adnimerge, file = 'C:/Users/qysdu/Desktop/adnimerge.csv')
  df1 = df[,c("RID","PTID","VISCODE","DX.bl","AGE","PTGENDER","PTEDUCAT",
              "PTMARRY","FDG","APOE4",'Hippocampus', "DX","Years.bl","ADAS11.bl")]
  #df1 = df1[which(df1$VISCODE == "bl"),]
  
  # indim = array(0,nrow(df1))
  # for(i in 1:nrow(df1)){  if(df1$RID[i] %in% idima){indim[i] = 1}  }
  
  # indim = 
  df2 = df1# [which(indim>0),]
  
  if(1){
    DX.bl = as.character(df2[,"DX.bl"])
    bl.13 = sort(unique(c(which(DX.bl == "AD"))))#which(DX.bl == "CN"),
    df2   = df2[-bl.13,]
  }
  
  DX    = as.character(df2[,"DX"])
  DX[which(DX == "CN")]       = 0
  DX[which(DX == "MCI")]      = 0
  DX[which(DX == "Dementia")] = 1
  DX    = as.numeric(DX)
  ind.DX= which(!is.na(DX))
  DX    = DX[ind.DX]
  df3   = df2[ind.DX,]
  RID   = df3[,"RID"]
  Years.bl = df3[,"Years.bl"]
  
  ID    = unique(RID)
  ID    = cbind(ID,rep(0,length(ID)),rep(0,length(ID)))
  colnames(ID) = c("ID","CN/MCI","AD")
  temp  = cbind(RID,DX)
  
  for(i in 1:nrow(ID)){
    sele = which(temp[,1] == ID[i,1])
    tem1 = unique(temp[sele,2])
    if(sum(tem1)==0){ID[i,2] = 1}
    if(sum(tem1)==1){ID[i,3] = 1}
  }
  ID = cbind(ID,1+ID[,3],rep(0,nrow(ID)),rep(0,nrow(ID)),rep(0,nrow(ID)))
  colnames(ID) = c("ID","CN/MCI","AD","event","delta","X.i","Final")
  
  for(i in 1:nrow(ID)){
    sele = which(temp[,1] == ID[i,1])
    if(ID[i,4]==1){
      ID[i,5] = 0
      ID[i,6] = max(Years.bl[sele])
      ID[i,7] = min(sele)
    }
    if(ID[i,4]==2){
      ID[i,5] = 1
      sele1   = min(which(DX[sele]==1))
      ID[i,6] = c(Years.bl[sele])[sele1]
      ID[i,7] = min(sele)
    }
  }
  final.ind = ID[which(ID[,"X.i"]>0),"Final"]
  
  selec = c("RID","PTID","AGE","PTGENDER","PTEDUCAT","PTMARRY","FDG","APOE4","ADAS11.bl",'Hippocampus')
  Z  = df3[final.ind, selec] 
  Z1 = Z[,-2]
  #Z  = as.matrix(Z)
  X  = ID[which(ID[,"X.i"]>0),"X.i"]
  delta = ID[which(ID[,"X.i"]>0),"delta"]
  
  data  = cbind(Z1,X,delta)#ID[which(ID[,"X.i"]>0),1],
  #colnames(data) = c("ID",colnames(data)[-1])
  ind1  = which(data[,"PTMARRY"] == "Unknown")
  data  = data[-ind1,]
  Z     = Z[-ind1,]
  data[which(data[,"PTGENDER"] == "Male"  ), "PTGENDER"] = 1
  data[which(data[,"PTGENDER"] == "Female"), "PTGENDER"] = 0
  marind = which(data[,"PTMARRY"]  == "Married")
  data[marind,"PTMARRY"]  = 0
  data[c(1:nrow(data))[-marind],"PTMARRY"] = 1
  #data  = as.matrix(data)
  #class(data) = "numeric"
  rownames(data) = 1:nrow(data)
  #indi.NA = sort(unique(which(is.na(data[,"V"]))))
  #data    = as.array(data[-indi.NA,])
  #temp = data.frame(data[,1]])
  if(0){ write.csv(Z[,c(1,2)], file = "C:/Users/YANG QI/Desktop/RID.csv") }
  
  for(j in 1:ncol(data)){ data[,j] = as.numeric(data[,j]) }
  data = na.omit(data)
  data = data[,-which(names(data) == "PTMARRY")]
  return(data)
}
data  = f()

fstan = function(data){
  if(is.vector(data)){
    data = (data - mean(data))/sd(data)
  }else if(is.matrix(data) |is.data.frame(data)){
    for(j in 1:ncol(data)){
      data[,j] = (data[,j] - mean(data[,j]))/sd(data[,j])
    }
  }
  return(data)
}
data[,c('AGE','FDG','Hippocampus')] = fstan(data[,c('AGE','FDG','Hippocampus')])
colnames(data) = c('RID','Age','Gender','Education','FDG','APOE4','ADAS11','Hippocampus','obstime','delta')

dta = as.matrix(data)
dta = dta[which(dta[,'obstime']<11),]

## npcox result
res = npcox(cva = dta[,2:8], delta = dta[,10], obstime = dta[,9], SE = T)

par(mfrow = c(2,4))
npplot(res, xrange = c(1.8,7))


## spcox result 

## Semi PH model:
if(0){
  cva_cons = dta[,c(2:8)]; cva_time = NULL; resamp = 100
  delta = dta[,10]; obstime = dta[,9]; SE = F; bandwidth = 1.46
  cva_cons = dta[,c(2,4:8)]; cva_time = dta[,3]; resamp = 100
  delta = dta[,10]; obstime = dta[,9]; SE = F; bandwidth = 1.46
}

print(Sys.time())
res_semi  = spcox(cva_cons = dta[,c(2:8)], cva_time = NULL,
                  delta = dta[,10], obstime = dta[,9], SE = F, bandwidth = 2)	
print(Sys.time())
semi = cbind(floor(res_semi$constant_coef*1000)/1000, 
             floor(res_semi$constant_coef_SEE*1000)/1000)
colnames(semi) = c('esti','see')
rownames(semi) = colnames(dta)[c(2:8)]
semi



se<-c(2:8)
nl<-length(se)
res<-lapply(1:nl,function(i) combn(se,i))

# s = 7
# options(scipen = 100)
# data1 = as.matrix(data)[,2:ncol(data)]
# phres = fcox(data1[,1:s],data1[,s+2],data1[,s+1]); phres
# write.table(phres, file = 'C:/Users/qysdu/Desktop/npcox/phres.txt', sep = ' & ')

library(survival)
coxres = coxph(Surv(obstime, delta) ~ Age + Gender + Education + FDG + APOE4 + ADAS11 + Hippocampus, data.frame(dta))
coxres



## Setting 1: bandwidth = 500, resamp = 100
# print(Sys.time())
# res15 = npcox(cva = dta[,2:8], delta = dta[,10], obstime = dta[,9], SE = F, bandwidth = 1.5)
# res2  = npcox(cva = dta[,2:8], delta = dta[,10], obstime = dta[,9], SE = F, bandwidth = 2)
# res3  = npcox(cva = dta[,2:8], delta = dta[,10], obstime = dta[,9], SE = F, bandwidth = 3)
# res   = npcox(cva = dta[,2:8], delta = dta[,10], obstime = dta[,9], SE = F)
# par(mfrow = c(3,3))
# npplot(res15, xrange = c(1.5,9.5))
# par(mfrow = c(3,3))
# npplot(res2, xrange = c(2,9))
# par(mfrow = c(3,3))
# npplot(res3, xrange = c(2,9))
# print(Sys.time())
# if(0){
#   cva = dta[,2:8]; delta = dta[,10];
#   obstime = dta[,9]; SE = F; bandwidth = 30
# }
# if(0){
#   par(mfrow = c(3,3))
#   fsummary(dta[,2:10])
#   par(mfrow = c(3,3))
#   fhist(dta[,2:10])
# }
# , bandwidth = 1.46
# if(0){
#   aa = readClipboard()
#   length(unique(aa))
#   sum(data[,'delta'])
# }

res   = npcox(cva = dta[,2:8], delta = dta[,10], obstime = dta[,9], SE = T, bandwidth = 1.46)
colnames(res$temporal_coef) = c('Age','Gender','Education','FDG','APOE4','ADAS11','Hippocampus')

par(mfrow = c(2,4))
npplot(res, xrange = c(1.8,7))


print(Sys.time())
write.table(data.frame(res$temporal_coef,res$temporal_coef_SEE, res$data),
            file = paste("C:/Users/qysdu/Desktop/npcox/real-adni.txt", sep = "" ))

if(0){
  cva_cons = dta[,3:4]; cva_time = dta[,5:7];
  delta = dta$status; obstime = dta$time; SE = F; bandwidth = 783; resamp = 100
}




# par(mfrow = c(3,3))
# datsum = function(data){
#   name = colnames(data)
#   for(j in 1:ncol(data)){
#     plot(sort(data[,j]),main = name[j])
#   }
# }
# datsum(data)
# summary(data)





# write.table(data, file = "C:/Users/qysdu/Desktop/npcox/adni.txt")





















# RID  = as.numeric(data$RID)
# ind  = array()
# for(i in 1:length(idima)){
#   if(idima[i] %in% RID){
#     ind[i] = which(RID == idima[i])
#   }else{
#     ind[i] = NA
#   }
# }

# write.csv(RID,file = "C:/Users/YANG QI/Desktop/RID.csv")
# setwd("C:/Users/Administrator/Desktop")
# install.packages("Hmisc_4.6-0.zip", repos = NULL)
# install.packages("C:/Users/YANG QI/Desktop/ADNIMERGE_0.0.1.tar.gz", repos=NULL, type="source")

# take the name from masked file
# f1 = function(){
#   setwd("D:/masked")
#   name  = list.files()
#   name1 = array()
#   ndir  = length(name)
#   for(i in 1:ndir){
#     if(grepl("m00", name[i])){
#       name1 = c(name1,name[i])
#     }
#   }
#   name1 = name1[-1] # unique(name) # 
#   IDima = array(length(name1))
#   temp  = 0
#   for(i in 1:length(name1)){
#     if(grepl("ADsim", name1[i])){
#       IDima[i] = c(substring(name1[i],17-temp,20))
#     }else if(grepl("MCIsim", name1[i])){
#       IDima[i] = c(substring(name1[i],18-temp,21))
#     }else if(grepl("NORMALsim", name1[i])){
#       IDima[i] = c(substring(name1[i],21-temp,24))
#     }
#   }#IDima1 = unique(IDima)
#   
#   return(IDima)
# }
# idima = f1()
# class(idima) = "numeric"


# if(0){
  # ID1 = data$ID
  # 
  # ind = array()
  # for(i in 1:length(idima)){
  #   if(idima[i] %in% ID1){
  #     ind[i] = which(ID1 == idima[i])
  #   }else{
  #     ind[i] = NA
  #   }
  # }
# }

setwd('C:/Users/qysdu/Desktop/npcox')
source('npcox funcs.R')
source('yqfuc.R')
# library(NPCox)

# re = 1
repe = 500
name = "n1c1"
for(re in 1:repe){
  
  print(paste("re = ",re, ", Time = ", Sys.time(), sep = ""))
  # data generation
  data = npsimu(500,cenpara = 3) # ;sum(data$delta)/300
  
  # Para estimation through the function pht
  for(h in c(0.1,0.2)){
    res  = npcox(data$cva, delta = data$delta, obstime = data$obstime, SE = T, bandwidth = h)
    res_rec = cbind(res$temporal_coef,res$temporal_coef_SEE,res$data[,3])
    res_rec[res$inconverged, ] = NA
    write.table(res_rec, file = paste("C:/Users/qysdu/Desktop/npcox/simurec/",name,"re",re,"h",h,".txt", sep = "" ))
  }
}



library(NPCox)
data = npsimu(n = 300,cenpara = 6)
# Para estimation through the function npcox
res  = npcox(data$cva, delta = data$delta, obstime = data$obstime, SE = T, bandwidth = 0.1)
res_rec[res$inconverged, ] = NA

rm(list = ls())

repe = 190
name = "n1c1"
nam1 = c("n1c1","n1c2")#,"n2c1","n2c2"
re   = 1
s    = 2
h    = 0.1
hr   = c(0.1,0.2)
v    = function(u){View(cbind(u))}
na   = function(u){which(is.na(u))}
su   = function(u){summary(u)}
de   = function(u){dev.off()}
num  = 200
rec  = array(NA,dim = c(2,repe,num,2*s+1))

fsmooth = function(ori){
  ori = ori[which(!is.na(ori[,1])),]
  n   = nrow(ori)
  Xs  = ori[,ncol(ori)]
  est = as.matrix(ori[,1:(ncol(ori)-1)])
  Xn  = seq(h,1, by = (1-h)/(num-1))
  new = array(0, dim = c(length(Xn),ncol(ori) -1))
  for(i in 1:length(Xn)){
    ind = min(which(sort(c(Xn[i],Xs)) == Xn[i]))
    tem = function(lam){ lam*Xs[ind-1] + (1-lam)*Xs[ind] - Xn[i]}
    lam = uniroot(tem,c(0,1))$root
    new[i,] = c(lam*est[ind-1,] + (1-lam)*est[ind,])
  }
  write.table(cbind(new,Xn),paste("C:/Users/qysdu/Desktop/test ubeta see/",name,"re",re,"h",h,".txt", sep = "" ))
  return(cbind(new,Xn))
}

smtf = function(name,h){
  print(paste('name = ', name, ', h = ', h, ', ', Sys.time(), sep = '' ))
  for(re in 1:repe){
    ori = read.table(paste("C:/Users/qysdu/Desktop/test ubeta see/",name,"re",re,"h",h,".txt", sep = "" ))
    if(ori[1,5] >= 0.1){ori[1:2,5] = c(0.099,0.0998)}
    if(ori[nrow(ori),5] <= 1){ori[nrow(ori):(nrow(ori)-1),5] = c(1.01,1)}
    m1 = min(which(ori[,5] >= 1))
    if(sum(is.na(ori[m1,])) > 0){ori[m1,1:4] = ori[m1-1,1:4]}
    
    rec[1,re,,] = fsmooth(ori)
    if(re%%100 == 0) print(paste('re = ', re, ', ', Sys.time()))
  }
  h1mean = array(0,dim = c(num,3*s))
  Xn     = seq(h,1, by = (1-h)/(num-1))
  for(j in 1:2){
    h1mean[,j]   = colMeans(rec[1,,,j])
    h1mean[,j+2] = colMeans(rec[1,,,j+2])
    h1mean[,j+4] = apply(rec[1,,,j], 2, sd)
  }
  b1 = function(t){1/(1+t)}
  b2 = function(t){t}
  true = cbind(b1(Xn),b2(Xn))
  res = data.frame(cbind(h1mean,Xn,true))
  colnames(res) = c('b1','b2','b1see','b2see','b1sd','b2sd','Xn','true1','true2')
  return(res)
}


## result summary and storage in txt files
name = NULL
for(h in hr){
  res = smtf(name,h)
  write.table(res, file = paste("C:/Users/qysdu/Desktop/test ubeta see/",name,"h",h,".txt", sep = "" ))
}
# for(name in nam1){
#   
# }


## For plot comparison
npcoxplot = function(res){
  Xn   = res[,'Xn']
  true = res[,c('true1','true2')]
  cc  = c(1:max(which(Xn <= 0.8)))
  for(j in 1:2){# "beta",j,"(t)",sep = ""
    ynam = c(expression(paste(beta[1](t), "= 1/(1+t)")),expression(paste(beta[2](t), "= t")))
    plot(Xn[cc],  res[cc,j], ylab = ynam[j], main = paste('beta',j,'(t)',sep = ''), cex.lab = 1.2, cex.main = 1,
         xlab = "Observed Time (t)", ylim = c(-3,4)*0.7, xlim = c(h,0.8), type = "l")
    lines(Xn[cc], true[cc,j],   col  = "red")
    lines(Xn[cc], res[cc,j] + 1.96*res[cc,j+2], col = "green")
    lines(Xn[cc], res[cc,j] - 1.96*res[cc,j+2], col = "green")
    lines(Xn[cc], res[cc,j] + 1.96*res[cc,j+4], col = "blue")
    lines(Xn[cc], res[cc,j] - 1.96*res[cc,j+4], col = "blue")
  }
}

par(mfrow = c(2, 2), mar = c(4,5,2,2))
res1 = read.table(paste("C:/Users/qysdu/Desktop/test ubeta see/",name,"h",0.1,".txt", sep = "" ))
res2 = read.table(paste("C:/Users/qysdu/Desktop/test ubeta see/",name,"h",0.2,".txt", sep = "" ))
npcoxplot(res1)
npcoxplot(res2)

madrec = c()
stdrec = c()
for(name in nam1){
  for(h in hr){
    temp = read.table(paste("C:/Users/qysdu/Desktop/npcox/simurec/",name,"h",h,".txt", sep = "" ))
    h1mean = temp[,1:6]
    Xn     = temp[,7]
    true   = temp[,8:9]
    
    n  = ifelse(grepl('n1', name), 300, 500)
    cr = ifelse(grepl('c1', name), 15, 30)
    namt   = paste('n=', n, ', CR=', cr,'%', ', h=', h, sep = "")
    
    MAD = array(0, s)
    std = array(0, s)
    sel = c() # seq(1,length(Xn), by = 10)
    for(i in 1:20){ sel[i] = max( which(Xn <= 0.21 + (i-1)*0.03 )) }
    for(j in 1:2){ MAD[j] = mean(abs(h1mean[sel,j] - true[sel,j])) }
    for(j in 1:2){ std[j] = mean(h1mean[sel,j+2]) }
    
    npcoxplot(h1mean)
    madrec = c(madrec, MAD)
    stdrec = c(stdrec, std)
  }
}
madrec = matrix(madrec, byrow = T, nrow = 2)
stdrec = matrix(stdrec, byrow = T, nrow = 2)






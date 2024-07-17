## Comparison with Cox(log-transform)
## PAM with log-transform, PAM with penalized spline.
## Package: pammtools;

## The nppbc2.txt required in plot is available in https://github.com/yqml/NPCox

# install.packages('pammtools')
rm(list = ls())
library(tidyr)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())
library(survival)
library(mgcv)
library(pammtools)
Set1    = RColorBrewer::brewer.pal(9, "Set1")
Greens  = RColorBrewer::brewer.pal(9, "Greens")
Purples = RColorBrewer::brewer.pal(9, "Purples")

find = function(ori, sub){
  inter = intersect(ori,sub)
  n     = length(inter)
  ind   = array(n)
  for(i in 1:n){
    ind[i] = min(which(ori == inter[i])) 
  }
  return(ind)
}


### data reading from our npcox result
data(pbc)
dta  = na.omit(pbc[,c('time', 'status', 'age', "albumin", "bili", "protime", "edema")])
dta[,'status'] = sign(dta[,'status'])
dta[,c("bili", "protime", "albumin")] = log(dta[,c("bili", "protime", "albumin")])
colnames(dta) = c('time', 'status', 'age', "logal", "logbi", "logpro", "edema")
### ped transformation
dta1 = dta %>% as_ped(Surv(time, status)~., id = "id") %>%
  mutate(logt25 = log(tstart + (tstart - tend) / 2 + 25))

real2 = read.table(file = paste("C:/Users/qysdu/Desktop/R journal submission/npcox/nppbc2.txt", sep = "" ))
namm  = c('age',"log(albumin)", "log(bili)", "log(protime)", "edema")
colnames(real2) = c( namm, paste(namm, '_SEE', sep = ''), 
                     paste('dat_',namm, sep = ""), 'obstime', 'delta')

##################################### age #####################################
## cox with log-transform
vfit1 = coxph(
  formula = Surv(time, status) ~ age+logal+logbi+logpro+edema + tt(age),
  data    = dta,
  tt      = function(x, t, ...) x * log(t + 25))
# coef(vfit1)

## PAM with log-transform
pamlog1 = gam(ped_status ~ s(tend) +age+logal+logbi+logpro+edema+age:logt25,
          data = dta1, offset = offset, family = poisson())


pena1 = gam(ped_status ~ s(tend)+age+logal+logbi+logpro+edema
         + s(tend, by = age), data = dta1, offset = offset, family = poisson())

term.df = dta1 %>% ped_info() %>% add_term(pena1, term = "age") %>%
  mutate_at(c("fit", "ci_lower", "ci_upper"), funs(. / .data$age)) %>%
  mutate(
    cox1 = coef(vfit1)["age"] + coef(vfit1)["tt(age)"] * log(tend + 25),
    pamlog1 = coef(pamlog1)["age"] + coef(pamlog1)["age:logt25"] * log(tend + 25))

ind = find(real2$obstime, term.df$tend)
subreal2 = real2[ind,c(1:10)]
term.df$npage    = subreal2$age
term.df$npagese  = subreal2$age_SEE

gg_tv_age = ggplot(term.df, aes(x = tend, y = fit)) +
  geom_step(aes(y = npage, col = "NPCox estimation")) +
  geom_stepribbon(aes(ymin = npage - 1.96*npagese, 
                      ymax = npage + 1.96*npagese), fill = "green", alpha = 0.2) +
  geom_step(aes(col = "PAM with penalized spline")) +
  geom_stepribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2) +
  geom_line(aes(y = cox1, col = "Cox with log-transform")) +
  geom_step(aes(y = pamlog1, col = "PAM with log-transform")) +
  scale_color_manual(name = "Method", values = c(Set1[1:3], "black")) +
  xlab("time") + ylab('beta(t) for age')
plot(gg_tv_age)


##################################### log(albumin) #####################################
## cox with log-transform
vfit2 = coxph(
  formula = Surv(time, status) ~ age+logal+logbi+logpro+edema + tt(logal),
  data    = dta,
  tt      = function(x, t, ...) x * log(t + 25))
# coef(vfit2)

## PAM with log-transform
pamlog2 = gam(ped_status ~ s(tend) +age+logal+logbi+logpro+edema+logal:logt25,
              data = dta1, offset = offset, family = poisson())


pena2 = gam(ped_status ~ s(tend)+age+logal+logbi+logpro+edema
            + s(tend, by = logal), data = dta1, offset = offset, family = poisson())

term.df = dta1 %>% ped_info() %>% add_term(pena2, term = "logal") %>%
  mutate_at(c("fit", "ci_lower", "ci_upper"), funs(. / .data$logal)) %>%
  mutate(
    cox2 = coef(vfit2)["logal"] + coef(vfit2)["tt(logal)"] * log(tend + 25),
    pamlog2 = coef(pamlog2)["logal"] + coef(pamlog2)["logal:logt25"] * log(tend + 25))

ind = find(real2$obstime, term.df$tend)
subreal2 = real2[ind,c(1:10)]
term.df$nplogal  = subreal2$`log(albumin)`
term.df$nplogalse  = subreal2$`log(albumin)_SEE`

gg_tv_logal = ggplot(term.df, aes(x = tend, y = fit)) +
  geom_step(aes(y = nplogal, col = "NPCox estimation")) +
  geom_stepribbon(aes(ymin = nplogal - 1.96*nplogalse, 
                      ymax = nplogal + 1.96*nplogalse), fill = "green", alpha = 0.2) +
  geom_step(aes(col = "PAM with penalized spline")) +
  geom_stepribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2) +
  geom_line(aes(y = cox2, col = "Cox with log-transform")) +
  geom_step(aes(y = pamlog2, col = "PAM with log-transform")) +
  scale_color_manual(name = "Method", values = c(Set1[1:3], "black")) +
  xlab("time") + ylab('beta(t) for log(albumin)')
plot(gg_tv_logal)


##################################### log(bili) #####################################
## cox with log-transform
vfit3 = coxph(
  formula = Surv(time, status) ~ age+logal+logbi+logpro+edema + tt(logbi),
  data    = dta,
  tt      = function(x, t, ...) x * log(t + 25))
# coef(vfit2)

## PAM with log-transform
pamlog3 = gam(ped_status ~ s(tend) +age+logal+logbi+logpro+edema+logbi:logt25,
              data = dta1, offset = offset, family = poisson())


pena3 = gam(ped_status ~ s(tend)+age+logal+logbi+logpro+edema
            + s(tend, by = logbi), data = dta1, offset = offset, family = poisson())

term.df = dta1 %>% ped_info() %>% add_term(pena3, term = "logbi") %>%
  mutate_at(c("fit", "ci_lower", "ci_upper"), funs(. / .data$logbi)) %>%
  mutate(
    cox3 = coef(vfit3)["logbi"] + coef(vfit3)["tt(logbi)"] * log(tend + 25),
    pamlog3 = coef(pamlog3)["logbi"] + coef(pamlog3)["logbi:logt25"] * log(tend + 25))

ind = find(real2$obstime, term.df$tend)
subreal2 = real2[ind,c(1:10)]
term.df$nplogbi  = subreal2$`log(bili)`
term.df$nplogbise  = subreal2$`log(bili)_SEE`

gg_tv_logbi = ggplot(term.df, aes(x = tend, y = fit)) +
  geom_step(aes(y = nplogbi, col = "NPCox estimation")) +
  geom_stepribbon(aes(ymin = nplogbi - 1.96*nplogbise, 
                      ymax = nplogbi + 1.96*nplogbise), fill = "green", alpha = 0.2) +
  geom_step(aes(col = "PAM with penalized spline")) +
  geom_stepribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2) +
  geom_line(aes(y = cox3, col = "Cox with log-transform")) +
  geom_step(aes(y = pamlog3, col = "PAM with log-transform")) +
  scale_color_manual(name = "Method", values = c(Set1[1:3], "black")) +
  xlab("time") + ylab('beta(t) for log(bili)')
plot(gg_tv_logbi)



##################################### log(protime) #####################################
## cox with log-transform
vfit4 = coxph(
  formula = Surv(time, status) ~ age+logal+logbi+logpro+edema + tt(logpro),
  data    = dta,
  tt      = function(x, t, ...) x * log(t + 25))
# coef(vfit2)

## PAM with log-transform
pamlog4 = gam(ped_status ~ s(tend) +age+logal+logbi+logpro+edema+logpro:logt25,
              data = dta1, offset = offset, family = poisson())


pena4 = gam(ped_status ~ s(tend)+age+logal+logbi+logpro+edema
            + s(tend, by = logpro), data = dta1, offset = offset, family = poisson())

term.df = dta1 %>% ped_info() %>% add_term(pena4, term = "logpro") %>%
  mutate_at(c("fit", "ci_lower", "ci_upper"), funs(. / .data$logpro)) %>%
  mutate(
    cox4 = coef(vfit4)["logpro"] + coef(vfit4)["tt(logpro)"] * log(tend + 25),
    pamlog4 = coef(pamlog4)["logpro"] + coef(pamlog4)["logpro:logt25"] * log(tend + 25))

ind = find(real2$obstime, term.df$tend)
subreal2 = real2[ind,c(1:10)]
term.df$nplogpro  = subreal2$`log(protime)`
term.df$nplogprose  = subreal2$`log(protime)_SEE`

gg_tv_logpro = ggplot(term.df, aes(x = tend, y = fit)) +
  geom_step(aes(y = nplogpro, col = "NPCox estimation")) +
  geom_stepribbon(aes(ymin = nplogpro - 1.96*nplogprose, 
                      ymax = nplogpro + 1.96*nplogprose), fill = "green", alpha = 0.2) +
  geom_step(aes(col = "PAM with penalized spline")) +
  geom_stepribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2) +
  geom_line(aes(y = cox4, col = "Cox with log-transform")) +
  geom_step(aes(y = pamlog4, col = "PAM with log-transform")) +
  scale_color_manual(name = "Method", values = c(Set1[1:3], "black")) +
  xlab("time") + ylab('beta(t) for log(protime)')
plot(gg_tv_logpro)


##################################### edema #####################################
## cox with log-transform
vfit5 = coxph(
  formula = Surv(time, status) ~ age+logal+logbi+logpro+edema + tt(edema),
  data    = dta,
  tt      = function(x, t, ...) x * log(t + 25))
# coef(vfit2)

## PAM with log-transform
pamlog5 = gam(ped_status ~ s(tend) +age+logal+logbi+logpro+edema+edema:logt25,
              data = dta1, offset = offset, family = poisson())

pena5 = gam(ped_status ~ s(tend)+age+logal+logbi+logpro+edema
            + s(tend, by = edema), data = dta1, offset = offset, family = poisson())

term.df = dta1 %>% ped_info() %>% add_term(pena5, term = "edema") %>%
  mutate_at(c("fit", "ci_lower", "ci_upper"), funs(. / .data$edema)) %>%
  mutate(
    cox5 = coef(vfit5)["edema"] + coef(vfit5)["tt(edema)"] * log(tend + 25),
    pamlog5 = coef(pamlog5)["edema"] + coef(pamlog5)["edema:logt25"] * log(tend + 25))

ind = find(real2$obstime, term.df$tend)
subreal2 = real2[ind,c(1:10)]
term.df$npedema  = subreal2$`edema`
term.df$npedemase  = subreal2$`edema_SEE`

gg_tv_edema = ggplot(term.df, aes(x = tend, y = fit)) +
  geom_step(aes(y = npedema, col = "NPCox estimation")) +
  geom_stepribbon(aes(ymin = npedema - 1.96*npedemase, 
                      ymax = npedema + 1.96*npedemase), fill = "green", alpha = 0.2) +
  geom_step(aes(col = "PAM with penalized spline")) +
  geom_stepribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2) +
  geom_line(aes(y = cox5, col = "Cox with log-transform")) +
  geom_step(aes(y = pamlog5, col = "PAM with log-transform")) +
  scale_color_manual(name = "Method", values = c(Set1[1:3], "black")) +
  xlab("time") + ylab('beta(t) for edema')
plot(gg_tv_edema)



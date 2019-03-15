illustrate_calibration <- function(data,gridcell,timeframe)

d_a_r=fitdistr(data$mod[1,,timeframe,gridcell[1],gridcell[2]],"normal")
d_n_r=fitdistr(data$mod[2,,timeframe,gridcell[1],gridcell[2]],"normal")
d_a_c=fitdistr(data$mod_cal[1,1,,timeframe,gridcell[1],gridcell[2]],"normal")
d_n_c=fitdistr(data$mod_cal[1,2,,timeframe,gridcell[1],gridcell[2]],"normal")

range.values <- range(data$mod[1,,timeframe,gridcell[1],gridcell[2]])

anom_vect<-seq(range.values[1]-2,range.values[2]+2,length=10000)

tau<-quantile(data$mod[2,,,gridcell[1],gridcell[2]],p=0.9,na.rm=TRUE)

p_a_r<-dnorm(anom_vect,mean=d_a_r$estimate[1],sd=d_a_r$estimate[2])
p_a_c<-dnorm(anom_vect,mean=d_a_c$estimate[1],sd=d_a_c$estimate[2])
p_n_r<-dnorm(anom_vect,mean=d_n_r$estimate[1],sd=d_n_r$estimate[2])
p_n_c<-dnorm(anom_vect,mean=d_n_c$estimate[1],sd=d_n_c$estimate[2])

e_a_r<-pnorm(tau,mean=d_a_r$estimate[1],sd=d_a_r$estimate[2],lower.tail=FALSE)
e_n_r<-pnorm(tau,mean=d_n_r$estimate[1],sd=d_n_r$estimate[2],lower.tail=FALSE)
e_a_c<-pnorm(tau,mean=d_a_c$estimate[1],sd=d_a_c$estimate[2],lower.tail=FALSE)
e_n_c<-pnorm(tau,mean=d_n_c$estimate[1],sd=d_n_c$estimate[2],lower.tail=FALSE)

br.palette<-brewer.pal(11,"RdYlGn")
c_all=br.palette[1]
br.palette<-brewer.pal(9,"GnBu")
c_nat=br.palette[5]

postscript('distributions.ps',width=6,height=4)
mp <- par(mfrow=c(2,1),mar=c(2, 4, 2, 2) + 0.1)

for (i in 1:2) { # Loop over raw and calibrated data

  if (i == 1) {
    p_a <- p_a_r
    p_n <- p_n_r
    e_a <- e_a_r
    e_n <- e_n_r
    title <- c('Raw HadGEM3-A')
  } else {
    p_a <- p_a_c
    p_n <- p_n_c
    e_a <- e_a_c
    e_n <- e_n_c
    title <- c('Calibrated HadGEM3-A')
  }

  plot(c(tau,tau),c(0,max(p_a)*0.85),'l',col='black',ylab='Temperature Anomaly',
  'xlab'='',main=title, ylim=c(0,max(cbind(p_a,p_a))*1.1),
  xlim=c(min(anom_vect),max(anom_vect)),yaxt='n',font.main=1)
  lines(anom_vect,p_n,col=c_nat,lwd=2)
  lines(anom_vect,p_a,col=c_all,lwd=2)
  text(tau,max(p_a)*0.95,expression('x'[EX]),cex=0.8)
  FAR <- 1-e_n/e_a
  RR <- e_a/e_n

  text(range.values[2],max(p_a)*0.95,paste("FAR=",round(FAR,2),sep=""),cex=0.8,font=2)
  text(range.values[2],max(p_a)*0.80,paste("RR=",round(RR),sep=""),cex=0.8,font=2)
}

par(mp)
dev.off()

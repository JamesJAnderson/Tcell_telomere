## ======================================================================== 
##                           Supplementary Information_C
## ========================================================================
##    R code to produce Figures 2,3,4,S2,S3
##
##    Short Telomeres and a T-Cell Shortfall in Covid-19: The Aging Effect
##    April 10 2021
##
##    R Core Team (2017). R: A language and environment for
##    statistical computing. R Foundation for Statistical Computing,
##    Vienna, Austria. URL https://www.R-project.org/
## ========================================================================
##   If transferring from pdf document copy without formatting and past 
##   into R Script.
## ========================================================================
# Functions =============================================================
  
    Cx.fun = function(TL_ref,x){
    # Clone size with age for given telomere length (TL) at age 20
    # n values <= 1  indicates no clone expansion
      Cmax = 2^20
      clone = rep(Cmax,length(x)) 
      n = rep(20,length(x))
      TL_x         = TL_x.fun(TL_ref, x, x_ref = 20)
      n            = (TL_x - TL_b)/r
      n            = replace(n,n > 20, 20)
      n            = replace(n,n < 1,  1)
      clone        = 2^n  
      return(clone) 
   }
  
  TL_x.fun <- function(TL_ref, x, x_ref = 20){ 
    # calculates TT for given age x 
    TL_x =  TL_ref - g * (x - x_ref)
    return(TL_x)
  }
  
  Get_TL_ref =function(n,TL_ref_mean = 7.5, TL_ref_sd = 0.6 ){
    # calculates distribution of TL at ref_age  20
    lower_TL = 5.5
    upper_TL = 9.5
    TL_ref = seq(n) + NA
    i = 0
    while (i < n ){
      o = rnorm(1,TL_ref_mean, TL_ref_sd)
      test = o > lower_TL & o < upper_TL
      if(test) {
        i = i + 1
        TL_ref[i] = o
      }
    }
    return(TL_ref)
  }

  XO_fun = function(TL_20) 20 + (TL_20 -TL_b -del_max)/g

# Parameters ============================================================
    x_ref       = 20    # reference year 
    TL_b        = 5     # TL brink in kb
    TL_ref_mean = 7.3   # TL mean for population at x_ref = 20
    TL_ref_sd   = 0.6   # TL standard deviation for population
    g           = 0.03 # TL loss per year (kb/yr)
    r           = 0.07  # TL loss in clone cell replication (kb/replication)
    Cmax        = 2^20  # maximum clone size (MCS)
    del_max     = 1.4   # TL to produce MCS (kb)
    c.scale     = 2^20 # clone size scale
    x = seq(x_ref,90, by = 1) # age vector
    TL_O        = 6.4
# Data =======================================================================
  #	CDC. 2021. Deaths involving coronavirus disease 2019 (COVID-19) 
  # reported to NCHS by sex and age group and week ending date. 
  # https://data.cdc.gov/NCHS/
  # Provisional-COVID-19-Death-Counts-by-Sex-Age-and-W/vsak-wrfu/data
  
  # CB. 2019. 2019: ACS 1-year estimates subject tables. 
  # Survey/Program: American Community Survey TableID:S0101. U.S.Census Bureau 
  # https://data.census.gov/cedsci/all?q=S0101
  #

    age.groups = c(  5,  10,  20,  30,  40,  50,  60,  70,  80,  90)
    
    covid.mort = c(0.00000129,0.00000187, 0.00001653, 
                   0.00006922, 0.00019807, 0.00056254, 
                   0.00138240,  0.00339544,  0.00854008,  0.02440538)
    
    total.mort = c(0.000194694,0.00014776,0.00091595,
                   0.001780577,0.002778419,0.005258215,
                   0.011727834,0.024365616,0.058232908,0.181834281)

    nonCV.mort = total.mort - covid.mort 
    plot(total.mort ,  nonCV.mort)
    out = lsfit(total.mort ,  nonCV.mort)
    abline(out)
    ls.print(out)
    plot(age.groups,total.mort)
    points(age.groups,nonCV.mort,col="red")
    
    ratio = covid.mort[3:10]/covid.mort[3]   
    c.mort = covid.mort[3:10]
    t.mort = nonCV.mort[3:10] # noncovid mortality
    t.ratio = t.mort/t.mort[1]
    c.ratio = c.mort/c.mort[1]
    c_t_ratio = c.mort/ t.mort/(c.mort[3]/ t.mort[3])

# Calculations ==========================================================

Nsamples =1000000 # large value assures smooth curves
C        = matrix(NA, nrow = Nsamples, ncol = length(x))  
TL_ref_random = Get_TL_ref(Nsamples,TL_ref_mean, TL_ref_sd) 
for(i in seq(Nsamples)) C[i,]  = Cx.fun(TL_ref_random[i] ,x)

# C[random sample,age}
# TK_ref_random = random TL length at age 20 
# x = age 20:90

# Figure 2 ABC ==============================================================
  Col   = c("#FF0000","#0070C0","#2C9044") # Col  = c("red","blue","green")
  TL.lab = expression(paste("TL "[20], " (kb)"))
  Plate.line = 3.5
  layout(matrix(c(1,2,3,1,1,1), 3, 1, byrow = TRUE),
         widths = c(1,1,1), heights =c(.8,1,1))
  layout.show(3)
  Cex = 1
  Cex.lab = 1.1
  Cex.axis = 1
  Cex.plate = 1.7
  Lwd = 2
  mx = 2.1
  Label = expression(paste('Clone size (10'^" 6"," T cells)"))
  
  par(cex = 1, mar=c(3.5,4.7,2.5,5),las = 1 ,cex = 1, cex.lab=1.05,  xpd = F)
  # Figure 2A ....................................................  
  { Lwd.A  = 3  
    omean= mean(TL_ref_random)
    olow = omean - var(TL_ref_random)^0.5
    ohi  = omean + var(TL_ref_random)^0.5
    O    = hist(TL_ref_random,breaks = 50, plot=F) 
    height.omean = which(O$density == max(O$density))
    ilow = min(which(O$mids > olow))
    ihi  = max(which(O$mids < ohi))
    
    plot(O$mids,O$density,xlab = "",ylab="Density",
         cex.lab =Cex.lab,lwd = Lwd,type ="l",cex.axis = Cex.axis)
    lines(c(olow,olow),c(0,O$density[ilow])*.93, lwd =  Lwd.A, col=Col[1] )
    lines(c(omean,omean),c(0,O$density[height.omean]),lwd= Lwd.A,col=Col[2])
    lines(c(ohi,ohi),c(0,O$density[ihi])*.92, lwd =  Lwd.A, col = Col[3])
    mtext("a",2,line = Plate.line, at = .8,cex=Cex.plate)
    mtext(TL.lab,1,line = mx,cex=Cex.lab)
    
  }
  
  # Figure 2B ....................................................
  {
    X = c(XO_fun(olow),XO_fun(omean),XO_fun(ohi))
    Xp = print(X , digits = 0)  
    TLO.lab = expression(paste("TL"[O]," = 6.4"))
    TLB.lab = expression(paste("TL"[B]," = 5"))
    XO1.lab = expression(paste("X"[O],"  = 30"))
    XO2.lab = expression(paste("X"[O],"  = 50"))
    XO3.lab = expression(paste("X"[O],"  = 70"))
    
    x1 = 15:90 ;   x2 = 15:95
    plot(x2,TL_x.fun(omean,x2),type="l",ylim = c(4.5,8),xlim = c(20,90),
         ylab = TL.lab, xlab = "", col=Col[2],
         cex.lab = Cex.lab, cex.axis = Cex.axis, lwd = Lwd,lty = 1)
    lines(x1,TL_x.fun(olow ,x1) , col=Col[1],lwd = Lwd,   lty = 1)
    lines(x2,TL_x.fun(ohi  ,x2) , col=Col[3],lwd = Lwd,   lty = 1)
    rect(17.1,4.38,103.2,5,col="lightgray",border = NA)
    
    axis(3, at = XO_fun(c(olow,omean,ohi)),
         labels = F,col.ticks = "black",lwd = 0,lwd.ticks=1,tcl = -.3) 
    axis(4, at = c(6.4,5), labels = F,col.ticks = "black",
         lwd = 0, lwd.ticks = 1,tcl = -.3) 
    
    mtext(TLO.lab,side= 4, line= .5, at = TL_O, col="black",cex=Cex)
    mtext(TLB.lab,side= 4, line= .5, at = 5,    col="black",cex=Cex)
    mtext("b",2,line = Plate.line,at = 8.7,cex =Cex.plate)
    mtext("Age (years)",1,line = mx,cex =Cex.lab)
    
    abline(h=6.4,lty = 3,lwd = Lwd)
    
    lines(x = c(X[1],X[1]),y =c(8.1,6.4),col = Col[1],lty = 3,lwd = Lwd)
    lines(x = c(X[2],X[2]),y =c(8.1,6.4),col = Col[2],lty = 3,lwd = Lwd)
    lines(x = c(X[3],X[3]),y =c(8.1,6.4),col = Col[3],lty = 3,lwd = Lwd)
    
    par(xpd = NA)
      yy =8.4
      text(X[1], yy, XO1.lab, cex = Cex, col=Col[1],lwd = Lwd)
      text(X[2], yy, XO2.lab, cex = Cex, col=Col[2],lwd = Lwd)
      text(X[3], yy, XO3.lab, cex = Cex, col=Col[3],lwd = Lwd)
    par(xpd = F)
  }
  
  # Figure 2C ....................................................
  {
    XO.lab = expression("X "[O])
    jitter = .015
    
    plot(x2,Cx.fun(omean,x2)/1e6,type="l",
         ylab = Label, xlab = "", col=Col[2],lwd = Lwd,lty = 1
         ,xlim = c(20,90),cex.lab = Cex.lab,cex.axis = Cex.axis)
    lines(x2,Cx.fun(olow ,x2)/1e6 - jitter, col=Col[1],lwd = Lwd,   lty = 1)
    lines(x2,Cx.fun(ohi  ,x2)/1e6 + jitter, col=Col[3],lwd = Lwd,   lty = 1)
    axis(3, at = XO_fun(c(olow,omean,ohi)),
         labels = F,col.ticks = "black", lwd = 0, lwd.ticks = 1,tcl = -.3) 
    abline(v=c(XO_fun(olow),XO_fun(omean),XO_fun(ohi)),col=Col,lty=3,lwd = Lwd)
    mtext("c",2,line = Plate.line,at = 1.2,cex =Cex.plate)
    mtext("Age (years)",1,line = mx,cex =Cex.lab)
    
    par(xpd = NA)
      yy =1.16
      text(X[1], yy, XO1.lab, cex = Cex, col=Col[1])
      text(X[2], yy, XO2.lab, cex = Cex, col=Col[2])
      text(X[3], yy, XO3.lab, cex = Cex, col=Col[3])
    par(xpd = F)
  }

# Calculating age at % drop ==================================================
CS_drop =  cbind(Cx.fun(olow,x2)/max(Cx.fun(olow,x2)),
                 Cx.fun(omean,x2)/max(Cx.fun(omean,x2)),
                 Cx.fun(ohi,x2)/max(Cx.fun(ohi,x2)))
print(cbind(x2,floor(CS_drop*10000)/10),digits = 3)

# Figure 3 AB =================================================================

 library("colorspace")

## panel A histograms of TL frequencies   
  layout(matrix(c(1,2,3,4,5,6,7,8),nrow = 2, ncol = 4, byrow = TRUE),
         widths = c(1,1,1,1), heights =c(1,1,1,1))
  layout.show(8)
  
  breaksTL = seq(4,9.6,by = 0.2) 
  Color.fig3 = rep(NA,length(breaksTL)-1) 
  Cex = .9
  Cex.lab = 1 
  Cex.axis = 1
  Lwd = 2 
  TL_O = 6.4             # age onset TL
  Age1 = c(20,30,50,70)  # ages for histograms
  
  par(cex = .9, cex.lab=Cex.lab , cex.axis = Cex.axis) 
  par(   mar=c(3.5,2,2,0), oma = c(1,2,1,0.5))

  for(i in 1:4){
        TL_x_random = TL_x.fun (TL_ref_random, Age1[i], x_ref = 20)
        OO = hist(TL_x_random, breaks = breaksTL, plot = F)
        
        for(ix in seq(length(Color.fig3))){
          if(OO$mids[ix] <= TL_O)Color.fig3[ix] = "#7F7FC1" # blue 
          if(OO$mids[ix] >  TL_O)Color.fig3[ix] = "#BF7F7F" # brown
        }
        
        if(i == 1){ 
          hist(TL_x_random, xlim=c(4. ,10),ylim = c(0,.7), main = "", 
               freq=F, breaks = breaksTL, xlab = "",axes = T,
               ylab = "",  col=Color.fig3, cex.axis=Cex.axis,
               cex.lab = Cex.lab)
          mtext(paste("Density"), side = 2, line = 2.7,
                las = 3, cex = Cex.lab)
          mtext("a", side = 2, line = 2.7, las = 1,
                cex = Cex.plate, at = .8)
        }
        
        if(i >1){
          hist(TL_x_random, xlim=c(4.,10),ylim = c(0,.7), main = "", 
               freq=F, breaks = breaksTL, xlab = "",axes = T,
               ylab = "",  col=Color.fig3)
          rect(2.6,-.037, 3.3  ,.8, col= "white" ,border="NA", par(xpd=NA))
        }
        
        mtext(paste("Age = ", Age1[i]), side = 3, line =1)
        mtext("TL (kb)", side = 1, line = 2.7,    las = 1)        
  }

## panel B histogram of Clone size
  breaksC = seq(0, 1.0,by =.1) 
  Color.fig3 = rep(NA,length(breaksC)-1) 
  cut.off = .9
  Label = expression(paste('Clone size (10'^" 6"," T cells)"))

  for(i in c(1,11,31,51)){
      OO = hist(C[,i]/Cmax, breaks = breaksC, plot = F)
      
      for(ix in seq(length(Color.fig3))){
        if(OO$breaks[ix] <  cut.off ) Color.fig3[ix]="#7F7FC1" # blue 
        if(OO$breaks[ix] >= cut.off ) Color.fig3[ix]="#BF7F7F" # brown 
      }
      
      hist(C[,i]/Cmax,freq=F,breaks = breaksC, xlab = "",axes = T, 
           ylab = "", main= "",col = Color.fig3,
           xlim = c(0,1.05),ylim = c(0,10), 
           cex.axis = Cex.axis, cex.lab = Cex.lab)
      
      if(i > 1) rect(-.25,-.6, -.12 ,11 , col= "white",
                    border="NA", par(xpd=NA))
      
      if(i == 1) {
        mtext(paste("Relative frequency"), side = 2, line = 2.7,
              cex = Cex.lab,las = 3)
        mtext("b", side = 2, line = 2.7, las = 1, at = 11,
              cex = Cex.plate,)
      }
      
      mtext(Label, side = 1, line = 2.7, las = 1,cex=Cex.lab)        
 }  

 # rbind(x,mLCS)      

# Figure 4 ABC ===============================================================

## Get clone size < MCS ---    
  { c.mean =  x + NA
a.mean = x + NA # all population clone size

cutoff = (2^20)*1        # this removes all Cmax individuals
for(i in seq(x)) {
  c.mean[i] = mean(C[C[,i] < cutoff,i])
  a.mean[i] = mean(C[,i])
  
} 
mLCS = c.mean/c.scale 
aCS = a.mean/c.scale 

}

## Plot information---
  {  
    par(mfrow =c(1,3), mar=c(4,3.5,2,2.5), oma = c(2,1,1,0),
        cex.lab=1., cex=1., las = 1, xpd = F)
    
    ylabAD  = expression("Hazards ratio  "[20] )
    xlabAB ="Age (years)"
    clone_lab.x = expression(paste('Mean LCS (10'^"6"," T cells)"))
    clone_lab.y = expression(paste('Mean  LCS (10 '^"6"," T cells)"))
    xx = c(1,11,21,31,41,51,61,71)
    
    c.scale   = 1e6
    line.col  = "#6666FF"  #   blueish
    line.col2 = "#CC3399" #   redish
    
    Cex = .9
    Cex.lab = 1.05 
    Cex.axis = 1
     
    li        = 3.2
    Lwd       = 3
    age       = c(20, 30, 40, 50, 60, 70, 80, 90)
  }

## Panel A ---
  {
    plot(age,c.ratio,xlab = xlabAB,ylab="", pch = 1, type = "n",
         cex.axis = Cex.axis,  cex.lab = Cex.lab)
    
    out = lsfit(age,log(c.ratio))
    ls.print(out)
    mortality.c = exp(coef(out)[1]+ coef(out)[2]*x)
    lines(x,mortality.c, col = line.col, lwd = Lwd)
    mtext(ylabAD,side = 2,line = li, las = 0,cex = Cex.lab)
    points(age,c.ratio, pch = 16, cex = 1.2, col= line.col)
    
    out = lsfit (age,log(t.ratio))    #; ls.print(out)
    mortality.t = exp(coef(out)[1]+ coef(out)[2]*x)
    lines(x,mortality.t, col = line.col2, lwd = Lwd)
    points(age,t.ratio, pch = 16, cex = 1.2, col= line.col2) 
  }

## Panel labels for all 3 graphs
  {
    par(xpd = NA)
    yy = 1850
    text( -5 , yy, "a", cex = Cex.plate)
    text( 105,yy, "b", cex = Cex.plate)
    text( 215,yy, "c", cex = Cex.plate)
    par(xpd = F)
  }

## Panel B ---   
  {
    plot(x, mLCS, type="l", 
         ylim = c(0,max(mLCS)), 
         xlim=c(20,90), ylab = "", xlab=xlabAB, 
         col = line.col, lwd = Lwd,cex.axis = Cex.axis,cex.lab = Cex.lab)
    mtext(clone_lab.y,side = 2,line = 3, las = 0, cex.lab = Cex.lab)
  }

## Panel C ---
  {
    plot(rev(mLCS), rev(mortality.c), type ="l", xlab=clone_lab.x, 
         ylab="", xlim=c(max(mLCS),0), ylim=c(0,1500),
         col = line.col, lwd = Lwd,cex.axis = Cex.axis, cex.lab = Cex.lab)
    mtext(ylabAD,side = 2,line = li, las = 0,cex = Cex.lab)
    lines(rev(mLCS), rev(mortality.t), col=line.col2, lwd = Lwd) 
    
    x.lines  = c(1,11,21,31,41,51,61,71,81)
    vertical.col="darkgoldenrod"
    x.at.LCS = rev(c(20,30,40,50,60,70,80,"","" ))
    x.LCS = rev (mLCS[x.lines])
    
    abline(v=mLCS[31],lty= 2, lwd=2, col= vertical.col)
    axis(side = 3,col="black",at=x.LCS,labels=F,line = 0)
    mtext("Age (years)",side = 3,line = 2., las = 0,cex = Cex.lab)
    mtext(x.at.LCS,side = 3,line = .7, las = 0,at = x.LCS ,cex = Cex.lab)
  }

#rbind(x, mLCS)

# Figure S2 ==========================================================
# calculations     
{ Cutoff = c(1,.5,.15)  
fraction.LCS = matrix(NA,length(Cutoff),length(x))
C.mean       = matrix(NA,length(Cutoff),length(x))
mLCS         = matrix(NA,length(Cutoff),length(x))
col3 = c("black","red","blue")
Ldw = 2 
yoff  = 1.1
for(ii in seq(Cutoff)){
  cutoff.test = (2^20)*Cutoff[ii]#  removes Cmax individuals
  for(i in seq(x)) {
    test =   C[,i] < cutoff.test 
    fraction.LCS[ii,i] = length(which(test == T))/Nsamples
  }
  for(i in seq(x)) C.mean[ii,i] = mean(C[C[,i] < cutoff.test,i]) 
  mLCS[ii,] = C.mean[ii,]/c.scale 
}}

# plot information     
{ par(mfrow =c(1,3), mar=c(4,4.3,1.5,1), oma=c(1,1,1 ,1),
      cex.lab=1., cex=1., las = 1, xpd = F)
  
  Cex = .9
  Cex.lab = 1.05 
  Cex.axis = 1
  Cex.plate = 1.7
  li        = 3.2
  Lwd       = 2
  age       = c(20, 30, 40, 50, 60, 70, 80, 90)
  plate.x = c(-5,-5,.37)
}

# Panel A  
{  plot(x, fraction.LCS[1,]*100,type="l", ylim=c(0,100),
        xlab="Age (years)",ylab ="% LCS",lwd = Lwd,
        cex.axis = Cex.axis, cex.lab = Cex.lab)
  lines(x, fraction.LCS[2,]*100,col="red",lwd = Lwd) 
  lines(x, fraction.LCS[3,]*100,col="blue", lty = 1,lwd = Lwd)
  mtext("a",3,line =  .5, at= plate.x[1], cex = Cex.plate)
}

# Panel B  
{   plot(x,mLCS[1,],type="l",
         xlim=c(20,90),ylim=c(0,max(mLCS[1,])),
         xlab ="Age (years)",ylab =clone_lab.y,lwd = Lwd,
         cex.axis = Cex.axis, cex.lab = Cex.lab) 
  lines(x,mLCS[2,],col="red",lwd = Lwd) 
  lines(x,mLCS[3,],col="blue", lty = 1,lwd = Lwd) 
  mtext("b",3,line =  .5, at= plate.x[2], cex = Cex.plate)
  legend("topright",c("LCS < 1 MCS","LCS < 0.5 MCS",
                      "LCS < 0.15 MCS"),col = col3, lty=c(1,1,1),
         text.col = col3,bty = "n",cex = 0.8)
  #  abline(v=50); abline( h=0.13119190)
}   
# Panel C   
{ plot(rev(mLCS[1,]) ,rev(mortality.c),type="l",
       ylim = c(0,1500), xlim = rev(c(0,max(mLCS[1,]))),
       ylab = ylabAD,xlab = clone_lab.x,lwd = Lwd,
       cex.axis = Cex.axis, cex.lab = Cex.lab) 
  lines(rev(mLCS[2,]) ,rev(mortality.c),col="red",lwd = Lwd) 
  lines(rev(mLCS[3,]) ,rev(mortality.c),col="blue",lty=1,lwd = Lwd)
  mtext("c",3,line =  .5, at= plate.x[3],cex = Cex.plate)
}

# Figure S3 ==========================================================

# Functions   
{
  sigma_xo = function(TL20 ) 
    ((sigma.B^2 
      + sigma.delMax^2
      + ((TL20 - TL_b - del_max)/g*sigma.g)^2)/g^2
    )^0.5
  
  sigma_TLO  = function()(sigma.B^2+sigma.delMax^2)^0.5
  
  sigma_N = function(age,TL_20=7.3)  
    (
      (
        sigma.B^2 + 
          ((age - 20)*sigma.g)^2 + 
          (N_fun(age,TL_20)*sigma.r)^2
      )/r^2
    )^0.5
  
  
  N_fun = function(age,TL_20)(TL_20 - TL_b - g * (age-20))/r 
  
  sigma_LCS = function(age,TL_20=7.3) log(2)*2^N_fun(age,TL_20)*
    sigma_N(age,TL_20) 
  
  sigma_MCS  = function(XO = 50) log(2)*2^20*sigma_N(XO)
}
# Data 
{ 
  sigma.B = 0.017/5   # kb CV = 1.74% mean LT 6.2
  sigma.delMax = 0.1  # kb  
  sigma.r = 0.02      # kb/replication
  sigma.g = 6.5e-04     # kb/year 0.02*30/1000  30*.03/1000
  TL_20 = TL_ref_mean # 7.3 kb
  xLCS = 50:90        # age range of LCS
  xMCS = 20:50        # age range of MCS
  # from main figures
  # del_max = 1.4 # kb 
  # TL_b = 5      # kb
  # r = 0.07      # kb/doubling
  # g = 0.03      # kb/year
  # x = 20:90     # year
  # Cmax = 2^20   # maximum clone size (kb)
  
  # for calculating g and sigma.6. Steenstrup 2013
  g.obs =c(31.3,40.7,31.4,33.5,32.2,25.2,30.8)/1000
  g.nobs=c(70,635,271,271,271,609,80)
  g.cv = c(1.5,1.4,2.4,2.4,2.4,2.2,2.8)/100 
  
  meang = sum(g.obs*g.nobs)/sum(g.nobs)
  sigma.g = sum(g.cv*g.obs*g.nobs)/sum(g.nobs)
  sigma_TLO()
}

# Plot information     
{ par(mfrow =c(1,3), mar=c(4.2,4.3,1.1,1), oma=c(1,1,1 ,0),
      cex.lab=1., cex=1., las = 1, xpd = F)
  Cex = 1
  Cex.lab = 1.2 
  Cex.axis = 1
  Cex.plate = 1.7
  li        = 3.2
  Lwd       = 2
  age       = c(20, 30, 40, 50, 60, 70, 80, 90)
  plate.x = c(3.5,-.7,1)
}

# Panel S3a  
{  
  Nmax = 20
  ylabA = expression(sigma[XO])
  ylabB = "(years)"
  xlabA = expression(TL[20])
  tl20 = seq(5,10,length.out = 50)
  sigmaxo = sigma_xo(tl20)
  plot(tl20, sigmaxo ,type="l",
       xlab=xlabA,ylab ="",lwd = Lwd,
       cex.axis = Cex.axis, cex.lab = Cex.lab, col="blue")
  mtext("a",3,line =  .5, at= plate.x[1], cex = Cex.plate)
  mtext(ylabA,side = 2,line = 3., cex = 1.5, las = 0, at = 3.65)
  mtext(ylabB,side = 2,line = 3.2, cex = 1 , las = 0, at = 3.9)
}

# Panel S3b  
{  
  Nmax = 20
  ylabA = expression(sigma[N/Nmax])
  sigmaN = sigma_N(x,TL_20)
  sigmaN[sigmaN > sigmaN[x == 50]] = sigmaN[x == 50]
  plot(x, sigmaN/Nmax ,type="l", ylim=c(0,max(sigmaN/Nmax)),
       xlab="Age (years)",ylab ="",lwd = Lwd,
       cex.axis = Cex.axis, cex.lab = Cex.lab, col="blue")
  mtext("b",3,line =  .5, at= plate.x[2], cex = Cex.plate)
  mtext(ylabA,side = 2,line = 3.1, cex = 1.5, las = 0)
}

# Panel S3C  
{  
  Cmax = 2^20
  ylabA = expression(sigma[CS/MCS])
  sigmaCS = sigma_LCS(x,TL_20)
  sigmaMCS = sigma_MCS()
  sigmaCS[sigmaCS > sigmaMCS] = sigmaMCS
  plot(x, sigmaCS/Cmax ,type="l", ylim=c(0,sigmaMCS/Cmax),
       xlab="Age (years)",ylab ="",lwd = Lwd,
       cex.axis = Cex.axis, cex.lab = Cex.lab, col="blue")
  mtext("c",3,line =  .5, at= plate.x[3], cex = Cex.plate)
  mtext(ylabA,side = 2,line = 2.4, cex = 1.5, las = 0)
}


##================================= END R CODE=============================
##
##

 

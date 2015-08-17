#### FUNCTION: Generate sequencies for naming folders.

numSeqGen=function(strTempl=NULL,first,last,algMax=5,sepStrNum=""){
  numSeq=sprintf(paste("%0",algMax,"d",sep=""),first:last)
  fullSeqName=paste(strTempl,sepStrNum,numSeq,sep="")
  return(fullSeqName)
}

mkLands=function(){
  
}
#### FUNCTION : run multiple replicates of a simulation given by a specific parameters combination.

multiReps=function(PATH,nrep=10,b0,d0,m0,inc.b,inc.d,step,radius,dens.t,config,N0,lands=land,tm=100,
                   ini.seq=1,death.mat=1){
  oldpath <- getwd()
  dir.create(PATH, recursive=TRUE)
  setwd(PATH)
  for(i in 1:nrep){
    TWoLife(taxa.basal=b0,
            taxa.morte=d0,
            move=m0,
            incl.birth=inc.b,
            incl.death=inc.d,
            passo=step,            
            raio=radius,
            density.type=dens.t,
            ini.config=config,
            #####################
            N=N0,
            AngVis=360,
            death.mat=death.mat,
            landscape=lands,
            tempo=tm,            
            out.code=ini.seq-1+i)
  }
  setwd(oldpath)
}

#### FUNCTION : run multiple paramters combinations, each one with an specific number of replicates.

multiRun=function(params,combs,pathsList=NULL,pathBase="~/Desktop/",algMax=3,strTempl="Comb-",
                  nrep=20,tmax=100,land){
  if(is.null(pathsList)){
    PATH=paste(pathBase,numSeqGen(strTempl=strTempl,min(combs),max(combs),algMax=algMax),sep="")
  } 
  else PATH=pathsList
  
  for (i in combs){
    multiReps(PATH=PATH[i],
             nrep=nrep,
             b0=params$b0[i],
             d0=params$d0[i],
             m0=params$movem.rate[i],
             inc.b=params$incl.birth[i],
             inc.d=params$incl.death[i],
             step=params$step.length[i],
             radius=params$radius[i],
             dens.t=params$dens.type[i],
             config=params$initial.config[i],
             N0=params$N0[i],
             lands=land,
             tm=tmax)
  }
}
#@@@@@@ End of functions for generating raw data 

##################################
#  Get N (t) and/or Speed (t) data
##################################

# For a set of parameters combination (i.e. parameter space sample). FULL DATASET
getDataPS=function(DIR,tmax,nrep,FUNCS=c("getNt","getSpeed"),ind=10){# DIR is the folder containing the outputs of each parameter combination of the evaluated parameter space
  old.path=getwd()
  setwd(DIR)
  combsOutFolders=dir()
  sizePS=length(combsOutFolders)
  arrNt=arrSpeed=array(0,c(tmax+1,nrep+1,sizePS),dimnames=list(NULL,c("Time",paste("R",1:nrep,sep="")),combsOutFolders))
  matExtTimes=matrix(0,nrep,sizePS,dimnames=list(paste("R",1:nrep,sep=""),combsOutFolders)) # toDo: look for better column names
  vecHasExt=numeric(sizePS)
  
  for(i in 1:sizePS){
    out1=getNSp(wdir=paste("./",combsOutFolders[i],sep=""),tmax=tmax,nrep=nrep,ind=10,FUNCS=FUNCS)
    if("getNt" %in% FUNCS) arrNt[,,i]=out1$Nt
    if("getSpeed" %in% FUNCS) arrSpeed[,,i]=out1$speed
    matExtTimes[,i]=out1$extincTimes
    vecHasExt[i]=out1$hasExt
  }
  setwd(old.path)
  return(list(Nt=ifelse("getNt" %in% FUNCS,arrNt,NA),speed=ifelse("getSpeed" %in% FUNCS,arrSpeed,NA),extincTimes=matExtTimes,hasExt=as.logical(vecHasExt)))
}
# For multiple replicates of a given parameter combination
getNSp=function(wdir=getwd(),tmax=50,nrep=20,ind=10,FUNCS=c("getNt","getSpeed")){ # BUG 1: Não funciona se tmax < que o tempo maximo da simulação.
  # Folder's paths management
  oldpath <- getwd()
  setwd(wdir)
  files=dir()
  # variables
  tExt=NA
  ### matrices template for Nt and for front speed data of each replicate 
  matNt=matSpeed=matrix(0,tmax+1,(nrep+1),
                        dimnames=list(paste(1:(tmax+1)),c("Time",numSeqGen("R",1,nrep,nchar(nrep)))))
  matNt[,1]=matSpeed[,1]=0:tmax
  # getting data
  for(i in 1:nrep){
    output2L=read.table(files[i],sep=" ") 
    if(is.na(output2L[dim(output2L)[1],2])){ # Is there extinction?
      if(output2L[dim(output2L)[1],1]<=tmax) 
        tExt[i]=output2L[dim(output2L)[1],1]
      output2L=output2L[-(dim(output2L)[1]),] # Remove last row, which indicates population extinction. 
    }
    else tExt[i]=NA
    # do calculations
    if("getNt"%in%FUNCS){ 
      cat(paste(i," ")) #### Remove this line after TESTS
      Nt=getNt(output2L=output2L,tEx=tExt[i],tmax=tmax)
      matNt[1:dim(Nt)[1],i+1]=Nt[,2]
    } 
    else matNt=NA 
    if("getSpeed"%in%FUNCS){
      speed=getSpeed(output2L=output2L,ind=ind,tEx=tExt[i],tmax=tmax)
      matSpeed[2:(length(speed)+1),i+1]=speed
    }
    else matSpeed=NA
  }
  cat("\n") #### Remove this line after TESTS
  if(sum(!is.na(tExt))>0) EXT=1 else EXT=0
  setwd(oldpath)
  return(list(Nt=matNt,speed=matSpeed,extincTimes=tExt,hasExt=as.logical(EXT)))   
}

### Aditional functions = Get data (N(t) or speed(t)). Used in getNSp, but can be used outside this functio if proper
### "output2L" is given.
getNt=function(output2L,tEx,tmax){
  Nt=as.data.frame(table(output2L[,1]))
  Nt[,1]=as.integer(levels(Nt[,1]))
  colnames(Nt)=c("times","Nt")
  MAXt=max(Nt[,1])
  
  Nt=rmTGaps(d.f=Nt,MAXt=MAXt,tEx=tEx)
  if(MAXt>tmax) Nt=Nt[1:(tmax+1),]
  
  return(Nt)
}

#
getSpeed=function(output2L,ind=10,tEx,tmax){
  displacements=displacement(output2L) # Individual displacement
  new.output2L=data.frame(output2L[,1],displacements) # Two columns data.frame (1: time, 2: displacement)
  front=aggregate(new.output2L[,2],by=list(new.output2L[,1]),FUN=frontDist,ind=ind) # Calculates the front distance based on the displacement vector
  colnames(front)=c("times","speed")
  MAXt=max(front[,1])
  
  front=rmTGaps(d.f=front,MAXt=MAXt,tEx=tEx)
  speed=diff(front[,2])
  if(MAXt>tmax) speed=speed[1:tmax]
  return(speed)
}

# Function to edit the data table in order to fill time gaps if they exist.
rmTGaps=function(d.f,MAXt=MAXt,tEx){# d.f = data.frame ==> Equivalent of front or Nt on previous version of the code.
  nROW=dim(d.f)[1]
  #if(nROW!=(MAXt+1)){# There are time gaps from begining until the last recoreded time before extinction
  chkMissing=data.frame(times=0:MAXt)
  d.f=merge(chkMissing,d.f,all.x=T,all.y=T)
  missing=sort(which(is.na(d.f[,2])),decreasing=T)
  for(j in missing){
    d.f[j,2]=d.f[j+1,2]
  }
  #}
  colnames(d.f)=c("times","Var")
  if(!is.na(tEx)&floor(tEx)>MAXt) { # Is there time gaps from the last recorded time to the extinction time?
    d.f=merge(data.frame(times=(MAXt+1):floor(tEx),Var=d.f[MAXt+1,2]),d.f,all.x=T,all.y=T) # Fill gaps
    #MAXt=floor(tEx)
  }
  return(d.f)
}

#@@@@@ End of functions for getting analytical datasets

#######################################
#  Simulation Summary and Ancillary functions
#######################################

# TO DO: ampliar para receber DATA = um objeto da função getDataPS. Muda porque os elementos $Nt ou $speed do objeto getDataPS 
# é um array3D e não uma matriz; o elemento $extincTimes é uma matriz e não um vetor e o elemento $hasExt é um vetor de comprimento > 1.
getMeanNSp=function(DATA,Type=c("n","v"),exclude=c("none","zeros","ifExt")){# DATA has a specific format: an object returned by getNSp
  if(Type=="n") DATA.=DATA$Nt
  else DATA.=DATA$speed
  dimData=dim(DATA.)
  isExt=which(DATA$extincTimes>0) # if there is no extinction, returns an empty integer vector
  
  if(exclude=="none"){
    Means=rowMeans(DATA.[,2:dimData[2]])
    return(Means)
  } 
  else {
    if(exclude=="zeros"){
      if(Type=="n"){
        DATA.[,2:dimData[2]][DATA.[,2:dimData[2]]==0]=NA # Convert 0's into NAs
        Means=rowMeans(DATA.[,2:dimData[2]],na.rm=T) 
        return(Means)
      } else {
        for(i in 2:dimData[2]){
          if(!is.na(isExt[i-1])){
            DATA.[ceiling(isExt):dimData[1],i]=NA
          }
        }
        Means=rowMeans(DATA.[,2:dimData[2]],na.rm=T)
        return(Means)
      }
    }
    else { # ifExt
      #mNt=rep(NA,dimData[1]-1)
      Means=rowMeans(DATA.[,-c(1,isExt)])
      #mNt[-(isExt-1)]=mNtShort      
      return(Means)
    }
  }       
}

# 
fitApred=function(DATA,anPred,p.crit=0.05,tInterval=NULL,modelType=NULL,method=c("lastF","lastMLH","lm","MS-AIC")){
  # BUG 2: methods lastMLH and MS-AIC need corrections
  # DATA: a data.frame with the independent variable on 
  # the first column and replicates (of the response variable) on the other columns. (e.g. getNSp.output$Nt, getNSp.output$speed)
  # TO DO: ampliar para receber DATA = um objeto da função getDataPS. Muda porque os elementos $Nt ou $speed do objeto getDataPS 
  # é um array3D e não uma matriz; o elemento $extincTimes é uma matriz e não um vetor e o elemento $hasExt é um vetor de comprimento > 1.
  last=DATA[dim(DATA)[1],-1]
  # method 1: Frequentist ==> last point
  if("lastF" %in% method){
   
    tTest=t.test(last,mu=anPred)
    #if(tTest$p.value<=p.crit){} # meanDiffMet1=tTest$estimate-anPred else meanDiffMet1=0  
    return(tTest)
  }
  # method 2: logLH test ==> last point
  if("lastMLH" %in% method){
    M1=mle2(last~dnorm(m=anPred,sd=sd(last)),start=list(m=anPred),data=list(last))
    M2=mle2(last~dnorm(m=mn,sd=sd(last)),start=list(mn=mean(last)),data=list(last))
    estMean=coef(M2)
    if(anPred<confint(M2)[1] & anPred>confint(M2)[2]){} # meanDiffMet2=tTest$estimate-anPred else meanDiffMet2=0  
    LR=-2*(logLik(M1)-logLik(M2))
  }
  
  if("lm" %in% method){
    DATA.=data.frame(IV=rep(tInterval,dim(DATA[,-1])[2]),RV=as.vector(as.matrix(DATA[tInterval,-1])))
    M1=lm(DATA.$RV~DATA.$IV)
    return(list(summary(aov(M1)),M1))
  }
  # method 4: model selection / AIC ==> time interval for estimating parameters.
  if("MS-AIC" %in% method){
    DATA.=data.frame(IV=rep(tInterval,dim(DATA[,-1])[2]),RV=as.vector(as.matrix(DATA[tInterval,-1])))
  
    
    M1=mle2(RV~dnorm(m=300,sd=sd(RV)),start=list(m=anPred),data=DATA.)
    M2=mle2(RV~dnorm(m=mn,sd=sd(RV)),start=list(mn=mean(DATA.$RV)),data=DATA.)
    M3=mle2(RV~dnorm(m=b0+b1*IV,sd=sd(RV)),start=list(b0=anPred,b1=0.0001),data=DATA.) # linear
    nObs=dim(DATA.)[1]
    tableAIC=ICtab(M1,M2,M3,type=ifelse(nObs/2<40,"AICc","AIC"),nobs=nObs,delta=T,weights=T,sort=T)
    return(list(tableAIC,M1,M2,M3))
    
  }
}

# Calculates de extinction probability of a population undergoing a TWoLife dynamics with a given parameters combination  
# This probability is given by the proportion of replicates in which extinction occured during the time inteval of the simulation
# (i.e. tmax)
getExtP=function(DATA){ # DATA = getNSp output
  nrep=dim(DATA$Nt[,-1])[2]
  pExt=ifelse(!DATA$hasExt,0,sum(DATA$extincTimes/DATA$extincTimes,na.rm=T)/nrep)
  return(pExt)
}

# This function get the probability of extinction data from many simulations run according to 
# a set of parameters combinations (i.e. sample of a parameters space)
getExtPPS=function(DATA){# DATA: an object from getDataPS
  dims=dim(DATA$Nt)
  pExt=numeric(dims[3])
  for(i in 1:dims[3]){
    pExt[i]=ifelse(!DATA$hasExt[i],0,sum(DATA$extincTimes[,i]/DATA$extincTimes[,i],na.rm=T)/(dims[2]-1))
    # getExtP(DATA$Nt[,,i])
  }
  names(pExt)=colnames(DATA$extincTimes)
  cat("\n *** Extinction Probabilities *** \n \n")
  return(pExt)
}

##### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##### Make table for results analysis
##### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

eqFrom2L=function(DATA,modelType){# DATA columns must have proper names (see .... for allowed names) according to the TWoLife parameters, initial conditions and analytical predictions
  newTab=switch(modelType,
         SG=data.frame(r=rCalc(DATA),Nt=NKCalc(DATA,F)),
         DDG=data.frame(r=rCalc(DATA),sumSlopes=DATA$bSlope+DATA$dSlope,K=NKCalc(DATA)),
         RW=data.frame(D=DCalc.2D(DATA),Speed=speedCalc(DATA,T)), 
         RDsk=data.frame(r=rCalc(DATA),Nt=NKCalc(DATA,F),D=DCalc.2D(DATA),Speed=speedCalc(DATA)),
         RDfk=data.frame(r=rCalc(DATA),sumSlopes=DATA$bSlope+DATA$dSlope,densK=NKCalc(DATA)$densK,
                         locArea=NKCalc(DATA)$locArea,Kl=NKCalc(DATA)$Kl,Kg=NKCalc(DATA)$Kg,D=DCalc.2D(DATA),
                         Speed=speedCalc(DATA)),
         stop(paste("modelType",modelType,"not included",sep=" ")))
  return(newTab)
  # In cases SG and RD, Nt and speed is a function of time. Thus, to summarize these variable in a single value, 
  # we used the mean among replicates at the end of the simulation. 
}

rCalc=function(DATA){
  r=DATA$b0-DATA$d0
  return(r)
}
DCalc.2D=function(DATA){
  D=DATA$w*DATA$step^2/4
  return(D)
}

NKCalc=function(DATA,is.DD=T){
  if(!is.DD){
    Nt=DATA$N0*exp(rCalc(DATA)*as.integer(DATA$tmax))
    return(Nt)
  }
  else {
    if("R" %in% colnames(DATA)){
      densK=rCalc(DATA)/(DATA$bSlope-DATA$dSlope)
      percepArea=pi*DATA$R*DATA$R
      Kl=percepArea*densK
      Kg=densK*DATA$landSize
      return(list(locArea=percepArea,densK=densK,Kl=Kl,Kg=Kg))
    }
    else {
      K=rCalc(DATA)/(DATA$bSlope-DATA$dSlope)
      return(K)
    }
  }  
}

speedCalc=function(DATA,is.RW=F){
  if(!is.RW) speed=2*sqrt(DCalc.2D(DATA)*rCalc(DATA))
  else speed=2*sqrt(DCalc.2D/as.integer(DATA$tmax))
  return(speed)
}

# Moving window mean

#
displacement=function(obj,XYref=c(0,0))
{
  disp=sqrt((obj[,3]-XYref[1])^2+(obj[,4]-XYref[2])^2)
  return(disp)
}

### toDo: implementar calculo da frente por percentil.
frontDist=function(data,ind=10) # Specific format required for "data". It is a vector containing the displacements of each individual at a specific time
{
  vec=sort(data,decreasing=T) # Sort data in decreasing order
  if(length(data)>=ind) front=vec[ind] 
  else front=vec[length(data)]
  return(front)
}

#########
#  Plots 
#########

# plot Nt data
plotNt=function(DATA,r,sum.incl=0,land.area=10^8,R=0,D=0,
                main="N(t)",ylab="N (population size)",xlab="t",
                ylim=c(0,max(c(max(DATA$Nt[,-1])+10,r*land.area/sum.incl))),xlim=NULL,
                bg="white",foreg="black",lineColors=c("gray50",1,2),cex=1.5,LINES=T,lwPred=2,pch=NULL,
                plotReps=T,plotEst=T,plotPred=T,meanType=c("none","zeros","ifExt")){ # DATA = list returned by getNSp()
  dimData=dim(DATA$Nt)
  old.par=par(bg=bg,fg=foreg,col.lab=foreg,col.axis=foreg,col.main=foreg)
  plot(DATA$Nt[,1],DATA$Nt[,2],type="n",xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,cex=cex,main=main) # plot template
  
  if(plotReps){ # Plot replicates?
    if(DATA$hasExt){ # Is there extinction in this parameter combination?
      for(i in 2:dimData[2]){    
        if(is.na(DATA$extincTimes[i-1])){# Is there extinction in this replicate?
          if(LINES) lines(DATA$Nt[,1],DATA$Nt[,i],col=lineColors[1])
          else points(DATA$Nt[,1],DATA$Nt[,i],col=lineColors[1],cex=0.2)
        } 
        else {
          lastPos=floor(DATA$extincTimes[i-1])
          if(LINES) lines(c(DATA$Nt[1:(lastPos+1),1],DATA$extincTimes[i-1],DATA$Nt[(lastPos+2):dim(DATA$Nt)[1],1]),
                          c(DATA$Nt[1:(lastPos+1),i],0,DATA$Nt[(lastPos+2):dim(DATA$Nt)[1],i]),col=lineColors[1])
          else points(c(DATA$Nt[1:(lastPos+1),1],DATA$extincTimes[i-1],DATA$Nt[(lastPos+2):dim(DATA$Nt)[1],1]),
                      c(DATA$Nt[1:(lastPos+1),i],0,DATA$Nt[(lastPos+2):dim(DATA$Nt)[1],i]),col=lineColors[1],cex=0.2)
        }    
      }
    }
    else {
      for(i in 2:dimData[2]){
        if(LINES) lines(DATA$Nt[,1],DATA$Nt[,i],col=lineColors[1]) 
        else points(DATA$Nt[,1],DATA$Nt[,i],col=lineColors[1],cex=0.2)
      }
    }
  }
  
  if(plotPred){# Plot prediction?
    if(sum.incl==0){
      curve(DATA$Nt[1,2]*exp(r*x),add=T,col=lineColors[3],lwd=lwPred,n=1001)
    }
    else {
      K=r*land.area/sum.incl
      if(D>0 & R>0) {
        #Kl=r*pi*R*R/sum.incl
        #Kg=Kl*land.area/(pi*R*R)
        curve(K+0*x,add=T,col=lineColors[3],lwd=lwPred,n=1001)        
      }
      else {
        curve(K/(1+((K/DATA$Nt[1,2])-1)*exp(-r*x)),add=T,col=lineColors[3],lwd=lwPred,n=1001)
      }
    }
  }
  
  if(plotEst){ # Plot estimated population sizes (mean(Nt))? 
    #lgth=length(meanType)
    lt=0    
    for(i in meanType){
      lt=lt+1
      meanN=getMeanNSp(DATA,Type="n",i)
      meanN[is.nan(meanN)]=0
      if(LINES) lines(DATA$Nt[,1],meanN,col=lineColors[2],lwd=2,lty=lt)
      else points(DATA$Nt[,1],meanN,col=lineColors[2],cex=cex,pch=pch)
    }
  }
  par(old.par)  
  # Legends need to be plotted manually.
}

#  Plot speed data
plotSpeed=function(DATA,r,D,
                   main="Speed of the expansion front",ylab="Speed (distance/time)",xlab="t",
                   ylim=NULL,xlim=NULL,
                   bg="white",foreg="black",lineColors=c("gray50",1,2),cex=1.5,LINES=T,pch=NULL,lwPred=NULL,
                   plotReps=T,plotEst=T,plotPred=T,meanType=c("none","zeros","ifExt")){ # DATA = list returned by getNSp()
  dimData=dim(DATA$speed)
  old.par=par(bg=bg,fg=foreg,col.lab=foreg,col.axis=foreg,col.main=foreg)
  plot(DATA$speed[,1],DATA$speed[,2],type="n",xlab=xlab,ylab=ylab,ylim=ylim,cex=cex,main=main) # plot template
  
  if(plotReps){ # Plot replicates?
    if(DATA$hasExt){ # Is there extinction in this parameter combination?
      for(i in 2:dimData[2]){    
        if(is.na(DATA$extincTimes[i-1])){# Is there extinction in this replicate?
          if(LINES) lines(DATA$speed[,1],DATA$speed[,i],col=lineColors[1])
          else points(DATA$speed[,1],DATA$speed[,i],col=lineColors[1],cex=0.2)
        } 
        else {
          lastPos=floor(DATA$extincTimes[i-1])
          if(LINES) lines(c(DATA$speed[1:(lastPos+1),1],DATA$extincTimes[i-1],DATA$speed[(lastPos+2):dim(DATA$speed)[1],1]),
                          c(DATA$speed[1:(lastPos+1),i],0,DATA$speed[(lastPos+2):dim(DATA$speed)[1],i]),col=lineColors[1])
          else points(c(DATA$speed[1:(lastPos+1),1],DATA$extincTimes[i-1],DATA$speed[(lastPos+2):dim(DATA$speed)[1],1]),
                      c(DATA$speed[1:(lastPos+1),i],0,DATA$speed[(lastPos+2):dim(DATA$speed)[1],i]),col=lineColors[1],cex=0.2)
        }    
      }
    }
    else {
      for(i in 2:dimData[2]){
        if(LINES) lines(DATA$speed[,1],DATA$speed[,i],col=lineColors[1]) 
        else points(DATA$speed[,1],DATA$speed[,i],col=lineColors[1],cex=0.2)
      }
    }
  }
  
  if(plotPred){# Plot prediction?
    if(r==0) curve(2*sqrt(D/x),lineColors[3],add=T,lwd=lwPred,n=1001)
    else curve(2*sqrt(r*D)+0*x,col=lineColors[3],add=T,lwd=lwPred,n=1001)
    }
    
  if(plotEst){ # Plot estimated population sizes (mean(Nt))? 
    #lgth=length(meanType)
    lt=0    
    for(i in meanType){
      lt=lt+1
      meanS=getMeanNSp(DATA,Type="v",i)
      meanS[is.nan(meanS)]=0
      if(LINES) lines(DATA$speed[,1],meanS,col=lineColors[2],lwd=2,lty=lt)
      else points(DATA$speed[,1],meanS,col=lineColors[2],cex=cex,pch=pch)
    }
  }
  par(old.par)  
  # Legends need to be plotted manually.
}














######## IMPLEMENTING #################


plotDenDis=function(output2L,times,plot3D=F){
  disp=displacement(output2L)
  for(i in 1:length(times)){
    par(mfrow=c(1,1))
    index=which(output2L[,1]==times[i])
    plot(disp[index],output2L[index,5],main=paste("t = ",times[i],sep="")) #### Ajustar parametros gráficos
    
  }

}

Anima=function(DATA,landscape,framePsec=10,Colors=c(indiv="white",Habitat="green",Matrix="black")){# landscape must have "landscape" class 
  library(animation)
  LatLong=landscape$numb.cells
  opts=ani.options(outdir=getwd(),ani.width=LatLong^2,ani.height=LatLong^2)
  lands=matrix(landscape$scape,landscape$numb.cells,landscape$numb.cells)
  saveVideo({mkPlots(DATA,unique(DATA$V1),land=lands,LatLong=LatLong,Colors=Colors)},video.name = "./movement.mp4",
            interval=1/framePsec,img.name="T")
  ani.options(opts)
}

mkPlots=function(DATA,times=unique(DATA$V1),maxX=max(DATA$V3),maxY=max(DATA$V4),LatLong,
                 land,Colors){  
  for (t in 1:length(times)){
    index<-which(DATA$V1==times[t]);
    #mycolors<-allcolors[DATA$V2[index]];
    image(1:LatLong,1:LatLong,land,col=Colors[2:3],xlab="Longitude",ylab="Latitude")
    points(DATA$V3[index],DATA$V4[index],main=paste("Time = ",times[t],sep=""),col=Colors[1],pch=19);
  }
}











####################################
#  Plot XY density distribution data
####################################
plot.XYdistrib=function(FILE="output-00001.txt",Dcoef,n.tests=5,tmax=50)
{
  times=sort(sample(1:tmax,n.tests,replace=F))
  t3=read.table(FILE,sep=" ")
  
  par(mfrow=c(n.tests,2))
  for(i in times)
  {
  #x
  hist(t3[which(t3[,1]==i),3],freq=F,breaks=13)
  curve(exp(-(x^2)/(4*Dcoef*i))/sqrt(4*pi*Dcoef*i),col=2,add=T)
  #y
  hist(t3[which(t3[,1]==i),4],freq=F,breaks=13)
  curve(exp(-(x^2)/(4*Dcoef*i))/sqrt(4*pi*Dcoef*i),col=2,add=T)
  }
}

############################
#  Plot sd(x) & sd(y) 
############################
plot.sdXY=function(FILE="output-00001.txt",Dcoef,tmax=50,name="Standard Deviation")
{
  ## Teste desvio padrao
  t3=read.table(FILE,sep=" ")
  #par(mfrow=c(2,1))  
  sd.x=aggregate(t3[,3],list(t3[,1]),sd)
  sd.y=aggregate(t3[,4],list(t3[,1]),sd)
  # sd x
  plot(sd.x[,1],sd.x[,2],xlim=c(0,tmax),type="n",xlab="Time",ylab= "Std. Deviation of X positions",main=name)
  lines(sd.x[,1],sd.x[,2],lwd=3)
  curve(sqrt(2*Dcoef*x),col=2,add=T,lwd=3)
  #sd y
  plot(sd.y[,1],sd.y[,2],xlim=c(0,tmax),type="n",xlab="Time",ylab= "Std. Deviation of Y positions",main=name)
  lines(sd.y[,1],sd.y[,2],lwd=3)
  curve(sqrt(2*Dcoef*x),col=2,add=T,lwd=3)
  #par(mfrow=c(1,1))
}

# Dcoef=10000*0.8/4
# plot.XYdistrib(Dcoef=Dcoef)

################
#  Plot XY data
################

# t3=read.table("output-00001.txt",sep=" ")
# for(i in 1:20){
# plot(t3[which(t3[,1]==i),3],t3[which(t3[,1]==i),4],pch=16)}
# 
# head(t3)
# t3
# TWoPlot <- function(pop, land, col1="gray20", col2="gray70") {
#   n = land$numb.cells
#   s <- seq(-n*land$cell.size/2, n*land$cell.size/2, length=n) # creates the x- and y- sequences for image
#   if (sum(land$scape) == n*n) { 
#     color = col1
#   } else {
#     color = c(col2, col1)
#   }
#   image(s, s, matrix(land$scape,ncol=n), col=color)
#   points(pop[,],pop[,], pch=4, col=2)
# }


###########################
#  Critical Patch Size Test
###########################

#params=c(0.007,0.004,0.8,0,0,500,0,0,0,20,0,1000)
multi.run.cP=function(PATH="~/Desktop/TWoTests/Critical_Patch/Comb-01",nrep=20,b0=0.007,d0=0.004,m0=0.8,inc.b=0,inc.d=0,
                      step=500,radius=0,dens.t=0,config=0,N0=20,tm=1000,lands.range=c(-5000,5000,1000))
{
  oldpath <- getwd()  
  dir.create(PATH)
  setwd(PATH)
  Dcoef=m0*(step^2)/4
  A.crit=2*(pi^2)*Dcoef/(b0-d0)
  land.sides=ceiling(sqrt(A.crit))+seq(lands.range[1],lands.range[2],by=lands.range[3])
  
  for(j in 1:length(land.sides))
  {
    lands=Landscape(cover=1,type="b",cell.size=land.sides[j]/10,numb.cells = 10)
    dir.create(paste(j,land.sides[j],sep="_"))
    setwd(paste(PATH,land.sides[j],sep="/"))          
    
    for(i in 1:nrep)
    {
      TWoLife(taxa.basal=b0,
              taxa.morte=d0,
              move=m0,
              incl.birth=inc.b,
              incl.death=inc.d,
              passo=step,            
              raio=radius,
              density.type=dens.t,
              ini.config=config,
              #####################
              N=N0,
              AngVis=360,
              death.mat=1,
              landscape=lands,
              tempo=tm,            
              out.code=i)
    }          
    setwd(PATH)
  }
  return(land.sides)
  setwd(oldpath)
}

##multi.run.cP(nrep=50)
##multi.run.cP(PATH="TWoTests/Critical_Patch/Comb-03",lands.range=c(2000,17000,1000),nrep=50)

data.Acrit=function(PATH="~/Desktop/TWoTests/Critical_Patch/Comb-01",nreps=50,tmaxi=1000)
{
  directories=paste(PATH,dir(PATH),sep="/")
  Nt.end=matrix(0,nreps,length(directories))
  for(i in 1:length(directories))
  {
    teste=Nt.data(wdir=directories[i],tmax=tmaxi)
    Nt.end[,i]=teste[tmaxi+1,2,]  
  }
  return(list(Nt.end,teste))
}
#testes=data.Acrit("~/Desktop/TWoTests/Critical_Patch/Comb-02")

plot.Acrit=function(object,move=0.8,step=500,growth=0.003,N0=20,range=c(-17000,17000,1000),tmaxi=1000,Ncrit=0)
{
  Dcoef=move*(step^2)/4
  A.crit=2*(pi^2)*Dcoef/growth
  sides=ceiling(sqrt(A.crit)+seq(range[1],range[2],by=range[3]))
  areas=sides^2
  
  par(mfrow=c(1,2))
  plot(rep(areas,each=dim(object)[1]),as.numeric(object),pch=1,xlab="Patch Size",ylab="Population size (Nt =1000)",
       main=paste("Reaction-diffusion: Exponential Growth \n growth = ",growth,", D = ", Dcoef, ", N0 = ", N0,sep=""),
       cex.main=0.8)
  lines(areas,apply(Nt.end,2,mean),col=2)
  abline(v=areas[ceiling(length(areas)/2)],h=N0*exp(growth*tmaxi),lty=c(1,3),col=c(4,1))
  
  
  plot(rep(areas,each=dim(object)[1]),as.numeric(object),type="n",ylim=c(0,1.1),xlab="Patch Size",
       ylab="Extinction Probability for t = 1000",
       main=paste("Reaction-diffusion: Exponential Growth \n growth = ",growth,", D = ", Dcoef, ", N0 = ", N0,sep=""),
       cex.main=0.8)
  for(i in Ncrit)
  {
    indexs=which(object<=i)
    occurence=object
    occurence[indexs]=0
    occurence[-indexs]=1
    points(areas,1-apply(occurence,2,sum)/dim(occurence)[1],col=i+1)
  }
  abline(v=areas[ceiling(length(areas)/2)],h=1,lty=c(1,3),col=c(2,1))
  legend(x =areas[ceiling(length(areas)/2)+2] ,y =0.6 ,"Critical Patch Size",bty="n",lty=3,col=1,cex=0.8)
}

#plot.Acrit(teste[[1]])

#################################
###### Fragmented Landscape Tests
#################################

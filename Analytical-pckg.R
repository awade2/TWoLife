#### FUNCTION: Generate sequencies for naming folders.

numSeqGen=function(strTempl=NULL,first,last,algMax=5,sepStrNum=""){
  numSeq=sprintf(paste("%0",algMax,"d",sep=""),first:last)
  fullSeqName=paste(strTempl,sepStrNum,numSeq,sep="")
  return(fullSeqName)
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

####################
#  Get N(t) data
####################

##### NOTE: This function is specific for intervals = 1. Other intervals were not implemented. 

getNt=function(wdir,intervals=1,tmax){
  oldpath <- getwd()
  setwd(wdir)
  files=dir()
  matNt=matrix(0,tmax+1,length(files)+1) # Matrix template for Nt data of each replicate 
  matNt[,1]=0:tmax
  tExt=NA
  
  for(i in 1:length(files)){
    output2L=read.table(files[i],sep=" ")
    lastLine=length(output2L[,1])
    if(is.na(output2L[lastLine,2])){ # Is there populaition extinction?
      EXT=1
      tExt[i]=output2L[lastLine,1] # Time to extinction
      output2L=output2L[-lastLine,] # For further purposes, remove last line (it has NAs).        
    } else EXT=0
    
    Nt=as.data.frame(table(output2L[,1]))
    Nt[,1]=as.integer(levels(Nt[,1]))
    colnames(Nt)[1]="times"
    MAXt=max(Nt[,1])
    
    if(length(Nt[,1])==(MAXt+1)) # There are no time gaps in output2L and Nt
      matNt[1:(MAXt+1),i+1]=Nt$Freq
    else{ # There are time gaps in output2L and Nt
      chkMissing=data.frame(times=0:MAXt)
      Nt=merge(chkMissing,Nt,all.x=T,all.y=T)
      missing=sort(which(is.na(Nt[,2])),decreasing=T)
      for(j in missing){ # Fill the gaps
        Nt$Freq[j]=Nt$Freq[j+1]
      }
      matNt[1:(MAXt+1),i+1]=Nt$Freq
    }
  }
  ### Outputs
  setwd(oldpath)
  if(EXT) 
    return(list(Nt=matNt,extincTimes=tExt,hasExt=as.logical(EXT)))
  else 
    return(list(Nt=matNt,extincTimes=NA,hasExt=as.logical(EXT)))    
}

################
#  Data Summary 
################
meanNt=function(DATA,exclude=c("none","zeros","ifExt")){
  dimData=dim(DATA$Nt)
  isExt=which(DATA$Nt[dimData[1],]==0) # if there is no extinction, returns an empty integer vector
  if(exclude=="none"){
    mNt=rowMeans(DATA$Nt[,2:dimData[2]])
  } 
  else{
    if(exclude=="zeros"){
      DATA$Nt[,2:dimData[2]][DATA$Nt[,2:dimData[2]]==0]=NA
      mNt=rowMeans(DATA$Nt[,2:dimData[2]],na.rm=T) 
    } 
    else { # ifExt
      #mNt=rep(NA,dimData[1]-1)
      mNt=rowMeans(DATA$Nt[,-c(1,isExt)])
      #mNt[-(isExt-1)]=mNtShort      
    }
  }       
  return(mNt)
}

####################
#  Plot N(t) data
####################
plotNt=function(DATA,r,sum.incl=0,land.area=10^8,R=0,D=0,
                main="N(t)",ylab="N (population size)",xlab="t",
                ylim=c(0,max(c(max(DATA$Nt[,-1])+10,r*land.area/sum.incl))),xlim=NULL,
                bg="white",foreg="black",lineColors=c("gray50",1,2),cex=1.5,LINES=T,pch=NULL,
                plotReps=T,plotEst=T,plotPred=T,meanType=c("none","zeros","ifExt")){
  
  dimData=dim(DATA$Nt)
  old.par=par(bg=bg,fg=foreg,col.lab=foreg,col.axis=foreg,col.main=foreg)
  plot(DATA$Nt[,1],DATA$Nt[,2],type="n",xlab=xlab,ylab=ylab,ylim=ylim,cex=cex,main=main) # plot template
  
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
      curve(DATA$Nt[1,2]*exp(r*x),add=T,col=lineColors[3],lwd=2,n=1001)
    }
    else {
      K=r*land.area/sum.incl
      if(D>0 & R>0) {
        #Kl=r*pi*R*R/sum.incl
        #Kg=Kl*land.area/(pi*R*R)
        curve(K+0*x,add=T,col=lineColors[3],lwd=2,n=1001)        
      }
      else {
        curve(K/(1+((K/DATA$Nt[1,2])-1)*exp(-r*x)),add=T,col=lineColors[3],lwd=2,n=1001)
      }
    }
  }
  
  if(plotEst){ # Plot estimated population sizes (mean(Nt))? 
    lgth=length(meanType)
    lt=0    
    for(i in meanType){
      lt=lt+1
      meanN=meanNt(DATA,i)
      meanN[is.nan(meanN)]=0
      if(LINES) lines(DATA$Nt[,1],meanN,col=lineColors[2],lwd=2,lty=lt)
      else points(DATA$Nt[,1],meanN,col=lineColors[2],cex=0.2)
    }
  }
par(old.par)  
# Legends need to be plotted manually.
}
#legend(dim(dataSim)[1]/4,max(dataSim[,2,]),legend=c("Observed","Predicted","Replicates"),
#lty=1,lwd=2,col=c(1,2,"gray50"),bty="n")


########################
#  Ancillary functions 
########################
displacement=function(obj)
{
  disp=sqrt(obj[,3]^2+obj[,4]^2)
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

#########################
#  Get velocity(t) data
#########################
getSpeed=function(wdir=getwd(),tmax=50,nrep=20,intervals=1,ind=10) # Analisar possibilidade de trocar wdir por ouput2L
{
  oldpath <- getwd()
  setwd(wdir)
  files=dir()
  matSpeed=matrix(0,tmax+1,(nrep+1),dimnames=list(paste(1:(tmax+1)),c("Time",numSeqGen("R",1,nrep,nchar(nrep)))))
  matSpeed[,1]=0:tmax
  tExt=NA
  #
  for(i in 1:nrep){
    output2L=read.table(files[i],sep=" ") ##### Analisar possibilidade de juntar com getNt. Assim cada output só teria que ser lido uma vez.
    if(is.na(output2L[dim(output2L)[1],2])){ # Is there extinction?
      tExt[i]=output2L[dim(output2L)[1],1]
      output2L=output2L[-(dim(output2L)[1]),] # Remove last row, which indicates population extinction. A similar procedure is implemented on function getNt  
    }
   #     
    displacements=displacement(output2L) # Individual displacement
    new.output2L=data.frame(output2L[,1],displacements) # Two columns data.frame (1: time, 2: displacement)
    front=aggregate(new.output2L[,2],by=list(new.output2L[,1]),FUN=frontDist,ind=ind) # Calculates the front distance based on the displacement vector
    colnames(front)=c("times","speed")
    MAXt=max(front$times)
    nROW=dim(front)[1]
    
    if(nROW==(MAXt+1)){
      if(!is.na(tExt[i])&floor(tExt[i])>MAXt) {
        front=merge(data.frame(times=(MAXt+1):floor(tExt[i]),speed=front$speed[MAXt+1]),front,all.x=T,all.y=T)
        MAXt=floor(tExt[i])
      }
      matSpeed[2:(MAXt+1),i+1]=diff(front$speed)
    }
    else{
      chkMissing=data.frame(times=0:MAXt)
      front=merge(chkMissing,front,all.x=T,all.y=T)
      missing=sort(which(is.na(front$speed)),decreasing=T)
      for(j in missing){
        front$speed[j]=front$speed[j+1]
      }
      if(!is.na(tExt[i])&floor(tExt[i])>MAXt){
        front=merge(data.frame(times=(MAXt+1):floor(tExt[i]),speed=front$speed[MAXt+1]),front,all.x=T,all.y=T)
        MAXt=floor(tExt[i])
      } 
      matSpeed[2:(MAXt+1),i+1]=diff(front$speed)
    }
  }
  setwd(oldpath)
  return(matSpeed)
}

####################
#  Plot velocity(t)
####################
plot.vel=function(object,nrep=20,Dcoef=D,b0=0,d0=0,name="V(t)",yrange=NULL)
{
  growth=b0-d0
  plot(object[,1],object[,2],type="n",main=name,xlab="t",ylab="velocity (distance/time)",ylim=yrange,cex.main=0.7)
  for(i in 1:nrep)
  {
    lines(object[,1],object[,i+1],col="gray50") # replicates
  }
  if(nrep>1){lines(object[,1],apply(object[,-1],1,mean),lwd=4)} # mean replicates (observed)
  if(b0==0 && d0==0){curve(2*sqrt(Dcoef/x),col=2,add=T,lwd=3)} else {abline(h=2*sqrt(growth*Dcoef),col=2,lwd=3)} # expected
}

# teste=velocity("~/Desktop/Comb-0014",nrep=100,tmax=50)
# Dcoef=10000*0.8/4
# plot.vel(teste,Dcoef=Dcoef,nrep=20,growth=0.03)

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

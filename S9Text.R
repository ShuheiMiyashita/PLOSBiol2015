#S9 Text. An R script used to generate Fig. 6A and 6B.
#R script for Fig. 6A and 6B.
#Seed value is set in order to link the simulation results of this simulation and a simulation for Fig. 6C. 

set.seed(3)

for (cell in 1:10){

# parameter settings
E <- 5*10^3
R <- 3*10^4
p <- 3*10^-10
d <- 1*10^-2
rep <- c(0.5,1)

# initial settings
table <- matrix(rep(0,6),nrow=2)
table[,1] <- c(1,2)      # ID
table[,2] <- c(E/2,E/2) # number of genomic RNA 
table[,3] <- c(0,0) # number of RC
RCO <- R
nsum <- E
t <- 0

# graphics settings
nresult <- matrix(rep("NA",2000),nrow=2)
rresult <- matrix(rep("NA",2000),nrow=2)
colorn <- c("green","magenta")

# main body of simulation
while (RCO > 0) {
if (nsum == 0) break

## showing graphics by every 100 unit time
if (t%%100 == 0){
nsum <- sum(table[,2])

idn <- c("variant","wild type")
idr <- c("","")

if (table[1,2]==0 && table[1,3]==0){
idn[1] <- ""
}else{
idn[1] <- "variant"
}
if (table[2,2]==0 && table[2,3]==0){
idn[2] <- ""
}else{
idn[2] <- "wild type"
}

if (table[1,3]==0){
}else{
idr[1] <- "variant"
}
if (table[2,3]==0){
}else{
idr[2] <- "wild type"
}


resultn <- table[,2]
matrixn1 <- rbind(resultn,idn,colorn)
matrixn2 <- matrix(c(as.character(matrixn1)),nrow=3)

if (RCO > 0){
idr <- c("open sites",idr)
resultr <- c(RCO,table[,3])
colorr <- c("#FFFFFF",colorn)
}else{
idr <- idn
resultr <- table[,3]
colorr <- colorn
}
matrixr1 <- rbind(resultr,idr,colorr)
matrixr2 <- matrix(c(as.character(matrixr1)),nrow=3)

if (RCO < R){
resultrr <- table[,3]
matrixrr1 <- rbind(resultrr,idn,colorn)
matrixrr2 <- matrix(c(as.character(matrixrr1)),nrow=3)
}else{
matrixrr2 <- matrix(c(1,"no RC","#FFFFFF"),nrow=3)
}

t1 <- t/100+1
nresult[1,t1] <- t
nresult[2,t1] <- nsum
rresult[1,t1] <- t
rresult[2,t1] <- R-RCO

x1 <- na.omit(nresult[1,])
y1 <- na.omit(nresult[2,])
x2 <- rresult[1,]
y2 <- rresult[2,]

par(mfrow=c(2,3))
pie(as.numeric(matrixn2[1,]),labels=matrixn2[2,],col=matrixn2[3,],clockwise=TRUE,main=paste("variant RNA: ",table[1,2],"\n wild-type RNA: ",table[2,2],sep=""),radius=sqrt(nsum/R*d))
pie(as.numeric(matrixr2[1,]),labels=matrixr2[2,],col=matrixr2[3,],clockwise=TRUE,main=paste("variant RC: ",table[1,3],"\n wild-type RC: ",table[2,3],"\n open site: ",RCO,sep=""))
pie(as.numeric(matrixrr2[1,]),labels=matrixrr2[2,],col=matrixrr2[3,],clockwise=TRUE,main=paste("RC proportion",sep=""))
plot(x1,y1,xlim=c(0,30000),ylim=c(0,R/d*6/5),col="red",xlab="t",ylab="genomic RNA")
plot(x2,y2,xlim=c(0,30000),ylim=c(0,R*6/5),col="blue",xlab="t",ylab="RC")
plot(0,0,xlim=c(0,1),ylim=c(0,1),col="white",xlab="",ylab="",axes=F,bty="n")
text(0.5,0.5,labels=paste("t = ",t),cex=2)

}else{
}

## genomic RNA degradation and synthesis
D <- NULL
G <- NULL
for (j in 1:2) {
D[j] <- rbinom(1,table[j,2],d)
G[j] <- rbinom(1,table[j,3],rep[j])
table[j,2] <- table[j,2]-D[j]+G[j]
}

nsum <- sum(table[,2])

## RC formation
if (nsum >0) {
irc <- rbinom(1,RCO,min(c(1,nsum*p)))
RCO <- RCO-irc
ircy <- rbinom(1,irc,table[1,2]/nsum)
ircc <- irc-ircy
table[1,3] <- table[1,3]+ircy
table[2,3] <- table[2,3]+ircc
} else {
}


t <- t+1

}

# showing graphics at the end of the simulation

if (table[1,2]==0 && table[1,3]==0){
idn[1] <- ""
}else{
idn[1] <- "variant"
}
if (table[2,2]==0 && table[2,3]==0){
idn[2] <- ""
}else{
idn[2] <- "wild type"
}

resultn <- table[,2]
matrixn1 <- rbind(resultn,idn,colorn)
matrixn2 <- matrix(c(as.character(matrixn1)),nrow=3)

idr <- idn
resultr <- table[,3]
colorr <- colorn

matrixr1 <- rbind(resultr,idr,colorr)
matrixr2 <- matrix(c(as.character(matrixr1)),nrow=3)

resultrr <- table[,3]
matrixrr1 <- rbind(resultrr,idn,colorn)
matrixrr2 <- matrix(c(as.character(matrixrr1)),nrow=3)

par(mfrow=c(2,3))
pie(as.numeric(matrixn2[1,]),labels=matrixn2[2,],col=matrixn2[3,],clockwise=TRUE,main=paste("variant RNA: ",table[1,2],"\n wild-type RNA: ",table[2,2],sep=""),radius=sqrt(nsum/R*d))
pie(as.numeric(matrixr2[1,]),labels=matrixr2[2,],col=matrixr2[3,],clockwise=TRUE,main=paste("variant RC: ",table[1,3],"\n wild-type RC: ",table[2,3],"\n open site: ",RCO,sep=""))
pie(as.numeric(matrixrr2[1,]),labels=matrixrr2[2,],col=matrixrr2[3,],clockwise=TRUE,main=paste("RC proportion",sep=""))
plot(x1,y1,xlim=c(0,30000),ylim=c(0,R/d*6/5),col="red",xlab="t",ylab="genomic RNA")
plot(x2,y2,xlim=c(0,30000),ylim=c(0,R*6/5),col="blue",xlab="t",ylab="RC")
plot(0,0,xlim=c(0,1),ylim=c(0,1),col="white",xlab="",ylab="",axes=F,bty="n")
text(0.5,0.5,labels=paste("t = ",t),cex=2)

}


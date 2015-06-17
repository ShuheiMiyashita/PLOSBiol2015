#S12 Text. Another R script used to obtain the data for Fig. 6D.
#R script for obtaining the simulation results summarized in Fig. 6D.
#A parameter value for "repr", which determine the relative efficiency of vRNA synthesis by a variant, need to be provided. This script works only when repr >= 1, i.e., only for simulating co-inoculation of wild-type and an equally or more-fit variant.
#The simulation results obtained by the authors are shown in S4 Data.

repr <- 1.5

cell <- 1000
Y <- c(rep(0,cell))
C <- c(rep(0,cell))

# parameter settings

E <- 5*10^3
d <- 1*10^-2
R <- 3*10^4
p <- 3*10^-10
Eyr <- 0.5	#proportion of wild-type (Y) in the introduced vRNAs.	

for (c in 1:cell){

# initial settings

t <- 1
ta <- 1000
Ey <- round(E*Eyr)
Ec <- E-Ey
table <- matrix(rep(0,6),nrow=2)
table[,1] <- c(1,2)
table[,2] <- c(Ey,Ec)
table[,3] <- c(0,0)
RCO <- R

# main body of simulation

nsum <- E

while (RCO > 0 || ta >0) {
if (nsum == 0) break

D <- rbinom(c(1,1),table[,2],c(d,d))
G <- c(table[1,3]+rbinom(1,table[1,3],repr-1),table[2,3])
table[,2] <- table[,2]-D+G

nsum <- sum(table[,2])

if (nsum > 0) {
irc <- rbinom(1,RCO,min(c(1,nsum*p)))
RCO <- RCO-irc
ircy <- rbinom(1,irc,table[1,2]/nsum)
ircc <- irc-ircy
table[1,3] <- table[1,3]+ircy
table[2,3] <- table[2,3]+ircc
} else {
}

if (RCO == 0) {
ta <- ta-1
}else{
}

t <- t+1

}

#showing the progress of simulation
par(mfrow=c(1,1))
plot(c,c)

Y[c] <- table[1,2]
C[c] <- table[2,2]

gc()
gc()


}

#Boot-strap analysis
BSresult <- matrix(rep(0,cell*3),nrow=3)
for(b in 1:cell){
YT <- 0
CT <- 0
s <- sample(1:cell,cell,replace=T)
for (i in 1:cell){
YT <- YT+Y[s[i]]
CT <- CT+C[s[i]]
}
YC <- YT/CT
BSresult[1,b] <- YT
BSresult[2,b] <- CT
BSresult[3,b] <- YC
}

mean(BSresult[1,])	#number of the variant RNA in 1000 cells
sd(BSresult[1,])	#standard deviation caluculated by boot-strap analysis

mean(BSresult[2,])	#number of wild-type RNA in 1000 cells
sd(BSresult[2,])	#standard deviation caluculated by boot-strap analysis

mean(BSresult[3,])	#ratio of the variant to wild-type RNA
sd(BSresult[3,])	#standard deviation caluculated by boot-strap analysis


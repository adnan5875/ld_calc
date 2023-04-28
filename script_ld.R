# calculating LD for pathogen population ####3
# 27/04/2023
#################
#Data:
# a vcf file containing combined genotyping population of three species of Aschochyta Lentis
# separate population of three populations

# step 1:
# read the vcf file directly in tassel
# calculate diversity, LD using a sliding window of 50
# Extract table of the LD calculations to create figures in R.

# step 2:
# Draw Ld (r2)
# Draw Ld decay
#####################3

# -------------------------------
# import TASSEL LD output file
ld <- read.delim(file.choose(),stringsAsFactors = FALSE,header=TRUE, sep = "\t") 

##remove sites that have NaN for distance or r2
ld_sub <- ld[ld$R.2 != "NaN",]
ld_sub$dist <- as.numeric(ld_sub$Dist_bp)
ld_sub2 <- ld_sub[ld_sub$Dist_bp != "NaN",]
ld_sub2$rsq <- ld_sub2$R.2

file <- ld_sub2[,c(1,2,7,8,13:19)]
file <- file[!duplicated(file$R.2), ]
options(bitmapType = "cairo")
library(ggplot2)

Cstart <- c(C=0.1)

# fit a non linear model using the arbitrary C value, 
# N is the number of the genotypes that have the SNP site
modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), 
              data=file, start=Cstart, control=nls.control(maxiter=100))

# extract rho, the recombination parameter in 4Nr
rho <- summary(modelC)$parameters[1]

# feed in the new value of rho to obtain LD values adjusted for their distances along the chromosome/genome
newrsq <- (( (10+(rho*file$dist) )/ ( (2+(rho*file$dist)) * (11+(rho*file$dist)) ) ) *
  ( 1 + (( (3+(rho * file$dist)) * (12+(12*(rho*file$dist)) + ((rho*file$dist)^2) )) / 
      (2*file$N*(2+(rho*file$dist)) * (11+(rho*file$dist) ) ))))

newfile <- data.frame(file$dist, newrsq)

maxld <- max(newfile$newrsq,na.rm=TRUE) #using max LD value from adjusted data
halfdecay = maxld*0.5
halfdecaydist <- newfile$file.dist[which.min(abs(newfile$newrsq-halfdecay))]
newfile <- newfile[order(newfile$file.dist),]

# plotting the values
pdf("LD_decay_allinone.pdf", height=5, width = 9)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 
plot(file_com$dist, file_com$rsq, pch=".", cex=2, xlab="Distance (bp)", ylab=expression(LD ~ (r^2)), col="grey")
lines(newfile_com$file.dist, newfile_com$newrsq, col=2, lwd=2)
lines(newfile_pk$file.dist, newfile_pk$newrsq, col=3, lwd=2)
lines(newfile_pp$file.dist, newfile_pp$newrsq, col=4, lwd=2)
lines(newfile_ppll$file.dist, newfile_ppll$newrsq, col=6, lwd=2)
legend(x = "topright", 
       legend = c("combined", "Pk", "Pp", "Ppll"),
       lty = c(1),
       lwd = c(2),
       col = c(2,3,4,6))

abline(h=0.1, col="blue") # if you need to add horizental line
#abline(v=halfdecaydist, col="green")
#mtext(round(halfdecaydist,2), side=1, line=0.05, at=halfdecaydist, cex=0.75, col="green")
dev.off()
#########################

file_com <- file
newfile_com <- newfile
newfile_com$specie <-"Combined"
# plotting the values
pdf("LD_decay_combined.pdf", height=5, width = 9)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 
plot(file_com$dist, file_com$rsq, pch=".", cex=2, xlab="Distance (bp)", ylab=expression(LD ~ (r^2)), col="grey")
lines(newfile_com$file.dist, newfile_com$newrsq, col=2, lwd=2)
#lines(newfile_pk$file.dist, newfile_pk$newrsq, col="green", lwd=2)
#lines(newfile_pp$file.dist, newfile_pp$newrsq, col="orange", lwd=2)
#lines(newfile_ppll$file.dist, newfile_ppll$newrsq, col="purple", lwd=2)
abline(h=0.1, col="blue") # if you need to add horizental line
#abline(v=halfdecaydist, col="green")
#mtext(round(halfdecaydist,2), side=1, line=0.05, at=halfdecaydist, cex=0.75, col="green")
dev.off()

file_pp <- file
newfile_pp <- newfile
newfile_pp$specie <-"Pp"
pdf("LD_decay_Pk.pdf", height=5, width = 9)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 
plot(file_pk$dist, file_pk$rsq, pch=".", cex=2, xlab="Distance (bp)", ylab=expression(LD ~ (r^2)), col="grey")
#lines(newfile_com$file.dist, newfile_com$newrsq, col=2, lwd=2)
lines(newfile_pk$file.dist, newfile_pk$newrsq, col=3, lwd=2)
#lines(newfile_pp$file.dist, newfile_pp$newrsq, col="orange", lwd=2)
#lines(newfile_ppll$file.dist, newfile_ppll$newrsq, col="purple", lwd=2)
abline(h=0.1, col="blue") # if you need to add horizental line
dev.off()

file_ppll <- file
newfile_ppll <- newfile
newfile_ppll$specie <-"Ppll"

pdf("LD_decay_Ppll.pdf", height=5, width = 9)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 
plot(file_ppll$dist, file_ppll$rsq, pch=".", cex=2, xlab="Distance (bp)", ylab=expression(LD ~ (r^2)), col="grey")
#lines(newfile_com$file.dist, newfile_com$newrsq, col=2, lwd=2)
#lines(newfile_pk$file.dist, newfile_pk$newrsq, col=3, lwd=2)
#lines(newfile_pp$file.dist, newfile_pp$newrsq, col="orange", lwd=2)
lines(newfile_ppll$file.dist, newfile_ppll$newrsq, col=5, lwd=2)
abline(h=0.1, col="blue") # if you need to add horizental line
dev.off()


file_pk <- file
newfile_pk <- newfile
newfile_pk$specie <-"Pk"

pdf("LD_decay_Pp.pdf", height=5, width = 9)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 
plot(file_pp$dist, file_pp$rsq, pch=".", cex=2, xlab="Distance (bp)", ylab=expression(LD ~ (r^2)), col="grey")
#lines(newfile_com$file.dist, newfile_com$newrsq, col=2, lwd=2)
#lines(newfile_pk$file.dist, newfile_pk$newrsq, col=3, lwd=2)
lines(newfile_pp$file.dist, newfile_pp$newrsq, col="orange", lwd=2)
#lines(newfile_ppll$file.dist, newfile_ppll$newrsq, col=5, lwd=2)

abline(h=0.1, col="blue") # if you need to add horizental line
dev.off()


#other method for plotting

newfile_all <- rbind(newfile_com, newfile_pp, newfile_ppll, newfile_pk)
ggplot(data = file_com, aes(x=dist, y=rsq)) +
  geom_point(alpha=0.2 ) +
  geom_line(data=newfile_all, aes(x=file.dist, y=newrsq, color =  specie),  size=0.5)

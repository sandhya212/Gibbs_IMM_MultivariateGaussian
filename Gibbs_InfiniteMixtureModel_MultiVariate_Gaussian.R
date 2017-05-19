## Infinite Mixture Model for multivariate Gaussian mixture, X(nxd)
## Likelihood; multivariate Gaussian distribution(mean,covariance)
## mean(1xd): multivariate Gaussian
## covariance(dxd): inverse Wishart
## 
## May 12th 2017
## Code author SP


library(lattice)
library(coda)
library(MASS)
library(MCMCpack)
library(Matrix)
library(bayesm)
library(robustbase)
library(mvtnorm)
library(chron)
library(mvtnorm)
library(mnormt)
library(RColorBrewer)
library(Rtsne)



###########################################
###########################################
####Change things within this block #######
###########################################
###########################################

##working path name
## change this to the folder location where code resides

Rprof(filename = paste(working_path,"Rprof.out",sep=""), append = FALSE, interval = 0.02, memory.profiling=FALSE)

#choose the dimension of the model (# of genes)
d <- 6

#fix number of Gibbs iterations
steps <- 20

## Simulate data from multivariate Gaussian mixture
## currently there are 6 clusters each with 500 observations

sim_matrix <- matrix(0,d,d)
for(i in 1:d){
    sim_matrix[i,i]<-2
}

n <- 500 # (number of rows)

x1 <- rmnorm(n, rep(0, d),sim_matrix)
x2  <- rmnorm(n, rep(7,d),sim_matrix)
x3  <- rmnorm(n, rep(-7,d),sim_matrix)
x4  <- rmnorm(n, rep(20,d),sim_matrix)
x5 <- rmnorm(n, rep(40, d),sim_matrix)
x6 <- rmnorm(n, rep(-20, d),sim_matrix)

x <- rbind(x1,x2,x3,x4,x5,x6)

N <- dim(x)[1] #total number of observations in all classes

z_true <- c(rep(1,n),rep(2,n),rep(3,n),rep(4,n),rep(5,n),rep(6,n))



###########################################
###########################################
#### Do NOT change anything below #########
###########################################
###########################################

if(! dir.exists(paste0(getwd(),"/output"))){
    dir.create(paste0(getwd(),"/output/"))
    dir.create(paste0(getwd(),"/output/plots/"))
    dir.create(paste0(getwd(),"/output/plots/Inferred_labels/"))
    dir.create(paste0(getwd(),"/output/plots/Inferred_labels_per_step/"))
}

#colorlist for plots
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) #len = 74
num_col <- 20
pie(rep(1,num_col), col=(col_vector[1:num_col]))
col_palette <- col_vector[1:num_col]; # or sample if you wish



X_tsne <- Rtsne(x,check_duplicates = FALSE,perplexity=10);

startTime <- Sys.time()


# Set initial values for the hyperparameters:
# nu, Delta, mu_0, kappa, alpha

nu <- d+1
Delta <- matrix(0,d,d)
for(i in 1:d){Delta[i,i]<-3}
mu0 <- rep(1.5,d)
kappa <- 1
alpha <- 1

# Choose initial values for:
# K, mu, Sigma, pi, C , N_k

K <- 2  
Kit <- K #counts how many new classes are created without taking into account take class can be left after some time
mu <- rep( rep(0,d),K+1)
dim(mu) <- c(d,K+1)
dimnames(mu)<-list( 1:d,c("c1","c2","dummy")) #dummy is here to insure that mu and sigma have at least two classes otherwise dim(mu) becomes null 
cnames <- dimnames(mu)[[2]]
cnames <- cnames[!cnames == "dummy"]
Sigma <- rep( matrix(1,d,d),K+1)
dim(Sigma)<- c(d,d,K+1)
dimnames(Sigma)<-list( 1:d,1:d,c("c1","c2","dummy"))
diag22<- c(rep( c(1, rep(0,d)),d-1),1)
dim(diag22) <- c(d,d)
for (k in 1:K){
  Sigma[,,k] <- diag22
}
Pi<- rep( 1/K, K)
Class <- rep(c("c1","c2"),floor(N/K)+1)
length(Class) <- N #cuts Class to the right size
Class <- as.factor(Class)
Nk<-table(Class)

#production of the data frame
mydata <- data.frame(observ=x,class=Class)

mydata <- data.frame(cbind(x,Class))
cloud(V1 ~ V2 * V3, data = mydata,
groups = Class,
perspective = TRUE,
border = TRUE)

# Update the parameters of the NIW distribution after observing data
####
updateNIW <- function(fctx,fctNc,fctnu,fctkappa,fctmu0,fctDelta)
  {
   fctnup <- fctnu + fctNc
   fctkappap <- fctkappa + fctNc
   sumxk <- rep(0,d)
   sumxkmatrix <- matrix(0,d,d)
   if (is.null(dim(fctx))){
     sumxk <- fctx
     sumxkmatrix <- fctx %o% fctx
   }else{  
     for(j in 1:length(fctx[,1])){ 
       sumxk <- sumxk + fctx[j,]
       sumxkmatrix <- sumxkmatrix + fctx[j,] %o% fctx[j,]
     }
   }
   fctmu0p <- (fctkappa*fctmu0 + sumxk) / fctkappap
   fctDeltap <-(fctDelta + sumxkmatrix + fctkappa * fctmu0 %o% fctmu0 - fctkappap * fctmu0p %o% fctmu0p )
   c(fctnup,fctkappap,fctmu0p,fctDeltap) 
  }  
Nc <- Nk[1]
test <- updateNIW(x,Nc,nu,kappa,mu0,Delta)

inferred_parm2 <- rep(0,N*(1+1+d+d*d))
dim(inferred_parm2) <- c(N,1+1+d+d*d)
for(i in 1:N){
inferred_parm2[i,] <- updateNIW(x[i,],1,nu,kappa,mu0,Delta)
}


timedmnorm <- 0
timedmnorm2 <- 0
mydata$class <- Class

###########################################
##### Gibbs sampling ######################
###########################################
step <- 0

while(step <steps){ #begin of sampling loop
step <- step +1

print(paste("step is ", step));

#evaluate the dmnorm grouped by classes
timedmnormstart2 <- Sys.time()
Q <- rep(0,N*K)
dim(Q) <- c(N,K)
dimnames(Q)<- list(1:N,cnames)
Q <- data.frame(Q,Class)
for(j in 1:K){
  for(k in 1:K){
      print(paste("cluster is", k));
      
      diag(Sigma[,,cnames[k]]) <- diag(Sigma[,,cnames[k]]) + 0.001;
      temp_Sigma <- matrix(forceSymmetric(Sigma[,,cnames[k]]),d,d);
      
    Q[Q$Class==cnames[j],cnames[k]]<-dmnorm(as.matrix(mydata[mydata$class==cnames[j],1:d]),mu[,cnames[k]],temp_Sigma)
  }
}
timedmnormend2 <- Sys.time()
timedmnorm2 <- timedmnorm2 + difftime(timedmnormstart2,timedmnormend2,tz="",units="mins")
cnamespreviousstep <- cnames

#go through the observations and reassign classes
for (i in 1:N){
    print(paste("cell is ", i));
  qi <- rep(0,K)
  q0 <- 0 
  names(qi)<- cnames
  inferred_parm <- inferred_parm2[i,]
  inferred_parm5 <- inferred_parm[-(1:(d+2))]
  dim(inferred_parm5) <- c(d,d)  
  tdeg <- inferred_parm[[1]]-d +1
  tmu <- inferred_parm[3:(3+d-1)]
  tsigma <- solve(inferred_parm5) * (inferred_parm[[2]]+1) / (inferred_parm[[2]]*tdeg)
  diag(tsigma) <- diag(tsigma)+ 0.001;
  
  #make symmetric psd for d >=10
  tsigma[lower.tri(tsigma)] <- 0;
  diag_tsigma <- diag(tsigma);
  te <- tsigma + t(tsigma);
  diag(te) <- diag(te) - diag_tsigma;

  
  q0<- alpha * mnormt::dmt(x[i,],tmu,te,tdeg)
  
  
  timedmnormstart <- Sys.time()
  for (k in 1:K){
      print(paste('cluster is ', k));
    if( cnames[k] %in% cnamespreviousstep){
      qi[cnames[k]] <- Nk[cnames[k]]*Q[i,cnames[k]]
    }else{
        diag(Sigma[,,cnames[k]]) <- diag(Sigma[,,cnames[k]]) + 0.001;
        temp_Sigma <- matrix(forceSymmetric(Sigma[,,cnames[k]]),d,d);
        qi[cnames[k]] <- Nk[cnames[k]]* (dmnorm(x[i,],mu[,cnames[k]],temp_Sigma)[1])
    }
  }
  timedmnormend <- Sys.time()
  timedmnorm <- timedmnorm + difftime(timedmnormstart,timedmnormend,tz="",units="mins")
  c <- sum(qi) + q0
  q0 <- q0/c
  qi <- qi /c
  ClassOld <- as.character(Class)
  NkOld <- Nk
  
  Classtemp <- sample(c("new",cnames),1, replace = TRUE,c(q0,qi) ) #resample the class indicator parameters, also sum(qi)+q0 = 1
 
  
  if( Classtemp=="new") #0 means i starts a new class
  {
    if (NkOld[ClassOld[i]]==1)#if the observation was alone in one class and changes to a new class then this is not really a new class
    {
      Kit <- Kit +1
      Class <- factor(Class, levels= c( levels(Class), paste("c",as.character(Kit),sep="")))
      Class[i] <- paste("c",as.character(Kit),sep="")
      Sigma[,,ClassOld[i]] <- riwish(inferred_parm[[1]],inferred_parm5)
      mu[,ClassOld[i]] <- rmnorm(1,inferred_parm[3:(3+d-1)], Sigma[,,ClassOld[i]]/inferred_parm[[2]])
      for(s in 1:length(dimnames(mu)[[2]]))
      {
        if(dimnames(mu)[[2]][s]==ClassOld[i]) {dimnames(mu)[[2]][s] <- paste("c",as.character(Kit),sep="")}
      }
      dimnames(Sigma)[[3]] <- dimnames(mu)[[2]]
      Nk <- table(Class)
      cnames <- dimnames(mu)[[2]]
      cnames <- cnames[!cnames == "dummy"] 
      
    } else
    {
      K<-K+1
      Kit <- Kit +1
      Class <- factor(Class, levels= c( levels(Class), paste("c",as.character(Kit),sep="")))
      Class[i] <- paste("c",as.character(Kit),sep="")
      Sigmatemp <- riwish(inferred_parm[[1]],inferred_parm5)
      
      mutemp <-  rmnorm(1,inferred_parm[3:(3+d-1)], Sigma[,,ClassOld[i]]/inferred_parm[[2]])
      dimnamesmu <- dimnames(mu)
      mu <- cbind(mu,mutemp);
      dimnames(mu) <- list( dimnamesmu[[1]],c( dimnamesmu[[2]],as.character(Class[i])))
      dimnamessig <- dimnames(Sigma)
      
      Sigma <- c(Sigma, Sigmatemp)
      dim(Sigma)<- c(d,d,K+1)
      dimnames(Sigma) <- list( dimnamessig[[1]], dimnamessig[[2]], c(dimnamessig[[3]],as.character(Class[i])))
      
      Nk <- table(Class)
      cnames <- dimnames(mu)[[2]]
      cnames <- cnames[!cnames == "dummy"] 
    }   
  } else
  {
    Class[i] <- Classtemp
    Nk <- table(Class)
    if ((NkOld[ClassOld[i]]==1) && (as.character(Class[i])!=ClassOld[i]) )
    {
       mu <- mu[ , !(colnames(mu) %in% ClassOld[i])]
       Sigma <- Sigma[,,!(dimnames(Sigma)[[3]] %in% ClassOld[i])]
       cnames <- dimnames(mu)[[2]]
       cnames <- cnames[!cnames == "dummy"] 
       K<- K-1 #correct the number of classes and remove its level from the factor(Class)
       Class <- factor(Class)
       Nk <- table(Class)
    }
    
  }   
}#end of the loop on the observations

ClassChar <- as.character(Class)
mydata$class <- Class
Nk <- table(Class)



#for each class resample the parameters
 for(k in 1:K){
   inferred_parm <-  updateNIW(as.matrix(mydata[mydata$class==cnames[k],1:d]),Nk[cnames[k]],nu,kappa,mu0,Delta)
   inferred_parm5 <- inferred_parm[-(1:(d+2))]
   dim(inferred_parm5) <- c(d,d) 
   Sigma[,,cnames[k]] <-  riwish(inferred_parm[[1]],inferred_parm5)
   mu[,cnames[k]] <- rmnorm(1,inferred_parm[3:(3+d-1)], Sigma[,,cnames[k]]/inferred_parm[[2]])
 }



f <- paste0(getwd(),"/output/plots/Inferred_labels_per_step/inferred_labels_(step)_",step,".pdf")
pdf(file=f)
plot(X_tsne$Y[,1],X_tsne$Y[,2],col = col_palette[1*(as.numeric(factor(Class)))],  main=paste("t-SNE of X (inferred labels) step", step));
dev.off();


} #end of sampling loop


endTime <- Sys.time()
diffTime <-  difftime(endTime, startTime,tz="",units="mins")
print(diffTime)
print(table(Class))


############

f <- paste0(getwd(),"/output/plots/Inferred_labels/final_inferred_labels_1.pdf")
pdf(file=f)
par(mfrow=c(2,1));
plot(X_tsne$Y[,1],X_tsne$Y[,2],col = col_palette[1*(z_true)],  main="t-SNE of X (true labels)");
plot(X_tsne$Y[,1],X_tsne$Y[,2],col = col_palette[1*(as.numeric(factor(Class)))],  main="t-SNE of X (inferred labels)");
dev.off();





















# things that have changed since the "true" original version
# 1. Correct mupdate
# 2. Set seed below
# 3. Take the l2 norm instead of the spectral norm
# 4. Save the ELBO
# 5. turn the alpha matrices into symmetric inclusion matrices
# 6. save the alpha matrices, inclusion probability matrices, and ELBO
# 7. Commented out E = rnorm(n, 1, 0)
# 8. Set.seed(resp_index)

set.seed(1)
#rm(list=ls())
rm(list=ls())
source("cov_vsvb.R")
source("ELBO_calculator.R")
library(Matrix)
library(magic)
library(psych)
library(ggplot2)
library(reshape2)
library(MASS)
library(pracma)
library(varbvs)
library(ks)
###################################################################

logit <- function(x) {if ((x==0)|(x==1)) return(0) else return(log(x/(1-x))) }
Afunc <- function(x) {if (x==0) return(-.125) else return(-tanh(x/2)/(4*x))}
Cfunc <- function(x)  x/2 - log(1+exp(x)) + x*tanh(x/2)/4

######### Data generation ############
n <- 100
p <- 10

# for saving the ELBO
ELBOs <- rep(NA, p + 1)

# sensitivity_1=matrix(0,MAXITER,1)
# specificity_1=sensitivity_1
# sensitivity_100=sensitivity_1
# specificity_100=sensitivity_1
# sensitivity_cs1=sensitivity_1
# specificity_cs1=sensitivity_cs1
# sensitivity_cs100=sensitivity_1
# specificity_cs100=sensitivity_cs1
# setdiff1=sensitivity_1
# setdiff100=setdiff1

############################# generating the precision matrix.:Assume two discrete covariate levels
Lam1=matrix(0,p+1,1)
Lam2=Lam1
Lam1=c(3,3,3,3,rep(0,p-3))*5#For Z[i]=-0.1
Lam2=Lam1 #Same lambda for both covariate levels, corresponds to covariate independent levels
#Lam2=c(3,3,3,3,3,3,3,0,0,0,0)*5#For Z[i]= 0.1
# Lam2=c(rep(0,p-3),3,3,3,3)*5#For Z[i]= 0.1 #corresponds to covariate dependent model, uncomment to try this out.
Var1=solve(Lam1%*%t(Lam1) + diag(rep(10,p+1))) #covariance matrix for covariate level 1
Var2=solve(Lam2%*%t(Lam2) + diag(rep(10,p+1))) #covariance matrix for covariate level 2

X1=mvrnorm(n/2,rep(0,p+1),Var1)
X2=mvrnorm(n/2,rep(0,p+1),Var2)
######### Generating the covariates ##########


Z=matrix(-1,n,p) #Initializing the covariate matrix

beta=matrix(0,n,p) # Ground truth of the dependence structure
resp_index=1;# The index we consider as response
mylist <- rep(list(beta),p+1) #The variable specific inclusion probability matrix:ith row corresponds to the dependence structure for the i th subject, j th matrix corresponds to
                       # the j th variable as response and the remaining as predictors.
data_mat=rbind(X1,X2)




Adj_Mat_vb <- array(0,dim=c(p+1,p+1))
###############################################
for(resp_index in 1:(p+1)){ #This loops over the p+1 variables
  for(i in 1:n){
    beta[i,]=(t(Lam1[-resp_index])>0)*(i<=n/2) + (t(Lam2[-resp_index])>0)*(i>n/2) #Ground truth

    for(j in 1:p){
      # Z[i,j]=rnorm(1,-2,.1)*(i<50) +rnorm(1,2,0.1)*(i>=50) #Uncomment to lay with other covariate values
      Z[i,j]=-.1*(i<=n/2)  + .1*(i>n/2)
    }
  }
  ######################

  ##############
  y=data_mat[,resp_index]; #Set variable number `resp_index` as the response

  X_mat=data_mat[,-resp_index] #Set the remaining p variables as predictor.
  X_vec <- matrix(0,n*p,1)
  #X<- matrix(rep(0,n^2*p),nrow=n,ncol=n*p)

  X<- matrix(rep(0,n^2*p),nrow=n,ncol=n*p)

  for(i in 1:n){
    for(j in 1:p){
      k=p*(i-1)+1
      X[i,k+j-1]=X_mat[i,j]
      X_vec[k+j-1]=X[i,k+j-1]
    }
  }
  ELBO_LBit=rep(0,10000)
  Big_diag_mat <- matrix(rep(0,n^2*p),nrow=n,ncol=n*p)
  for(i in 1:n){
    k=p*(i-1)
    for(j in 1:p){
      Big_diag_mat[i,k+j]=1
    }
  }

  q=matrix(2,n,1)


  sigmasq=1 #Initialization of the hyperparameter value
  #E <- rnorm(n,0,sigmasq)


  XtX=t(X)%*%X

  DXtX=diag(XtX)
  DXtX_rep=rep(DXtX,p); dim(as.matrix(DXtX_rep))
  DXtX_mat=matrix(DXtX_rep,n*p,p,byrow=FALSE)
  Diff_mat=XtX-diag(DXtX)


  D=matrix(1,n,n)
  for(i in 1:n){
    for(j in 1:n){
      diff_vect <- Z[i,]-Z[j,]
      diff_norm <- as.numeric(sqrt(crossprod(diff_vect)))
      D[i,j]= dnorm(diff_norm,0,.1)
    }
  }
  for(i in 1:n){
    D[,i]=n*(D[,i]/sum(D[,i])) #Scaling the weights so that they add up to n
    #      D[,i]=1 # When there is no covariate information, set the weights to be 1 throughout.
  }


  alpha= rep(0.2,n*p) #Initialization of the inclusion probability matrix for a fixed variable, with i th row corresponding to i th subject.
  sigmabeta_sq=3 #Initialization for hyperparameter
  mu=rep(0,p) # Variational parameter
  true_pi=0.5 #Hyperparameter


  ###########



  y_long_vec=as.vector(t(y%*%matrix(1,1,p)))
  Xty=t(X)%*%as.vector(y)
  beta_mat=matrix(beta,n,p,byrow=TRUE)
  mu_mat=beta_mat


  D_long=matrix(0,n*p,n)
  for( i in 1:n){
    D_long[,i]=matrix(t(D[,i]%*%matrix(1,1,p)),n*p,1)
  }


  S_sq=matrix(sigmasq*(DXtX + 1/sigmabeta_sq)^(-1),n,p) #Initialization

  iter=1
  ###############################################




  ind_vec=seq(0,(n-1)*p,by=p)
  Ind_mat=matrix(0,n,p)
  for(j in 1:p){
    Ind_mat[,j]=ind_vec+j
  }
  Big_ind=matrix(0,n*p,p)
  Big_ind_1=matrix(0,n*p,p)
  for(j in 1:p){
    Big_ind[Ind_mat[,j],j]=1
    Big_ind_1[Ind_mat[,j],j]=0
  }

  DXtX_Big_ind=DXtX_mat*Big_ind
  ###################################################
  candL=seq(0.1,0.9,.2)# Different values of hyperparameter true_pi
  #candL=0.5
  like=rep(0,length(candL))
  elb=like

  est_pi=rep(0,n)
  est_q=est_pi
  beta_matr=matrix(0,n,p)

  ####################tuning hyperparameters##################################
  set.seed(resp_index)
  idmod=varbvs(X_mat,y,Z=Z[,1],verbose=FALSE)#Setting hyperparameter value as in Carbonetto Stephens model
  inprob=idmod$pip
  rest_index_set=setdiff(c(1:(p+1)),resp_index)

  sigmasq=mean(idmod$sigma)
  pi_est=mean(1/(1+exp(-idmod$logodds)))
  sigmavec=c(0.01,0.05,0.1,0.5,1,3,7,10)
  elb1=matrix(0,length(sigmavec),1)
  for(j in 1:length(sigmavec)){
    res=cov_vsvb(y,X,Z,XtX,DXtX,Diff_mat,Xty,sigmasq,sigmavec[j],pi_est)
    elb1[j]=res$var.elbo

  }
  sigmabeta_sq=sigmavec[which.max(elb1)] #Choosing hyperparameter based on ELBO maximization
  ELBOs[resp_index] <- max(elb1)
  result=cov_vsvb(y,X,Z,XtX,DXtX,Diff_mat,Xty,sigmasq, sigmabeta_sq,pi_est)
  incl_prob=result$var.alpha
  mu0_val=result$var.mu0_lambda

  #



  heat_alpha=matrix(incl_prob,n,p,byrow=TRUE)
  mylist[[resp_index]]=heat_alpha
}

alpha_matrices <- mylist

# Creating the graphs:
# transform p + 1 n by n matrices to n p + 1 by p + 1 matrices using alpha_matrices
# the j, k entry in the l-th matrix is the probability of inclusion of an edge
# between the j, k variables for the l-th individual
incl_probs <- replicate(n, matrix(0, p + 1, p + 1), simplify = F)

# iterate over the p matrices
for (j in 1:(p + 1)){

  # fix the j-th alpha matrix
  alpha_mat_j <- alpha_matrices[[j]]

  # iterate over the rows of alpha_mat_j
  for (l in 1:n){

    # the j-th row of the l-th individual's graph is the l-th row of
    # alpha_mat_j with a 0 in the j-th position
    incl_probs[[l]][j, -j] <- alpha_mat_j[l,]
  }
}

original_alpha_matrices <- incl_probs
original_incl_probs <- lapply(incl_probs, function(mat) (mat + t(mat)) / 2)
out_original <- list(original_alpha_matrices = original_alpha_matrices,
                     original_incl_probs = original_incl_probs,
                     original_ELBO = ELBOs)
save(out_original, file = "out_original_discrete.Rdata")

#################BELOW IS FOR VISUALIZATION#####################


original_graphs <- vector("list", n)

for (SUBJECT in 1:n){
  alph=matrix(0,p+1,p+1)
  for(i in 1:(p+1)){
    alph[i,-i]=mylist[[i]][SUBJECT,];
  }
  heat_alpha=alph


  a=heat_alpha
  for(i in 1:(p+1)){
    for(j in i:(p+1)){
      #    a[i,j]=max(heat_alpha[i,j],heat_alpha[j,i])
      a[i,j]=0.5*(heat_alpha[i,j]+heat_alpha[j,i])
      a[j,i]=a[i,j]
    }
  }
  original_graphs[[SUBJECT]] <- a
}
#save(original_graphs, file = "original_discrete_graphs.Rdata")
#heat_alpha=a



alph=matrix(0,p+1,p+1)

alph=matrix(0,p+1,p+1)

SUBJECT=1
for(i in 1:(p+1)){
  alph[i,-i]=mylist[[i]][SUBJECT,]; #Individual specific inclusion probability matrix
}
beta=matrix(0,p+1,p+1)
for(i in 1:(p+1)){
  for(j in 1:(p+1)){
    beta[i,j]=(Lam1[i]!=0 & Lam1[j]!=0)
  }}
diag(beta)=0

heat_alpha=alph

a=heat_alpha
for(i in 1:(p+1)){
  for(j in i:(p+1)){
    #  a[i,j]=max(heat_alpha[i,j],heat_alpha[j,i])
    a[i,j]=mean(c(heat_alpha[i,j],heat_alpha[j,i]))
    a[j,i]=a[i,j]
  }
}

heat_alpha=a



alphvec=sort(as.vector(heat_alpha[which(heat_alpha!=0)]))

selection1=1*(heat_alpha>0.5)



SUBJECT=100

for(i in 1:(p+1)){
  alph[i,-i]=mylist[[i]][SUBJECT,];
}
beta=matrix(0,p+1,p+1)
for(i in 1:(p+1)){
  for(j in 1:(p+1)){
    beta[i,j]=(Lam2[i]!=0 & Lam2[j]!=0)
  }}
diag(beta)=0

heat_alpha=alph


a=heat_alpha
for(i in 1:(p+1)){
  for(j in i:(p+1)){
    #  a[i,j]=max(heat_alpha[i,j],heat_alpha[j,i])
    a[i,j]=mean(c(heat_alpha[i,j],heat_alpha[j,i]))
    a[j,i]=a[i,j]
  }
}

heat_alpha=a


alphvec=sort(as.vector(heat_alpha[which(heat_alpha!=0)]))

selection1=1*(heat_alpha>0.5)



SUBJECT=1
for(i in 1:(p+1)){
  alph[i,-i]=mylist[[i]][SUBJECT,];
}
beta=matrix(0,p+1,p+1)
for(i in 1:(p+1)){
  for(j in 1:(p+1)){
    beta[i,j]=(Lam1[i]!=0 & Lam1[j]!=0)
  }}
diag(beta)=0

data = melt(t(beta))
fig = ggplot(data, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "brown") +
  scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0),guide = "legend" )

fig = fig + scale_x_continuous(expand = c(0, 0))
fig = fig + scale_y_continuous( expand = c(0, 0))


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("True Dependence")) )

fig = fig + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(),

                   panel.background = element_blank())
fig=fig + theme(plot.title = element_text(hjust = 0.5))


fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"),

                  axis.title = element_text(size=30, face = "bold"))

fig = fig + theme( legend.title = element_text(face = "bold",size = 25),

                   legend.text = element_text(face="bold",size = 25),

                   legend.key.size = unit(2,'lines'))
fig=fig+ coord_equal()

fig


heat_alpha=alph


a=heat_alpha
for(i in 1:(p+1)){
  for(j in i:(p+1)){
  #  a[i,j]=max(heat_alpha[i,j],heat_alpha[j,i])
    a[i,j]=0.5*(heat_alpha[i,j]+heat_alpha[j,i])
    a[j,i]=a[i,j]
  }
}

#heat_alpha=a

data = melt(t(a))
fig = ggplot(data, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "brown") +
  scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0),guide = "colorbar" )

fig = fig + scale_x_continuous(expand = c(0, 0))
fig = fig + scale_y_continuous( expand = c(0, 0))


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("Inclusion Probability for Covariate Level 1")) )

fig = fig + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(),

                   panel.background = element_blank())
fig=fig + theme(plot.title = element_text(hjust = 0.5))


fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"),

                  axis.title = element_text(size=30, face = "bold"))

fig = fig + theme( legend.title = element_text(face = "bold",size = 25),

                   legend.text = element_text(face="bold",size = 25),

                   legend.key.size = unit(2,'lines'))
fig=fig+ coord_equal()

fig

alphvec=setdiff(as.vector(heat_alpha),diag(heat_alpha))

selection1=1*(a>0.5)
data = melt(t(selection1))
fig = ggplot(data, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "brown") +
  scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0),guide = "legend" )

fig = fig + scale_x_continuous(expand = c(0, 0))
fig = fig + scale_y_continuous( expand = c(0, 0))


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("Graph Estimate For Covariate Level 1")) )

fig = fig + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(),

                   panel.background = element_blank())
fig=fig + theme(plot.title = element_text(hjust = 0.5))


fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"),

                  axis.title = element_text(size=30, face = "bold"))

fig = fig + theme( legend.title = element_text(face = "bold",size = 25),

                   legend.text = element_text(face="bold",size = 25),

                   legend.key.size = unit(2,'lines'))
fig=fig+ coord_equal()

fig

SUBJECT=100

for(i in 1:(p+1)){
  alph[i,-i]=mylist[[i]][SUBJECT,];
}
beta=matrix(0,p+1,p+1)
for(i in 1:(p+1)){
  for(j in 1:(p+1)){
    beta[i,j]=(Lam2[i]!=0 & Lam2[j]!=0)
  }}
diag(beta)=0
data = melt(t(beta))
fig = ggplot(data, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "brown") +
  scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0),guide = "legend" )

fig = fig + scale_x_continuous(expand = c(0, 0))
fig = fig + scale_y_continuous( expand = c(0, 0))


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("True Dependence")) )

fig = fig + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(),

                   panel.background = element_blank())
fig=fig + theme(plot.title = element_text(hjust = 0.5))


fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"),

                  axis.title = element_text(size=30, face = "bold"))

fig = fig + theme( legend.title = element_text(face = "bold",size = 25),

                   legend.text = element_text(face="bold",size = 25),

                   legend.key.size = unit(2,'lines'))
fig=fig+ coord_equal()

fig

heat_alpha=alph


a=heat_alpha
for(i in 1:(p+1)){
  for(j in i:(p+1)){
  #  a[i,j]=max(heat_alpha[i,j],heat_alpha[j,i])
    a[i,j]=0.5*(heat_alpha[i,j]+heat_alpha[j,i])
    a[j,i]=a[i,j]
  }
}

#heat_alpha=a
data = melt(t(a))
fig = ggplot(data, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "brown") +
  scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0),guide = "colorbar" )

fig = fig + scale_x_continuous(expand = c(0, 0))
fig = fig + scale_y_continuous( expand = c(0, 0))


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("Inclusion Probability For Covariate Level 2")) )

fig = fig + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(),

                   panel.background = element_blank())
fig=fig + theme(plot.title = element_text(hjust = 0.5))


fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"),

                  axis.title = element_text(size=30, face = "bold"))

fig = fig + theme( legend.title = element_text(face = "bold",size = 25),

                   legend.text = element_text(face="bold",size = 25),

                   legend.key.size = unit(2,'lines'))
fig=fig+ coord_equal()
fig

alphvec=setdiff(as.vector(heat_alpha),diag(heat_alpha))


selection1=1*(a>0.5)
data = melt(t(selection1))
fig = ggplot(data, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "brown") +
  scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0),guide = "legend" )

fig = fig + scale_x_continuous(expand = c(0, 0))
fig = fig + scale_y_continuous( expand = c(0, 0))


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("Graph Estimate For Covariate Level 2")) )

fig = fig + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(),

                   panel.background = element_blank())

fig=fig+ coord_equal()
fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"),

                  axis.title = element_text(size=30, face = "bold"))
fig=fig + theme(plot.title = element_text(hjust = 0.5))
fig = fig + theme( legend.title = element_text(face = "bold",size = 25),

                   legend.text = element_text(face="bold",size = 25),

                   legend.key.size = unit(2,'lines'))

fig


# things that have changed since the "true" original version
# 1. Correct mupdate
# 2. Set seed below
# 3. Take the l2 norm instead of the spectral norm
# 4. Save the ELBO
# 5. turn the alpha matrices into symmetric inclusion matrices
# 6. save the alpha matrices, inclusion probability matrices, and ELBO

set.seed(1)
rm(list=ls())
setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/dev/Archive/original_code_Covariate-dependent-Graph-Estimation-main")
source("cov_vsvb.R")
source("ELBO_calculator.R")
library(Matrix)
library(magic)
library(psych)
library(ggplot2)
library(reshape2)
library(MASS)
library(pracma)
library(latex2exp)
library(ks)
library(varbvs)
###################################################################

logit <- function(x) {if ((x==0)|(x==1)) return(0) else return(log(x/(1-x))) }
Afunc <- function(x) {if (x==0) return(-.125) else return(-tanh(x/2)/(4*x))}
Cfunc <- function(x)  x/2 - log(1+exp(x)) + x*tanh(x/2)/4
#####FUNCTION DEFINITION#####
#############################
#############################

#############################
#############################
#############################
######### Data generation ############
n <- 180
p <- 4

# for saving the ELBO
ELBOs <- rep(NA, p + 1)

MAXITER=1
STR=1
in_pr_13=matrix(0,MAXITER,n)
in_pr_12=in_pr_13
Var_cont=function(z){
  pr=matrix(0,p+1,p+1)
  diag(pr)=2
  #  pr[1,2]=STR*((z>0) && (z< .33)) + (STR - STR*((z-.33)/.33))*((z>0.33) && (z<0.66)) + (0)*((z>0.66) && (z<1))
  #  pr[1,3]=0*((z>0) && (z< .33)) + (STR*((z-.33)/.33))*((z>0.33) && (z<0.66)) + (STR)*((z>0.66) && (z<1))
  pr[2,3]=STR
  pr[1,2]=STR*((z>-1) && (z< -.33)) + (STR - STR*((z+.23)/.56))*((z>-0.23) && (z<0.33)) + (0)*((z>0.43) && (z<1))
  pr[1,3]=0*((z>-1) && (z< -.33)) + (STR*((z+.23)/.56))*((z>-0.23) && (z<0.33)) + (STR)*((z>0.43) && (z<1))

  pr[2,1]=pr[1,2]
  pr[3,1]=pr[1,3]
  pr[3,2]=pr[2,3]


  Var=solve(pr)
  return(Var)
}
sensitivity_20=matrix(0,MAXITER,1)
specificity_20=sensitivity_20
sensitivity_90=sensitivity_20
specificity_90=sensitivity_20
sensitivity_160=sensitivity_20
specificity_160=sensitivity_20
for(overiter in 1:MAXITER){
  # Z=c(runif(n/3,-1,-0.33),runif(n/3, -0.23,0.33),runif(n/3,0.43,1))
  Z=c(seq(-0.99,-0.331,(-.331+.99)/59),seq(-0.229,0.329,(.329+.229)/59),seq(0.431,.99,(.99-.431)/59))
  # Z=seq(0.01,.99,.98/(n-1))
  Z=matrix(Z,n,1)
  X=matrix(0,n,p+1)
  for(i in 1:n){
    X[i,]=  mvrnorm(1,rep(0,p+1),Var_cont(Z[i]))
  }

  ######### Generating the covariates ##########


  beta=matrix(0,n,p)
  resp_index=1;# The index we consider as response
  mylist <- rep(list(beta),p+1)
  data_mat=X





  ###############################################
  for(resp_index in 1:(p+1)){


    ##############
    y=data_mat[,resp_index];

    X_mat=data_mat[,-resp_index]
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
    A_xi=rep(1,n)
    #X <- matrix(rnorm(n*p,0,1), ncol=p)
    q=matrix(2,n,1)


    # for(i in 1:n){
    #   q[i]=floor()
    #   for(j in 1:q[i]){
    #     beta[(i-1)*p+j] <- 4 # runif(1,3,3.8) ###c(rep(3.5,q),rep(0,p-q))
    #   }
    # }
    #beta ### Done
    sigmasq=1
    E <- rnorm(n,0,sigmasq)
    #snr=sd(X%*%as.vector(beta0))/sd(sigmasq*E)

    #y

    XtX=t(X)%*%X

    DXtX=diag(XtX)
    DXtX_rep=rep(DXtX,p); dim(as.matrix(DXtX_rep))
    DXtX_mat=matrix(DXtX_rep,n*p,p,byrow=FALSE)
    Diff_mat=XtX-diag(DXtX)


    tt1=kde(Z[,1],eval.points=Z[,1])
    h1=tt1$h
    # tt2=kde(Z[,2],eval.points=Z[,2])
    # h2=tt2$h

    D=matrix(1,n,n)
    for(i in 1:n){
      for(j in 1:n){
        # D[j,i]= dnorm(sqrt(norm((Z[i,1]-Z[j,1])/h1,"2")^2 + norm((Z[i,2]-Z[j,2])/h2,"2")^2),0,1)
        diff_vect <- Z[i,]-Z[j,]
        diff_norm <- as.numeric(sqrt(crossprod(diff_vect)))
        D[j,i]= dnorm(diff_norm,0,0.56)
      }
    }
    for(i in 1:n){
      D[,i]=n*(D[,i]/sum(D[,i]))
      #    D[,i]=1
    }
    #true_lambda=rep(-log((p-q)/q)/p,p)
    true_lambda=0.5*rep(1,p)
    L0=0.5
    lambda_mean=true_lambda##rep(0,p) ###rnorm(p,0,4)
    lambda_var=.001*diag(p)
    mu0_lambda<-L0*rep(1,p)## rep(0,p)
    Sigma0_lambda=lambda_var###diag(p)
    alpha= rep(0.2,n*p)
    sigmabeta_sq=3
    mu=rep(0,p)
    true_pi=0.5


    ###########



    y_long_vec=as.vector(t(y%*%matrix(1,1,p)))
    Xty=t(X)%*%as.vector(y)
    beta_mat=matrix(beta,n,p,byrow=TRUE)
    mu_mat=beta_mat


    D_long=matrix(0,n*p,n)
    for( i in 1:n){
      D_long[,i]=matrix(t(D[,i]%*%matrix(1,1,p)),n*p,1)
    }


    S_sq=matrix(sigmasq*(DXtX + 1/sigmabeta_sq)^(-1),n,p)
    Sigma_xi=diag(p)
    S0=solve(lambda_var)
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

    ####################tuning hyperparameters using Carbonetto Stephens##################################
    idmod=varbvs(X_mat,y,Z=Z[,1],verbose=FALSE)
    sigmasq=mean(idmod$sigma)
    pi_est=mean(1/(1+exp(-idmod$logodds)))
    sigmavec=c(0.01,0.05,0.1,0.5,1,3,7,10)
    elb1=matrix(0,length(sigmavec),1)
    for(j in 1:length(sigmavec)){
      res=cov_vsvb(y,X,Z,XtX,DXtX,Diff_mat,Xty,sigmasq,sigmavec[j],pi_est)
      elb1[j]=res$var.elbo

    }
    sigmabeta_sq=sigmavec[which.max(elb1)]
    ELBOs[resp_index] <- max(elb1)
    result=cov_vsvb(y,X,Z,XtX,DXtX,Diff_mat,Xty,sigmasq, sigmabeta_sq,pi_est)
    incl_prob=result$var.alpha
    mu0_val=result$var.mu0_lambda

    #



    heat_alpha=matrix(incl_prob,n,p,byrow=TRUE)
    mylist[[resp_index]]=heat_alpha
  }
  alph=matrix(0,p+1,p+1)
  #For first subject
  SUBJECT=20
  for(i in 1:(p+1)){
    alph[i,-i]=mylist[[i]][SUBJECT,];
  }
  beta=1*(solve(Var_cont(Z[SUBJECT]))>0.01)
  diag(beta)=0

  heat_alpha=alph


  a=heat_alpha
  for(i in 1:(p+1)){
    for(j in i:(p+1)){
  #    a[i,j]=max(heat_alpha[i,j],heat_alpha[j,i])
      a[i,j]=0.5*(heat_alpha[i,j]+heat_alpha[j,i])
      a[j,i]=a[i,j]
    }
  }

  heat_alpha=a


  alphvec=setdiff(as.vector(heat_alpha),diag(heat_alpha))

  selection1=1*(heat_alpha>0.5)


  sensitivity_20[overiter]=sum(selection1==1 & beta==1)/sum(beta==1)
  specificity_20[overiter]=sum(selection1==0 & beta==0)/sum(beta==0)

  SUBJECT=90

  for(i in 1:(p+1)){
    alph[i,-i]=mylist[[i]][SUBJECT,];
  }
  beta=1*(solve(Var_cont(Z[SUBJECT]))>0.01)
  diag(beta)=0

  heat_alpha=alph


  a=heat_alpha
  for(i in 1:(p+1)){
    for(j in i:(p+1)){
  #    a[i,j]=max(heat_alpha[i,j],heat_alpha[j,i])
      a[i,j]=0.5*(heat_alpha[i,j]+heat_alpha[j,i])
      a[j,i]=a[i,j]
    }
  }

  heat_alpha=a


  alphvec=setdiff(as.vector(heat_alpha),diag(heat_alpha))

  selection1=1*(heat_alpha>0.5)


  sensitivity_90[overiter]=sum(selection1==1 & beta==1)/sum(beta==1)
  specificity_90[overiter]=sum(selection1==0 & beta==0)/sum(beta==0)

  SUBJECT=160

  for(i in 1:(p+1)){
    alph[i,-i]=mylist[[i]][SUBJECT,];
  }
  beta=1*(solve(Var_cont(Z[SUBJECT]))>0.01)
  diag(beta)=0

  heat_alpha=alph


  a=heat_alpha
  for(i in 1:(p+1)){
    for(j in i:(p+1)){
  #    a[i,j]=max(heat_alpha[i,j],heat_alpha[j,i])
      a[i,j]=0.5*(heat_alpha[i,j]+heat_alpha[j,i])
      a[j,i]=a[i,j]
    }
  }

  heat_alpha=a


  alphvec=setdiff(as.vector(heat_alpha),diag(heat_alpha))

  selection1=1*(heat_alpha>0.5)


  sensitivity_160[overiter]=sum(selection1==1 & beta==1)/sum(beta==1)
  specificity_160[overiter]=sum(selection1==0 & beta==0)/sum(beta==0)
  for(SUBJECT in 1:n){
    alph=matrix(0,p+1,p+1)


    for(i in 1:(p+1)){
      alph[i,-i]=mylist[[i]][SUBJECT,];
    }
    in_pr_13[overiter,SUBJECT]=alph[1,3]
    in_pr_12[overiter,SUBJECT]=alph[1,2]
  }

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
save(out_original, file = "out_original_continuous.Rdata")
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
#save(original_graphs, file = "original_continuous_graphs.Rdata")
#heat_alpha=a

data = melt(t(a))
fig = ggplot(data, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "brown") +
  scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0),guide = "colorbar" )

fig = fig + scale_x_continuous(expand = c(0, 0))
fig = fig + scale_y_continuous( expand = c(0, 0))


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("Inclusion Probability for Subject 20")) )

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


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("Graph Estimate for Subject 20")) )

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
###########
SUBJECT=90
for(i in 1:(p+1)){
  alph[i,-i]=mylist[[i]][SUBJECT,];
}
beta=1*(solve(Var_cont(Z[SUBJECT]))>0.01)
diag(beta)=0

data = melt(t(beta))
fig = ggplot(data, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "brown") +
  scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0),guide = "legend" )

fig = fig + scale_x_continuous(expand = c(0, 0))
fig = fig + scale_y_continuous( expand = c(0, 0))


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("True Dependence Structure")) )

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
 #   a[i,j]=max(heat_alpha[i,j],heat_alpha[j,i])
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


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("Inclusion Probability For Subject 90")) )

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


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("Graph Estimate for Subject 90")) )

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


###########
SUBJECT=160

for(i in 1:(p+1)){
  alph[i,-i]=mylist[[i]][SUBJECT,];
}
beta=1*(solve(Var_cont(Z[SUBJECT]))>0.01)
diag(beta)=0
data = melt(t(beta))
fig = ggplot(data, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "brown") +
  scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0),guide = "legend" )

fig = fig + scale_x_continuous(expand = c(0, 0))
fig = fig + scale_y_continuous( expand = c(0, 0))


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("True Dependence Structure")) )

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


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("Inclusion Probability for Subject 160")) )

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

#selection1=1*(heat_alpha>alphvec[kmeans(alphvec,c(min(alphvec),max(alphvec)))$size[1]])
selection1=1*(a>0.5)
data = melt(t(selection1))
fig = ggplot(data, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "brown") +
  scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0),guide = "legend" )

fig = fig + scale_x_continuous(expand = c(0, 0))
fig = fig + scale_y_continuous( expand = c(0, 0))


fig = fig + labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("Graph Estimate for Subject 160")) )

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

data=melt(t(sensitivity_20))
fig=ggplot(data, aes(x=sensitivity_20)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")
fig = fig + labs( x=expression(bold("Sensitivity For Subject 20")), y=expression(bold("density")) )
fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"),
                  axis.title = element_text(size=30, face = "bold"))
fig


data=melt(t(specificity_20))
fig=ggplot(data, aes(x=specificity_20)) +
  geom_histogram(aes(y=..density..), colour="black", bins=10,  fill="white")+
  geom_density(alpha=.2, fill="#FF6666")
fig = fig + labs( x=expression(bold("Specificity For Subject 20")), y=expression(bold("density")) )
fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"),
                  axis.title = element_text(size=30, face = "bold"))
fig

data=melt(t(sensitivity_90))
fig=ggplot(data, aes(x=sensitivity_90)) +
  geom_histogram(aes(y=..density..), colour="black",bins=10, fill="white")+
  geom_density(alpha=.2, fill="#FF6666")
fig = fig + labs( x=expression(bold("Sensitivity For Subject 90")), y=expression(bold("density")) )
fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"),
                  axis.title = element_text(size=30, face = "bold"))
fig


data=melt(t(specificity_90))
fig=ggplot(data, aes(x=specificity_90)) +
  geom_histogram(aes(y=..density..), colour="black", bins=10,  fill="white")+
  geom_density(alpha=.2, fill="#FF6666")
fig = fig + labs( x=expression(bold("Specificity For Subject 90")), y=expression(bold("density")) )
fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"),
                  axis.title = element_text(size=30, face = "bold"))
fig
data=melt(t(sensitivity_160))
fig=ggplot(data, aes(x=sensitivity_160)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")
fig = fig + labs( x=expression(bold("Sensitivity For Subject 160")), y=expression(bold("density")) )
fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"),
                  axis.title = element_text(size=30, face = "bold"))
fig


data=melt(t(specificity_160))
fig=ggplot(data, aes(x=specificity_160)) +
  geom_histogram(aes(y=..density..), colour="black", bins=10,  fill="white")+
  geom_density(alpha=.2, fill="#FF6666")
fig = fig + labs( x=expression(bold("Specificity For Subject 160")), y=expression(bold("density")) )
fig = fig + theme(axis.text = element_text(size=15, face = "bold", colour = "black"),
                  axis.title = element_text(size=30, face = "bold"))
fig

al=matrix(0,1,n)

for(SUBJECT in 1:n){
  alph=matrix(0,p+1,p+1)


  for(i in 1:(p+1)){
    alph[i,-i]=mylist[[i]][SUBJECT,];
  }
  al[SUBJECT]=alph[1,3]
}

grd=seq(1,n,1)
grd=matrix(grd,n,1)
al=matrix(al,n,1)
data=as.data.frame(cbind(grd,al))
fig=ggplot(data, aes(x=grd, y=al)) +
  geom_line()+
  geom_point()
fig=fig+labs( x=expression(bold("Subject Index")), y=expression(bold("Inclusion Probability")) )
fig=fig +  theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA),
  plot.title = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
  axis.title.x = element_text(size=16, face="bold", colour = "black"),
  axis.title.y = element_text(size=16, face="bold", colour = "black"),
  axis.text.x = element_text(size=16, face="bold", colour = "black"),
  axis.text.y = element_text(size=16, face="bold", colour = "black"),
  strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
  strip.text.y = element_text(size = 14, face="bold", colour = "black"),
  axis.line.x = element_line(color="black", size = 0.3),
  axis.line.y = element_line(color="black", size = 0.3),
  panel.border = element_rect(colour = "black", fill=NA, size=0.3)
)
fig=fig+ ggtitle(TeX('Inclusion Probability of edge between $x_1$ and $x_3$ versus covariate'))
fig
#save(mylist,file="converged_alpha_per_regression_continuouscov_p4.RData")

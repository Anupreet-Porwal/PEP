gvs.crpep<-function(y,X,iter,discard,family,n.star=nrow(y),model='search',model.prob='beta',hyper=FALSE,hyper.type='hyper-g/n',hyper.multiplier=1,ini.b='none',ini.b0='none',ini.y='none',variables=rep(1,ncol(X))) {
Begin<-Sys.time()

if(model=='single') {
if(length(variables)!=dim(X)[2]) stop('Length of variables is not equal to the number of columns of design matrix')
ini.gamma<-c(1,variables)
include<-ini.gamma==1
exclude<-ini.gamma==0
}

############################### BINOMIAL

if(family=='binomial') {

trials<-y[,1]
success<-y[,2]
fail<-trials-success
Y<-cbind(success,fail)
n<-length(trials)

 
if(ini.b0=='none') {
glm.ini.b0<-glm(Y~1,family=binomial)
ini.b0<-glm.ini.b0$coefficients
}

ini.p0<-exp(ini.b0)/(1+exp(ini.b0))

if(class(ini.y)=='numeric') {
n.star<-length(ini.y)
trials.star<-rep(trials,length.out=n.star)
}
if(ini.y=='none'){
trials.star<-rep(trials,length.out=n.star)
ini.y<-rbinom(n.star,trials.star,ini.p0)
}
if(n.star<n) stop('Size of imaginary data must be at least equal to real sample size')

delta<-n.star
inv.delta<-1/n.star

ini.fail<-trials.star-ini.y        
Y.star<-cbind(ini.y,ini.fail)
Y.all<-rbind(Y,Y.star)

X_1<-cbind(1,X)
if(is.null(colnames(X))) {colnames(X)<-paste(rep('X',ncol(X)),c(1:ncol(X)),sep='')}
colnames(X_1)<-c('intercept',colnames(X))
X.star<-X[rep(1:n,length.out=n.star),]
X.star_1<-cbind(1,X.star)
X.all<-rbind(X,X.star)
	
weights.star<-rep(inv.delta,n.star)
weights.all<-c(rep(1,n),rep(inv.delta,n.star))

p<-dim(X_1)[2]

if(model=='search'){
ini.gamma<-c(1,rbinom((p-1),1,0.5))     ##yparxei provlima gia p=0 exei sxesh me thn log posterior twn beta pou einai se comment
include<-ini.gamma==1
exclude<-ini.gamma==0
act.model<-rep(NA,iter)
act.predictors<-rep(NA,iter)
}

if(model=='single') {
if(ini.b=='none')   {
glm.ini.b<-glm(Y.all~X.all[,include[-1]],family=binomial)
ini.b<-glm.ini.b$coefficients
}
ini.b.gamma<-rep(0,p)
ini.b.gamma[include]<-ini.b
ini.p.gamma<-exp(X_1%*%ini.b.gamma)/(1+exp(X_1%*%ini.b.gamma))
}


if(model=='search') {
if(ini.b=='none')   {
glm.ini.b<-glm(Y.all~X.all,family=binomial)
ini.b<-glm.ini.b$coefficients
}
ini.b.gamma<-ini.b*ini.gamma
ini.p.gamma<-exp(X_1%*%ini.b.gamma)/(1+exp(X_1%*%ini.b.gamma))
}

log.posterior.b<-function(beta,Pi,success,fail,success_star,fail_star,betaML,SigmaML) {

logPi<-log(Pi)
log1_Pi<-log(1-Pi)
logPi[logPi==-Inf]<- -exp(709.78)
log1_Pi[log1_Pi==-Inf]<- -exp(709.78)
logPi_star<-rep(logPi,length.out=n.star)
log1_Pi_star<-rep(log1_Pi,length.out=n.star)

if(sum(include)==1){
prob<-sum(success*logPi+fail*log1_Pi)+
      inv.delta*sum(success_star*logPi_star+fail_star*log1_Pi_star)-
      dnorm(beta,betaML,SigmaML,log=TRUE)
} else {
prob<-sum(success*logPi+fail*log1_Pi)+
      inv.delta*sum(success_star*logPi_star+fail_star*log1_Pi_star)-
      dmvnorm(beta,betaML,SigmaML,log=TRUE)
}

 if(prob>= exp(710))  {prob<- exp(709.78)}
 if(prob<= -exp(710)) {prob<- -exp(709.78)}
prob

}

log.posterior.b0<-function(beta0,success_star,fail_star,Pi0,beta0ML,Sigma0ML) {

Pi0<-rep(Pi0,length.out=n.star)
prob<-sum(success_star*log(Pi0)+fail_star*log(1-Pi0))-
      dnorm(beta0,beta0ML,Sigma0ML,log=TRUE)
 if(prob>= exp(710))  {prob<- exp(709.78)}
 if(prob<= -exp(710)) {prob<- -exp(709.78)}
prob

}

log.posterior.y<-function(y_star,trials_star,Pi_star,Pi0_star,PiML_star,DetML_star,PiProp_star) {

Pi_star<-rep(Pi_star,length.out=n.star)
Pi0_star<-rep(Pi0_star,length.out=n.star)

prob<-inv.delta*sum(dbinom(y_star,trials_star,Pi_star,log=TRUE))+
      sum(dbinom(y_star,trials_star,Pi0_star,log=TRUE))-
      0.5*log(DetML_star)-inv.delta*sum(dbinom(y_star,trials_star,PiML_star,log=TRUE))-
      sum(dbinom(y_star,trials_star,PiProp_star,log=TRUE))
 if(prob>= exp(710))  {prob<- exp(709.78)}
 if(prob<= -exp(710)) {prob<- -exp(709.78)}
prob

}

if(hyper==TRUE) {

log.posterior.delta<-function(cand.delta,y_star,trials_star,Pi_star,PiML_star,DetML_star,propA,propB) {

Pi_star<-rep(Pi_star,length.out=n.star)
cand.inv.delta<-1/cand.delta

if(hyper.type=='gamma')     {prior.d<-dgamma(cand.delta,shape=n.star*0.001,rate=0.001,log=TRUE)}
if(hyper.type=='ZS')        {prior.d<- -3/2*log(cand.delta)-n/(2*cand.delta)}
if(hyper.type=='hyper-g')   {prior.d<- -(3/2)*log(1+cand.delta)}
if(hyper.type=='hyper-g/n') {prior.d<- -(3/2)*log(1+cand.delta/n.star)}
if(hyper.type=='maruyama')  {prior.d<- ((n-sum(include)+1)/2-7/4)*log(cand.delta)-((n-sum(include))/2)*log(1+cand.delta)}

prob<-prior.d+
      cand.inv.delta*sum(dbinom(y_star,trials_star,Pi_star,log=TRUE))-
      0.5*log(DetML_star)-cand.inv.delta*sum(dbinom(y_star,trials_star,PiML_star,log=TRUE))-
      dgamma(cand.delta,shape=propA,rate=propB,log=TRUE)
 if(prob>= exp(710))  {prob<- exp(709.78)}
 if(prob<= -exp(710)) {prob<- -exp(709.78)}
prob

}

}

cand.b<-rep(NA,p)
cand.b.gamma<-matrix(0,iter,p)
cand.y<-matrix(NA,iter,n.star)
if(model=='search'){
cand.gamma<-matrix(NA,iter+1,p-1)
cand.gamma[1,]<-ini.gamma[-1]
}

if(hyper==TRUE) {
cand.delta<-rep(NA,iter)
ratio<-rep(0,4)
names(ratio)<-c('betas','beta0','y*','delta')
} else {
ratio<-rep(0,3)
names(ratio)<-c('betas','beta0','y*')
}

if(model=='single'){
glm.b.gamma<-glm(Y.all~X.all[,include[-1]],weights=weights.all,family=binomial)
BetaML<-glm.b.gamma$coefficients
CovML<-vcov(glm.b.gamma)
}


if(model=='search'){
glm.b.gamma<-glm(Y~X,family=binomial)
BetaML<-glm.b.gamma$coefficients
CovML<-vcov(glm.b.gamma)
CovML.ind<-diag(diag(CovML))
}

for(i in 1:iter) {

##### Generation and MH for beta's gamma=1


u.b<-log(runif(1))

if(model=='single'){
cand.b<-rmvnorm(1,BetaML,CovML)
cand.b.gamma[i,include]<-cand.b
cand.p.gamma<-exp(X_1%*%cand.b.gamma[i,])/(1+exp(X_1%*%cand.b.gamma[i,]))

mh.ratio.b<-log.posterior.b(cand.b.gamma[i,include],cand.p.gamma,success,fail,ini.y,ini.fail,BetaML,CovML)-
            log.posterior.b(ini.b.gamma[include],ini.p.gamma,success,fail,ini.y,ini.fail,BetaML,CovML)

if(u.b <= mh.ratio.b) {                  
ini.b.gamma<-cand.b.gamma[i,]
ini.p.gamma<-cand.p.gamma
ratio[1]<-ratio[1]+1
} else {
cand.b.gamma[i,]<-ini.b.gamma
cand.p.gamma<-ini.p.gamma
cand.b<-ini.b                    
}
}


if(model=='search'){

if(sum(include)==1){
glm.b.gamma<-glm(Y.all~1,weights=weights.all,family=binomial)
BetaML.INC<-glm.b.gamma$coefficients
CovML.INC<-sqrt(vcov(glm.b.gamma))
cand.b<-rnorm(1,BetaML.INC,CovML.INC)
}
if(sum(include)!=1){
glm.b.gamma<-glm(Y.all~X.all[,include[-1]],weights=weights.all,family=binomial)
BetaML.INC<-glm.b.gamma$coefficients
CovML.INC<-vcov(glm.b.gamma)
cand.b<-rmvnorm(1,BetaML.INC,CovML.INC)
}

cand.b.gamma[i,include]<-cand.b
cand.p.gamma<-exp(X_1%*%cand.b.gamma[i,])/(1+exp(X_1%*%cand.b.gamma[i,]))

mh.ratio.b<-log.posterior.b(cand.b,cand.p.gamma,success,fail,ini.y,ini.fail,BetaML.INC,CovML.INC)-
            log.posterior.b(ini.b.gamma[include],ini.p.gamma,success,fail,ini.y,ini.fail,BetaML.INC,CovML.INC)

if(u.b <= mh.ratio.b) {                     
ini.b.gamma<-cand.b.gamma[i,]
ini.p.gamma<-cand.p.gamma
ratio[1]<-ratio[1]+1
} else {                 
cand.b.gamma[i,]<-ini.b.gamma
cand.p.gamma<-ini.p.gamma
cand.b<-ini.b[include]                    
}
}

##### Generation beta's gamma=0

if(model=='search') {
cand.b.gamma.exl<-rmvnorm(1,BetaML,CovML.ind)
cand.b.gamma.exl<-cand.b.gamma.exl[exclude]
}

##### Generation and MH for beta0 - NULL model

glm.b0<-glm(Y.star~1,family=binomial)
BetaML0<-glm.b0$coefficients
CovML0<-sqrt(vcov(glm.b0))

cand.b0<-rnorm(1,BetaML0,CovML0)
cand.p0<-exp(cand.b0)/(1+exp(cand.b0))
u.b0<-log(runif(1))

mh.ratio.b0<-log.posterior.b0(cand.b0,ini.y,ini.fail,cand.p0,BetaML0,CovML0)-
             log.posterior.b0(ini.b0,ini.y,ini.fail,ini.p0,BetaML0,CovML0)

if(u.b0 <= mh.ratio.b0)  {
ini.b0<-cand.b0
ini.p0<-cand.p0
ratio[2]<-ratio[2]+1
} else {
cand.b0<-ini.b0
cand.p0<-ini.p0
}

##### Generation and MH for y*

p.STAR<-(cand.p0*cand.p.gamma^inv.delta)/(cand.p0*cand.p.gamma^inv.delta+(1-cand.p0)*(1-cand.p.gamma)^inv.delta)
p.STAR<-rep(p.STAR,length.out=n.star)

cand.y[i,]<-rbinom(n.star,trials.star,p.STAR)
cand.fail<-trials.star-cand.y[i,]       
Y.star.cand<-cbind(cand.y[i,],cand.fail)


# Laplace approximation MLE at iteration t-1

if(sum(include)==1){glm.prev<-glm(Y.star~1,family=binomial)
betaML.prev<-glm.prev$coefficients
piML.prev<-exp(betaML.prev)
if(piML.prev==0) {
piML.prev<-10^-1
} else {
piML.prev[piML.prev==Inf]<-exp(709.78)
piML.prev<-piML.prev/(1+piML.prev)
}
H.prev<-diag(c(trials.star*piML.prev*(1-piML.prev)))
Det.prev<-delta/sum(H.prev)
}
if(sum(include)!=1){glm.prev<-glm(Y.star~X.star[,include[-1]],family=binomial)
betaML.prev<-glm.prev$coefficients
piML.prev<-exp(X.star_1[,include]%*%betaML.prev)
if(sum(piML.prev)==0) {
piML.prev<-rep(10^-1,n.star)
} else {
piML.prev[piML.prev==Inf]<-exp(709.78)
piML.prev<-piML.prev/(1+piML.prev)
}
H.prev<-diag(c(trials.star*piML.prev*(1-piML.prev)))
Det.prev<-det(delta*solve(t(X.star_1[,include])%*%H.prev%*%X.star_1[,include]))
}

# Laplace approximation MLE at iteration t

if(sum(include)==1){glm.curr<-glm(Y.star.cand~1,family=binomial)
betaML.curr<-glm.curr$coefficients
piML.curr<-exp(betaML.curr)
if(piML.curr==0) {
piML.curr<-10^-1
} else {
piML.curr[piML.curr==Inf]<-exp(709.78)
piML.curr<-piML.curr/(1+piML.curr)
}
H.curr<-diag(c(trials.star*piML.curr*(1-piML.curr)))
Det.curr<-delta/sum(H.curr)
}
if(sum(include)!=1){glm.curr<-glm(Y.star.cand~X.star[,include[-1]],family=binomial)
betaML.curr<-glm.curr$coefficients
piML.curr<-exp(X.star_1[,include]%*%betaML.curr)
if(sum(piML.curr)==0) {
piML.curr<-rep(10^-1,n.star)
} else {
piML.curr[piML.curr==Inf]<-exp(709.78)
piML.curr<-piML.curr/(1+piML.curr)
}
H.curr<-diag(c(trials.star*piML.curr*(1-piML.curr)))
Det.curr<-det(delta*solve(t(X.star_1[,include])%*%H.curr%*%X.star_1[,include]))
}

u.y<-log(runif(1))

mh.ratio.y<-log.posterior.y(cand.y[i,],trials.star,cand.p.gamma,cand.p0,piML.curr,Det.curr,p.STAR)-
            log.posterior.y(ini.y,trials.star,cand.p.gamma,cand.p0,piML.prev,Det.prev,p.STAR)
                 
if(u.y <= mh.ratio.y) {
ini.y<-cand.y[i,]
ini.fail<-cand.fail
Y.star<-Y.star.cand
ratio[3]<-ratio[3]+1
} else {
cand.y[i,]<-ini.y
cand.fail<-ini.fail
Y.star.cand<-Y.star
}

Y.all<-rbind(Y,Y.star)

##### Generation and MH for delta

if(hyper==TRUE) {

A.prev<-delta*hyper.multiplier
B<-hyper.multiplier
cand.delta[i]<-rgamma(1,shape=A.prev,rate=B)
cand.inv.delta<-1/cand.delta[i]       
A.curr<-cand.delta[i]*hyper.multiplier

# Laplace approximation 

if(sum(include)==1){glm.delta<-glm(Y.star~1,family=binomial)
betaML.delta<-glm.delta$coefficients
piML.delta<-exp(betaML.delta)/(1+exp(betaML.delta))
H.delta<-diag(c(trials.star*piML.delta*(1-piML.delta)))
Det.prev<-delta/sum(H.delta)
Det.curr<-cand.delta[i]/sum(H.delta)
}
if(sum(include)!=1){glm.delta<-glm(Y.star~X.star[,include[-1]],family=binomial)
betaML.delta<-glm.delta$coefficients
piML.delta<-exp(X.star_1[,include]%*%betaML.delta)/(1+exp(X.star_1[,include]%*%betaML.delta))
H.delta<-diag(c(trials.star*piML.delta*(1-piML.delta)))
Det.prev<-det(delta*solve(t(X.star_1[,include])%*%H.delta%*%X.star_1[,include]))
Det.curr<-det(cand.delta[i]*solve(t(X.star_1[,include])%*%H.delta%*%X.star_1[,include]))
}


u.delta<-log(runif(1))

mh.ratio.delta<-log.posterior.delta(cand.delta[i],ini.y,trials.star,cand.p.gamma,piML.delta,Det.curr,A.prev,B)-
                log.posterior.delta(delta,ini.y,trials.star,cand.p.gamma,piML.delta,Det.prev,A.curr,B)
                 
if(u.delta <= mh.ratio.delta) {
inv.delta<-cand.inv.delta
delta<-cand.delta[i]
ratio[4]<-ratio[4]+1
} else {
cand.inv.delta<-inv.delta
cand.delta[i]<-delta
}

weights.star<-rep(inv.delta,n.star)
weights.all<-c(rep(1,n),rep(inv.delta,n.star))

}

##### Generation of gamma's

if(model=='search'){

cand.b.all<-cand.b.gamma[i,]
cand.b.all[exclude]<-cand.b.gamma.exl

for(j in 2:p) {

inc.state<-rbind(include,include)

inc.state[1,j]<-1
inc.state[2,j]<-0

exc.state<-inc.state
exc.state<-inc.state+1
exc.state[exc.state==2]<-0

cand.b.inc1<-c(cand.b.all*inc.state[1,])
cand.b.inc0<-c(cand.b.all*inc.state[2,])

cand.b.exc1<-c(cand.b.all*exc.state[1,])
cand.b.exc0<-c(cand.b.all*exc.state[2,])

Pi_y1<-exp(X_1%*%cand.b.inc1)
Pi_y1[Pi_y1==Inf]<-exp(709.78)
Pi_y1<-Pi_y1/(1+Pi_y1)
Pi_y0<-exp(X_1%*%cand.b.inc0)
Pi_y0[Pi_y0==Inf]<-exp(709.78)
Pi_y0<-Pi_y0/(1+Pi_y0)
Lik_y1<-sum(dbinom(success,trials,Pi_y1,log=TRUE))
Lik_y0<-sum(dbinom(success,trials,Pi_y0,log=TRUE))

Lik_y.star1<-inv.delta*sum(dbinom(cand.y[i,],trials.star,Pi_y1,log=TRUE))
Lik_y.star0<-inv.delta*sum(dbinom(cand.y[i,],trials.star,Pi_y0,log=TRUE))

inc.ind1<-inc.state[1,]==1
inc.ind0<-inc.state[2,]==1

glm.inc1<-glm(Y.star.cand~X.star[,inc.ind1[-1]],family=binomial)
betaML.inc1<-glm.inc1$coefficients
piML.inc1<-exp(X.star_1[,inc.ind1]%*%betaML.inc1)
piML.inc1[piML.inc1==Inf]<-exp(709.78)
piML.inc1<-piML.inc1/(1+piML.inc1)
piML.inc1<-(piML.inc1^inv.delta)/(piML.inc1^inv.delta+(1-piML.inc1)^inv.delta)
H.inc1<-diag(c(trials.star*piML.inc1*(1-piML.inc1)))
Det.inc1<-det(delta*solve(t(X.star_1[,inc.ind1])%*%H.inc1%*%X.star_1[,inc.ind1]))

if(sum(inc.ind0)==1){
glm.inc0<-glm(Y.star.cand~1,family=binomial)
betaML.inc0<-glm.inc0$coefficients
piML.inc0<-exp(betaML.inc0)
piML.inc0[piML.inc0==Inf]<-exp(709.78)
piML.inc0<-piML.inc0/(1+piML.inc0)
piML.inc0<-(piML.inc0^inv.delta)/(piML.inc0^inv.delta+(1-piML.inc0)^inv.delta)
H.inc0<-diag(c(trials.star*piML.inc0*(1-piML.inc0)))
Det.inc0<-delta/sum(H.inc0)
}

if(sum(inc.ind0)!=1){
glm.inc0<-glm(Y.star.cand~X.star[,inc.ind0[-1]],family=binomial)
betaML.inc0<-glm.inc0$coefficients
piML.inc0<-exp(X.star_1[,inc.ind0]%*%betaML.inc0)
piML.inc0[piML.inc0==Inf]<-exp(709.78)
piML.inc0<-piML.inc0/(1+piML.inc0)
piML.inc0<-(piML.inc0^inv.delta)/(piML.inc0^inv.delta+(1-piML.inc0)^inv.delta)
H.inc0<-diag(c(trials.star*piML.inc0*(1-piML.inc0)))
Det.inc0<-det(delta*solve(t(X.star_1[,inc.ind0])%*%H.inc0%*%X.star_1[,inc.ind0]))
}

Marg1<-sum(inc.state[1,])/2*log(2*pi)+0.5*log(Det.inc1)+inv.delta*sum(dbinom(cand.y[i,],trials.star,piML.inc1,log=TRUE))
Marg0<-sum(inc.state[2,])/2*log(2*pi)+0.5*log(Det.inc0)+inv.delta*sum(dbinom(cand.y[i,],trials.star,piML.inc0,log=TRUE))

exc.ind1<-exc.state[1,]==1
exc.ind0<-exc.state[2,]==1

if(sum(exc.ind1)==0){Pseudo1<-0}
if(sum(exc.ind1)==1){Pseudo1<-dnorm(cand.b.exc1[exc.ind1],BetaML[exc.ind1],sqrt(CovML.ind[exc.ind1,exc.ind1]),log=TRUE)}
if(sum(exc.ind1)>1){Pseudo1<-dmvnorm(cand.b.exc1[exc.ind1],BetaML[exc.ind1],CovML.ind[exc.ind1,exc.ind1],log=TRUE)}

if(sum(exc.ind0)==1){Pseudo0<-dnorm(cand.b.exc0[exc.ind0],BetaML[exc.ind0],sqrt(CovML.ind[exc.ind0,exc.ind0]),log=TRUE)}
if(sum(exc.ind0)>1){Pseudo0<-dmvnorm(cand.b.exc0[exc.ind0],BetaML[exc.ind0],CovML.ind[exc.ind0,exc.ind0],log=TRUE)}

d<-sum(inc.state[1,-1])-1
if(model.prob=='uniform') {Oj<-exp(Lik_y1+Lik_y.star1+Pseudo1-Marg1-Lik_y0-Lik_y.star0-Pseudo0+Marg0)}
if(model.prob=='beta') {Oj<-exp(Lik_y1+Lik_y.star1+Pseudo1-Marg1-Lik_y0-Lik_y.star0-Pseudo0+Marg0+log(d+1)-log(p-1-d))}
if(Lik_y1 == - Inf | Lik_y.star1 == - Inf) {Oj<-exp(-exp(709.78)-Lik_y0-Lik_y.star0-Pseudo0+Marg0)}
if(Oj==Inf) {Oj<-exp(709.78)}
prop.gamma.j<-Oj/(1+Oj)

ini.gamma[j]<-rbinom(1,1,prop.gamma.j)
include<-ini.gamma==1
exclude<-ini.gamma==0

} # end of inner for loop for the gamma's

cand.gamma[i+1,]<-ini.gamma[2:p]

ini.b.gamma<-cand.b.all*ini.gamma
ini.p.gamma<-exp(X_1%*%ini.b.gamma)/(1+exp(X_1%*%ini.b.gamma))

act.model[i]<-sum(cand.gamma[i+1,]*2^(0:(p-2)))
act.predictors[i]<-c(paste0("X", c(1:(p-1))[include[-1]],collapse='+'))
} # end of if model=='search'

} # end of iterations

} # end of binomial

############################  POISSON

if(family=='poisson') {

Y<-y
n<-length(y)

 
if(ini.b0=='none') {
glm.ini.b0<-glm(Y~1,family=poisson)
ini.b0<-glm.ini.b0$coefficients
}

ini.lambda0<-exp(ini.b0)

if(class(ini.y)=='numeric') {
n.star<-length(ini.y)
}
if(ini.y=='none'){
ini.y<-rpois(n.star,ini.lambda0)
}
if(n.star<n) stop('Size of imaginary data must be at least equal to real sample size')

delta<-n.star
inv.delta<-1/n.star
    
Y.star<-ini.y 
Y.all<-c(Y,Y.star)

X_1<-cbind(1,X)
if(is.null(colnames(X))) {colnames(X)<-paste(rep('X',ncol(X)),c(1:ncol(X)),sep='')}
colnames(X_1)<-c('intercept',colnames(X))
X.star<-X[rep(1:n,length.out=n.star),]
X.star_1<-cbind(1,X.star)
X.all<-rbind(X,X.star)

weights.star<-rep(inv.delta,n.star)
weights.all<-c(rep(1,n),rep(inv.delta,n.star))

p<-dim(X_1)[2]

if(model=='search'){
ini.gamma<-c(1,rbinom((p-1),1,0.5))     
include<-ini.gamma==1
exclude<-ini.gamma==0
act.model<-rep(NA,iter)
act.predictors<-rep(NA,iter)
}

if(model=='single') {
if(ini.b=='none')   {
glm.ini.b<-glm(Y.all~X.all[,include[-1]],family=poisson)
ini.b<-glm.ini.b$coefficients
}
ini.b.gamma<-rep(0,p)
ini.b.gamma[include]<-ini.b
ini.lambda.gamma<-exp(X_1%*%ini.b.gamma)
}


if(model=='search') {
if(ini.b=='none')   {
glm.ini.b<-glm(Y.all~X.all,family=poisson)
ini.b<-glm.ini.b$coefficients
}
ini.b.gamma<-ini.b*ini.gamma
ini.lambda.gamma<-exp(X_1%*%ini.b.gamma)
}

log.posterior.b<-function(beta,Lambda,real,star,betaML,SigmaML) {

Lambda_star<-rep(Lambda,length.out=n.star)
if(sum(include)==1){
prob<-sum(dpois(real,Lambda,log=TRUE))+
      inv.delta*sum(dpois(star,Lambda,log=TRUE))-
      dnorm(beta,betaML,SigmaML,log=TRUE)
} else {
prob<-sum(dpois(real,Lambda,log=TRUE))+
      inv.delta*sum(dpois(star,Lambda,log=TRUE))-
      dmvnorm(beta,betaML,SigmaML,log=TRUE)
}
 if(prob>= exp(710))  {prob<- exp(709.78)}
 if(prob<= -exp(710)) {prob<- -exp(709.78)}
prob

}

log.posterior.b0<-function(beta0,star,Lambda0,beta0ML,Sigma0ML) {

Lambda0<-rep(Lambda0,length.out=n.star)
prob<-sum(dpois(star,Lambda0,log=TRUE))-
      dnorm(beta0,beta0ML,Sigma0ML,log=TRUE)
 if(prob>= exp(710))  {prob<- exp(709.78)}
 if(prob<= -exp(710)) {prob<- -exp(709.78)}
prob

}

log.posterior.y<-function(y_star,Lambda_star,Lambda0_star,LambdaML_star,DetML_star,LambdaProp_star) {

Lambda_star<-rep(Lambda_star,length.out=n.star)
Lambda0_star<-rep(Lambda0_star,length.out=n.star)

prob<-inv.delta*sum(dpois(y_star,Lambda_star,log=TRUE))+
      sum(dpois(y_star,Lambda0_star,log=TRUE))-
      0.5*log(DetML_star)-inv.delta*sum(dpois(y_star,LambdaML_star,log=TRUE))-
      sum(dpois(y_star,LambdaProp_star,log=TRUE))
 if(prob>= exp(710))  {prob<- exp(709.78)}
 if(prob<= -exp(710)) {prob<- -exp(709.78)}
prob

}



if(hyper==TRUE) {

log.posterior.delta<-function(cand.delta,y_star,Lambda_star,LambdaML_star,DetML_star,propA,propB) {

Lambda_star<-rep(Lambda_star,length.out=n.star)
cand.inv.delta<-1/cand.delta

if(hyper.type=='gamma')     {prior.d<-dgamma(cand.delta,shape=n.star*0.001,rate=0.001,log=TRUE)}
if(hyper.type=='ZS')        {prior.d<- -3/2*log(cand.delta)-n/(2*cand.delta)}
if(hyper.type=='hyper-g')   {prior.d<- -(3/2)*log(1+cand.delta)}
if(hyper.type=='hyper-g/n') {prior.d<- -(3/2)*log(1+cand.delta/n.star)}
if(hyper.type=='maruyama')  {prior.d<- ((n-sum(include)+1)/2-7/4)*log(cand.delta)-((n-sum(include))/2)*log(1+cand.delta)}

prob<-prior.d+
      cand.inv.delta*sum(dpois(y_star,Lambda_star,log=TRUE))-
      0.5*log(DetML_star)-cand.inv.delta*sum(dpois(y_star,LambdaML_star,log=TRUE))-
      dgamma(cand.delta,shape=propA,rate=propB,log=TRUE)
 if(prob>= exp(710))  {prob<- exp(709.78)}
 if(prob<= -exp(710)) {prob<- -exp(709.78)}
prob

}

}


cand.b<-rep(NA,p)
cand.b.gamma<-matrix(0,iter,p)
cand.y<-matrix(NA,iter,n.star)
if(model=='search'){
cand.gamma<-matrix(NA,iter+1,p-1)
cand.gamma[1,]<-ini.gamma[-1]
}

if(hyper==TRUE) {
cand.delta<-rep(NA,iter)
ratio<-rep(0,4)
names(ratio)<-c('betas','beta0','y*','delta')
} else {
ratio<-rep(0,3)
names(ratio)<-c('betas','beta0','y*')
}

if(model=='single'){
glm.b.gamma<-glm(Y.all~X.all[,include[-1]],weights=weights.all,family=poisson)
BetaML<-glm.b.gamma$coefficients
CovML<-vcov(glm.b.gamma)
}


if(model=='search'){
glm.b.gamma<-glm(Y~X,family=poisson)
BetaML<-glm.b.gamma$coefficients
CovML<-vcov(glm.b.gamma)
CovML.ind<-diag(diag(CovML))
}

for(i in 1:iter) {

##### Generation and MH for beta's gamma=1


u.b<-log(runif(1))

if(model=='single'){
cand.b<-rmvnorm(1,BetaML,CovML)
cand.b.gamma[i,include]<-cand.b
cand.lambda.gamma<-exp(X_1%*%cand.b.gamma[i,])

mh.ratio.b<-log.posterior.b(cand.b.gamma[i,include],cand.lambda.gamma,Y,ini.y,BetaML,CovML)-
            log.posterior.b(ini.b.gamma[include],ini.lambda.gamma,Y,ini.y,BetaML,CovML)

if(u.b <= mh.ratio.b) {                   
ini.b.gamma<-cand.b.gamma[i,]
ini.lambda.gamma<-cand.lambda.gamma
ratio[1]<-ratio[1]+1
} else {
cand.b.gamma[i,]<-ini.b.gamma
cand.lambda.gamma<-ini.lambda.gamma
cand.b<-ini.b                    
}
}


if(model=='search'){

if(sum(include)==1){
glm.b.gamma<-glm(Y.all~1,weights=weights.all,family=poisson)
BetaML.INC<-glm.b.gamma$coefficients
CovML.INC<-sqrt(vcov(glm.b.gamma))
cand.b<-rnorm(1,BetaML.INC,CovML.INC)
}
if(sum(include)!=1){
glm.b.gamma<-glm(Y.all~X.all[,include[-1]],weights=weights.all,family=poisson)
BetaML.INC<-glm.b.gamma$coefficients
CovML.INC<-vcov(glm.b.gamma)
cand.b<-rmvnorm(1,BetaML.INC,CovML.INC)
}

cand.b.gamma[i,include]<-cand.b
cand.lambda.gamma<-exp(X_1%*%cand.b.gamma[i,])

mh.ratio.b<-log.posterior.b(cand.b,cand.lambda.gamma,Y,ini.y,BetaML.INC,CovML.INC)-
            log.posterior.b(ini.b.gamma[include],ini.lambda.gamma,Y,ini.y,BetaML.INC,CovML.INC)

if(u.b <= mh.ratio.b) {                   
ini.b.gamma<-cand.b.gamma[i,]
ini.lambda.gamma<-cand.lambda.gamma
ratio[1]<-ratio[1]+1
} else {                 
cand.b.gamma[i,]<-ini.b.gamma
cand.lambda.gamma<-ini.lambda.gamma
cand.b<-ini.b[include]                    
}
}

##### Generation beta's gamma=0

if(model=='search') {
cand.b.gamma.exl<-rmvnorm(1,BetaML,CovML.ind)
cand.b.gamma.exl<-cand.b.gamma.exl[exclude]
}

##### Generation and MH for beta0 - NULL model

glm.b0<-glm(Y.star~1,family=poisson)
BetaML0<-glm.b0$coefficients
CovML0<-sqrt(vcov(glm.b0))

cand.b0<-rnorm(1,BetaML0,CovML0)
cand.lambda0<-exp(cand.b0)
u.b0<-log(runif(1))

mh.ratio.b0<-log.posterior.b0(cand.b0,ini.y,cand.lambda0,BetaML0,CovML0)-
             log.posterior.b0(ini.b0,ini.y,ini.lambda0,BetaML0,CovML0)

if(u.b0 <= mh.ratio.b0)  {
ini.b0<-cand.b0
ini.lambda0<-cand.lambda0
ratio[2]<-ratio[2]+1
} else {
cand.b0<-ini.b0
cand.lambda0<-ini.lambda0
}

##### Generation and MH for y*

lambda.STAR<-(cand.lambda0*cand.lambda.gamma^inv.delta)
lambda.STAR<-rep(lambda.STAR,length.out=n.star)

cand.y[i,]<-rpois(n.star,lambda.STAR)     
Y.star.cand<-cand.y[i,]


# Laplace approximation MLE at iteration t-1

if(sum(include)==1){glm.prev<-glm(Y.star~1,family=poisson)
betaML.prev<-glm.prev$coefficients
lambdaML.prev<-exp(betaML.prev)
H.prev<-diag(lambdaML.prev,n.star)
Det.prev<-delta/sum(H.prev)
}
if(sum(include)!=1){glm.prev<-glm(Y.star~X.star[,include[-1]],family=poisson)
betaML.prev<-glm.prev$coefficients
lambdaML.prev<-exp(X.star_1[,include]%*%betaML.prev)
H.prev<-diag(c(lambdaML.prev))
Det.prev<-det(delta*solve(t(X.star_1[,include])%*%H.prev%*%X.star_1[,include]))
}

# Laplace approximation MLE at iteration t

if(sum(include)==1){glm.curr<-glm(Y.star.cand~1,family=poisson)
betaML.curr<-glm.curr$coefficients
lambdaML.curr<-exp(betaML.curr)
H.curr<-diag(lambdaML.curr,n.star)
Det.curr<-delta/sum(H.curr)
}
if(sum(include)!=1){glm.curr<-glm(Y.star.cand~X.star[,include[-1]],family=poisson)
betaML.curr<-glm.curr$coefficients
lambdaML.curr<-exp(X.star_1[,include]%*%betaML.curr)
H.curr<-diag(c(lambdaML.curr))
Det.curr<-det(delta*solve(t(X.star_1[,include])%*%H.curr%*%X.star_1[,include]))
}

u.y<-log(runif(1))

mh.ratio.y<-log.posterior.y(cand.y[i,],cand.lambda.gamma,cand.lambda0,lambdaML.curr,Det.curr,lambda.STAR)-
            log.posterior.y(ini.y,cand.lambda.gamma,cand.lambda0,lambdaML.prev,Det.prev,lambda.STAR)
                 
if(u.y <= mh.ratio.y) {
ini.y<-cand.y[i,]
Y.star<-Y.star.cand
ratio[3]<-ratio[3]+1
} else {
cand.y[i,]<-ini.y
Y.star.cand<-Y.star
}

Y.all<-c(Y,Y.star)

##### Generation and MH for delta

if(hyper==TRUE) {

A.prev<-delta*hyper.multiplier
B<-hyper.multiplier
cand.delta[i]<-rgamma(1,shape=A.prev,rate=B)
cand.inv.delta<-1/cand.delta[i]       
A.curr<-cand.delta[i]*hyper.multiplier     

# Laplace approximation 

if(sum(include)==1){glm.delta<-glm(Y.star~1,family=poisson)
betaML.delta<-glm.delta$coefficients
lambdaML.delta<-exp(betaML.delta)
H.delta<-diag(lambdaML.delta,n.star)
Det.prev<-delta/sum(H.delta)
Det.curr<-cand.delta[i]/sum(H.delta)
}
if(sum(include)!=1){glm.delta<-glm(Y.star~X.star[,include[-1]],family=poisson)
betaML.delta<-glm.delta$coefficients
lambdaML.delta<-exp(X.star_1[,include]%*%betaML.delta)
H.delta<-diag(c(lambdaML.delta))
Det.prev<-det(delta*solve(t(X.star_1[,include])%*%H.delta%*%X.star_1[,include]))
Det.curr<-det(cand.delta[i]*solve(t(X.star_1[,include])%*%H.delta%*%X.star_1[,include]))
}

u.delta<-log(runif(1))

mh.ratio.delta<-log.posterior.delta(cand.delta[i],ini.y,cand.lambda.gamma,lambdaML.delta,Det.curr,A.prev,B)-
                log.posterior.delta(delta,ini.y,cand.lambda.gamma,lambdaML.delta,Det.prev,A.curr,B)
                 
if(u.delta <= mh.ratio.delta) {
inv.delta<-cand.inv.delta
delta<-cand.delta[i]
ratio[4]<-ratio[4]+1
} else {
cand.inv.delta<-inv.delta
cand.delta[i]<-delta
}

weights.star<-rep(inv.delta,n.star)
weights.all<-c(rep(1,n),rep(inv.delta,n.star))

}

##### Generation of gamma's

if(model=='search'){

cand.b.all<-cand.b.gamma[i,]
cand.b.all[exclude]<-cand.b.gamma.exl

for(j in 2:p) {

inc.state<-rbind(include,include)

inc.state[1,j]<-1
inc.state[2,j]<-0

exc.state<-inc.state
exc.state<-inc.state+1
exc.state[exc.state==2]<-0

cand.b.inc1<-c(cand.b.all*inc.state[1,])
cand.b.inc0<-c(cand.b.all*inc.state[2,])

cand.b.exc1<-c(cand.b.all*exc.state[1,])
cand.b.exc0<-c(cand.b.all*exc.state[2,])

Lambda_y1<-exp(X_1%*%cand.b.inc1)
Lambda_y0<-exp(X_1%*%cand.b.inc0)
Lik_y1<-sum(dpois(Y,Lambda_y1,log=TRUE))
Lik_y0<-sum(dpois(Y,Lambda_y0,log=TRUE))

Lik_y.star1<-inv.delta*sum(dpois(cand.y[i,],Lambda_y1,log=TRUE))
Lik_y.star0<-inv.delta*sum(dpois(cand.y[i,],Lambda_y0,log=TRUE))

inc.ind1<-inc.state[1,]==1
inc.ind0<-inc.state[2,]==1

glm.inc1<-glm(Y.star.cand~X.star[,inc.ind1[-1]],family=poisson)
betaML.inc1<-glm.inc1$coefficients
LambdaML.inc1<-exp(X.star_1[,inc.ind1]%*%betaML.inc1)
H.inc1<-diag(c(LambdaML.inc1))
Det.inc1<-det(delta*solve(t(X.star_1[,inc.ind1])%*%H.inc1%*%X.star_1[,inc.ind1]))

if(sum(inc.ind0)==1){
glm.inc0<-glm(Y.star.cand~1,family=poisson)
betaML.inc0<-glm.inc0$coefficients
LambdaML.inc0<-exp(betaML.inc0)
H.inc0<-diag(LambdaML.inc0,n.star)
Det.inc0<-delta/sum(H.inc0)
}

if(sum(inc.ind0)!=1){
glm.inc0<-glm(Y.star.cand~X.star[,inc.ind0[-1]],family=poisson)
betaML.inc0<-glm.inc0$coefficients
LambdaML.inc0<-exp(X.star_1[,inc.ind0]%*%betaML.inc0)
H.inc0<-diag(c(LambdaML.inc0))
Det.inc0<-det(delta*solve(t(X.star_1[,inc.ind0])%*%H.inc0%*%X.star_1[,inc.ind0]))
}

Marg1<-sum(inc.state[1,])/2*log(2*pi)+0.5*log(Det.inc1)+inv.delta*sum(dpois(cand.y[i,],LambdaML.inc1,log=TRUE))
Marg0<-sum(inc.state[2,])/2*log(2*pi)+0.5*log(Det.inc0)+inv.delta*sum(dpois(cand.y[i,],LambdaML.inc0,log=TRUE))

exc.ind1<-exc.state[1,]==1
exc.ind0<-exc.state[2,]==1

if(sum(exc.ind1)==0){Pseudo1<-0}
if(sum(exc.ind1)==1){Pseudo1<-dnorm(cand.b.exc1[exc.ind1],BetaML[exc.ind1],sqrt(CovML.ind[exc.ind1,exc.ind1]),log=TRUE)}
if(sum(exc.ind1)>1){Pseudo1<-dmvnorm(cand.b.exc1[exc.ind1],BetaML[exc.ind1],CovML.ind[exc.ind1,exc.ind1],log=TRUE)}

if(sum(exc.ind0)==1){Pseudo0<-dnorm(cand.b.exc0[exc.ind0],BetaML[exc.ind0],sqrt(CovML.ind[exc.ind0,exc.ind0]),log=TRUE)}
if(sum(exc.ind0)>1){Pseudo0<-dmvnorm(cand.b.exc0[exc.ind0],BetaML[exc.ind0],CovML.ind[exc.ind0,exc.ind0],log=TRUE)}

d<-sum(inc.state[1,-1])-1
if(model.prob=='uniform') {Oj<-exp(Lik_y1+Lik_y.star1+Pseudo1-Marg1-Lik_y0-Lik_y.star0-Pseudo0+Marg0)}
if(model.prob=='beta') {Oj<-exp(Lik_y1+Lik_y.star1+Pseudo1-Marg1-Lik_y0-Lik_y.star0-Pseudo0+Marg0+log(d+1)-log(p-1-d))}
if(Oj==Inf) {Oj<-exp(709.78)}
prop.gamma.j<-Oj/(1+Oj)

ini.gamma[j]<-rbinom(1,1,prop.gamma.j)
include<-ini.gamma==1
exclude<-ini.gamma==0

} # end of inner for loop for the gamma's

cand.gamma[i+1,]<-ini.gamma[2:p]

ini.b.gamma<-cand.b.all*ini.gamma
ini.lambda.gamma<-exp(X_1%*%ini.b.gamma)

act.model[i]<-sum(cand.gamma[i+1,]*2^(0:(p-2)))
act.predictors[i]<-c(paste0("X", c(1:(p-1))[include[-1]],collapse='+'))
} # end of if model=='search'

} # end of iterations

} # end of poisson

End<-Sys.time()
runtime<-difftime(End,Begin)
if(model=='single'){
final.b<-cand.b.gamma[(discard+1):iter,include]
colnames(final.b)<-colnames(X_1[,include])
if(hyper==FALSE) {results<-list(betas=final.b,acceptance.ratios=ratio/iter,runtime=runtime)}
if(hyper==TRUE) {
final.delta<-cand.delta[(discard+1):iter]
results<-list(betas=final.b,delta=final.delta,acceptance.ratios=ratio/iter,runtime=runtime)
}
}
if(model=='search'){
final.b<-cand.b.gamma[(discard+1):iter,]
colnames(final.b)<-colnames(X_1)
final.gamma<-cand.gamma[(discard+1):iter,]
colnames(final.gamma)<-colnames(X)
PIPs<-round(apply(final.gamma,2,sum)/(iter-discard),3)
act.model<-act.model[(discard+1):iter]
act.predictors<-act.predictors[(discard+1):iter]
crosstab<-table(act.model)
rows<-c()
for(i in 1:length(crosstab)){
rows[i]<-which(act.model==names(crosstab)[i])[1]
}
summary<-data.frame(crosstab,act.predictors[rows])
summary<-summary[order(summary[,2],decreasing=TRUE),]
rownames(summary)<-1:dim(summary)[1]
summary<-data.frame(summary[,1:2],summary[,2]/(iter-discard),summary[,3])
names(summary)<-c('Model','Frequency','PostProb','Variables')
if(hyper==FALSE) {results<-list(betas=final.b,gammas=final.gamma,search=summary,PIPs=PIPs,acceptance.ratios=ratio/iter,runtime=runtime)}
if(hyper==TRUE) {
final.delta<-cand.delta[(discard+1):iter]
results<-list(betas=final.b,gammas=final.gamma,delta=final.delta,search=summary,PIPs=PIPs,acceptance.ratios=ratio/iter,runtime=runtime)
}
}
results
}



#################
# Sample R code and Simulations to illustrate estimation and inference
#	for the sample average treatment effect in trials with and without pair-matching.
#	Demonstrates the unadjusted estimator, 
#	TMLE with logistic regression for outcome regression E_0(Y|A,W)=Qbar_0(A,W), 
#	and TMLE with SuperLearner for Qbar_0(A,W)
#
# Programmer: Laura Balzer (lbbalzer@hsph.harvard.edu)
#	Please email with questions, concerns or requests
#
# R version 3.2.1
#
# Last update: Dec 8, 2015
####################

#-------------------------------
# simulate.data.and.run: function to generate the simulated data
#	and run the estimators
#-------------------------------
simulate.data.and.run<- function(){
	
	# directly simulate the full data (covariates and counterfactual outcomes)
	X.full<- generateData(n)
	
	# Sample average counterfactual outcome under A=a 
	R1<- mean(X.full$Y.1)
	R0<- mean(X.full$Y.0)

	# SATE is the sample average difference in the counterfactuals
	SATE= R1-R0
	
	if(PAIRED){
		# if pair-matching trial, match units and randomize the treatment
		X.all <- doPairMatching(matchData=X.full[, matchOn], fullData=X.full)	

	} else{	
		# Otherwise, assign the treatment - guarantee that n/2 are treated 
		A<- rbinom(n/2, 1, 0.5)
		A.2<- ifelse(A==1, 0, 1)
		A <- sample( c(A,A.2))
		X.all <- cbind( X.full, paired=rep(NA, n), A)	
	}
	
	# we observe the counterfactual outcome corresponding to the observed exp
	Y<- ifelse(X.all$A, X.all$Y.1, X.all$Y.0)
	X.all<- cbind(X.all, Y)

	#-----------------
	# Estimation and inference
	#-----------------
	
	unadj<- doTMLE(SATE, data=X.all, Qadj='U', family='binomial')
	
	adj.AW1<- doTMLE(SATE, data=X.all, Qadj='W.1', Qform=as.formula(Y~A*W.1), family='binomial')
	
	Qadj<- c('W.1','W.2', 'W.3', 'W.4', 'W.5')
	
	adj.SL<- doTMLE(SATE, data=X.all, Qadj=Qadj, family='binomial', Do.SL=T)

	RETURN<- list(unadj=unlist(unadj),  adj.AW1=unlist(adj.AW1),
		adj.SL =unlist(adj.SL) )

	RETURN
}

#---------------------------
# generateData: function to generate the full data
#	including baseline covariates and the counterfactual outcomes
#---------------------------
generateData<- function(n){
	
	U.Y<- generateU.Y(n)
	W<- generateW(n)
	Y.0<- generateY(W=W, A=0, U.Y=U.Y)
	if(EFFECT){
		Y.1<- generateY(W=W, A=1, U.Y=U.Y)
	} else{
		Y.1 <- Y.0
	}
	data.frame(W, Y.0,Y.1)
}

#-------------------
# additional functions to generate the simulated data
#--------------------

# generate unmeasured U.Y
generateU.Y<- function(n){
	rnorm(n, 0, SD)  
}

# generate the baseline covariates W
generateW<- function(n) {
	
	Sigma<- matrix(CORR.W*SD*SD, nrow=3, ncol=3)
	diag(Sigma)<- SD^2
	
	W<- cbind(rnorm(n,0,1), rnorm(n,0,1), mvrnorm(n, rep(0,3), Sigma))
	
	data.frame(U=1, W=W )
}	

# generate the outcome Y
generateY<- function(W, A, U.Y) {
	.2*plogis(1*A + .75*W$W.1 + .75*W$W.2 + 1.25*W$W.3 + U.Y + .75*W$W.1*A - .5*W$W.2*A - A*U.Y ) 
}


#-----------------------
# get.PATE: function to calculate the true value of the PATE
#	over a population of 500,000 units
#-----------------------

get.PATE<- function(pop= 500000){
	
	X.full<- generateData(pop)
	
	# average counterfactual outcome under A=a for the population
	R1<- mean( X.full$Y.1)
	R0<- mean( X.full$Y.0)

	RD= R1-R0
	
	c(R1, R0, RD)
}

#--------------------	
# doPairMatching - function to pair-match units (fullData) 
#	based on the matching covariates (matchData)
#	and assign the treatment within the resulting matched pairs
# Requires nbpMatching package
#-----------------------

doPairMatching<- function(matchData, fullData){
	
	dist<- distancematrix(gendistance(data.frame(matchData)))
	matches<- nonbimatch(dist)
	
	# matches contains ids for the pair as well as the distance measure
	grpA<- as.numeric(matches$halves[,'Group1.Row'])
	grpB<- as.numeric(matches$halves[,'Group2.Row'])	
	
	npairs<- length(grpA)
	X1<- data.frame(fullData[grpA, ], pair=1:npairs, A= rbinom(npairs, 1, .5))
	X2<- data.frame(fullData[grpB, ], pair=1:npairs, A= ifelse(X1$A==1, 0, 1 ))
	
	Xpaired<- NULL
	for(i in 1:npairs){		
		Xpaired<- rbind(Xpaired, X1[i,], X2[i,]) 
	}
	Xpaired
}

#--------------------------------
# doTMLE: function to run full TMLE and get inference
#	input: SATE (sample ATE for that study), data, Qadj (candidate adjustment variables for Qbar_0(A,W))
#		Qform (the form of the outcome regression), family (binomial for logistic regression), 
#		Do.SL (whether or not do SuperLearner)
#	output: estimation and inference for the population and sample effect
#
#	Requires the SuperLearner package
#
#	For further information about coding TMLE and calling SuperLearner, 
#	please see http://www.ucbbiostat.com/
#-----------------------------------------------

doTMLE<- function(SATE, data, Qadj, Qform=as.formula(Y~.), family='binomial', Do.SL=F){
	
	if(!Do.SL){ # if not doing SuperLearner
		
		X1 = X0 = X= data[,c(Qadj, 'A','Y')]
		X1$A<-1; X0$A<- 0	
		glm.out<- suppressWarnings(  glm(Qform, family=family, data=X) )
	
		# get predicted outcomes under obs exp, txt and control 
		QbarAW<- suppressWarnings( predict(glm.out, newdata=X, type="response"))
		Qbar1W<- suppressWarnings( predict(glm.out, newdata=X1, type='response'))
		Qbar0W<- suppressWarnings( predict(glm.out, newdata=X0, type='response') )
	
	} else{  	# do super learner
		
		X1 = X0 = X= data[,c(Qadj, 'A')]
		X1$A<-1; X0$A<- 0	
		newX<- rbind(X,X1, X0)
		
		# call SuperLearner
		if(PAIRED){	
			# for the cross-validation step, we need to respect the unit of (conditional) independence 
			Qinit<-SuperLearner(Y=data$Y, X=X, newX=newX, SL.library= QSL.LIBRARY, family="binomial", 
				cvControl=list(V=n/2), id=data$pairs )
		} else{
			Qinit<-SuperLearner(Y=data$Y, X=X, newX=newX, SL.library= QSL.LIBRARY, family="binomial", 
				cvControl=list(V=n/2) )
		}
				
		QbarAW<-Qinit$SL.predict[1:n]
		Qbar1W<-Qinit$SL.predict[(n+1): (2*n)]
		Qbar0W<-Qinit$SL.predict[(2*n+1): (3*n)]
	
	}	

	# We're not estimating the known exposure mechanism, 
	# 	but we could for greater efficiency
	# For further details, email lbbalzer@hsph.harvard.edu
	pscore = rep(0.5,n)
	
	# Calculating the clever covariate
	H.1W<- 1/ (pscore)
	H.0W<-  -1/ (1-pscore)
	H.AW<- rep(NA, n)
	H.AW[data$A==1]<- H.1W[data$A==1]
	H.AW[data$A==0]<- H.0W[data$A==0]
	
	# updating step
	logitUpdate<- suppressWarnings( glm(data$Y ~ -1 +offset(qlogis(QbarAW)) + H.AW, family="binomial"))
	
	# estimated coefficient on the clever covariate
	eps<-logitUpdate$coef
	
	# targeted estimates of the outcome regression
	QbarAW<-plogis( qlogis(QbarAW)+eps*H.AW)
	Qbar0W<-plogis( qlogis(Qbar0W)+eps*H.0W)
	Qbar1W<-plogis( qlogis(Qbar1W)+eps*H.1W)
	
	# risk estimates under txt, under control and risk difference
	R1<- mean(Qbar1W)
	R0<- mean(Qbar0W)
	RD<- mean(Qbar1W- Qbar0W)
	
	#--------------------
	# get inference via the influence curve
	#--------------------
	
	# the relevant components of the influence curve
	DY<- H.AW*(data$Y- QbarAW)
	DW<- Qbar1W - Qbar0W - RD
		
	if(!PAIRED){
		
		var.PATE<- var(DY+DW)/n			
		var.SATE<- var(DY)/n
		df=(n-2)
		
		est.PATE<- get.inference(truth=PATE[3], RD=RD, var=var.PATE, df=df)
		est.SATE<- get.inference(truth=SATE, RD=RD, var=var.SATE, df=df)
	
		RETURN<- data.frame(R1=R1, R0=R0, RD=RD, PATE=est.PATE, SATE=est.SATE)  

	}else{
		
		pairs<- data$pair
		temp<- unique(pairs) 
		n.pairs<- length(temp)
		
		# DbarY= 1/2 sum_{i in pairs} HAW_i*(Y_i -Qbar_i)
		# Serves as the upper bound on the IC	
		DY.paired<- rep(NA, n.pairs)
	
		for(i in 1:n.pairs){		
			DY.paired[i]<- 0.5*sum(DY[ pairs== temp[i]] ) 
		}
		
		var.SATE<- var(DY.paired)/n.pairs
		
		df= (n.pairs -1)
		est.SATE<- get.inference(truth=SATE, RD=RD, var=var.SATE, df=df)
		
		# for estimation of the PATE in an Adaptive Pair-Matched Trial, see van der Laan et al. 2012
		RETURN<-  data.frame(R1=R1, R0=R0, RD=RD, SATE=est.SATE)  

	}

	RETURN
}

#---------------------
# get.inference: 
#	input: true value of target parameter, estimate (RD), variance estimate 
#		and df for t-dist
#	output: variance est, test statistic, indicator of 95% CI contained the truth
#		and indicator that rejected the null at the alpha=0.05 level
#-------------------------

get.inference<- function(truth, RD, var, df){
	
	se<- sqrt(var)	
	cutoff <- qt(0.05/2, df=df, lower.tail=F)
	
	cov<- (RD - cutoff*se) <= truth & truth <= (RD + cutoff*se)
	tstat <- RD/se	
	reject <- abs(tstat) > cutoff

	data.frame(truth=truth, var=var, tstat=tstat, cov=cov, reject=reject)
}

#==========================================================
# The remaining are helper functions to run SuperLearner
#===========================================================

SL.glmAW1int<- function (Y, X, newX, family, obsWeights, ...) {
    fit.glm <- glm(Y ~  A + W.1 + A:W.1, data = X, family = family, weights = obsWeights)
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glmAW1int"
    out <- list(pred = pred, fit = fit)
    return(out)
}
SL.glmAW2int<- function (Y, X, newX, family, obsWeights, ...) {
    fit.glm <- glm(Y ~  A + W.2 + A:W.2, data = X, family = family, weights = obsWeights)
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glmAW2int"
    out <- list(pred = pred, fit = fit)
    return(out)
}
SL.glmAW3int<- function (Y, X, newX, family, obsWeights, ...) {
    fit.glm <- glm(Y ~  A + W.3 + A:W.3, data = X, family = family, weights = obsWeights)
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glmAW3int"
    out <- list(pred = pred, fit = fit)
    return(out)
}
SL.glmAW4int<- function (Y, X, newX, family, obsWeights, ...) {
    fit.glm <- glm(Y ~  A + W.4 + A:W.4, data = X, family = family, weights = obsWeights)
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glmAW4int"
    out <- list(pred = pred, fit = fit)
    return(out)
}
SL.glmAW5int<- function (Y, X, newX, family, obsWeights, ...) {
    fit.glm <- glm(Y ~  A + W.5 + A:W.5, data = X, family = family, weights = obsWeights)
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glmAW5int"
    out <- list(pred = pred, fit = fit)
    return(out)
}

#===============================================================

set.seed(1)

library(MASS)
library(SuperLearner)
library(nbpMatching)

# the following are global variables - specified by the user
n<<- 30

SD<<- 1 #std deviation of baseline covariate
CORR.W<<- .65	# correlation in (W3,W4,W5)

# SuperLearner library for Qbar_0(A,W)
QSL.LIBRARY<<- c('SL.glmAW1int', 'SL.glmAW2int', 'SL.glmAW3int', 'SL.glmAW4int', 'SL.glmAW5int')

PAIRED<<- T
matchOn<<- c('W.1', 'W.4','W.5')

EFFECT<<- T

PATE<<- get.PATE()

out<- simulate.data.and.run()

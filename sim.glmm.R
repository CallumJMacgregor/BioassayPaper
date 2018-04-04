sim.glmm_CJM<-
  function(mer.fit=NULL, design.data=NULL, fixed.eff=NULL, rand.V=NULL,
           distribution=c("gaussian","poisson","binomial","negbinomial"),
           SD=NULL, theta=NULL, drop.effects=NULL)
  {
    
    # if mer.fit is supplied, extract all inputs from it
    
    if(!is.null(mer.fit))
    {
      design.data <- model.frame(mer.fit)
      fixed.eff <- lme4::fixef(mer.fit)
      rand.V <- lme4::VarCorr(mer.fit)
      if(grepl("Negative Binomial", as.character(family(mer.fit))[1])){
        distribution <- "negbinomial"
      } else {
        distribution <- as.character(family(mer.fit))[1]
      }
      
      if(distribution == "binomial") design.data$n <- get("n",envir=slot(slot(mer.fit,"resp"),".xData"))
      SD <- attr(rand.V,"sc")
      if(distribution == "negbinomial") theta <- summary(mer.fit)$AICtab[4] / summary(mer.fit)$AICtab[5]
    }
    
    # add column of 1s. this is the "data value" for the intercept
    
    design.data$"(Intercept)"<-1
    
    
    # add the fixed effects to the data frame, first changing "intercept" to "(Intercept)"
    
    if(is.null(mer.fit))
    {
      names(fixed.eff)[names(fixed.eff)=="intercept"] <- "(Intercept)"
      fix.names<-names(fixed.eff)
      design.data[,paste(fix.names,".fixed",sep="")]<-
        sapply(fix.names,
               function(fe)
               {
                 if(is.factor(design.data[,fe])) fixed.eff[[fe]][as.character(design.data[,fe])]
                 else fixed.eff[[fe]]*design.data[,fe]
               })
    } else
    {
      design.data$combinedeffects.fixed <- model.matrix(mer.fit) %*% fixed.eff
      fix.names <- "combinedeffects"
    }
    eff.names <- paste(fix.names,".fixed",sep="")
    
    
    # if random effects are specified, simulate them
    
    if(!is.null(rand.V))
    {
      
      # if rand.V is a vector, convert it to a list
      
      if(!is.list(rand.V)) rand.V <- as.list(rand.V)
      rand.names<-names(rand.V)
      
      # simulate the random effects
      
      a.list<-
        lapply(rand.names,
               function(rn)
               {
                 design.data[,rn] <- factor(design.data[,rn])
                 nt <- nlevels(design.data[,rn])
                 vcv <- rand.V[[rn]]
                 if(is.null(dim(vcv))) vcv <- structure(vcv, dim=c(1,1), dimnames=list("(Intercept)","(Intercept)"))
                 rownames(vcv)[rownames(vcv)=="intercept"] <- "(Intercept)"
                 colnames(vcv)[colnames(vcv)=="intercept"] <- "(Intercept)"
                 a <- MASS::mvrnorm(nt,rep(0,nrow(vcv)),vcv)
                 rownames(a) <- levels(design.data[,rn])
                 ax <-
                   if(is.null(mer.fit)) {
                     a[design.data[,rn],colnames(vcv)] * design.data[,colnames(vcv)] } else {
                       a[design.data[,rn],colnames(vcv)] * model.matrix(mer.fit)[,colnames(vcv)] }
                 if(is.null(dim(ax))) dim(ax) <- c(nrow(design.data),nrow(vcv))
                 apply(ax,1,sum)
               })
      names(a.list) <- rand.names
      
      
      # add the random effect residuals to the data frame
      
      design.data[,paste(rand.names,".random",sep="")] <- do.call("cbind",a.list[rand.names])
      eff.names <- c(eff.names,paste(rand.names,".random",sep=""))
    }
    
    
    # sum the fixed effects and random effect residuals to give the linear predictor
    # of the GLMM
    
    design.data$linear.predictor<-apply(design.data[,eff.names],1,sum)
    
    # if design.data includes an offset, add it to the linear predictor
    
    if(!is.null(design.data$offset)) {
      design.data$linear.predictor <- design.data$linear.predictor + design.data$offset
    }
    
    
    # simulate the response values from the linear predictor
    
    if(distribution[1]=="gaussian")
      design.data$response<-rnorm(n=1:nrow(design.data),mean=design.data$linear.predictor,sd=SD)
    if(distribution[1]=="poisson")
      design.data$response<-rpois(n=1:nrow(design.data),lambda=exp(design.data$linear.predictor))
    if(distribution[1]=="binomial")
      design.data$response<-rbinom(n=1:nrow(design.data),size=design.data$n,prob=plogis(design.data$linear.predictor))
    if(distribution[1]=="negbinomial")
      design.data$response<-rnbinom(n=1:nrow(design.data),size=theta,mu=exp(design.data$linear.predictor))
    
    
    # drop the intercept, fixed effects and random effect residuals from the data frame
    
    design.data<-design.data[,!(names(design.data) %in% c("(Intercept)","linear.predictor",eff.names))]
    
    
    # return the data frame including the simulated response vector
    
    return(design.data)
  }

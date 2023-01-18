library(lubridate) #convertVars()


################ FUNCTIONS TO CONVERT CHARACTER FACTORS TO NUMBERS, FACTORS OR DATE,
####  WRITTEN BY VIVEK PHILLIP
####  RELIES ON A DATA DICTIONARY, THAT HAS INFORMATION ABOUT WHAT THE FACTORS SHOULD BE
#### 'Dictionary_KOMP.csv' has 4 columns samples below
####  Phenotype,Parameter,Field,Sample
####  HB,Parameters:Total,Hole Pokes,Numeric,37
####  Dependency: 'lubridate' package
ConvertVars <- function(data,dict,pheno, verbose=FALSE) {
  if(verbose){cat("\n### ConvertVars Start ###\n")}
  ### this will select OFA OR ALL, important for common factors
  dict.sub <- dict[which(dict$Phenotype == pheno | dict$Phenotype == "ALL"),]
  data.colnames <- colnames(data)
  for(i in 1:length(data.colnames)) {
    idx <- which(dict.sub$Parameter.fixed == data.colnames[i])
    if (length(idx) > 0) {
      if(dict.sub[idx,3] == "Numeric") {
        data[,i] <- as.numeric(data[,i])
      }
      if(dict.sub[idx,3] == "Factor") {
        data[,i] <- as.factor(data[,i])
      }
      if(dict.sub[idx,3] == "Date_ymd_hms") {
        data[,i] <- ymd_hms(data[,i])
      }
      if(dict.sub[idx,3] == "Date_mdy") {
        data[,i] <- mdy(data[,i])
      }
      if(dict.sub[idx,3] == "Date_mdy_hm") {
        data[,i] <- mdy_hm(data[,i])
      }
    }
  }
  if(verbose){cat("### ConvertVars End ###\n")}
  return(data)
}

ConvertVars2 <- function(data,dict,pheno, verbose=FALSE) {
  if(verbose){cat("\n### ConvertVars2 Start ###\n")}
  ## this will select OFA OR ALL, important for common factors
  dict.sub <- dict[which(dict$Phenotype == pheno | dict$Phenotype == "ALL"),]
  data.colnames <- colnames(data)
  for(i in 1:length(data.colnames)) {
    idx <- which(dict.sub$Parameter.fixed == data.colnames[i])
    if (length(idx) > 0) {
      if(dict.sub[idx,3] == "Numeric") {
        data[,i] <- as.numeric(data[,i])
      }
      if(dict.sub[idx,3] == "Factor") {
        data[,i] <- as.factor(data[,i])
      }
    }
  }
  if(verbose){cat("### ConvertVars2 End ###\n")}
  return(data)
}

#### DON: make date format same.
MakeDateFormatSame <- function(date.list, verbose=FALSE){
  if(verbose){cat("\n### MakeDateFormatSame Start ###\n")}
  date.list <- as.Date(parse_date_time(x = date.list, orders=c("Ymd", "Ymd H:M:S", "Y/m/d", "m/d/Y", "mdy", "mdy H:M")))
  if(verbose){cat("### MakeDateFormatSame End ###\n")}
  return(date.list)
}


################ Remove Outliers
##### will remove, on a per genotypes basis, outliers
##### will only work on genotypes with 5 animals
##### tested on Hematology, there are graphs there for with and without outliers
##### this function removes outliers

## DON: I used this function to exclude all outliers in phenotypes
## DON: ,which lie beyond plus/minus 25qnt*stdv.
## DON: stdv is the stdev of the outlier - this is a parameter that is flexible, normally it is 1.5
RemoveOutliers <- function(dt, list, stdv, verbose=FALSE, na.rm = TRUE, ...) {
  if(verbose){
    cat("\n### RemoveOutliers Start ###\n")
    cat("BEFORE: \n")
    for(i in list){
      print(paste0(i, " :",sum(complete.cases(dt[,i]))))
    }
  }
  col.ids <- which(colnames(dt) %in% list)
  for(i in 1:length(col.ids)) {
    x <- dt[,col.ids[i]]
    qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
    H <- stdv * IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    dt[,col.ids[i]] <- y
  }

  if(verbose){
    cat("AFTER: \n")
    for(i in list){
      print(paste0(i, " :",sum(complete.cases(dt[,i]))))
    }
    cat("### RemoveOutliers End ###\n")
  }
  return(dt)
}

RemoveSampleWithNAPredictor <- function(data, pred.list, verbose=FALSE){
  if(verbose){
    cat("\n### RemoveSampleWithNAPredictor Start ###\n")
    for(i in pred.list){
      print(paste("BEFORE: num of NAs in", i, ":", sum(is.na(data[,i]))))
    }
    cat("\nDim of data:")
    print(dim(data))
  }
  data <- na.omit(data)

  if(verbose){
    for(i in pred.list){
      print(paste("After: num of NAs in", i, ":", sum(is.na(data[,i]))))
    }
    cat("\nDim of data:")
    print(dim(data))
    cat("### RemoveSampleWithNAPredictor End ###\n")
  }
  return(data)
}

##Check variance of GenotypeSymbol*Sex interaction and correlation btw interaction and main effects
##need to test wether interaction term has a certain level of variance. If not, drop it.
##If the interaction term is higly correlated with main effects, then drop it.
CheckVarCorGenoSexInter <- function(data, pl, verbose=FALSE){
  if(verbose){cat("\n### CheckVarCorGenoSexInter Start ###\n")}
  keep.geno.sex.inter <- TRUE
  formula <- as.formula(paste0(pl,"~GenotypeSymbol*Sex"))
  d.mat <- model.matrix(formula, data=data)[,-1]
  cov.mat <- var(d.mat)
  cor.mat <- cor(d.mat)

  ## check variance of Genotype*Sex
  if(cov.mat[3,3]<0.03){
    keep.geno.sex.inter <- FALSE
    if(verbose){
      cat("Genotype * Sex term zero variance \n")
    }
  } else {
    ## check collinearity
    if((sum(cor.mat[,3]>0.7)-1) > 0){
      keep.geno.sex.inter <- FALSE
      if(verbose){
        cat("Genotype * Sex term collinear \n")
      }
    }
  }
  if(verbose){
    print(cor.mat)
    print(cov.mat)
    cat(paste0("\nkeep.geno.sex.inter: ", keep.geno.sex.inter, "\n"))
    cat("### CheckVarCorGenoSexInter End ###\n")
  }
  return(keep.geno.sex.inter)
}

SubStrRight <- function(x, n){
  return(substr(x, nchar(x)-n+1, nchar(x)))
}

RmvObsWthNA <- function(data, verbose=FALSE){
  if(verbose){
    cat("\n### Function RmvObsWthNA()  Start ###\n")
    print(dim(data))
  }
  data <- na.omit(data)
  if(verbose){
    cat("\n After removing obs with NA \n")
    print(dim(data))
    cat("\n### Function RmvObsWthNA()  End ###\n")
  }
  return(data)
}

DropEmptyLevel <- function(data, verbose=FALSE){
  if(verbose){
    cat("\n### Function DropEmptyLevel()  Start ###\n")
    print(summary(data))
  }
  data <- droplevels(data)
  if(verbose){
    cat("\n After dropping empty levels \n")
    print(summary(data))
    cat("\n### Function DropEmptyLevel()  End ###\n")
  }
  return(data)
}

if(FALSE){
  MakeTestData <- function(test.mutant, test.control, pl, pred.list, verbose=FALSE){
    if(verbose){cat("\n### MakeTestData Start ###\n")}

    col.list <- c(pl, pred.list)
    test.data <- rbind(test.mutant[,col.list], test.control[,col.list])

    if(verbose){
      cat("\n### Dim of test.data ###\n")
      print(dim(test.data))
    }

    test.data <- na.omit(test.data)  # remove samples with NA
    test.data <- droplevels(test.data) # remove unused levels of factors in test.data

    if(verbose){
      cat("\n### Dim of test.data after na.omit&droplevels ###\n")
      print(dim(test.data))
      cat("### MakeTestData End ###\n")
    }
    return(test.data)
  }
}

MakeTestData <- function(test.data, pl, pred.list, verbose=FALSE){
  if(verbose){cat("\n### Function MakeTestData() Start ###\n")}

  col.list <- c(pl, pred.list)
  test.data <- test.data[,col.list]

  if(verbose){
    print(dim(test.data))
    print(summary(test.data))
  }
  test.data <- na.omit(test.data)  # remove observations with NA
  test.data <- droplevels(test.data) # remove unused levels of factors in test.data

  if(verbose){
    cat("\nAfter na.omit&droplevels\n")
    print(dim(test.data))
    print(summary(test.data))
    cat("### Function MakeTestData() End ###\n")
  }
  return(test.data)
}


DropOneLevelLowVarPred <- function(test.data, pred.list){
  if(verbose){cat("\n### DropOneLevelLowVarPred Start ###\n")}
  v <- NULL
  for(pred in pred.list){
    v <- c(v, length(table(test.data[,pred]))>1 & var(as.numeric(test.data[,pred]))>0.03 )
  }
  if(verbose){cat("### DropOneLevelLowVarPred End ###\n")}
  return(pred.list[v])
}


DropPredWthOneLevel <- function(test.data, pred.list, verbose=FALSE){
  if(verbose){
    cat("\n### Function DropPredWthOneLevel() Start ###\n")
    print(pred.list)
  }
  v <- NULL
  for(pred in pred.list){
    v <- c(v, length(table(test.data[,pred]))>1)
  }
  pred.list <- pred.list[v]
  if(verbose){
    cat("\nAfter removing predictors with one level \n")
    print(pred.list)
    cat("\n### Function DropOneLevelLowVarPred End ###\n")
  }
  return(pred.list)
}

readUrl <- function(url) {
  out <- tryCatch(
    {
      ## Just to highlight: if you want to use more than one
      ## R expression in the "try" part then you'll have to
      ## use curly brackets.
      ## 'tryCatch()' will return the last evaluated expression
      ## in case the "try" part was completed successfully

      message("This is the 'try' part")

      readLines(con=url, warn=FALSE)
      ## The return value of `readLines()` is the actual value
      ## that will be returned in case there is no condition
      ## (e.g. warning or error).
      ## You don't need to state the return value via `return()` as code
      ## in the "try" part is not wrapped insided a function (unlike that
      ## for the condition handlers for warnings and error below)
    },
    error=function(cond) {
      message(paste("URL does not seem to exist:", url))
      message("Here's the original error message:")
      message(cond)
      ## Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(paste("URL caused a warning:", url))
      message("Here's the original warning message:")
      message(cond)
      ## Choose a return value in case of warning
      return(NULL)
    },
    finally={
      ## NOTE:
      ## Here goes everything that should be executed at the end,
      ## regardless of success or error.
      ## If you want more than one expression to be executed, then you
      ## need to wrap them in curly brackets ({...}); otherwise you could
      ## just have written 'finally=<expression>'
      message(paste("Processed URL:", url))
      message("Some other message at the end")
    }
  )
  return(out)
}

GetLassoResid <- function(test.data, pheno, pred.list, s='lambda.min', seed=100, nfolds=5, verbose=FALSE){
  if(verbose){
    cat("\n### Function GetLassoResid() Start ###\n")
  }
  out <- tryCatch({
    formula <- as.formula(paste(pheno,"~",paste(c(pred.list), collapse="+")))
    covars <- model.matrix(formula, test.data)[,-1]
    resid <- NULL
    pred.list.lasso <- NULL
    code <- NULL
    if(ncol(as.matrix(covars))>1){
      fit <- cv.glmnet(covars, test.data[,1], alpha=1, nfolds=nfolds, family='gaussian')
      fitted.pheno <- predict(fit, newx=covars, type="response", s=s)
      resid <- as.vector(test.data[,1] - fitted.pheno)
      if(verbose){
        cat("\nLasso coefficients \n")
        print(coef(fit, s=s))
        cat("\n### Function GetLassoResid() End ###\n")
      }
      coef.mat <- as.matrix(coef(fit, s=s))
      pred.list.lasso <- rownames(coef.mat)
      coef.vals <- as.vector(coef.mat)
      pred.list.lasso <- pred.list.lasso[coef.vals!=0]
      code <- "lasso_fitted"
    } else {
      fit <- lm(test.data[,1]~covars)
      resid <- resid(fit)
      pred.list.lasso <- names(coef(fit))
      code <- "simp_reg_fitted"
    }
    list("resid"=resid, "pred.list"=pred.list.lasso, "code"=code) ##return
  },
  error=function(cond) {
    message("Here's the original error message:")
    message(cond)
    return(NULL)
  })
  return (out)
}

GetResid <- function(test.data, pheno, pred.list, verbose=FALSE){
  if(verbose){
    cat("\n### Function GetResid() Start ###\n")
  }
  out <- tryCatch({
    formula <- as.formula(paste(pheno,"~",paste(c(pred.list), collapse="+")))
    covars <- model.matrix(formula, test.data)[,-1]
    resid <- NULL
    pred.list.lasso <- NULL
    code <- NULL
    fit <- lm(test.data[,1]~covars)
    resid <- resid(fit)
    code <- "simp_reg_fitted"
    list("resid"=resid, "code"=code) ##return
  },
  error=function(cond) {
    message("Here's the original error message:")
    message(cond)
    return(NULL)
  })
  return (out)
}

LassoVarSelect <- function(test.data, pred.list.lasso, pheno, seed=1000, nfolds=10, verbose=FALSE){
  if(verbose){cat("\n### LassoVarSelect Start ###\n")}
  ##Lasso Variable Selection
  formula.lasso <- as.formula(paste(pheno,"~",paste(c(pred.list.lasso), collapse="+")))
  covs <- model.matrix(formula.lasso, test.data)[,-1]

  ##set.seed(seed)
  lasso.res <- cv.glmnet(covs, test.data[,pheno], alpha=1, nfolds=nfolds, family='gaussian')
  lasso.res$lambda.min
  lasso.res$lambda.1se
  lambda.mean <- (lasso.res$lambda.min + lasso.res$lambda.1se)/2
  lambda.mean

  ##coefs <- coef(lasso.res, s=lambda.mean)
  coefs <- coef(lasso.res, s='lambda.1se')
  ##coefs <- coef(lasso.res, s='lambda.min')
  inds <- which(coefs!=0)
  ##inds <- which(abs(coefs) > 0.2)
  vars <- row.names(coefs)[inds]
  vars <- vars[!(vars %in% '(Intercept)')]

  pred.list.all <- colnames(test.data)[-1]

  p.vec <- rep(0, length(pred.list.all))
  for(i in 1:length(pred.list.all)){
    p.vec[i] <- sum(grepl(pred.list.all[i], vars))>0
  }
  if(verbose){
    cat("#### lasso coefs ###\n")
    print(coefs)
    cat("### LassoVarSelect End ###\n")
  }
  return(pred.list.all[as.logical(p.vec)])
}

TestBatch <- function(test.data, formula, verbose=FALSE){
  if(verbose){cat("\n### TestBatch Start ###\n")}
  model.batch   <- do.call("lme", args=list(formula, random=~1|Date.of.test.New,
                                            test.data, na.action="na.omit",method="REML"))
  model.nobatch <- do.call("gls", args=list(formula, test.data, na.action="na.omit"))

  ### ANOVA test to compare modelwithbatch and modelwithoutbatch
  keep.batch <- anova(model.batch, model.nobatch)$p[2]/2 < 0.05 ## model comparison

  if(verbose){
    cat("\n keep.batch: ",keep.batch,"\n")
    if(keep.batch){
      cat("\nRandom batch effect significant: include in the model (lme)\n")
    } else {
      cat("\nRandom batch effect non-significant: exclude from the model (gls)\n")
    }
    cat("### TestBatch End ###\n")
  }
  return(keep.batch)
}

GetPredStats <- function(model, pred.name, sd.ratio, verbose=FALSE){
  if(verbose){cat("\n### GetPredStats Start ###\n")}

  beta <- summary(model)$tTable[,"Value"][pred.name]
  tstat <- summary(model)$tTable[,"t-value"][pred.name]
  pval <- summary(model)$tTable[,"p-value"][pred.name]

  stat.vec <- c(sd.ratio*beta, beta, tstat, pval)
  names(stat.vec) <- c("std.beta", "beta", "tstat", "pval")
  if(verbose){
    cat("#### Pred: ", pred.name, "\n")
    print(stat.vec)
    cat("### GetPredStats End ###\n")
  }
  return(stat.vec)
}

SaveResultTable <- function(mat, pheno.group, mutant.list, pheno.list, res.dir, filename.str, verbose=FALSE){
  if(verbose){cat("\n### SaveResultTable Start ###\n")}
  df <- as.data.frame(mat)
  row.names(df) <- NULL
  df <- cbind(mutant.list, df)
  colnames(df) <- c("gene",pheno.list)
  write.table(df, paste0(res.dir, pheno.group,".",filename.str,".",date.str,".csv"),
              row.names=FALSE, col.names=TRUE, quote=FALSE, sep=",", append=FALSE)
  if(verbose){cat("### SaveResultTable End ###\n")}
}


## remove collinear predictors using chisq-test, length should be >= 2
RmvColliPred <- function(test.data, pred.list, verbose=FALSE){
  if(verbose){cat("\n### RmvColliPred Start ###\n")}
  pred.df <- test.data[, pred.list]
  pd.vec <-rep(TRUE,ncol(pred.df))
  for(i in ncol(pred.df):2){
    for(j in (i-1):1){
      if(chisq.test(as.matrix(table(pred.df[,c(i,j)])))$p.value < 1e-05){
        pd.vec[i] <- FALSE
        break
      }
    }
  }
  if(verbose){cat("### RmvColliPred End ###\n")}
  return(pred.list[pd.vec])
}

RmvColliPredCramerV <- function(test.data, pred.list0, cutoff=0.5, verbose=FALSE){
  if(verbose){cat("\n### RmvColliPredCramerV Start ###\n")}
  pred.df <- test.data[, pred.list0]
  pd.vec <-rep(TRUE,ncol(pred.df))
  for(i in (ncol(pred.df)):2){
    for(j in (i-1):1){
      crV <- CramerV(pred.df[,i], pred.df[,j])
      if(verbose){cat(colnames(pred.df)[i],"&",colnames(pred.df)[j], ": Cramer's V=",crV,"\n")}
      if(crV > cutoff){
        pd.vec[i] <- FALSE
        break
      }
    }
  }
  if(verbose){cat("### RmvColliPredCramerV End ###\n")}
  return(pred.list0[pd.vec])
}


## IMPC model
RunIMPC <- function(test.data, pred.list, verbose=FALSE){
  if(verbose){cat("\n### RunIMPC Start ###\n")}
  pl <- colnames(test.data)[1]
  pred.list <- pred.list[grep("GenotypeSymbol|Sex|BW", pred.list)]

  keep.geno.sex.inter <- CheckVarCorGenoSexInter(test.data, pl, verbose=verbose)
  pred.str <- NULL
  if(keep.geno.sex.inter){
    pred.str <- paste(c(pred.list,"GenotypeSymbol*Sex"), collapse="+")
  } else {
    pred.str <- paste(pred.list, collapse="+")
  }

  formula <- as.formula(paste(pl,"~",pred.str))

  ## Test batch
  keep.batch <- TestBatch(test.data, formula, verbose=verbose)
  model <- NULL
  if(keep.batch){
    model <- lme(formula, random=~1|Date.of.test.New, test.data, na.action="na.omit", method="REML")
  } else {
    model <- gls(formula, test.data, na.action="na.omit")
  }

  ## Test fixed effects
  anova.mar.res <- anova(model, type="marginal") # ANOVA Type III test.
  if(verbose){
    cat("\n### test fixed effects using ANOVA Type III test ###\n")
    print(anova.mar.res)
  }

  keep.fixed.effects <- anova.mar.res$"p-value" < 0.05
  names(keep.fixed.effects) <- rownames(anova.mar.res)
  keep.fixed.effects <- keep.fixed.effects[c(-1,-2)] # Don't test Intercept, Geno

  fixed.effects <- NULL
  for(name in names(keep.fixed.effects)){
    if(keep.fixed.effects[name]){
      fixed.effects <- c(fixed.effects, name)
    }
  }

  fixed.effects.str <- paste(fixed.effects, collapse="+")
  if(fixed.effects.str==""){
    pred.str <- paste("GenotypeSymbol", sep="+")
  } else {
    pred.str <- paste("GenotypeSymbol", fixed.effects.str, sep="+")
  }
  final.formula <- as.formula(paste(pl,"~", pred.str))

  final.model <- NULL
  aic.val <- NULL
  if(keep.batch){
    final.model <- lme(final.formula, random=~1|Date.of.test.New, test.data,
                       na.action="na.omit", method="REML")

    formula.lmer <- as.formula(paste(pl,"~", paste(c(pred.str,"(1|Date.of.test.New)"), collapse="+")))
    aic.model <- lmer(formula.lmer, test.data)
    aic.val <- cAIC(aic.model)$caic

  } else {
    final.model <- gls(final.formula, test.data, na.action="na.omit")
    aic.val <- AIC(final.model)
  }

  if(verbose){cat("### RunIMPC End ###\n")}
  return(list(model=final.model, aic=aic.val, batch=keep.batch))
}


## Regression of residuals
RunREGRES <- function(test.data, pred.list, verbose=FALSE){
  if(verbose){cat("\n### RunREGRES Start###\n")}
  pl <- colnames(test.data)[1]
  pred.list.ols <- pred.list[grep("GenotypeSymbol|Sex|BW", pred.list)]
  pred.list.res <- pred.list[!pred.list%in%c(pred.list.ols,"Date.of.test.New")]

  pred.list.res <- DropOneLevelLowVarPred(test.data, pred.list.res)
  if(length(pred.list.res)>1){
    pred.list.res <- RmvColliPredCramerV(test.data, pred.list.res)
  }

  ##Get residuals
  if(length(pred.list.res)){
    formula.res <- as.formula(paste(pl,"~",paste(pred.list.res, collapse="+")))
    model.res <- gls(formula.res, test.data, na.action="na.exclude")
    test.data$res <- residuals(model.res)
  } else {
    test.data$res <- test.data[,pl]
  }

  keep.geno.sex.inter <- CheckVarCorGenoSexInter(test.data, pl, verbose=verbose)

  pred.str <- NULL
  if(keep.geno.sex.inter){
    pred.str <- paste(c(pred.list.ols,"GenotypeSymbol*Sex"), collapse="+")
  } else {
    pred.str <- paste(pred.list.ols, collapse="+")
  }

  formula <- as.formula(paste("res ~",pred.str))

  ## Test batch
  keep.batch <- TestBatch(test.data, formula, verbose=verbose)
  model <- NULL
  if(keep.batch){
    model <- lme(formula, random=~1|Date.of.test.New, test.data, na.action="na.omit", method="REML")
  } else {
    model <- gls(formula, test.data, na.action="na.omit")
  }
  if(verbose){
    cat("\n### summary of model ###\n")
    print(summary(model))
  }

  ## Test fixed effects
  anova.mar.res <- anova(model, type="marginal") # ANOVA Type III test.

  keep.fixed.effects <- anova.mar.res$"p-value" < 0.05
  names(keep.fixed.effects) <- rownames(anova.mar.res)
  keep.fixed.effects <- keep.fixed.effects[c(-1,-2)] # Don't test Intercept, Geno

  fixed.effects <- NULL
  for(name in names(keep.fixed.effects)){
    if(keep.fixed.effects[name]){
      fixed.effects <- c(fixed.effects, name)
    }
  }

  fixed.effects.str <- paste(fixed.effects, collapse="+")
  if(fixed.effects.str==""){
    pred.str <- paste("GenotypeSymbol", sep="+")
  } else {
    pred.str <- paste("GenotypeSymbol", fixed.effects.str, sep="+")
  }
  final.formula <- as.formula(paste("res ~", pred.str))

  final.model <- NULL
  aic.val <- NULL
  if(keep.batch){
    final.model <- lme(final.formula, random=~1|Date.of.test.New, test.data,
                       na.action="na.omit", method="REML")

    formula.lmer <- as.formula(paste("res ~", paste(c(pred.str,"(1|Date.of.test.New)"), collapse="+")))
    aic.model <- lmer(formula.lmer, test.data)
    aic.val <- cAIC(aic.model)$caic
  } else {
    final.model <- gls(final.formula, test.data, na.action="na.omit")
    aic.val <- AIC(final.model)
  }

  if(verbose){cat("### RunREGRES End ###\n")}
  return(list(model=final.model, aic=aic.val, batch=keep.batch))
}

## OLS post Lasso variable selection
RunLASSO <- function(test.data, pred.list, verbose=FALSE){
  if(verbose){cat("\n### RunLasso Start ###\n")}
  pl <- colnames(test.data)[1]
  pred.list.lasso <- DropOneLevelLowVarPred(test.data, pred.list)
  pred.list.ols <- LassoVarSelect(test.data, pred.list.lasso, pl, seed=1200, verbose=verbose)

  pred.indicator <- as.integer(pred.list%in%pred.list.ols)
  names(pred.indicator) <- pred.list

  ##Remove highly correlated categorical predictors
  bw.in <- "BW"%in%pred.list.ols
  clist <- pred.list.ols[!pred.list.ols%in%"BW"]
  geno.in <- sum(pred.list.ols%in%"GenotypeSymbol")
  if(!geno.in){
    clist <- c("GenotypeSymbol", clist)
  }
  if(length(clist)>1){
    clist <- RmvColliPredCramerV(test.data, clist,verbose=verbose)
    ##clist <- RmvColliPred(test.data, clist)
    if(bw.in){
      clist <- c(clist, "BW")
    }
    pred.list.ols <- clist
  }

  if(FALSE){ ## for debugging
    formula <- as.formula(paste(pl,"~", paste(pred.list.ols, collapse="+")))
    model <- gls(formula, test.data, na.action="na.omit")
  }

  test.batch <- FALSE
  if("Date.of.test.New"%in%pred.list.ols){
    test.batch <- TRUE
    pred.list.ols <- pred.list.ols[!(pred.list.ols %in% 'Date.of.test.New')]
  }

  if("GenotypeSymbol"%in%pred.list.ols){
    pred.list.ols <- pred.list.ols[!(pred.list.ols %in% 'GenotypeSymbol')]
  }

  keep.geno.sex.inter <- FALSE
  if("Sex"%in%pred.list.ols){
    keep.geno.sex.inter <- CheckVarCorGenoSexInter(test.data, pl, verbose=verbose)
    if(keep.geno.sex.inter)
      pred.list.ols <- c("GenotypeSymbol*Sex",pred.list.ols)
  }

  pred.list.ols <- c("GenotypeSymbol",pred.list.ols)
  formula <- as.formula(paste(pl,"~", paste(pred.list.ols, collapse="+")))


  keep.batch <- FALSE
  if(test.batch){
    keep.batch <- TestBatch(test.data, formula, verbose=verbose)
  }


  final.model <- NULL
  aic.val <- NULL
  if(keep.batch){
    final.model <- lme(formula, random=~1|Date.of.test.New, test.data, na.action="na.omit", method="REML")

    formula.lmer <- as.formula(paste(pl,"~", paste(c(pred.list.ols,"(1|Date.of.test.New)"), collapse="+")))
    aic.model <- lmer(formula.lmer, test.data)
    aic.val <- cAIC(aic.model)$caic
  } else {
    final.model <- gls(formula, test.data, na.action="na.omit")
    aic.val <- AIC(final.model)
  }
  if(verbose){cat("### RunLasso End ###\n")}
  return(list(model=final.model, pred.list=pred.indicator, aic=aic.val, batch=keep.batch))
}

CramerV <- function(x,y){
  CV = sqrt(chisq.test(x, y, correct=FALSE)$statistic /
              (length(x) * (min(length(unique(x)),length(unique(y))) - 1)))
  ##print.noquote("Cramér V / Phi:")
  return(as.numeric(CV))
}


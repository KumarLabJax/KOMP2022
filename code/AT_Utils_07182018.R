## Association Testing, IMPC model
RunIMPC.AT <- function(test.data, pheno, verbose=FALSE){
    if(verbose){
        cat("\n### Function RunIMPC() Start ###\n")
    }
    ##out <- tryCatch({

        pred.list <- pred.list[grep("GenotypeSymbol|Sex|BW", )]
    
    #keep.geno.sex.inter <- CheckVarCorGenoSexInter(test.data, pl, verbose=verbose)
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
RunREGRES.AT <- function(test.data, pred.list, verbose=FALSE){
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
RunLASSO.AT <- function(test.data, pheno, pred.list, varhetero.pval.cutoff=0.05, verbose=FALSE){

    if(verbose){
        cat("\n### Function RunLASSO() Start ###\n")
    }
    
    out <- tryCatch({        
        test.data$Geno.Sex.Inter <- model.matrix(~GenotypeSymbol*Sex,data=test.data)[,4]
        pred.list.wth.gs <- c("Geno.Sex.Inter", pred.list)
        pred.list.lasso <- DropOneLevelLowVarPred(test.data, pred.list.wth.gs)
        pred.list.ols <- LassoVarSelect(test.data, pred.list.lasso, pl, seed=1200, verbose=verbose)
        
        pred.indicator <- as.integer(pred.list.wth.gs%in%pred.list.ols)
        names(pred.indicator) <- pred.list.wth.gs

        surv.env.pred.list <- pred.list.ols[!pred.list.ols%in%c("Geno.Sex.Inter","GenotypeSymbol","Sex","Date.of.test.New","BW")]
    
        if(!"GenotypeSymbol"%in%pred.list.ols){
            ## add GenotypeSymbol from pred.list.ols
            pred.list.ols <- c("GenotypeSymbol",pred.list.ols)    
        }

        keep.batch <- FALSE
        if("Date.of.test.New"%in%pred.list.ols){
            keep.batch <- TRUE
            ## remove Date.of.test.New from pred.list.ols to model it as random effect
            pred.list.ols <- pred.list.ols[!(pred.list.ols %in% 'Date.of.test.New')] 
        }
        
        formula <- as.formula(paste(pheno,"~", paste(pred.list.ols, collapse="+")))

        
        ## ----------------------------------------
        ## Testing Residual Varaince Heterogeneity
        ## 
        ## We adopted codes (lines 120-133) from MMFramework.R of the PhenStat R package (PhenStat_2.17.0.tar.gz)
        ## and slightly modified for our analysis.
        if(keep.batch){            
            model.varhomo <- do.call("lme", args=list(formula, random=~1|Date.of.test.New, test.data, na.action="na.omit", method="REML"))
            model.varhetero <- do.call("lme", args=list(formula, random=~1|Date.of.test.New, test.data, weights=varIdent(form=~1|GenotypeSymbol), na.action="na.omit", method="REML"))

        } else {
            model.varhomo <- do.call("gls", args=list(formula, test.data, na.action="na.omit"))
            model.varhetero <- do.call("gls", args=list(formula, test.data, weights=varIdent(form=~1|GenotypeSymbol), na.action="na.omit"))
        }

        final.model <- NULL
        var.equal <- NULL
        if(anova(model.varhomo, model.varhetero)$"p-value"[2] >= varhetero.pval.cutoff){
            final.model <- model.varhomo
            var.equal <- TRUE
        } else {
            final.model <- model.varhetero
            var.equal <- FALSE
        }
        ##
        ##--------------------------------
        
        
        ## SD of phenotype
        sd.pheno <- sd(test.data[,pheno], na.rm=TRUE)
        ## Genotype phenotype SD ratio
        geno.pheno.sd.ratio <- sd(test.data[,"GenotypeSymbol"])/sd.pheno

        geno.std.beta <- summary(final.model)$tTable[,"Value"]["GenotypeSymbol"] * geno.pheno.sd.ratio
        geno.beta <- summary(final.model)$tTable[,"Value"]["GenotypeSymbol"]
        geno.stde  <- summary(final.model)$tTable[,"Std.Error"]["GenotypeSymbol"]
        geno.tstat <- summary(final.model)$tTable[,"t-value"]["GenotypeSymbol"]
        geno.pval <- summary(final.model)$tTable[,"p-value"]["GenotypeSymbol"]
        
        stat.vec <- c(geno.std.beta, geno.beta, geno.stde, geno.tstat, geno.pval)
        names(stat.vec) <- c("geno.stat.beta","geno.beta","geno.stde","geno.tstat","geno.pval")
        list(stat.vec=stat.vec, var.equal=var.equal, code="model_success", pred.ind=pred.indicator, surv.env.pred.list=paste(surv.env.pred.list, collapse=":"))
    },
    error=function(cond) {
        message("Here's the original error message:")
        message(cond)
        pred.indicator <- rep(NA, length(pred.list)+1)
        names(pred.indicator) <- c("Geno.Sex.Inter", pred.list)
        #names(pred.indicator) <- pred.list
        return(list(stat.vec=rep(NA,5), var.equal= NA , code="model_error", pred.ind=pred.indicator, surv.env.pred.list=NA))
    })
    return (out)
}

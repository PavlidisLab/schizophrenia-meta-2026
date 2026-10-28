# Primary model on several person tables. Usage: Rscript fit_primary.R outdir nboot table1.tsv table2.tsv ...
suppressMessages(library(lme4)); args <- commandArgs(TRUE); out <- args[1]; nb <- as.integer(args[2]); files <- args[-(1:2)]; set.seed(20261014)
ctrl <- lmerControl(optimizer='bobyqa', optCtrl=list(maxfun=2e5), check.conv.singular='ignore')
res <- list(); pcs <- list()
for(f in files){
  d <- read.delim(f, comment.char='#', stringsAsFactors=FALSE, na.strings=c('','NA')); m <- d[d$group %in% c('SZ','CTRL'),]
  m$sz <- as.numeric(m$group=='SZ'); m$cohort <- factor(m$cohort); m$sc <- paste(m$study,m$cohort,sep=':')
  fit <- lmer(ph ~ sz + (1+sz|cohort) + (1|sc), m, REML=TRUE, control=ctrl)
  pr <- tryCatch(confint(fit, method='profile', parm='sz', quiet=TRUE), error=function(e) c(NA,NA))
  pa <- tryCatch(confint(fit, method='profile', quiet=TRUE), error=function(e) NULL)
  bs <- tryCatch(bootMer(fit, function(x) c(fixef(x)['sz'], attr(VarCorr(x)$cohort,'stddev')['sz']), nsim=nb, type='parametric', use.u=FALSE, seed=20261014), error=function(e) NULL)
  bci <- if(!is.null(bs)) quantile(bs$t[,1], c(.025,.975), na.rm=TRUE) else c(NA,NA); bsd <- if(!is.null(bs)) quantile(bs$t[,2], c(.025,.975), na.rm=TRUE) else c(NA,NA)
  w <- fixef(fit)['sz']; se <- sqrt(vcov(fit)['sz','sz'])
  sdv <- attr(VarCorr(fit)$cohort,'stddev')['sz']
  nm <- sub('.tsv','',basename(f))
  res[[nm]] <- data.frame(table=nm, n_persons=nrow(m), n_SZ=sum(m$sz==1), n_CTRL=sum(m$sz==0), n_cohorts=nlevels(droplevels(m$cohort)), n_studies=length(unique(m$sc)), estimate=w, profile_lo=pr[1], profile_hi=pr[2], bootstrap_lo=bci[1], bootstrap_hi=bci[2], wald_lo=w-1.96*se, wald_hi=w+1.96*se, between_cohort_sd=sdv, sd_profile_lo=if(!is.null(pa)&&'.sig03'%in%rownames(pa)) pa['.sig03',1] else NA, sd_profile_hi=if(!is.null(pa)&&'.sig03'%in%rownames(pa)) pa['.sig03',2] else NA, sd_boot_lo=bsd[1], sd_boot_hi=bsd[2], singular=isSingular(fit))
  re <- ranef(fit, condVar=TRUE)$cohort; pv <- attr(re,'postVar')
  for(c in levels(droplevels(m$cohort))){
    x <- m[m$cohort==c,]; x$st <- factor(x$study); mm <- if(nlevels(x$st)>1) lm(ph ~ sz + st, x) else lm(ph ~ sz, x); ci <- tryCatch(confint(mm)['sz',], error=function(e) c(NA,NA))
    k <- which(rownames(re)==c); sh <- w + re[k,'sz']; s2 <- sqrt(pv[2,2,k]+vcov(fit)['sz','sz'])
    pcs[[length(pcs)+1]] <- data.frame(table=nm, cohort=c, n_SZ=sum(x$sz==1), n_CTRL=sum(x$sz==0), raw=if(sum(x$sz==1)>0&sum(x$sz==0)>0) coef(mm)['sz'] else NA, raw_lo=ci[1], raw_hi=ci[2], shrunken=sh, shrunk_lo=sh-1.96*s2, shrunk_hi=sh+1.96*s2) }
}
write.table(do.call(rbind,res), file.path(out,'fit_primary.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
write.table(do.call(rbind,pcs), file.path(out,'fit_primary_per_cohort.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
cat('done\n')

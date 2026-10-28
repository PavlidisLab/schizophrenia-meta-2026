# Pooled donor-level pH models. R + lme4 + nlme. One thread. Usage: Rscript run_models.R <donor_level_table.tsv> <outdir>
suppressMessages({library(lme4);library(nlme)})
args <- commandArgs(TRUE); inp <- args[1]; out <- args[2]; dir.create(out, showWarnings=FALSE)
set.seed(20261014)
d <- read.delim(inp, comment.char='#', stringsAsFactors=FALSE, na.strings=c('','NA'))
d$sz <- as.numeric(d$group=='SZ'); d$bd <- as.numeric(d$group=='BD')
main <- d[d$group %in% c('SZ','CTRL'),]
main$cohort <- factor(main$cohort); main$sc <- paste(main$study, main$cohort, sep=':')
ctrl <- c(lmerControl(optimizer='bobyqa', optCtrl=list(maxfun=2e5)))
ctrl <- lmerControl(optimizer='bobyqa', optCtrl=list(maxfun=2e5), check.conv.singular='ignore')
rows <- list()
add <- function(model, label, fit=NULL, est=NA, lo=NA, hi=NA, method='', n=NA, ncoh=NA, nstu=NA, slope_sd=NA, sd_lo=NA, sd_hi=NA, note='', term='SZ minus control (pH units)')
  rows[[length(rows)+1]] <<- data.frame(model=model, label=label, term=term, estimate=est, lo=lo, hi=hi, interval=method, n_donors=n, n_cohorts=ncoh, n_studies=nstu, between_cohort_sd_of_effect=slope_sd, sd_lo=sd_lo, sd_hi=sd_hi, note=note, stringsAsFactors=FALSE)
vc_sd <- function(f, g='cohort', nm='sz'){ v <- VarCorr(f)[[g]]; if(is.null(v)) return(NA); attr(v,'stddev')[nm] }
profci <- function(f, parm='sz'){ p <- tryCatch(confint(f, method='profile', parm=parm, quiet=TRUE), error=function(e) NULL)
  if(is.null(p) || any(is.na(p))){ w <- wald(f, parm); list(lo=w[2], hi=w[3], method='Wald 95% (profile failed)') } else list(lo=p[1], hi=p[2], method='profile likelihood 95%') }
wald <- function(f, nm='sz'){ b <- fixef(f)[nm]; s <- sqrt(vcov(f)[nm,nm]); c(b, b-1.96*s, b+1.96*s) }
nn <- function(x) c(nrow(x), nlevels(droplevels(x$cohort)), length(unique(x$sc)))
fitm <- function(form, dat) lmer(form, data=dat, REML=TRUE, control=ctrl)
note_fit <- function(f) paste0(if(isSingular(f)) 'singular fit; ' else '', 'optimizer msgs: ', paste(f@optinfo$conv$lme4$messages, collapse='|'))
# ---- Model 1 (primary)
f1 <- fitm(ph ~ sz + (1+sz|cohort) + (1|sc), main)
pr <- tryCatch(confint(f1, method='profile', parm='sz', quiet=TRUE), error=function(e) NULL)
prall <- tryCatch(confint(f1, method='profile', quiet=TRUE), error=function(e) NULL)
w <- wald(f1)
sdp <- if(!is.null(prall)) prall['.sig03',] else c(NA,NA)
cat('profile ok:', !is.null(pr), ' sd rownames:', if(!is.null(prall)) paste(rownames(prall),collapse=',') else 'none','\n')
nb <- 1000
bs <- tryCatch(bootMer(f1, function(x) c(fixef(x)['sz'], attr(VarCorr(x)$cohort,'stddev')['sz']), nsim=nb, type='parametric', use.u=FALSE, seed=20261014), error=function(e) NULL)
bci <- if(!is.null(bs)) quantile(bs$t[,1], c(.025,.975), na.rm=TRUE) else c(NA,NA)
bsd <- if(!is.null(bs)) quantile(bs$t[,2], c(.025,.975), na.rm=TRUE) else c(NA,NA)
sd1 <- vc_sd(f1)
add('1','Primary: pH ~ SZ + (1+SZ|cohort) + (1|study:cohort), Gaussian, REML', est=w[1], lo=if(!is.null(pr)) pr[1] else NA, hi=if(!is.null(pr)) pr[2] else NA, method='profile likelihood 95%', n=nn(main)[1], ncoh=nn(main)[2], nstu=nn(main)[3], slope_sd=sd1, sd_lo=sdp[1], sd_hi=sdp[2], note=note_fit(f1))
add('1','Primary, parametric bootstrap interval (1000 sims, percentile)', est=w[1], lo=bci[1], hi=bci[2], method='parametric bootstrap percentile 95%', n=nn(main)[1], ncoh=nn(main)[2], nstu=nn(main)[3], slope_sd=sd1, sd_lo=bsd[1], sd_hi=bsd[2], note=paste0('bootstrap sims ',nb))
add('1','Primary, Wald interval', est=w[1], lo=w[2], hi=w[3], method='Wald 95%', n=nn(main)[1], ncoh=nn(main)[2], nstu=nn(main)[3], slope_sd=sd1)
vcs <- as.data.frame(VarCorr(f1)); write.table(vcs, file.path(out,'model1_varcorr.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
if(!is.null(prall)) write.table(data.frame(parm=rownames(prall), prall), file.path(out,'model1_profile_all.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
saveRDS(list(res=residuals(f1,scaled=TRUE), fitted=fitted(f1)), file.path(out,'model1_resid.rds'))
# ---- per-cohort estimates: raw (within-cohort lm with study fixed effects) and shrunken (BLUP)
re <- ranef(f1, condVar=TRUE)$cohort; pv <- attr(re,'postVar')
pc <- list()
for(c in levels(main$cohort)){
  x <- main[main$cohort==c,]; x$st <- factor(x$study)
  m <- if(nlevels(x$st)>1) lm(ph ~ sz + st, x) else lm(ph ~ sz, x)
  ci <- confint(m)['sz',]; k <- which(rownames(re)==c)
  sh <- fixef(f1)['sz'] + re[k,'sz']; se <- sqrt(pv[2,2,k] + vcov(f1)['sz','sz'])
  pc[[c]] <- data.frame(cohort=c, n_studies=nlevels(x$st), n_SZ=sum(x$sz==1), n_CTRL=sum(x$sz==0), raw=coef(m)['sz'], raw_lo=ci[1], raw_hi=ci[2], shrunken=sh, shrunk_lo=sh-1.96*se, shrunk_hi=sh+1.96*se,
    mean_SZ=mean(x$ph[x$sz==1]), mean_CTRL=mean(x$ph[x$sz==0]), median_SZ=median(x$ph[x$sz==1]), median_CTRL=median(x$ph[x$sz==0]))
}
pc <- do.call(rbind, pc); write.table(pc, file.path(out,'per_cohort_raw.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
# per study (raw Welch) for reference
ps <- do.call(rbind, lapply(split(main, main$sc), function(x){ if(sum(x$sz==1)<2||sum(x$sz==0)<2) return(NULL); t <- t.test(ph~sz, x); data.frame(study_cohort=x$sc[1], cohort=as.character(x$cohort[1]), n_SZ=sum(x$sz==1), n_CTRL=sum(x$sz==0), diff=diff(t$estimate), se=t$stderr, lo=-t$conf.int[2], hi=-t$conf.int[1]) }))
write.table(ps, file.path(out,'per_study_raw.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
dl <- function(y,s){ w <- 1/s^2; fe <- sum(w*y)/sum(w); Q <- sum(w*(y-fe)^2); k <- length(y); tau2 <- max(0,(Q-(k-1))/(sum(w)-sum(w^2)/sum(w))); ws <- 1/(s^2+tau2); m <- sum(ws*y)/sum(ws); se <- sqrt(1/sum(ws)); c(m, m-1.96*se, m+1.96*se, k, sqrt(tau2), sum(y<0)) }
mm <- dl(ps$diff, ps$se)
add('0','Study-level DerSimonian-Laird random-effects mean on the de-duplicated persons (studies with >=2 per group); for comparison', est=mm[1], lo=mm[2], hi=mm[3], method='DL Wald 95%', n=sum(ps$n_SZ+ps$n_CTRL), nstu=mm[4], slope_sd=mm[5], note=paste0(mm[6],' of ',mm[4],' studies lower; study-level tau'))
# ---- Model 2: covariates
m2 <- main[complete.cases(main[,c('age','sex','pmi_h')]),]; m2$age10 <- m2$age/10; m2$pmi10 <- m2$pmi_h/10; m2$male <- as.numeric(m2$sex=='M')
f2 <- fitm(ph ~ sz + age10 + male + pmi10 + (1+sz|cohort) + (1|sc), m2)
pp2 <- profci(f2); p2 <- c(pp2$lo,pp2$hi); w2 <- wald(f2)
add('2','Covariates: + age (per 10 y) + sex + PMI (per 10 h), complete cases', est=w2[1], lo=p2[1], hi=p2[2], method=pp2$method, n=nn(m2)[1], ncoh=nn(m2)[2], nstu=nn(m2)[3], slope_sd=vc_sd(f2), note=paste0(note_fit(f2),'; dropped for missing covariates: ',nrow(main)-nrow(m2)))
cf <- summary(f2)$coefficients; write.table(data.frame(term=rownames(cf),cf), file.path(out,'model2_fixed.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
# ---- Model 3: centred on control mean
cm <- tapply(main$ph[main$sz==0], main$cohort[main$sz==0], mean); main$ph_cc <- main$ph - cm[as.character(main$cohort)]
sm <- tapply(main$ph[main$sz==0], main$sc[main$sz==0], mean); main$ph_sc <- main$ph - sm[main$sc]
for(v in c('ph_cc','ph_sc')){
  f <- fitm(as.formula(paste(v,'~ sz + (1+sz|cohort) + (1|sc)')), main[!is.na(main[[v]]),]); pp <- profci(f); p <- c(pp$lo,pp$hi); w3 <- wald(f)
  add('3', if(v=='ph_cc') 'pH minus the cohort control mean' else 'pH minus the study control mean (variant)', est=w3[1], lo=p[1], hi=p[2], method=pp$method, n=sum(!is.na(main[[v]])), ncoh=nn(main)[2], nstu=nn(main)[3], slope_sd=vc_sd(f), note=note_fit(f))
}
# ---- Model 4: heteroscedastic residual variance by study (nlme)
tab <- table(main$sc); main$vc <- ifelse(tab[main$sc] >= 10, main$sc, 'small studies (n<10)'); main$vc <- factor(main$vc)
f4 <- tryCatch(lme(ph ~ sz, random=list(cohort=pdSymm(~sz), sc=~1), weights=varIdent(form=~1|vc), data=main, method='REML', control=lmeControl(maxIter=500, msMaxIter=500, opt='optim')), error=function(e) {message('nlme het failed: ',e$message); NULL})
if(!is.null(f4)){ iv <- tryCatch(intervals(f4, which='fixed')$fixed['sz',], error=function(e) c(NA,NA,NA)); vcm <- VarCorr(f4); 
  sdsl <- tryCatch(as.numeric(vcm[grep('sz',rownames(vcm))[1],'StdDev']), error=function(e) NA)
  add('4','Heteroscedastic: residual variance by study (varIdent; studies with n<10 share one variance), nlme', est=iv[2], lo=iv[1], hi=iv[3], method='nlme intervals() approximate 95%', n=nrow(main), ncoh=nn(main)[2], nstu=nn(main)[3], slope_sd=sdsl, note=paste0('variance classes: ',nlevels(main$vc))) }
f4b <- tryCatch(lme(ph ~ sz, random=list(cohort=pdSymm(~sz), sc=~1), data=main, method='REML', control=lmeControl(opt='optim')), error=function(e) NULL)
if(!is.null(f4b)){ iv <- intervals(f4b, which='fixed')$fixed['sz',]; add('4','Same nlme model, homoscedastic (reference for model 4)', est=iv[2], lo=iv[1], hi=iv[3], method='nlme intervals() approximate 95%', n=nrow(main), ncoh=nn(main)[2], nstu=nn(main)[3]); if(!is.null(f4)) { an <- anova(f4b, f4); cat('LRT het:', an$`p-value`[2], '\n'); add('4','LRT homoscedastic vs heteroscedastic (P in estimate column)', est=an$`p-value`[2], method='likelihood ratio', note='estimate column holds the P value') } }
# ---- Model 5: matched sets
tryCatch({
m5 <- main[!is.na(main$pair),]
f5a <- fitm(ph ~ sz + (1+sz|cohort) + (1|sc), m5); f5 <- fitm(ph ~ sz + (1+sz|cohort) + (1|sc) + (1|pair), m5)
for(nm in c('5a','5')){ f <- get(paste0('f',nm)); pp <- profci(f); p <- c(pp$lo,pp$hi); w5 <- wald(f)
  add('5', if(nm=='5a') 'Donors in matched-set designs (Arion 2010, 2015, 2017; Bowden; Weidenhofer), without pair term' else 'same donors, + (1|pair)', est=w5[1], lo=p[1], hi=p[2], method=pp$method, n=nrow(m5), ncoh=nn(m5)[2], nstu=nn(m5)[3], slope_sd=vc_sd(f), note=paste0(note_fit(f),'; pairs/sets: ',length(unique(m5$pair)))) }

}, error=function(e) message('model 5 skipped: ', e$message))
# ---- Model 6 sensitivities (Wald interval unless noted)
sens <- function(label, dat, form=ph ~ sz + (1+sz|cohort) + (1|sc), prof=FALSE, extra=NULL){
  dat <- droplevels(dat); f <- tryCatch(fitm(form, dat), error=function(e) NULL); if(is.null(f)){ add('6',label,note='fit failed'); return(invisible()) }
  w <- wald(f); pp <- if(prof) profci(f) else list(lo=w[2],hi=w[3],method='Wald 95%'); p <- c(pp$lo,pp$hi)
  add('6', label, est=w[1], lo=p[1], hi=p[2], method=pp$method, n=nrow(dat), ncoh=nlevels(dat$cohort), nstu=length(unique(dat$sc)), slope_sd=tryCatch(vc_sd(f),error=function(e) NA), note=paste0(note_fit(f), if(!is.null(extra)) extra(f) else '')) }
for(c in levels(main$cohort)) sens(paste0('leave out cohort: ',c), main[main$cohort!=c,], prof=TRUE)
sens('drop donors with text-derived pH (article-text tables)', main[main$text_derived==0,], prof=TRUE)
sens('drop donors that are samples of the 14 dataset-regions (added donors only)', main[main$in_14_dataset_regions==0,], prof=TRUE)
sens('only donors that are samples of the 14 dataset-regions (9 original studies)', main[main$in_14_dataset_regions==1,], prof=TRUE)
sens('only donors that are samples of any of the 18 dataset-regions with pH', main[main$in_expression_datasets==1,], prof=TRUE)
sens('drop donors that are samples of any of the 18 dataset-regions with pH', main[main$in_expression_datasets==0,], prof=TRUE)
sens('only donors tabulated in one record (no pH choice made)', main[main$n_records==1,], prof=TRUE)
m6 <- main; m6$ph <- m6$ph_mean_all; sens('pH = mean of all sources printing a pH for the person', m6, prof=TRUE)
sens('person-level study = random slope as well: + (0+SZ|study:cohort) -> (1+SZ|study:cohort)', main, form=ph ~ sz + (1+sz|cohort) + (1+sz|sc), prof=FALSE, extra=function(f) paste0('; study-level SD of the SZ effect = ',round(vc_sd(f,'sc'),4)))
x <- main; x$cohort <- as.character(x$cohort); x$cohort[x$cohort=='other single-study'] <- paste0('single:',x$study[x$cohort=='other single-study']); x$cohort <- factor(x$cohort)
sens('each single-study pool member as its own cohort', x, prof=FALSE)
# outliers and rank-normal
r <- residuals(f1, scaled=TRUE); isout <- abs(r)>3
sens(paste0('drop |scaled residual|>3 (',sum(isout),' donors) and refit'), main[!isout,], prof=TRUE)
main$z <- NA; for(s in unique(main$sc)){ k <- main$sc==s; main$z[k] <- qnorm((rank(main$ph[k])-0.5)/sum(k)) }
fz <- fitm(z ~ sz + (1+sz|cohort) + (1|sc), main); wz <- wald(fz); ppz <- profci(fz); pz <- c(ppz$lo,ppz$hi)
add('7','Within-study rank-based normal scores of pH (effect in SD-like units, not pH units)', est=wz[1], lo=pz[1], hi=pz[2], method=ppz$method, n=nrow(main), ncoh=nn(main)[2], nstu=nn(main)[3], slope_sd=vc_sd(fz), note=note_fit(fz), term='SZ minus control (normal-score units)')
# bipolar third level
bdat <- d[d$group %in% c('SZ','CTRL','BD'),]; bdat$cohort <- factor(bdat$cohort); bdat$sc <- paste(bdat$study, bdat$cohort, sep=':')
fb <- fitm(ph ~ sz + bd + (1+sz|cohort) + (1|sc), bdat); wb <- wald(fb); wbd <- wald(fb,'bd')
pb <- tryCatch(confint(fb, method='profile', parm=c('sz','bd'), quiet=TRUE), error=function(e) matrix(NA,2,2))
add('6','bipolar as third level: SZ minus control', est=wb[1], lo=pb[1,1], hi=pb[1,2], method='profile likelihood 95%', n=nrow(bdat), ncoh=nlevels(bdat$cohort), nstu=length(unique(bdat$sc)), slope_sd=vc_sd(fb), note=paste0(note_fit(fb),'; bipolar n=',sum(bdat$bd)))
add('6','bipolar as third level: bipolar minus control', est=wbd[1], lo=pb[2,1], hi=pb[2,2], method='profile likelihood 95%', n=nrow(bdat), term='bipolar minus control (pH units)', note='bipolar fixed effect only; no bipolar random slope; bipolar persons de-duplicated approximately (same cohort, sex, age, PMI)')
# cohort-level fixed-effect alternative (no random slope): pH ~ sz + cohort + (1|sc)
fa <- fitm(ph ~ sz + cohort + (1|sc), main); wa <- wald(fa); add('6','No random slope: pH ~ SZ + cohort (fixed) + (1|study:cohort)', est=wa[1], lo=wa[2], hi=wa[3], method='Wald 95%', n=nrow(main), ncoh=nn(main)[2], nstu=nn(main)[3], note=note_fit(fa))
res <- do.call(rbind, rows); write.table(res, file.path(out,'models.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
# ---- descriptive thresholds
th <- list()
for(c in c(levels(main$cohort),'all')){ x <- if(c=='all') main else main[main$cohort==c,]
  for(t in c(6.0,6.3,6.5)){ a <- x$ph[x$sz==1]; b <- x$ph[x$sz==0]
    th[[length(th)+1]] <- data.frame(cohort=c, cutoff=t, n_SZ=length(a), n_CTRL=length(b), SZ_below=sum(a<t), CTRL_below=sum(b<t), SZ_prop=mean(a<t), CTRL_prop=mean(b<t), diff_prop=mean(a<t)-mean(b<t)) } }
write.table(do.call(rbind,th), file.path(out,'thresholds.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
# ---- distribution checks
dig1 <- function(x) mean(abs(x*10-round(x*10))<1e-6); dig2 <- function(x) mean(abs(x*100-round(x*100))<1e-6)
sk <- function(x) mean((x-mean(x))^3)/sd(x)^3; ku <- function(x) mean((x-mean(x))^4)/sd(x)^4
em2 <- function(x){ x <- x; mu <- quantile(x,c(.25,.75)); s <- rep(sd(x)/2,2); p <- .5; for(i in 1:500){ d1 <- p*dnorm(x,mu[1],s[1]); d2 <- (1-p)*dnorm(x,mu[2],s[2]); g <- d1/(d1+d2); p <- mean(g); mu <- c(sum(g*x)/sum(g), sum((1-g)*x)/sum(1-g)); s <- sqrt(c(sum(g*(x-mu[1])^2)/sum(g), sum((1-g)*(x-mu[2])^2)/sum(1-g))); s <- pmax(s,.02)}; ll <- sum(log(p*dnorm(x,mu[1],s[1])+(1-p)*dnorm(x,mu[2],s[2]))); list(ll=ll, mu=mu, s=s, p=p) }
chk <- function(x, lab){ e <- em2(x); n <- length(x); bic1 <- -2*sum(dnorm(x,mean(x),sd(x),log=TRUE))+2*log(n); bic2 <- -2*e$ll+5*log(n)
  data.frame(set=lab, n=n, min=min(x), max=max(x), mean=mean(x), sd=sd(x), skew=sk(x), kurtosis=ku(x), prop_below_5.8=mean(x<5.8), prop_above_7.2=mean(x>7.2), prop_one_decimal=dig1(x), prop_two_decimals=dig2(x), at_min=sum(x==min(x)), at_max=sum(x==max(x)), bic_1comp=bic1, bic_2comp=bic2, mix_means=paste(round(e$mu,2),collapse='/'), mix_weight1=round(e$p,2)) }
dc <- chk(main$ph,'all donors'); dc <- rbind(dc, chk(main$ph[main$sz==1],'SZ'), chk(main$ph[main$sz==0],'control'))
for(c in levels(main$cohort)) dc <- rbind(dc, chk(main$ph[main$cohort==c], paste0('cohort: ',c)))
rr <- r; dc$resid_shapiro_P <- NA; dc$resid_skew <- NA; dc$resid_kurtosis <- NA; dc$resid_n_abs_gt3 <- NA
dc$resid_shapiro_P[1] <- shapiro.test(rr)$p.value; dc$resid_skew[1] <- sk(rr); dc$resid_kurtosis[1] <- ku(rr); dc$resid_n_abs_gt3[1] <- sum(abs(rr)>3)
write.table(dc, file.path(out,'distribution_checks.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
ds <- do.call(rbind, lapply(split(main, main$study), function(x) data.frame(study=x$study[1], cohort=as.character(x$cohort[1]), n=nrow(x), prop_one_decimal=dig1(x$ph), prop_two_decimals=dig2(x$ph), min=min(x$ph), max=max(x$ph), sd=sd(x$ph))))
write.table(ds, file.path(out,'ph_digits_by_study.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
main$resid1 <- r; main$fitted1 <- fitted(f1); write.table(main[,c('person','cohort','study','group','ph','resid1','fitted1')], file.path(out,'residuals_model1.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
sink(file.path(out,'sessionInfo.txt')); print(sessionInfo()); sink()
cat('done\n')

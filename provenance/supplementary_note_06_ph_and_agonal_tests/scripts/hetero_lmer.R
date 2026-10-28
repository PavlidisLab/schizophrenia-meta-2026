# Heterogeneity of the SZ-minus-control brain pH difference. Usage: Rscript hetero_lmer.R persons.tsv outdir nsim
# persons.tsv = 81_C.tsv (variant B2, bank rule; persons are the unit). Primary model: ph ~ sz + (1+sz|cohort) + (1|sc), sc = study:cohort.
suppressMessages(library(lme4)); args <- commandArgs(TRUE); inp <- args[1]; out <- args[2]; nsim <- as.integer(args[3]); set.seed(20261016)
ctrl <- lmerControl(optimizer='bobyqa', optCtrl=list(maxfun=2e5), check.conv.singular='ignore')
d <- read.delim(inp, comment.char='#', stringsAsFactors=FALSE, na.strings=c('','NA')); m <- d[d$group %in% c('SZ','CTRL'),]
m$sz <- as.numeric(m$group=='SZ'); m$cohort <- factor(m$cohort); m$sc <- factor(paste(m$study,m$cohort,sep=':')); m$e14 <- as.numeric(m$in_14_dataset_regions==1)
cat('n', nrow(m), 'SZ', sum(m$sz), 'CTRL', sum(1-m$sz), 'cohorts', nlevels(m$cohort), 'studies', nlevels(m$sc), 'expr14', sum(m$e14), '\n')
res <- list(); add <- function(...) res[[length(res)+1]] <<- data.frame(..., stringsAsFactors=FALSE)
sdof <- function(fit, g, v='sz') { vc <- VarCorr(fit); if(!g %in% names(vc)) return(NA); s <- attr(vc[[g]],'stddev'); if(v %in% names(s)) unname(s[v]) else NA }
qci <- function(x) quantile(x, c(.025,.975), na.rm=TRUE)

# ---- primary (REML) and ML versions
P  <- lmer(ph ~ sz + (1+sz|cohort) + (1|sc), m, REML=TRUE, control=ctrl)
PM <- update(P, REML=FALSE)
cat('primary: est', fixef(P)['sz'], 'SD cohort slope', sdof(P,'cohort'), '\n')
add(item='primary_fixed_sz', est=fixef(P)['sz'], se=sqrt(vcov(P)['sz','sz']), stat=NA, df=NA, P=NA, lo=NA, hi=NA, note='REML')

# ---- (1a) random diagnosis slope by cohort: LRT vs no slope (ML). Models: slope variance + intercept-slope covariance = 2 df; chibar mixture reference 0.5*chi2(1)+0.5*chi2(2)
N0 <- lmer(ph ~ sz + (1|cohort) + (1|sc), m, REML=FALSE, control=ctrl)
lr <- 2*(as.numeric(logLik(PM))-as.numeric(logLik(N0)))
p2 <- pchisq(lr, 2, lower.tail=FALSE); pmix <- 0.5*pchisq(lr,1,lower.tail=FALSE)+0.5*pchisq(lr,2,lower.tail=FALSE)
add(item='LRT_cohort_slope_2df', est=NA, se=NA, stat=lr, df=2, P=p2, lo=NA, hi=NA, note=paste0('ML; naive chi2(2) P; boundary mixture 0.5chi2(1)+0.5chi2(2) P=', signif(pmix,4)))
# uncorrelated-slope version, 1 df, boundary 0.5*chi2(1)
PU <- lmer(ph ~ sz + (1|cohort) + (0+sz|cohort) + (1|sc), m, REML=FALSE, control=ctrl)
lr1 <- 2*(as.numeric(logLik(PU))-as.numeric(logLik(N0)))
add(item='LRT_cohort_slope_uncorrelated_1df', est=NA, se=NA, stat=lr1, df=1, P=0.5*pchisq(lr1,1,lower.tail=FALSE), lo=NA, hi=NA, note='ML; boundary 0.5*chi2(1)')
# parametric bootstrap LRT for the primary-vs-no-slope comparison
pb <- replicate(nsim, { y <- simulate(N0)[[1]]; mm <- m; mm$ph <- y
  a <- tryCatch(lmer(ph ~ sz + (1|cohort) + (1|sc), mm, REML=FALSE, control=ctrl), error=function(e) NULL)
  b <- tryCatch(lmer(ph ~ sz + (1+sz|cohort) + (1|sc), mm, REML=FALSE, control=ctrl), error=function(e) NULL)
  if(is.null(a)||is.null(b)) NA else 2*(as.numeric(logLik(b))-as.numeric(logLik(a))) })
add(item='LRT_cohort_slope_param_bootstrap', est=NA, se=NA, stat=lr, df=2, P=(1+sum(pb>=lr, na.rm=TRUE))/(1+sum(!is.na(pb))), lo=NA, hi=NA, note=paste0('parametric bootstrap nsim=', sum(!is.na(pb))))

# ---- (1b) SD of diagnosis effect across cohorts (primary): profile and bootstrap
pa <- tryCatch(confint(P, method='profile', quiet=TRUE), error=function(e) NULL)
vn <- names(attr(VarCorr(P)$cohort,'stddev'))
sigrows <- if(!is.null(pa)) rownames(pa)[grepl('^\\.sig', rownames(pa))] else character(0)
print(pa)
bs <- bootMer(P, function(x) c(fixef(x)['sz'], sdof(x,'cohort')), nsim=nsim, type='parametric', use.u=FALSE, seed=20261016)
bq <- qci(bs$t[,2])
add(item='SD_diagnosis_across_cohorts', est=sdof(P,'cohort'), se=NA, stat=NA, df=NA, P=NA, lo=if(length(sigrows)>=3) pa[sigrows[3],1] else NA, hi=if(length(sigrows)>=3) pa[sigrows[3],2] else NA, note=paste0('profile (row ', paste(sigrows,collapse=','), '); bootstrap 95% ', signif(bq[1],3), ' to ', signif(bq[2],3), '; boot share SD<0.01=', mean(bs$t[,2]<0.01,na.rm=TRUE)))
add(item='primary_fixed_sz_boot', est=fixef(P)['sz'], se=NA, stat=NA, df=NA, P=NA, lo=qci(bs$t[,1])[1], hi=qci(bs$t[,1])[2], note='parametric bootstrap')
prof <- tryCatch(confint(P, method='profile', parm='sz', quiet=TRUE), error=function(e) c(NA,NA))
add(item='primary_fixed_sz_profile', est=fixef(P)['sz'], se=NA, stat=NA, df=NA, P=NA, lo=prof[1], hi=prof[2], note='profile')

# ---- (1c) random diagnosis slope by STUDY (within cohort): add (0+sz|sc) to primary; and study slope without cohort slope
S1 <- lmer(ph ~ sz + (1+sz|cohort) + (1|sc) + (0+sz|sc), m, REML=TRUE, control=ctrl)  # cohort slope + extra study-within-cohort slope
S1M <- update(S1, REML=FALSE)
l1 <- 2*(as.numeric(logLik(S1M))-as.numeric(logLik(PM)))
add(item='LRT_study_slope_within_cohort_1df', est=sdof(S1,'sc'), se=NA, stat=l1, df=1, P=0.5*pchisq(l1,1,lower.tail=FALSE), lo=NA, hi=NA, note='ML; adds (0+sz|study:cohort) to primary; boundary 0.5*chi2(1); est = REML SD of study slope within cohort')
pa1 <- tryCatch(confint(S1, method='profile', quiet=TRUE), error=function(e) NULL); print(pa1)
sr <- if(!is.null(pa1)) rownames(pa1)[grepl('^\\.sig', rownames(pa1))] else character(0)
bs1 <- tryCatch(bootMer(S1, function(x) c(sdof(x,'cohort'), sdof(x,'sc')), nsim=nsim, type='parametric', use.u=FALSE, seed=20261016), error=function(e) NULL)
add(item='SD_study_slope_within_cohort', est=sdof(S1,'sc'), se=NA, stat=NA, df=NA, P=NA, lo=NA, hi=NA, note=paste0('profile table rows: ', paste(sr,collapse=','), '; bootstrap SD(study) 95% ', if(!is.null(bs1)) paste(signif(qci(bs1$t[,2]),3),collapse=' to ') else 'NA', '; SD(cohort slope) in this model ', signif(sdof(S1,'cohort'),3), ' (boot ', if(!is.null(bs1)) paste(signif(qci(bs1$t[,1]),3),collapse=' to ') else 'NA', ')'))
write.table(pa1, file.path(out,'profile_study_model.tsv'), sep='\t', quote=FALSE)
write.table(pa,  file.path(out,'profile_primary.tsv'), sep='\t', quote=FALSE)
# study slope only (all heterogeneity as study level, no cohort level)
S2 <- lmer(ph ~ sz + (1|cohort) + (1|sc) + (0+sz|sc), m, REML=TRUE, control=ctrl); S2M <- update(S2, REML=FALSE)
l2 <- 2*(as.numeric(logLik(S2M))-as.numeric(logLik(N0)))
pa2 <- tryCatch(confint(S2, method='profile', quiet=TRUE), error=function(e) NULL); write.table(pa2, file.path(out,'profile_studyonly_model.tsv'), sep='\t', quote=FALSE)
sr2 <- if(!is.null(pa2)) rownames(pa2)[grepl('^\\.sig', rownames(pa2))] else character(0)
add(item='LRT_study_slope_only_1df', est=sdof(S2,'sc'), se=NA, stat=l2, df=1, P=0.5*pchisq(l2,1,lower.tail=FALSE), lo=NA, hi=NA, note=paste0('ML; (0+sz|study:cohort) with no cohort slope vs no slope; boundary; profile rows ', paste(sr2,collapse=',')))

# ---- (1d) fixed-effect diagnosis-by-cohort interaction (ML), study intercepts random
F0 <- lmer(ph ~ sz + cohort + (1|sc), m, REML=FALSE, control=ctrl); F1 <- lmer(ph ~ sz*cohort + (1|sc), m, REML=FALSE, control=ctrl)
a <- anova(F0,F1); add(item='LRT_fixed_sz_x_cohort', est=NA, se=NA, stat=a$Chisq[2], df=a$Df[2], P=a$`Pr(>Chisq)`[2], lo=NA, hi=NA, note='ML; ph ~ sz+cohort+(1|study:cohort) vs sz*cohort')
# fixed-effects LS, study as fixed effect: sz x cohort F test, and sz x study(within cohort) F test
l0 <- lm(ph ~ sc + sz, m); l1c <- lm(ph ~ sc + sz*cohort, m); l1s <- lm(ph ~ sc + sz:sc, m)
aa <- anova(l0,l1c); add(item='F_LS_sz_x_cohort', est=NA, se=NA, stat=aa$F[2], df=paste0(aa$Df[2],',',aa$Res.Df[2]), P=aa$`Pr(>F)`[2], lo=NA, hi=NA, note='ordinary LS, study:cohort fixed; sz x cohort (6 df)')
l1s <- lm(ph ~ sc + sz*cohort + sz:sc, m)  # sz:sc partly aliased with sz*cohort; anova gives the study-within-cohort remainder
ab <- anova(l1c, lm(ph ~ sc + sz:sc, m)); add(item='F_LS_sz_x_study_within_cohort', est=NA, se=NA, stat=ab$F[2], df=paste0(ab$Df[2],',',ab$Res.Df[2]), P=ab$`Pr(>F)`[2], lo=NA, hi=NA, note='ordinary LS; extra diagnosis-by-study terms beyond diagnosis-by-cohort (sc-specific slopes vs cohort-specific slopes)')

# ---- (3) expr14 vs other
E0 <- lmer(ph ~ sz + e14 + (1+sz|cohort) + (1|sc), m, REML=FALSE, control=ctrl); E1 <- lmer(ph ~ sz*e14 + (1+sz|cohort) + (1|sc), m, REML=FALSE, control=ctrl)
a <- anova(E0,E1); E1r <- update(E1, REML=TRUE)
pe <- tryCatch(confint(E1r, method='profile', parm='sz:e14', quiet=TRUE), error=function(e) c(NA,NA))
add(item='expr14_sz_x_e14_interaction_primary_structure', est=fixef(E1r)['sz:e14'], se=sqrt(vcov(E1r)['sz:e14','sz:e14']), stat=a$Chisq[2], df=a$Df[2], P=a$`Pr(>Chisq)`[2], lo=pe[1], hi=pe[2], note='difference (expr14 minus other) in SZ-minus-control; profile CI; LRT ML; random slope by cohort retained')
bse <- bootMer(E1r, function(x) fixef(x)['sz:e14'], nsim=nsim, type='parametric', use.u=FALSE, seed=20261016)
add(item='expr14_interaction_bootstrap', est=fixef(E1r)['sz:e14'], se=NA, stat=NA, df=NA, P=NA, lo=qci(bse$t[,1])[1], hi=qci(bse$t[,1])[2], note='parametric bootstrap')
add(item='expr14_sz_effect_in_other_(sz_main)', est=fixef(E1r)['sz'], se=sqrt(vcov(E1r)['sz','sz']), stat=NA, df=NA, P=NA, lo=NA, hi=NA, note='SZ effect when e14=0 in the interaction model')
# no-random-slope variant of the interaction (cohort intercept only)
G0 <- lmer(ph ~ sz + e14 + (1|cohort) + (1|sc), m, REML=FALSE, control=ctrl); G1 <- lmer(ph ~ sz*e14 + (1|cohort) + (1|sc), m, REML=FALSE, control=ctrl); ag <- anova(G0,G1)
add(item='expr14_interaction_no_random_slope', est=fixef(G1)['sz:e14'], se=sqrt(vcov(G1)['sz:e14','sz:e14']), stat=ag$Chisq[2], df=ag$Df[2], P=ag$`Pr(>Chisq)`[2], lo=NA, hi=NA, note='ML fit, Wald SE; cohort random intercept only')
# study-fixed LS: within-study comparison
w <- lm(ph ~ sc + sz*e14, m); cw <- confint(w)['sz:e14',]; sw <- summary(w)$coefficients['sz:e14',]
add(item='expr14_interaction_LS_study_fixed', est=sw[1], se=sw[2], stat=sw[3], df=w$df.residual, P=sw[4], lo=cw[1], hi=cw[2], note='ordinary LS, study:cohort fixed effect, sz*e14; t test')
# separate-fit two-group comparison (independent estimates, primary model refit on each set)
fit_sub <- function(x){ f <- lmer(ph ~ sz + (1+sz|cohort) + (1|sc), x, REML=TRUE, control=ctrl); ci <- tryCatch(confint(f, method='profile', parm='sz', quiet=TRUE), error=function(e) c(NA,NA)); c(est=unname(fixef(f)['sz']), se=sqrt(vcov(f)['sz','sz']), lo=ci[1], hi=ci[2], n=nrow(x), nSZ=sum(x$sz), nC=sum(1-x$sz), nco=nlevels(droplevels(x$cohort))) }
fe <- fit_sub(m[m$e14==1,]); fo <- fit_sub(m[m$e14==0,])
dd <- unname(fe['est']-fo['est']); sd <- sqrt(fe['se']^2+fo['se']^2)
add(item='expr14_separate_fit_expr14', est=fe['est'], se=fe['se'], stat=NA, df=NA, P=NA, lo=fe['lo'], hi=fe['hi'], note=paste0('n=',fe['n'],' SZ=',fe['nSZ'],' CTRL=',fe['nC'],' cohorts=',fe['nco']))
add(item='expr14_separate_fit_other', est=fo['est'], se=fo['se'], stat=NA, df=NA, P=NA, lo=fo['lo'], hi=fo['hi'], note=paste0('n=',fo['n'],' SZ=',fo['nSZ'],' CTRL=',fo['nC'],' cohorts=',fo['nco']))
add(item='expr14_separate_fit_difference', est=dd, se=sd, stat=dd/sd, df=NA, P=2*pnorm(-abs(dd/sd)), lo=dd-1.96*sd, hi=dd+1.96*sd, note='expr14 minus other; independent Wald z from the two refits (Wald SEs)')
# cohort-level slope SDs in each subset already come with fits; also the primary-model sz within cohort x e14 is not estimable in all cohorts (Harvard has no expr14 persons; multi-bank no other)
r <- do.call(rbind, res); write.table(r, file.path(out,'hetero_lmer_results.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
writeLines(capture.output(sessionInfo()), file.path(out,'sessionInfo.txt')); print(r); cat('done\n')

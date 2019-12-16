library("meta")

# dataframe df of wide format, e.g. contains the beta's and standard errors for the various studies in 1 row
#
# SNP EUR.BETA EUR.SE AFR.BETA AFR.SE
# rs1     0.30   0.03     0.25   0.04
# rs2     0.12   0.03     0.08   0.02
# rsx     0.03   0.01     0.08   0.01

# meta-analysis of EUR and AFR
for(i in 1:nrow(df)) {
  df$pooledEffect[i] = metagen(c(df$EUR.BETA[i], df$AFR.BETA[i]), c(df$EUR.SE[i], df$AFR.SE[i]), sm = "RD")$TE.fixed[[1]]
  df$pooledSE[i]     = metagen(c(df$EUR.BETA[i], df$AFR.BETA[i]), c(df$EUR.SE[i], df$AFR.SE[i]), sm = "RD")$seTE.fixed[[1]]
  df$pooledLower[i]  = metagen(c(df$EUR.BETA[i], df$AFR.BETA[i]), c(df$EUR.SE[i], df$AFR.SE[i]), sm = "RD")$lower.fixed[[1]]
  df$pooledUpper[i]  = metagen(c(df$EUR.BETA[i], df$AFR.BETA[i]), c(df$EUR.SE[i], df$AFR.SE[i]), sm = "RD")$upper.fixed[[1]]
  df$pooledLower[i]  = metagen(c(df$EUR.BETA[i], df$AFR.BETA[i]), c(df$EUR.SE[i], df$AFR.SE[i]), sm = "RD")$lower.fixed[[1]]
  df$pooledZ[i]      = metagen(c(df$EUR.BETA[i], df$AFR.BETA[i]), c(df$EUR.SE[i], df$AFR.SE[i]), sm = "RD")$zval.fixed[[1]]
  df$pooledP[i]      = metagen(c(df$EUR.BETA[i], df$AFR.BETA[i]), c(df$EUR.SE[i], df$AFR.SE[i]), sm = "RD")$pval.fixed[[1]]
  # optional 
  # if output is from logistic regression, not linear
  # df$pooledOR[i]    = exp(df$pooledEffect[i])
}

# make sure that you are meta-analyzing beta's and not OR's, because the standard error is for the beta not the OR.

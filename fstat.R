#CALCULATE F STATISTIC: https://www.thelancet.com/cms/10.1016/j.ebiom.2021.103397/attachment/8cf16507-a4e4-444b-9ab3-0e36866ebf3e/mmc1.pdf

df$EAF2 <- (1 - check_dat$eaf.exposure)
df$MAF <- pmin(check_dat$eaf.exposure, check_dat$EAF2)
PVEfx <- function(BETA, MAF, SE, N){
  pve <- (2*(BETA^2) * MAF * (1 - MAF))/ ((2*(BETA^2) * MAF * (1 - MAF)) + ((SE^2) * 2 * N * MAF * (1 - MAF)))
  return(pve)
}

df$PVE <- mapply(PVEfx,
                        df$beta.exposure,
                        df$MAF,
                        df$se.exposure,
                        df$samplesize.exposure)

check_dat$FSTAT <- ((df$samplesize.exposure - 1 - 1)/1)*(df$PVE/(1 - check_dat$PVE))

write.csv(df, file ="output.csv")

#Total instrument F statistic -> will need to assume the largest sample size in the GWAS
((184765- 6 - 1)/6)*(0.00186892/(1 - 0.00186892)) = 57.65741

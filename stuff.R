args = commandArgs(trailingOnly = TRUE)

# continous to ordinal
ordinize = function(x, n.cat) {
  cut(x, breaks = c(quantile(x, probs = seq(0, 1, 1/n.cat))), labels = seq(1:n.cat))
}

# stating the obvious
# bgzip -c file.vcf > file.vcf.gz
# tabix -p vcf file.vcf.gz                   

# don't restrict R rounding of p-value to 2.2E-16 (default: 53 binary digits accuracy)
.Machine$double.eps # 2.220446e-16
.Machine$double.xmin # 2.225074e-308
# soluation
.Machine$double.eps = .Machine$double.xmin

#  pvalue to standard error
df$SE = abs(df$BETA) / abs(qnorm(1 - df$P/2))

# formatting of regression model output (OR, CI, P)
glm.fit <- glm(OUT ~ SNP + AGE + SEX + PCA1 + PCA2 + PCA3, data = df, family = "binomial")
round(cbind(OR = exp(coef(glm.fit)), exp(confint(glm.fit)), P = summary(glm.fit)$coef[, "Pr(>|z|)"]), digits = 3)

# systematically rename column names
# remove everything after _ in column name
names(df) = gsub(pattern = "_*", replacement = "", x = names(df))

# remove everyting after dash (-) 
# e.g. 1:1520-1520 > 1:1520
df$CHRCBP = gsub(pattern = "(-).*", replacement = "", x = df$Location)

# keep position (e.g. 1520)
df$BP = gsub(pattern = "*.(-)", replacement = "", x = df$CHRCBP)

# keep chr (e.g. 1) 
df$CHR = gsub(pattern = "(:).*", replacement = "", x = df$CHRCBP)

# alternatively
tmp = as.data.frame(do.call(rbind, strsplit(df$variant, '\\:')))
colnames(tmp) = c("chr", "pos", "ref", "alt")
df = cbind(df, tmp)

# take multiple rows and merge them into one variable (long into single)
df.new = aggregate( . ~ grouping_variable, df, function(x) toString(unique(x)))                               

# adjust for genomic inflation
# remove outliers first (p < 0.01)
df$X    = qchisq(1 - in.data$P, 1)
LAMBDA  = median(in.data$stats) / 0.4549
df$Xadj = df$X / LAMBDA
#df$Padj = calculate p from new ChiSquare

# find best-matched control for a case based on age, gender and PC's (1x)
library('e1071')
myMatch   = matchControls(DISEASE ~ AGE + SEX + PCA1 + PCA2 + PCA3, caselabel = 1, contlabel = 0, data = df, replace = F)
dfMatched = rbind(df[myMatch$cases, ], df[myMatch$controls, ])
                   
# support vector machines and such (in train and test dataset)
library('caret')
fitControl = trainControl(method = "repeatedcv", number = 5, repeats = 5)
svm.fit    = train(OUTCOME ~ PCA1 + PCA2 + PCA3, method = 'svmRadial', trControl = fitControl, data = df.train)
svm.pred   = predict(svm.fit, df.test)
df.test    = cbind(df.test, svm.pred)

# collapse multiple rows into one field
paste(unique(unlist(df$variable_name)), collapse = ';')

# time
t1 = as.POSIXlt(strptime("06.21.2018 9:00", "%m.%d.%Y %H:%M"))
t2 = as.POSIXlt(strptime("06.21.2018 10:15", "%m.%d.%Y %H:%M"))
t1-t2
difftime(t1, t2, units = 'mins')
as.numeric(difftime(t1, t2, units = 'mins'))
                   
# split column by delimiter (unequal lenghts)
require(data.table)
dt <- data.table(df)
#            ids gene
# 1: 266;372;572  ABBA
# 2: 908;202;896  BETA
#
# expand the patients
dt[ , list( patient = unlist(strsplit(ids, ";" ))) , by = gene]
#    gene patient
# 1: ABBA     266
# 2: ABBA     372
# 3: ABBA     572
# 4: BETA     908
# 5: BETA     202
# collapse the genes 
paste(unlist(df[which(df$patient %in% random.list), "gene"]), collapse = ";")

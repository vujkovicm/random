# don't restrict R rounding of p-value to 2.2E-16 (default: 53 binary digits accuracy)
.Machine$double.eps # 2.220446e-16
.Machine$double.xmin # 2.225074e-308
# soluation
.Machine$double.eps = .Machine$double.xmin

# systematically rename column names
# remove everything after _ in column name
names(df) = gsub(pattern = "_*", replacement = "", x = names(df))


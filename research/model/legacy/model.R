#%%
library(fixest)

data <- read.csv("Yield_Climate/research/model/kanto.csv", header = TRUE)
head(data)

# yieldを数値に変換し、NAを除外
data$yield <- as.numeric(data$yield)

# paddy_ratioを固定
data <- data[!is.na(data$yield) & data$year >= 2009, ]

#%%
# 固定効果モデル (TWFE)
model <- fixest::feols(
  yield ~ APCP + GSR + heat + heat * paddy_ratio | year + city_id,
  data    = data,
  cluster = ~city_id
)

etable(model,
       title  = "TWFE",
       digits = 3)

#%%
#
model_log <- fixest::feols(
  log(yield) ~ APCP + GSR + heat + heat * paddy_ratio | year + city_id,
  data    = data,
  cluster = ~city_id
)

etable(model_log,
       title  = "TWFE(log)",
       digits = 3)

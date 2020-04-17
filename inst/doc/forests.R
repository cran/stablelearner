## ---- GermanCredit------------------------------------------------------------
data("GermanCredit", package = "evtree")

## ---- GermanCreditstr1--------------------------------------------------------
set.seed(2409)
dat <- droplevels(GermanCredit[sample(1000, size = 500), ])

## ---- GermanCreditstr2, eval = FALSE------------------------------------------
#  str(dat)

## ---- GermanCreditstr3, echo = FALSE------------------------------------------
str(dat, width = 80, strict.width = "cut")

## ---- ctree-------------------------------------------------------------------
set.seed(2906)
ct_partykit <- partykit::ctree(credit_risk ~ ., data = dat,
  control = partykit::ctree_control(mtry = 5, teststat = "quadratic",
    testtype = "Univariate", mincriterion = 0, saveinfo = FALSE))

## ---- stablelearner_cforest---------------------------------------------------
set.seed(2907)
cf_stablelearner <- stablelearner::stabletree(ct_partykit, sampler = stablelearner::subsampling,
  savetrees = TRUE, B = 100, v = 0.632)

## ---- stablelearner_methods---------------------------------------------------
summary(cf_stablelearner, original = FALSE)

## ---- stablelearner_barplot, fig.height = 4, fig.width = 6--------------------
barplot(cf_stablelearner, original = FALSE, cex.names = 0.7)

## ---- stablelearner_image, fig.height = 4, fig.width = 6----------------------
image(cf_stablelearner, original = FALSE, cex.names = 0.7)

## ---- stablelearner_plot, fig.height = 12, fig.width = 6----------------------
plot(cf_stablelearner, original = FALSE,
  select = c("status", "employment_duration", "duration"))

## ---- stablelearner_trees, eval = FALSE---------------------------------------
#  cf_stablelearner$tree[[1]]

## ---- partykit_cforest, eval = FALSE------------------------------------------
#  set.seed(2908)
#  cf_partykit <- partykit::cforest(credit_risk ~ ., data = GermanCredit,
#    ntree = 100, mtry = 5)

## ---- as.stabletree_partykit, eval = FALSE------------------------------------
#  cf_partykit_st <- stablelearner::as.stabletree(cf_partykit)
#  summary(cf_partykit_st, original = FALSE)
#  barplot(cf_partykit_st, cex.names = 0.7)
#  image(cf_partykit_st, cex.names = 0.7)
#  plot(cf_partykit_st, select = c("status", "employment_duration", "duration"))

## ---- party_cforest, eval = FALSE---------------------------------------------
#  set.seed(2909)
#  cf_party <- party::cforest(credit_risk ~ ., data = dat,
#    control = party::cforest_unbiased(ntree = 100, mtry = 5))

## ---- as.stabletree_party, eval = FALSE---------------------------------------
#  cf_party_st <- stablelearner::as.stabletree(cf_party)
#  summary(cf_party_st, original = FALSE)
#  barplot(cf_party_st, cex.names = 0.7)
#  image(cf_party_st, cex.names = 0.7)
#  plot(cf_party_st, select = c("status", "employment_duration", "duration"))

## ---- randomForest------------------------------------------------------------
set.seed(2910)
rf <- randomForest::randomForest(credit_risk ~ ., data = dat,
  ntree = 100, mtry = 5)

## ---- as.stabletree_randomForest----------------------------------------------
rf_st <- stablelearner::as.stabletree(rf)
summary(rf_st, original = FALSE)

## ---- rf_barplot, fig.height = 4, fig.width = 6-------------------------------
barplot(rf_st, cex.names = 0.7)

## ---- rf_image, fig.height = 4, fig.width = 6---------------------------------
image(rf_st, cex.names = 0.7)

## ---- rf_plot, fig.height = 12, fig.width = 6---------------------------------
plot(rf_st, select = c("status", "employment_duration", "duration"))

## ---- ranger, eval = FALSE----------------------------------------------------
#  set.seed(2911)
#  rf_ranger <- ranger::ranger(credit_risk ~ ., data = dat,
#    num.trees = 100, mtry = 5)

## ---- as.stabletree_ranger, eval = FALSE--------------------------------------
#  rf_ranger_st <- stablelearner::as.stabletree(rf_ranger)
#  summary(rf_ranger_st, original = FALSE)
#  barplot(rf_ranger_st, cex.names = 0.7)
#  image(rf_ranger_st, cex.names = 0.7)
#  plot(rf_ranger_st, select = c("status", "employment_duration", "duration"))


#         Author: Pedro Salas Rojo
#         Date: 03/04/2021
#         Dataset: PGM
#         Name of project: Select circumstances with forest / Estimate IOp in many ways

rm(list = ls(all.names = TRUE)) 
library(tidyverse)
library(haven)
library(stats)
library(spatstat)
library(dineq)
library(radiant.data)
library(Hmisc)
library(partykit)
library(party)
library(glmnet)
library(glmnetUtils)
library(caret)
library(mice)
library(xlsx)
set.seed(1)

# Functions ----

# Plot the importance of variables in forests
implot <- function(x) {
  cforest_importance <<- v <- varimp(x, OOB = TRUE)
  dotchart(v[order(v)])
}

# Load database and check there are no NA's----
iop <- read_csv("C:/Users/pedro/Google Drive/mispapers/missing/working/comp_inc_18.csv")

summary(is.na(iop))

# Types creation ----

iop <- iop %>%
  mutate(types = group_indices(iop, sesso, highoc, highed, region)) %>%
  select(redtot, coeff, sesso, highoc, highed, reg, region, types)

length(unique(iop$types))
summary(is.na(iop))

# IOp Checchi and Peragine non parametric ----

iop <- iop %>%
  group_by(types) %>%
  mutate(cyp_np = weighted.mean(redtot, coeff)) %>%
  ungroup(types)

cp <- iop %>%
  summarise(method = "Checchi and Peragine (2010)",
            gini = gini.wtd(redtot, coeff),
            abs_iop_gini = gini.wtd(cyp_np, coeff),
            rel_iop_gini = 100*abs_iop_gini/gini,
            mld = mld.wtd(redtot, coeff),
            abs_iop_mld = mld.wtd(cyp_np, coeff),
            rel_iop_mld = 100*abs_iop_mld/mld,
            types_used = length(unique(iop$cyp_np)))

tab <- rbind(cp)

# IOp Ferreira and Guignoux parametric ----

iop <- iop %>%
  mutate(gyf_p = predict(lm(redtot ~ factor(sesso) + factor(highed) +
                            factor(highoc) + factor(region), weights=coeff)))

cp <- iop %>%
  summarise(method = "Ferreira and Guignoux (2011)",
            gini = gini.wtd(redtot, coeff),
            abs_iop_gini = gini.wtd(gyf_p, coeff),
            rel_iop_gini = 100*abs_iop_gini/gini,
            mld = mld.wtd(redtot, coeff),
            abs_iop_mld = mld.wtd(gyf_p, coeff),
            rel_iop_mld = 100*abs_iop_mld/mld,
            types_used = length(unique(iop$gyf_p)))

tab <- rbind(tab, cp) 

# IOp Regularized Regression (ML parametric)  -----

mat <- iop %>%
  dplyr::select(sesso, highoc, highed, region)
vec <- data.matrix(mat)

# Train and obtain the tuned parameters in "a" (alpha and lambda)
cv5 <- caret::trainControl(method = "cv", number = 5)

el_train <- caret::train(redtot ~ factor(sesso) + factor(highed) +
                           factor(highoc) + factor(region), data = iop, method = "glmnet", trControl = cv5)
print(el_train)
plot(el_train)
a <- el_train$bestTune

# Get imputed vector and estimate IOp (note that vec requires 2+ columns)

iop <- iop %>%
  mutate(en_p = predict(glmnet(vec, iop$redtot, lambda = a$lambda, alpha=a$alpha), newx = vec))

cp <- iop %>%
  summarise(method = "Regularized Regression (Brunori et al 2021)",
            gini = gini.wtd(redtot, coeff),
            abs_iop_gini = gini.wtd(en_p, coeff),
            rel_iop_gini = 100*abs_iop_gini/gini,
            mld = mld.wtd(redtot, coeff),
            abs_iop_mld = mld.wtd(en_p, coeff),
            rel_iop_mld = 100*abs_iop_mld/mld,
            types_used = length(unique(iop$en_p)))

tab <- rbind(tab, cp) 

# IOp Trees (ML non parametric) ----

# Train and obtain the best parameters (1 - alpha = a). 
cv5 <- trainControl(method = "cv", number = 5)

tr_train <- caret::train(redtot ~ sesso + highed + highoc + region, data = iop,
                    method = "ctree", trControl = cv5)
print(tr_train)
plot(tr_train)
a <- tr_train$bestTune

# Run tree with tunned alpha

detach("package:party", unload = TRUE)
tree <- partykit::ctree(redtot ~ sesso + highed + highoc + region, 
                  data = iop, control = ctree_control(testtype = c("Bonferroni"), 
                                           teststat = c("quad"), mincriterion = as.numeric(a)))

iop$tr_np <- predict(tree, type="response")
iop$types <- predict(tree, type="node")
ct_node <- as.list(tree$node)

pred <- aggregate(predict(tree),list(predict(tree, type = "node")), FUN = mean)
qi<-pred

for (t in 1:length(qi[,1])){
  typ<-as.numeric(names(table(iop$types)[t]))
  qi[t,2]<-length(iop$types[iop$types==typ])/length(iop$types) 
}
qi

for(t in 1:nrow(pred)) {
  ct_node[[pred[t,1]]]$info$prediction <- as.numeric(paste(format(round(pred[t, -1], digits = 2), nsmall = 2)))
  ct_node[[pred[t,1]]]$info$nobs       <- as.numeric(paste(format(round(qi[t, -1]  , digits = 3), nsmall = 2)))
}

tree$node <- as.partynode(ct_node)
as.list(ct_node)

plot(tree,  terminal_panel=node_terminal, tp_args = list(FUN = function(node) 
  c("Exp. outcome",node$prediction, "Pop. Share", node$nobs)))

cp <- iop %>%
  summarise(method = "Trees (Brunori et al 2019)",
            gini = gini.wtd(redtot, coeff),
            abs_iop_gini = gini.wtd(tr_np, coeff),
            rel_iop_gini = 100*abs_iop_gini/gini,
            mld = mld.wtd(redtot, coeff),
            abs_iop_mld = mld.wtd(tr_np, coeff),
            rel_iop_mld = 100*abs_iop_mld/mld,
            types_used = length(unique(iop$tr_np)))

tab <- rbind(tab, cp)

# IOp Forests (ML non parametric) ----

# Tuning forests takes **A LOT** It delivers mtry=3
#cv5 <- trainControl(method = "cv", number = 5)
#
#rf_train <- train(redtot ~ sesso + highed + highoc + region, data = iop,
#                  method = "cforest", trControl = cv5)
#print(rf_train)
#plot(rf_train)
#a <- rf_train$bestTune
a <- 3

library(party)
forest <- cforest(redtot ~ sesso + highoc + highed + region,
                  data = iop, controls=cforest_control(ntree=50L, mtry = a, 
                                                       trace = TRUE, mincriterion = 0))
iop$rf_np <- predict(forest, type="response")

cp <- iop %>%
  summarise(method = "Forests (Brunori et al 2019)",
            gini = gini.wtd(redtot, coeff),
            abs_iop_gini = gini.wtd(rf_np, coeff),
            rel_iop_gini = 100*abs_iop_gini/gini,
            mld = mld.wtd(redtot, coeff),
            abs_iop_mld = mld.wtd(rf_np, coeff),
            rel_iop_mld = 100*abs_iop_mld/mld,
            types_used = NA)

tab <- rbind(tab, cp)

importance <- as.table(varimp(forest, conditional = TRUE))

importance
plot(importance)

# We can get a bit deeper into the variable importance (time consuming, it took more than an hour),
# but we get similar results

# cvcont <- trainControl(method="cv", number = 5,
#                        allowParallel=TRUE)
# 
# train <- train(redtot ~ sesso + highoc + highed + reg,
#                data=iop, method = "rf", trControl = cvcont,
#                importance = TRUE)
#
# plot(varImp(train))

# Summary of IOp ----
tab

write.xlsx(as.data.frame(tab), file="C:/Users/pedro/Google Drive/mispapers/missing/results/iop_2018.xlsx", sheetName = "results", append=TRUE, col.names=TRUE, row.names=TRUE)


###!!!!!!!!!!!######################!!!!!!!!!!!###
###!!!!!!!!!!!####protein.RData#####!!!!!!!!!!!###
###!!!!!!!!!!!######################!!!!!!!!!!!###

# loading protein data. Data
load("C:/Users/Igor/Desktop/Bioinformatyka I rok II stopnia/II sem/SAD1/projekt zaliczeniowy/Grupa5/protein.RData")

df.train.protein <- as.data.frame(data.train)
df.test.protein <- as.data.frame(data.test)


min(df.train.protein$Y)
max(df.train.protein$Y)

plot(df.train.protein$Y)

# using the histogram I am trying to detect if there are variables whose unique 
# values are rare (below 10), what would indicate the existence of qualitative variables
unique_values <- apply(df.train.protein, 2, function(x){length(unique(x))})
hist(unique_values)

# the histogram shows that there are variable, where numbers of unique variables less than 10 (exactly 2)

# from the above I replace these variables as qualitative variables (as.factor), 
# numerical variables, continuous scale, 
# except for the explanatory variable 'Y'


for(i in 1:dim(df.train.protein)[2]){
  if(length(unique(df.train.protein[,i])) < 3){
    df.train.protein[,i][df.train.protein[,i] < 0] <- 0
    df.train.protein[,i][df.train.protein[,i] > 0] <- 1
    #df.train.protein[,i] <- as.factor(df.train.protein[,i])
    #df.train.protein[,i] <- as.numeric(levels(df.train.protein[,i]))[df.train.protein[,i]]
    }
  else if(all(df.train.protein[,i] == df.train.protein$Y)){df.train.protein[,i] <- df.train.protein[,i]}
  else df.train.protein[,i] <- scale(df.train.protein[,i])
}

for(i in 1:dim(df.test.protein)[2]){
  if(length(unique(df.test.protein[,i])) < 3){
    df.test.protein[,i][df.test.protein[,i] < 0] <- 0
    df.test.protein[,i][df.test.protein[,i] > 0] <- 1
    #df.train.protein[,i] <- as.factor(df.train.protein[,i])
    #df.train.protein[,i] <- as.numeric(levels(df.train.protein[,i]))[df.train.protein[,i]]
  }
  else if(all(df.test.protein[,i] == df.test.protein$Y)){df.test.protein[,i] <- df.test.protein[,i]}
  else df.test.protein[,i] <- scale(df.test.protein[,i])
}

# Data are too large to make any meaningful model on them, therefore
# I will analyze the correlation of explanatory variables with the explained variable and
# I keep only these explanatory variables for further analysis,
# which have a high correlation with the explained variable

protein.correlation <- df.train.protein[1,1:ncol(df.train.protein)]

for(i in 1:ncol(df.train.protein)){
  protein.correlation[,i] <- cor(x = as.numeric(df.train.protein[,i]),y = df.train.protein$Y)
}

a <- seq(0.01,0.9, by = 0.01)
plot_of_significance <- as.data.frame(matrix(,nrow = 2,ncol = length(a)))

for(b in 1:length(a)){
  correlation_good <- protein.correlation[,which(abs(protein.correlation) >= a[b])]
  plot_of_significance[1,b] <- a[b]
  plot_of_significance[2,b] <- dim(correlation_good)[2]
}

regulator <- 0.05

df.train.protein.correlated <- df.train.protein[,which(abs(protein.correlation) >= regulator)]
df.train.protein <- df.train.protein.correlated

###VIF###
#install.packages('car')
library(car)

protein_model <- lm(Y~.,data = df.train.protein.correlated)
protein.vif <- vif(protein_model)
hist(protein.vif)


#############################################################################################

sample.rows <- sample(x = 1:nrow(df.train.protein), size = floor(nrow(df.train.protein) * 0.8))
train.protein <- df.train.protein[sample.rows,]
valid.protein <- df.train.protein[-sample.rows,]
test.protein <- as.data.frame(data.test)



####################################################################
#CHOOSING THE MODEL USING GREEDY ALGORITHM TO CHOOSE PREDICTORS#####
#ASSESMENT METHOD: TEST-F###########################################
####################################################################

# in the next iterations I add one predictor at the time from the pool of all predictors
# check if the model with the added predictor has all relevant predictors
# if yes then
# we check if any variable significantly improves the model:
# we remove the first item because it is NA
library(MASS)

null.model <- lm(Y ~ 1, train.protein)
full.formula <- terms(Y ~., data = train.protein)
full.formula <- formula(full.formula)
f.terms <- addterm(null.model, full.formula, test = 'F')
model.protein.F <- null.model

alpha <- 0.01

f.terms <- addterm(model.protein.F, full.formula, test = 'F')
while(any(f.terms$'Pr(F)'[-1] < alpha )){
  new.term <- rownames(f.terms)[which.max(f.terms$'F Value')]
  new.formula <- paste('.~.+', new.term)
  model.protein.F <- update(model.protein.F, new.formula)
  f.terms <-addterm(model.protein.F, full.formula, test = 'F')
}

summary(model.protein.F)
F.coeficients <- coef(model.protein.F)
plot(F.coeficients)

F.coeficients[abs(F.coeficients)>2]

################################################################
#CHOOSING THE MODEL USING GREEDY ALGORITHM TO CHOOSE PREDICTORS#
#ASSESMENT METHOD: INFORMATION CRITERION BIC####################
################################################################

full.model <- lm(Y~., data = train.protein)
null.model <- lm(Y~1, data = train.protein)
model.protein.BIC <- stepAIC(null.model, k=log(nrow(train.protein)), direction = "forward", scope=list(upper=full.model,lower=null.model))

summary(model.protein.BIC)
BIC.coeficient <- coef(model.protein.BIC)

BIC.coeficient[abs(BIC.coeficient)>2]

################################################################
#CHOOSING THE MODEL USING ELASTIC NET AND LASSO#################
################################################################
################################################################
#install.packages("glmnet")
library(glmnet)

data.protein.train.df.lasso <- as.data.frame(df.train.protein)
df.protein.train.lasso <- data.protein.train.df.lasso[sample.rows,]
df.protein.valid.lasso <- data.protein.train.df.lasso[-sample.rows,]


dims.protein.train <- dim(df.protein.train.lasso)
dims.protein.valid <- dim(df.protein.valid.lasso)

X.train.protein <- df.protein.train.lasso[,-dims.protein.train[2]]
Y.train.protein <- df.protein.train.lasso[,dims.protein.train[2]]
X.valid.protein <- df.protein.valid.lasso[,-dims.protein.valid[2]]
Y.valid.protein <- df.protein.valid.lasso[,dims.protein.valid[2]]

X.train.protein <- as.matrix(X.train.protein)
Y.train.protein <- as.matrix(Y.train.protein)
X.valid.protein <- as.matrix(X.valid.protein)
Y.valid.protein <- as.matrix(Y.valid.protein)

#################
###ELASTIC NET###
#################
model.protein.ELASTIC.NET <- cv.glmnet(X.train.protein, Y.train.protein, nfolds = 10, standardize = F, intercept = T, thresh = 1e-16)

optimal.betas.protein.lambda.min <- coef(model.protein.ELASTIC.NET, s = "lambda.min")
optimal.betas.protein.lambda.1se <- coef(model.protein.ELASTIC.NET, s = "lambda.1se")
optimal.betas.protein.lambda.min 
optimal.betas.protein.lambda.1se

betas.protein.names <- names(df.protein.train.lasso)
beta.protein.table <- as.data.frame(cbind("LAMBDA_MIN" = optimal.betas.protein.lambda.min[,1], "LAMBDA_1SE" = optimal.betas.protein.lambda.1se[,1]))
beta.protein.table <- data.frame("NAMES" = row.names(beta.protein.table), beta.protein.table )
beta.protein.table <- beta.protein.table[which(beta.protein.table$LAMBDA_MIN != 0 | beta.protein.table$LAMBDA_1SE != 0),]
betas.protein.1se <- as.data.frame(as.character(beta.protein.table$NAMES[which(beta.protein.table$LAMBDA_1SE != 0)]))
betas.protein.min <- as.data.frame(as.character(beta.protein.table$NAMES[which(beta.protein.table$LAMBDA_MIN != 0)]))

model.protein.string.1se = betas.protein.1se[2,1]
for(i in 3:nrow(betas.protein.1se)){
  model.protein.string.1se = paste(model.protein.string.1se, betas.protein.1se[i,1], sep = " + ")
}

model.protein.string.min = betas.protein.min[2,1]
for(i in 3:nrow(betas.protein.min)){
  model.protein.string.min = paste(model.protein.string.min, betas.protein.min[i,1], sep = " + ")
}

model.protein.ELASTIC.NET.string.1se <- paste("Y ~ ",model.protein.string.1se, sep = "")
model.protein.ELASTIC.NET.string.min <- paste("Y ~ ",model.protein.string.min, sep = "")

model.protein.ELASTIC.NET.lambda_1SE <- lm(as.formula(model.protein.ELASTIC.NET.string.1se), data = df.protein.train.lasso)
model.protein.ELASTIC.NET.lambda_min <- lm(as.formula(model.protein.ELASTIC.NET.string.min), data = df.protein.train.lasso)

model.protein.string.1se
model.protein.string.min

###########
###LASSO###
###########
model.protein.LASSO <- cv.glmnet(X.train.protein, Y.train.protein, nfolds = 10, standardize = F, intercept = T, alpha = 1)

optimal.protein.betas.lambda.min <- coef(model.protein.LASSO, s = "lambda.min")
optimal.protein.betas.lambda.1se <- coef(model.protein.LASSO, s = "lambda.1se")
optimal.protein.betas.lambda.min 
optimal.protein.betas.lambda.1se

betas.protein.names <- names(df.protein.train.lasso)
beta.protein.table <- as.data.frame(cbind("LAMBDA_MIN" = optimal.protein.betas.lambda.min[,1], "LAMBDA_1SE" = optimal.protein.betas.lambda.1se[,1]))
beta.protein.table <- data.frame("NAMES" = row.names(beta.protein.table), beta.protein.table )
beta.protein.table <- beta.protein.table[which(beta.protein.table$LAMBDA_MIN != 0 | beta.protein.table$LAMBDA_1SE != 0),]
betas.protein.1se <- as.data.frame(as.character(beta.protein.table$NAMES[which(beta.protein.table$LAMBDA_1SE != 0)]))
betas.protein.min <- as.data.frame(as.character(beta.protein.table$NAMES[which(beta.protein.table$LAMBDA_MIN != 0)]))

model.protein.string.1se = betas.protein.1se[2,1]
for(i in 3:nrow(betas.protein.1se)){
  model.protein.string.1se = paste(model.protein.string.1se, betas.protein.1se[i,1], sep = " + ")
}

model.protein.string.min = betas.protein.min[2,1]
for(i in 3:nrow(betas.protein.min)){
  model.protein.string.min = paste(model.protein.string.min, betas.protein.min[i,1], sep = " + ")
}

model.protein.LASSO.string.1se <- paste("Y ~ ",model.protein.string.1se, sep = "")
model.protein.LASSO.string.min <- paste("Y ~ ",model.protein.string.min, sep = "")
model.protein.string.1se
model.protein.string.min

model.protein.LASSO.lambda_1SE <- lm(as.formula(model.protein.LASSO.string.1se), data = df.protein.train.lasso)
model.protein.LASSO.lambda_min <- lm(as.formula(model.protein.LASSO.string.min), data = df.protein.train.lasso)


################################################################
#CHOOSING THE MODEL USING RANDOM FOREST#########################
################################################################
################################################################
#install.packages("randomForest")
library(randomForest)

'small.subset <- sample(x = 1:nrow(dt.train), size = floor(nrow(dt.train) * 0.8))

dt.train.small <- dt.train[small.subset,]
dt.valid.small <- dt.train[-small.subset,]
dt.valid.small.no.Y <- dt.valid.small[, names(dt.valid.small) != "Y"]'
valid.protein.no.Y = valid.protein[, names(valid.protein) != "Y"]

predictors.square.root = floor((ncol(train.protein))^0.5)
predictors.square.root
predictors.protein <- seq(210,230,by = 3)
predictors.protein
trenuj_i_daj_MSE = function(m) {
  las_losowy <- randomForest(Y~.,train.protein,mtry=m)
  predykcja_rf <- predict(las_losowy, newdata = valid.protein.no.Y)
  val <- sum((predykcja_rf - valid.protein$Y)^2)
  return(val/nrow(valid.protein))
}

mses.protein <- sapply(predictors.protein, trenuj_i_daj_MSE)
plot(predictors.protein, mses.protein)

optimal.protein.predictor <- predictors.protein[which(mses.protein == min(mses.protein))]
optimal.protein.predictor

model.protein.RandomForest <- randomForest(Y~.,train.protein,mtry=optimal.protein.predictor)

################################################################
#CHOOSING THE MODEL USING BOOSTINGU#############################
################################################################
################################################################
library(gbm)
library(ISLR)
valid.protein.no.Y = valid.protein[, names(valid.protein) != "Y"]

lambdas <- seq(0.05, 0.15, by = 0.001)
lambdas
trenuj_gbm_i_daj_mse <- function(lambda) {
  model <- gbm(Y ~ ., data = train.protein, n.trees = 300, shrinkage = lambda)
  prediction <- predict(model, newdata = valid.protein.no.Y, n.trees = 300)
  mse <- (sum((valid.protein$Y - prediction)^2))/nrow(valid.protein)
  return(mse)
}

mses.protein.gbm <- sapply(lambdas, trenuj_gbm_i_daj_mse)
plot(lambdas, mses.protein.gbm)

lambda_gbm <- lambdas[which.min(mses.protein.gbm)]
lambda_gbm
model.protein.GBMTree <- gbm(Y ~ ., data = train.protein, n.trees = 300, shrinkage = lambda_gbm)

###???????????#####################???????????###
####_______________VALIDATION________________####
###???????????#####################???????????###

prediction.protein.F <- predict(model.protein.F,valid.protein)

prediction.protein.BIC <- predict(model.protein.BIC, valid.protein)

prediction.protein.ELASTIC.NET.1se.FALSE <- predict(model.protein.ELASTIC.NET.lambda_1SE, valid.protein)
prediction.protein.ELASTIC.NET.min.FALSE <- predict(model.protein.ELASTIC.NET.lambda_min, valid.protein)
prediction.protein.LASSO.1se.FALSE <- predict(model.protein.LASSO.lambda_1SE, valid.protein)
prediction.protein.LASSO.min.FALSE <- predict(model.protein.LASSO.lambda_min, valid.protein)

prediction.protein.ELASTIC.NET.1se.TRUE <- predict.cv.glmnet(object = model.protein.ELASTIC.NET, newx = as.matrix(valid.protein[,-dim(valid.protein)[2]]), s = "lambda.1se")
prediction.protein.ELASTIC.NET.min.TRUE <- predict.cv.glmnet(object = model.protein.ELASTIC.NET, newx = as.matrix(valid.protein[,-dim(valid.protein)[2]]), s = "lambda.min")
prediction.protein.LASSO.1se.TRUE <- predict.cv.glmnet(object = model.protein.LASSO, newx = as.matrix(valid.protein[,-dim(valid.protein)[2]]), s = "lambda.1se")
prediction.protein.LASSO.min.TRUE <- predict.cv.glmnet(object = model.protein.LASSO, newx = as.matrix(valid.protein[,-dim(valid.protein)[2]]), s = "lambda.min")

prediction.protein.RandomForest <- predict(model.protein.RandomForest, valid.protein)

prediction.protein.GBMTree <- predict(model.protein.GBMTree, valid.protein,n.trees = 300)


MSE_test1 <- (sum((valid.protein$Y - prediction.protein.F)^2))/nrow(valid.protein) #MSE for the model with variables according to F-test
MSE_test2 <- (sum((valid.protein$Y - prediction.protein.BIC)^2))/nrow(valid.protein) #MSE for the model with variables according to BIC

MSE_test3 <- (sum((valid.protein$Y - prediction.protein.LASSO.1se.FALSE)^2))/nrow(valid.protein) #MSE for the model with variables according to LASSO with lambda 1se
MSE_test4 <- (sum((valid.protein$Y - prediction.protein.LASSO.min.FALSE)^2))/nrow(valid.protein) #MSE for the model with variables according to LASSO with lambda min
MSE_test5 <- (sum((valid.protein$Y - prediction.protein.ELASTIC.NET.1se.FALSE)^2))/nrow(valid.protein) #MSE for the model with variables according to LASSO with lambda 1se
MSE_test6 <- (sum((valid.protein$Y - prediction.protein.ELASTIC.NET.min.FALSE)^2))/nrow(valid.protein) #MSE for the model with variables according to LASSO with lambda min

MSE_test7 <- (sum((valid.protein$Y - prediction.protein.LASSO.1se.TRUE)^2))/nrow(valid.protein) #MSE for the model with variables according to LASSO with lambda 1se
MSE_test8 <- (sum((valid.protein$Y - prediction.protein.LASSO.min.TRUE)^2))/nrow(valid.protein) #MSE for the model with variables according to LASSO with lambda min
MSE_test9 <- (sum((valid.protein$Y - prediction.protein.ELASTIC.NET.1se.TRUE)^2))/nrow(valid.protein) #MSE for the model with variables according to LASSO with lambda 1se
MSE_test10 <- (sum((valid.protein$Y - prediction.protein.ELASTIC.NET.min.TRUE)^2))/nrow(valid.protein) #MSE for the model with variables according to LASSO with lambda min


MSE_test11 <- (sum((valid.protein$Y - prediction.protein.RandomForest)^2))/nrow(valid.protein) #MSE for the model with variables according to Random Forest
MSE_test12 <- (sum((valid.protein$Y - prediction.protein.GBMTree)^2))/nrow(valid.protein) #MSE for the model with variables according to modelu Boosting Tree


MSE_test1
MSE_test2
MSE_test3
MSE_test4
MSE_test5
MSE_test6
MSE_test7
MSE_test8
MSE_test9
MSE_test10
MSE_test11
MSE_test12

###???????????#####################???????????###
####____________CROSS-VALIDATION_____________####
###???????????#####################???????????###

###TRAINED FOR ALL OF THE VARIABLES###

chunks <- seq(1,nrow(df.train.protein), length.out = 6)
for(s in 1:length(chunks)){
  chunks[s] <- as.integer(chunks[s])
}

cross.validates.mses.protein <- as.data.frame(matrix(,nrow = 5, ncol = 12))
names(cross.validates.mses.protein) <- c("MODEL.F","BIC",
                                         "ELASTIC.NET.1se.FALSE", "ELASTIC.NET.min.FALSE", "LASSO.1se.FALSE", "LASSO.min.FALSE", 
                                         "ELASTIC.NET.1se.TRUE", "ELASTIC.NET.min.TRUE", "LASSO.1se.TRUE", "LASSO.min.TRUE",
                                         "RANDOM.FOREST", "BOOSTING")
chunks[2]
for(i in 2:(length(chunks))-1){
  
  rows.to.use <- seq(chunks[i],chunks[i+1], by = 1)
  training.set <- df.train.protein[-c(rows.to.use),]
  validation.set <- df.train.protein[c(rows.to.use),]
  
  print(i)
  
  cv.model.protein.ELASTIC.NET.lambda_1SE.FALSE  <- lm(as.formula(model.protein.ELASTIC.NET.string.1se), data = training.set)
  cv.model.protein.ELASTIC.NET.lambda_min.FALSE  <- lm(as.formula(model.protein.ELASTIC.NET.string.min), data = training.set)
  cv.model.protein.LASSO.lambda_1SE.FALSE        <- lm(as.formula(model.protein.LASSO.string.1se), data = training.set)
  cv.model.protein.LASSO.lambda_min.FALSE        <- lm(as.formula(model.protein.LASSO.string.min), data = training.set)
  
  cv.model.protein.ELASTIC.NET.TRUE  <- cv.glmnet(as.matrix(training.set[,-dim(training.set)[2]]), as.matrix(training.set[,dim(training.set)[2]]), nfolds = 10, standardize = F, intercept = T, thresh = 1e-16)
  cv.model.protein.LASSO.TRUE        <- cv.glmnet(as.matrix(training.set[,-dim(training.set)[2]]), as.matrix(training.set[,dim(training.set)[2]]), nfolds = 10, standardize = F, intercept = T, alpha = 1)
  
  cv.model.protein.RandomForest            <- randomForest(Y~.,training.set,mtry=optimal.protein.predictor)
  cv.model.protein.GBMTree                 <- gbm(Y ~ ., data = training.set, n.trees = 300, shrinkage = lambda_gbm)
  
  ##########################################
  cv.prediction.protein.MODEL.F         <- predict(model.protein.F, validation.set[-c(dim(validation.set)[2])])
  cv.prediction.protein.BIC             <- predict(model.protein.BIC, validation.set[-c(dim(validation.set)[2])])
  
  cv.prediction.protein.ELASTIC.NET.1se.FALSE <- predict(cv.model.protein.ELASTIC.NET.lambda_1SE.FALSE, validation.set[-c(dim(validation.set)[2])])
  cv.prediction.protein.ELASTIC.NET.min.FALSE <- predict(cv.model.protein.ELASTIC.NET.lambda_min.FALSE, validation.set[-c(dim(validation.set)[2])])
  cv.prediction.protein.LASSO.1se.FALSE       <- predict(cv.model.protein.LASSO.lambda_1SE.FALSE, validation.set[-c(dim(validation.set)[2])])
  cv.prediction.protein.LASSO.min.FALSE       <- predict(cv.model.protein.LASSO.lambda_min.FALSE, validation.set[-c(dim(validation.set)[2])])
  
  cv.prediction.protein.ELASTIC.NET.1se.TRUE <- predict.cv.glmnet(object = cv.model.protein.ELASTIC.NET.TRUE, newx = as.matrix(validation.set[,-dim(validation.set)[2]]), s = "lambda.1se")
  cv.prediction.protein.ELASTIC.NET.min.TRUE <- predict.cv.glmnet(object = cv.model.protein.ELASTIC.NET.TRUE, newx = as.matrix(validation.set[,-dim(validation.set)[2]]), s = "lambda.min")
  cv.prediction.protein.LASSO.1se.TRUE       <- predict.cv.glmnet(object = cv.model.protein.LASSO.TRUE, newx = as.matrix(validation.set[,-dim(validation.set)[2]]), s = "lambda.1se")
  cv.prediction.protein.LASSO.min.TRUE       <- predict.cv.glmnet(object = cv.model.protein.LASSO.TRUE, newx = as.matrix(validation.set[,-dim(validation.set)[2]]), s = "lambda.min")
  
  cv.prediction.protein.RandomForest    <- predict(cv.model.protein.RandomForest, validation.set)
  cv.prediction.protein.GBMTree         <- predict(cv.model.protein.GBMTree, validation.set,n.trees = 300)
  ##########################################
  
  cross.validates.mses.protein$MODEL.F[i]          <- (sum((validation.set$Y - cv.prediction.protein.MODEL.F)^2))/nrow(validation.set) #MSE for the model with variables according to LASSO with lambda 1se
  cross.validates.mses.protein$BIC[i]              <- (sum((validation.set$Y - cv.prediction.protein.BIC)^2))/nrow(validation.set) #MSE for the model with variables according to LASSO with lambda min
  
  cross.validates.mses.protein$ELASTIC.NET.1se.FALSE[i]  <- (sum((validation.set$Y - cv.prediction.protein.ELASTIC.NET.1se.FALSE)^2))/nrow(validation.set) #MSE for the model with variables according to LASSO with lambda 1se
  cross.validates.mses.protein$ELASTIC.NET.min.FALSE[i]  <- (sum((validation.set$Y - cv.prediction.protein.ELASTIC.NET.min.FALSE)^2))/nrow(validation.set) #MSE for the model with variables according to LASSO with lambda min
  cross.validates.mses.protein$LASSO.1se.FALSE[i]        <- (sum((validation.set$Y - cv.prediction.protein.LASSO.1se.FALSE)^2))/nrow(validation.set) #MSE for the model with variables according to LASSO with lambda 1se
  cross.validates.mses.protein$LASSO.min.FALSE[i]        <- (sum((validation.set$Y - cv.prediction.protein.LASSO.min.FALSE)^2))/nrow(validation.set) #MSE for the model with variables according to LASSO with lambda min
  
  cross.validates.mses.protein$ELASTIC.NET.1se.TRUE[i]  <- (sum((validation.set$Y - cv.prediction.protein.ELASTIC.NET.1se.TRUE)^2))/nrow(validation.set) #MSE for the model with variables according to LASSO with lambda 1se
  cross.validates.mses.protein$ELASTIC.NET.min.TRUE[i]  <- (sum((validation.set$Y - cv.prediction.protein.ELASTIC.NET.min.TRUE)^2))/nrow(validation.set) #MSE for the model with variables according to LASSO with lambda min
  cross.validates.mses.protein$LASSO.1se.TRUE[i]        <- (sum((validation.set$Y - cv.prediction.protein.LASSO.1se.TRUE)^2))/nrow(validation.set) #MSE for the model with variables according to LASSO with lambda 1se
  cross.validates.mses.protein$LASSO.min.TRUE[i]        <- (sum((validation.set$Y - cv.prediction.protein.LASSO.min.TRUE)^2))/nrow(validation.set) #MSE for the model with variables according to LASSO with lambda min
  
  cross.validates.mses.protein$RANDOM.FOREST[i]    <- (sum((validation.set$Y - cv.prediction.protein.RandomForest)^2))/nrow(validation.set) #MSE for the model with variables according to Random Forest
  cross.validates.mses.protein$BOOSTING[i]         <- (sum((validation.set$Y - cv.prediction.protein.GBMTree)^2))/nrow(validation.set) #MSE for the model with variables according to Boosting Tree
}

cross.validates.mses.summary.protein <- cross.validates.mses.protein[1,]
cross.validates.mses.summary.protein <- apply(cross.validates.mses.protein,2,mean)
cross.validates.mses.summary.protein <- as.data.frame(cross.validates.mses.summary.protein)

###PREDITCED Y###
#1. selection of variables essential for prediction using ELASTIC NET
protein.final.train <- as.data.frame(data.train)
protein.final.test <- df.test.protein

dims.protein.train.final <- dim(protein.final.train)

X.train.protein.final <- protein.final.train[,-dims.protein.train.final[2]]
Y.train.protein.final <- protein.final.train[,dims.protein.train.final[2]]

X.train.protein.final <- as.matrix(X.train.protein.final)
Y.train.protein.final <- as.matrix(Y.train.protein.final)
dim(Y.train.protein.final)

model.protein.ELASTIC.NET.final <- cv.glmnet(X.train.protein.final, Y.train.protein.final, nfolds = 10, standardize = F, intercept = T, thresh = 1e-16)

optimal.betas.protein.lambda.1se.final <- coef(model.protein.ELASTIC.NET.final, s = "lambda.1se")
optimal.betas.protein.lambda.1se.final

betas.protein.names.final <- names(protein.final.train)
beta.protein.table.final <- as.data.frame(cbind("LAMBDA_1SE" = optimal.betas.protein.lambda.1se.final[,1]))
beta.protein.table.final <- data.frame("NAMES" = row.names(beta.protein.table.final), beta.protein.table.final )
beta.protein.table.final <- beta.protein.table.final[which(beta.protein.table.final$LAMBDA_1SE != 0),]
betas.protein.1se.final <- as.data.frame(as.character(beta.protein.table.final$NAMES[which(beta.protein.table.final$LAMBDA_1SE != 0)]))

model.protein.string.1se.final = betas.protein.1se.final[2,1]
for(i in 3:nrow(betas.protein.1se.final)){
  model.protein.string.1se.final = paste(model.protein.string.1se.final, betas.protein.1se.final[i,1], sep = " + ")
}

#2. building linear regression model using previously selected variables
model.protein.ELASTIC.NET.string.1se.final <- paste("Y ~ ",model.protein.string.1se.final, sep = "")
model.protein.ELASTIC.NET.string.1se.final
model.protein.ELASTIC.NET.lambda_1SE.final <- lm(as.formula(model.protein.ELASTIC.NET.string.1se.final), data = protein.final.train)

model.protein.string.1se.final

ELASTIC.coeficient <- coef(model.protein.ELASTIC.NET.lambda_1SE.final)

ELASTIC.coeficient[abs(ELASTIC.coeficient)>1]
#3. predicting Y
pred.protein <- predict(model.protein.ELASTIC.NET.lambda_1SE.final, protein.final.test[-c(dim(protein.final.test)[2])])
pred.protein
predictors.protein.dummy <- ELASTIC.coeficient[abs(ELASTIC.coeficient)>1]
predictors.protein <- names(predictors.protein.dummy)[2:6]
#===============================================================================================================================#
#===============================================================================================================================#
#===============================================================================================================================#
#===============================================================================================================================#
#===============================================================================================================================#
#===============================================================================================================================#

###!!!!!!!!!!!#####################!!!!!!!!!!!###
###!!!!!!!!!!!####cancer.RData#####!!!!!!!!!!!###
###!!!!!!!!!!!#####################!!!!!!!!!!!###

#loading cancer
load("C:/Users/Igor/Desktop/Bioinformatyka I rok II stopnia/II sem/SAD1/projekt zaliczeniowy/Grupa5/cancer.RData")

df.train.cancer <- as.data.frame(data.train)
df.test.cancer <- as.data.frame(data.test)

#=======================#PRELIMINARY DATA ANALYSIS#=======================#

plot(df.train.cancer$Y)
max(df.train.cancer$Y)
min(df.train.cancer$Y)

# using the histogram I am trying to detect if there are variables whose unique 
# values are rare (below 10), what would indicate the existence of qualitative variables
cancer.data.unique.values <- apply(df.train.cancer, 2, function(x){length(unique(x))})
hist(cancer.data.unique.values)

# the histogram shows that there are no qualitative variables in the cancer set
# numeric variables, I scale continuously
# except the explanatory variable 'Y'
headers<-c("minimum","maximum")
histogram_data <- as.data.frame(matrix(,ncol=2,nrow=ncol(df.train.cancer)))
names(histogram_data)<-headers


for(i in 1:dim(df.train.cancer)[2]){
  histogram_data$minimum[i] <-min(df.train.cancer[,i])
  histogram_data$maximum[i] <-max(df.train.cancer[,i])
}

hist(histogram_data$minimum)
hist(histogram_data$maximum)

df.train.cancer.normalized <- df.train.cancer
df.test.cancer.normalized <- df.test.cancer
for(i in 1:dim(df.test.cancer)[2]){
  df.train.cancer.normalized[,i] <- scale(df.train.cancer[,i])
  df.test.cancer.normalized[,i] <- scale(df.test.cancer[,i])
}

for(i in 1:(dim(df.train.cancer.normalized)[2]-1)){
  histogram_data$minimum[i] <-min(df.train.cancer.normalized[,i])
  histogram_data$maximum[i] <-max(df.train.cancer.normalized[,i])
}

hist(histogram_data$minimum)
hist(histogram_data$maximum)

# Data are too large to make meaningful data on them, therefore
# I will analyze the correlation of explanatory variables with the explained variable and
# I keep only these explanatory variables for further analysis,
# which have a high correlation with the explained variable

cancer.correlation <- df.train.cancer[1,1:ncol(df.train.cancer.normalized)]

for(i in 1:ncol(df.train.cancer.normalized)){
  cancer.correlation[,i] <- cor(x = df.train.cancer.normalized[,i],y = df.train.cancer.normalized$Y)
}

a <- seq(0.01,0.45, by = 0.02)
plot_of_significance <- as.data.frame(matrix(,nrow = 2,ncol = length(a)))

for(b in 1:length(a)){
  correlation_good <- cancer.correlation[,which(abs(cancer.correlation) >= a[b])]
  plot_of_significance[1,b] <- a[b]
  plot_of_significance[2,b] <- dim(correlation_good)[2]
}

regulator <- 0.25

#df.train.cancer.correlated <- df.train.cancer.normalized[,which(abs(cancer.correlation) >= 0.14)]
df.train.cancer.correlated <- df.train.cancer.normalized[,which(abs(cancer.correlation) >= regulator)]


# From a fairly quick analysis of the covariance between predictors and the explained variable,
# it cannot be stated unequivocally that there is a small, specific group of variables on which it depends
# variability of the Y variable. Based on this I will use methods that predict Y well from many variables i.e.
# elastic net and tree methods

###VIF###
cancer_model <- lm(Y~.,data = df.train.cancer.correlated)
cancer.vif <- vif(cancer_model)
cancer.vif.5 <- cancer.vif[cancer.vif <= 5]
cancer.vif.10 <- cancer.vif[cancer.vif <= 10]
hist(cancer.vif)
hist(cancer.vif.5)
hist(cancer.vif.10)
#df.train.cancer.correlated <- df.train.cancer.correlated[cancer.vif <= 8]

#=======================#THE END OF PRELIMINARY ANALYSIS#=======================#

cancer.data.sample.rows <- sample(x = 1:nrow(df.train.cancer.correlated), size = floor(nrow(df.train.cancer.correlated) * 0.8))
train.cancer <- df.train.cancer.correlated[cancer.data.sample.rows,]
valid.cancer <- df.train.cancer.correlated[-cancer.data.sample.rows,]
test.cancer <- as.data.frame(df.test.cancer.normalized)

################################################################
#CHOOSING THE MODEL USING ELASTIC NET AND LASSO############
################################################################
################################################################
#install.packages("glmnet")
library(glmnet)

dims.train.cancer <- dim(train.cancer)
dims.valid.cancer <- dim(test.cancer)

X.train.cancer <- train.cancer[,-dims.train.cancer[2]]
Y.train.cancer <- train.cancer[,dims.train.cancer[2]]
X.valid.cancer <- test.cancer[,-dims.valid.cancer[2]]
Y.valid.cancer <- test.cancer[,dims.valid.cancer[2]]

X.train.cancer <- as.matrix(X.train.cancer)
Y.train.cancer <- as.matrix(Y.train.cancer)
X.valid.cancer <- as.matrix(X.valid.cancer)
Y.valid.cancer <- as.matrix(Y.valid.cancer)

#################
###ELASTIC NET###
#################
model.cancer.ELASTIC.NET <- cv.glmnet(X.train.cancer, Y.train.cancer, nfolds = 10, standardize = F, intercept = T, thresh = 1e-16)

optimal.cancer.betas.lambda.min <- coef(model.cancer.ELASTIC.NET, s = "lambda.min")
optimal.cancer.betas.lambda.1se <- coef(model.cancer.ELASTIC.NET, s = "lambda.1se")
optimal.cancer.betas.lambda.min 
optimal.cancer.betas.lambda.1se

betas.cancer.names <- names(df.train.cancer.correlated)
beta.cancer.table <- as.data.frame(cbind("LAMBDA_MIN" = optimal.cancer.betas.lambda.min[,1], "LAMBDA_1SE" = optimal.cancer.betas.lambda.1se[,1]))
beta.cancer.table <- data.frame("NAMES" = row.names(beta.cancer.table), beta.cancer.table )
beta.cancer.table <- beta.cancer.table[which(beta.cancer.table$LAMBDA_MIN != 0 | beta.cancer.table$LAMBDA_1SE != 0),]
betas.cancer.1se <- as.data.frame(as.character(beta.cancer.table$NAMES[which(beta.cancer.table$LAMBDA_1SE != 0)]))
betas.cancer.min <- as.data.frame(as.character(beta.cancer.table$NAMES[which(beta.cancer.table$LAMBDA_MIN != 0)]))

model.cancer.string.1se = betas.cancer.1se[2,1]
for(i in 3:nrow(betas.cancer.1se)){
  model.cancer.string.1se = paste(model.cancer.string.1se, betas.cancer.1se[i,1], sep = " + ")
}

model.cancer.string.min = betas.cancer.min[2,1]
for(i in 3:nrow(betas.cancer.min)){
  model.cancer.string.min = paste(model.cancer.string.min, betas.cancer.min[i,1], sep = " + ")
}

model.cancer.ELASTIC.NET.string.1se <- paste("Y ~ ", model.cancer.string.1se, sep = "")
model.cancer.ELASTIC.NET.string.min <- paste("Y ~ ", model.cancer.string.1se, sep = "")
model.cancer.ELASTIC.NET.string.1se
model.cancer.ELASTIC.NET.lambda_1SE <- lm(as.formula(model.cancer.ELASTIC.NET.string.1se), data = train.cancer)
model.cancer.ELASTIC.NET.lambda_min <- lm(as.formula(model.cancer.ELASTIC.NET.string.min), data = train.cancer)

###########
###LASSO###
###########
model.cancer.LASSO <- cv.glmnet(X.train.cancer, Y.train.cancer, nfolds = 10, standardize = F, intercept = T, alpha = 0)

optimal.cancer.betas.lambda.min <- coef(model.cancer.LASSO, s = "lambda.min")
optimal.cancer.betas.lambda.1se <- coef(model.cancer.LASSO, s = "lambda.1se")
optimal.cancer.betas.lambda.min 
optimal.cancer.betas.lambda.1se

betas.cancer.names <- names(df.train.cancer.correlated)
beta.cancer.table <- as.data.frame(cbind("LAMBDA_MIN" = optimal.cancer.betas.lambda.min[,1], "LAMBDA_1SE" = optimal.cancer.betas.lambda.1se[,1]))
beta.cancer.table <- data.frame("NAMES" = row.names(beta.cancer.table), beta.cancer.table )
beta.cancer.table <- beta.cancer.table[which(beta.cancer.table$LAMBDA_MIN != 0 | beta.cancer.table$LAMBDA_1SE != 0),]
betas.cancer.1se <- as.data.frame(as.character(beta.cancer.table$NAMES[which(beta.cancer.table$LAMBDA_1SE != 0)]))
betas.cancer.min <- as.data.frame(as.character(beta.cancer.table$NAMES[which(beta.cancer.table$LAMBDA_MIN != 0)]))

model.cancer.string.1se = betas.cancer.1se[2,1]
for(i in 3:nrow(betas.cancer.1se)){
  model.cancer.string.1se = paste(model.cancer.string.1se, betas.cancer.1se[i,1], sep = " + ")
}

model.cancer.string.min = betas.cancer.min[2,1]
for(i in 3:nrow(betas.cancer.min)){
  model.cancer.string.min = paste(model.cancer.string.min, betas.cancer.min[i,1], sep = " + ")
}

model.cancer.LASSO.string.1se <- paste("Y ~ ", model.cancer.string.1se, sep = "")
model.cancer.LASSO.string.min <- paste("Y ~ ", model.cancer.string.1se, sep = "")
model.cancer.LASSO.string.1se
model.cancer.LASSO.lambda_1SE <- lm(as.formula(model.cancer.LASSO.string.1se), data = train.cancer)
model.cancer.LASSO.lambda_min <- lm(as.formula(model.cancer.LASSO.string.min), data = train.cancer)

################################################################
#CHOOSING THE MODEL USING RANDOM FORESTS########################
################################################################
################################################################
#install.packages("randomForest")
library(randomForest)

'small.subset <- sample(x = 1:nrow(dt.train), size = floor(nrow(dt.train) * 0.8))

dt.train.small <- dt.train[small.subset,]
dt.valid.small <- dt.train[-small.subset,]
dt.valid.small.no.Y <- dt.valid.small[, names(dt.valid.small) != "Y"]'
valid.cancer.no.Y = valid.cancer[, names(valid.cancer) != "Y"]

predictors.cancer.square.root = floor((ncol(train.cancer))^0.5)
predictors.cancer.square.root

predictors.cancer <- seq(floor(4*predictors.cancer.square.root),220, by = 10)
predictors.cancer
trenuj_i_daj_MSE = function(m) {
  las_losowy <- randomForest(Y~.,train.cancer,mtry=m)
  predykcja_rf <- predict(las_losowy, newdata = valid.cancer.no.Y)
  val <- sum((predykcja_rf - valid.cancer$Y)^2)
  return(val/nrow(valid.cancer))
}

mse.cancer <- sapply(predictors.cancer, trenuj_i_daj_MSE)
plot(predictors.cancer, mse.cancer)
predictors.cancer

optimal.cancer.predictor <- predictors.cancer[which(mse.cancer == min(mse.cancer))]
optimal.cancer.predictor
predictors.cancer.square.root = floor((ncol(train.cancer))^0.5)
model.cancer.RandomForest <- randomForest(Y~.,train.cancer,mtry=optimal.cancer.predictor)
################################################################
#CHOOSING THE MODEL USING BOOSTINGU#############################
################################################################
################################################################
library(gbm)
library(ISLR)
valid.cancer.no.Y = valid.cancer[, names(valid.cancer) != "Y"]

lambdas.cancer <- seq(0.001, 0.04, by = 0.001)
lambdas.cancer
trenuj_gbm_i_daj_mse <- function(lambda) {
  model <- gbm(Y ~ ., data = train.cancer, n.trees = 400, shrinkage = lambda)
  prediction <- predict(model, newdata = valid.cancer.no.Y, n.trees = 400)
  mse <- (sum((valid.cancer$Y - prediction)^2))/nrow(valid.cancer)
  return(mse)
}

mses.cancer <- sapply(lambdas.cancer, trenuj_gbm_i_daj_mse)
plot(lambdas.cancer, mses.cancer)

lambda_gbm <- lambdas.cancer[which.min(mses.cancer)]
lambda_gbm
model.cancer.GBMTree <- gbm(Y ~ ., data = train.cancer, n.trees = 400, shrinkage = lambda_gbm)

###???????????#####################???????????###
####________________VALIDATION_______________####
###???????????#####################???????????###

#prediction.F <- predict(model.F,dt.valid)
#prediction.BIC <- predict(model.BIC, dt.valid)
prediction.cancer.ELASTIC.NET.1se <- predict(model.cancer.ELASTIC.NET.lambda_1SE, valid.cancer)
prediction.cancer.ELASTIC.NET.min <- predict(model.cancer.ELASTIC.NET.lambda_min, valid.cancer)
prediction.cancer.LASSO.1se <- predict(model.cancer.LASSO.lambda_1SE, valid.cancer)
prediction.cancer.LASSO.min <- predict(model.cancer.LASSO.lambda_min, valid.cancer)
prediction.cancer.RandomForest <- predict(model.cancer.RandomForest, valid.cancer)
prediction.cancer.GBMTree <- predict(model.cancer.GBMTree, valid.cancer,n.trees = 400)

#MSE_test1 <- (sum((dt.valid$Y - prediction.F)^2))/nrow(dt.valid) #MSE for the model with variables according to F-test
#MSE_test2 <- (sum((dt.valid$Y - prediction.BIC)^2))/nrow(dt.valid) #MSE for the model with variables according to BIC
MSE_test3.cancer <- (sum((valid.cancer$Y - prediction.cancer.LASSO.1se)^2))/nrow(valid.cancer) #MSE for the model with variables according to LASSO with lambda 1se
MSE_test4.cancer <- (sum((valid.cancer$Y - prediction.cancer.LASSO.min)^2))/nrow(valid.cancer) #MSE for the model with variables according to LASSO with lambda min
MSE_test5.cancer <- (sum((valid.cancer$Y - prediction.cancer.ELASTIC.NET.1se)^2))/nrow(valid.cancer) #MSE for the model with variables according to LASSO with lambda 1se
MSE_test6.cancer <- (sum((valid.cancer$Y - prediction.cancer.ELASTIC.NET.min)^2))/nrow(valid.cancer) #MSE for the model with variables according to LASSO with lambda min

MSE_test7.cancer <- (sum((valid.cancer$Y - prediction.cancer.RandomForest)^2))/nrow(valid.cancer) #MSE for the model with variables according to Random Forest
MSE_test8.cancer <- (sum((valid.cancer$Y - prediction.cancer.GBMTree)^2))/nrow(valid.cancer) #MSE for the model with variables according to Boosting Tree


#MSE_test1.cancer
#MSE_test2.cancer
MSE_test3.cancer
MSE_test4.cancer
MSE_test5.cancer
MSE_test6.cancer
MSE_test7.cancer
MSE_test8.cancer

###???????????#####################???????????###
####_____________CROSS-VALIDATION____________####
###???????????#####################???????????###

chunks <- seq(1,nrow(df.train.cancer.normalized), length.out = 6)
for(s in 1:length(chunks)){
  chunks[s] <- as.integer(chunks[s])
}

cross.validates.mses.cancer <- as.data.frame(matrix(,nrow = 5, ncol = 11))
names(cross.validates.mses.cancer) <- c("LINEAR.REGRESSION", "ELASTIC.NET.1se.FALSE", "ELASTIC.NET.min.FALSE", "RIDGE.1se.FALSE", "RIDGE.min.FALSE",
                                 "ELASTIC.NET.1se.TRUE", "ELASTIC.NET.min.TRUE", "RIDGE.1se.TRUE", "RIDGE.min.TRUE",
                                 "RANDOM.FOREST", "BOOSTING")
chunks[2]
for(i in 2:(length(chunks))-1){
  
  rows.to.use <- seq(chunks[i],chunks[i+1], by = 1)
  training.set <- df.train.cancer.correlated[-c(rows.to.use),]
  validation.set <- df.train.cancer.correlated[c(rows.to.use),]
  
  print(i)
  cv.model.cancer.linear.regression             <- lm(Y~., data = training.set)
  
  cv.model.cancer.ELASTIC.NET.lambda_1SE.FALSE  <- lm(as.formula(model.cancer.ELASTIC.NET.string.1se), data = training.set)
  cv.model.cancer.ELASTIC.NET.lambda_min.FALSE  <- lm(as.formula(model.cancer.ELASTIC.NET.string.min), data = training.set)
  cv.model.cancer.LASSO.lambda_1SE.FALSE        <- lm(as.formula(model.cancer.LASSO.string.1se), data = training.set)
  cv.model.cancer.LASSO.lambda_min.FALSE        <- lm(as.formula(model.cancer.LASSO.string.min), data = training.set)
  
  cv.model.cancer.ELASTIC.NET.TRUE  <- cv.glmnet(as.matrix(training.set[,-dim(training.set)[2]]), as.matrix(training.set[,dim(training.set)[2]]), nfolds = 5, standardize = F, intercept = T, thresh = 1e-16)
  cv.model.cancer.LASSO.TRUE        <- cv.glmnet(as.matrix(training.set[,-dim(training.set)[2]]), as.matrix(training.set[,dim(training.set)[2]]), nfolds = 5, standardize = F, intercept = T, alpha = 0)
  
  cv.model.cancer.RandomForest            <- randomForest(Y~.,training.set,mtry=optimal.cancer.predictor)
  cv.model.cancer.GBMTree                 <- gbm(Y ~ ., data = training.set, n.trees = 400, shrinkage = lambda_gbm)
  
  ######################################
  cv.prediction.cancer.linear.regression     <- predict(cv.model.cancer.linear.regression, validation.set[-c(dim(validation.set)[2])])
  
  cv.prediction.cancer.ELASTIC.NET.1se.FALSE <- predict(cv.model.cancer.ELASTIC.NET.lambda_1SE.FALSE, validation.set[-c(dim(validation.set)[2])])
  cv.prediction.cancer.ELASTIC.NET.min.FALSE <- predict(cv.model.cancer.ELASTIC.NET.lambda_min.FALSE, validation.set[-c(dim(validation.set)[2])])
  cv.prediction.cancer.LASSO.1se.FALSE       <- predict(cv.model.cancer.LASSO.lambda_1SE.FALSE, validation.set[-c(dim(validation.set)[2])])
  cv.prediction.cancer.LASSO.min.FALSE       <- predict(cv.model.cancer.LASSO.lambda_min.FALSE, validation.set[-c(dim(validation.set)[2])])
  
  cv.prediction.cancer.ELASTIC.NET.1se.TRUE <- predict.cv.glmnet(object = cv.model.cancer.ELASTIC.NET.TRUE, newx = as.matrix(validation.set[,-dim(validation.set)[2]]), s = "lambda.1se")
  cv.prediction.cancer.ELASTIC.NET.min.TRUE <- predict.cv.glmnet(object = cv.model.cancer.ELASTIC.NET.TRUE, newx = as.matrix(validation.set[,-dim(validation.set)[2]]), s = "lambda.min")
  cv.prediction.cancer.LASSO.1se.TRUE       <- predict.cv.glmnet(object = cv.model.cancer.LASSO.TRUE, newx = as.matrix(validation.set[,-dim(validation.set)[2]]), s = "lambda.1se")
  cv.prediction.cancer.LASSO.min.TRUE       <- predict.cv.glmnet(object = cv.model.cancer.LASSO.TRUE, newx = as.matrix(validation.set[,-dim(validation.set)[2]]), s = "lambda.min")
  
  
  cv.prediction.cancer.RandomForest    <- predict(cv.model.cancer.RandomForest, validation.set)
  cv.prediction.cancer.GBMTree         <- predict(cv.model.cancer.GBMTree, validation.set,n.trees = 400)
  
  #####################################
  cross.validates.mses.cancer$LINEAR.REGRESSION[i]      <- (sum((validation.set$Y - cv.prediction.cancer.linear.regression)^2))/nrow(validation.set)
  
  cross.validates.mses.cancer$ELASTIC.NET.1se.FALSE[i]  <- (sum((validation.set$Y - cv.prediction.cancer.ELASTIC.NET.1se.FALSE)^2))/nrow(validation.set) #MSE for the model with variables according to LASSO with lambda 1se
  cross.validates.mses.cancer$ELASTIC.NET.min.FALSE[i]  <- (sum((validation.set$Y - cv.prediction.cancer.ELASTIC.NET.min.FALSE)^2))/nrow(validation.set) #MSE for the model with variables according to LASSO with lambda min
  cross.validates.mses.cancer$RIDGE.1se.FALSE[i]        <- (sum((validation.set$Y - cv.prediction.cancer.LASSO.1se.FALSE)^2))/nrow(validation.set) #MSE for the model with variables according to LASSO with lambda 1se
  cross.validates.mses.cancer$RIDGE.min.FALSE[i]        <- (sum((validation.set$Y - cv.prediction.cancer.LASSO.min.FALSE)^2))/nrow(validation.set) #MSE for the model with variables according to LASSO with lambda min
  
  cross.validates.mses.cancer$ELASTIC.NET.1se.TRUE[i]  <- (sum((validation.set$Y - cv.prediction.cancer.ELASTIC.NET.1se.TRUE)^2))/nrow(validation.set) #MSE for the model with variables according to LASSO with lambda 1se
  cross.validates.mses.cancer$ELASTIC.NET.min.TRUE[i]  <- (sum((validation.set$Y - cv.prediction.cancer.ELASTIC.NET.min.TRUE)^2))/nrow(validation.set) #MSE for the model with variables according to LASSO with lambda min
  cross.validates.mses.cancer$RIDGE.1se.TRUE[i]        <- (sum((validation.set$Y - cv.prediction.cancer.LASSO.1se.TRUE)^2))/nrow(validation.set) #MSE for the model with variables according to LASSO with lambda 1se
  cross.validates.mses.cancer$RIDGE.min.TRUE[i]        <- (sum((validation.set$Y - cv.prediction.cancer.LASSO.min.TRUE)^2))/nrow(validation.set) #MSE for the model with variables according to LASSO with lambda min
  
  
  cross.validates.mses.cancer$RANDOM.FOREST[i]    <- (sum((validation.set$Y - cv.prediction.cancer.RandomForest)^2))/nrow(validation.set) #MSE for the model with variables according to Random Forest
  cross.validates.mses.cancer$BOOSTING[i]         <- (sum((validation.set$Y - cv.prediction.cancer.GBMTree)^2))/nrow(validation.set) #MSE for the model with variables according to Boosting Tree
}
warnings()
cross.validates.mses.summary.cancer <- cross.validates.mses.cancer[1,]
cross.validates.mses.summary.cancer <- apply(cross.validates.mses.cancer,2,mean)
cross.validates.mses.summary.cancer <- as.data.frame(cross.validates.mses.summary.cancer)

coeficient.LASSO.final <- coef(cv.model.cancer.LASSO.TRUE, s="lambda.min")
u <- as.data.frame(as.matrix(coeficient.LASSO.final))
u[,which(abs(u[,1])>1)]
coeficient.LASSO.final$[which(abs(coeficient.LASSO.final[,2]) >1)]
dim(coeficient.LASSO.final)
cancer.predictors <- row.names(as.data.frame(coeficient.LASSO.final[,1]))
cancer.predictors <- cancer.predictors[2:101]
cancer.predictors
cancer.predictors
###predicting Y###
#1. selection of variables essential for prediction using ELASTIC NET
cancer.final.train <- as.data.frame(df.train.cancer.correlated)
cancer.final.test <- as.data.frame(test.cancer)
cancer.final.test <- cancer.final.test[,cancer.predictors]

dims.cancer.train.final <- dim(cancer.final.train)

X.train.cancer.final <- cancer.final.train[,-dims.cancer.train.final[2]]
Y.train.cancer.final <- cancer.final.train[,dims.cancer.train.final[2]]

X.train.cancer.final <- as.matrix(X.train.cancer.final)
Y.train.cancer.final <- as.matrix(Y.train.cancer.final)
dim(cancer.final.test)

#2. building linear regression model using previously selected variables
model.cancer.RIDGE.final <- cv.glmnet(X.train.cancer.final, Y.train.cancer.final, nfolds = 10, standardize = F, intercept = T, alpha = 0)

#3. predicting Y
pred.cancer <- predict.cv.glmnet(object = model.cancer.RIDGE.final, newx = as.matrix(cancer.final.test), s = "lambda.min")
pred.cancer

coeficient.RIDGE.final <- coef(model.cancer.RIDGE.final, s="lambda.min")
dim(coeficient.RIDGE.final)
cancer.predictors <- row.names(as.data.frame(coeficient.RIDGE.final[,1]))
predictors.cancer <- cancer.predictors[1:100]
length(predictors.cancer)

save(pred.protein, pred.cancer, predictors.protein, predictors.cancer, file = "C:/Users/Igor/Desktop/Bioinformatyka I rok II stopnia/II sem/SAD1/projekt zaliczeniowy/IgorFilipiuk.RData")

library(MASS)
data("biopsy")

#K-fold cross-validation method: 
#1. We divide data into k (generally 10) equivocal, disjoint subsets, assigning 
#each observation to a set at random.
#2. Each subset will serve as a test set, a other observations as 
#the corresponding training set. 
#3. Instead of keeping observations in variables from the test set, 
#it's a good idea to keep the corresponding row numbers to get a training set easily.
#4. The cut function with the labels = FALSE attribute is useful here.

only.numeric <- sapply(biopsy, is.numeric)
dane_biopsy <- biopsy[, only.numeric]
dane_biopsy <- na.omit(dane_biopsy)

#K-fold cross-validation
podzial <- cut(1:nrow(dane_biopsy), 10, labels=F)
podzial <- sample(podzial)
bledy.kfold <- numeric(10) #keep estimated MSE
bledy.treningowe <- numeric(10) 

for(i in 1:10){
  train.set <- which(podzial != i)
  test.set <- which(podzial == i)
  cv.model <- lm(V2~., dane_biopsy[train.set,])
  cv.pred <- predict(cv.model, newdata=dane_biopsy[test.set,])
  
  bledy.kfold[i] <- mean((cv.pred-dane_biopsy[test.set,"V2"])^2)
  bledy.treningowe[i] <- mean((cv.model$residuals)^2)
}
bledy.kfold
bledy.treningowe

# Monte Carlo method cross-validation: 
#1. We create many (e.g. 1000) test sets, containing e.g. 1/10 of observations per set, 
#each time selecting a test set randomly (without returning).
#2. As in the previous method, each test set corresponds to a training set composed of the others
#observations.

bledy.MC <- numeric(1000) #keep estimated MSE
bledy.treningowe <- numeric(1000) 

for(i in 1:1000){
  podzial = sample(1:nrow(dane_biopsy), round(nrow(dane_biopsy)/10))
  train.set <- (1:nrow(dane_biopsy))[-podzial]
  test.set <- podzial
  cv.model <- lm(V2~., dane_biopsy, subset=train.set)
  cv.pred <- predict(cv.model, newdata=dane_biopsy[test.set,])
  
  bledy.MC[i] <- mean((cv.pred - dane_biopsy[test.set,"V2"])^2)
  bledy.treningowe[i] <- mean((cv.model$residuals)^2)
}

bledy.MC
bledy.treningowe

boxplot(bledy.kfold, bledy.MC)
par(mfrow=c(1,2))
density(bledy.kfold)
density(bledy.MC)


# Finally, compare the results of your k-fold method with the train function from the caret package. 
install.packages("caret")
library(caret)

train.control <- trainControl(method = "cv", number = 10)
kfold.train <- train(V2~., data = dane_biopsy, method = "lm", trControl = train.control)
print(kfold.train)
sqrt(mean(bledy.kfold))

kfold.train.rmse = kfold.train$resample$RMSE
my.rmse = sqrt(bledy.kfold)
my.rmse.MC = sqrt(bledy.MC)

boxplot(list(my.rmse, my.rmse.MC, kfold.train.rmse))
kfold.train$finalModel



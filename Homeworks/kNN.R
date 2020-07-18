# 0. Consider the example from the beginning of this section for data on class A and B.
# Consider a classifier that assigns class A to all observations, and
# such that assigns class B to everyone. What rA, rB, pA, pB values will both classifiers have?

simpleExample <- matrix(ncol = 3, nrow = 100)
simpleExample <- as.data.frame(simpleExample)
colnames(simpleExample) <- c("real_values", "only_As", "only_Bs")

x <- rep("A", 99)
x[100] <- "B"

simpleExample$real_values <- x
simpleExample$only_As <- rep("A", 100)
simpleExample$only_Bs <- rep("B", 100)

classes <- unique(simpleExample[,1])

#precision: TP/(TP + FP)
#recall:    TP/(TP + FN)

#classifier A
confusion_matrix_A = matrix(rep(0,length(classes)^2),nrow = length(classes))
colnames(confusion_matrix_A) <- c("Actual: A", "Actual: B")
rownames(confusion_matrix_A) <- c("Predicted: A", "Predicted: B")

confusion_matrix_A[1,1] <- nrow(simpleExample[which(simpleExample$real_values == simpleExample$only_As),])
confusion_matrix_A[1,2] <- nrow(simpleExample[which(simpleExample$real_values != simpleExample$only_As),])

p_A <- confusion_matrix_A[1,1]/sum(confusion_matrix_A[1,])
r_A <- confusion_matrix_A[1,1]/sum(confusion_matrix_A[,1])

#classifier B
confusion_matrix_B = matrix(rep(0,length(classes)^2),nrow = length(classes))
colnames(confusion_matrix_B) <- c("Actual: A", "Actual: B")
rownames(confusion_matrix_B) <- c("Predicted: A", "Predicted: B")

confusion_matrix_B[2,1] <- nrow(simpleExample[which(simpleExample$real_values != simpleExample$only_Bs),])
confusion_matrix_B[2,2] <- nrow(simpleExample[which(simpleExample$real_values == simpleExample$only_Bs),])

p_B <- confusion_matrix_B[2,2]/sum(confusion_matrix_B[2,])
r_B <- confusion_matrix_B[2,2]/sum(confusion_matrix_B[,2])

# 1. Create a confusion matrix for the results of the kNN classifier on
# test set from previous tasks. Use the table function.
wine = read.csv('C:/Users/Igor/Desktop/Bioinformatyka I rok II stopnia/II sem/SAD1/Prace Domowe/wine.csv')
wine.class = wine[,1]
wine.class <- factor(wine.class) 
nclass = length(levels(wine.class))

wine.data <- scale(wine[,-1])


# Select data at random
train.set <- sample(1:nrow(wine), 4400, rep=F)
wine.train.data <- wine.data[train.set, ]
wine.train.class <- wine.class[train.set]
wine.test.data <- wine.data[-train.set, ]
wine.test.class <- wine.class[-train.set]

# kNN classifier 
classifyKNN = function(x,K){
  distances = apply(wine.train.data, 1, function(y){sqrt(sum((x-y)^2))})
  kNNno = order(distances)[1:K]
  KNN = wine.train.class[kNNno]
  common.classes <- table(KNN)
  common.classes
  #most commong class index
  most.common.class = which.max(common.classes)
  #The value of the class over the most common index
  most.common.class = levels(wine.class)[most.common.class]
  return(most.common.class)
}

classification = apply(wine.test.data, 1, classifyKNN, K=4)

wine.test.predict = classification 
confusion <- table("Real" = wine.test.class, "Prediction" = wine.test.predict)

recall <- sapply(1:ncol(confusion), function(i){ confusion[i,i]/sum(confusion[i,])})
precision <- sapply(1:ncol(confusion), function(i){ confusion[i,i]/sum(confusion[,i])})

recall
precision

# 2. Calculate precision and recall for each class. Compare the results for 
# the 3 selected values of the parameter k using the knn function.
r.classification.3 = knn(wine.train.data, wine.test.data, wine.train.class, 3)
r.classification.5 = knn(wine.train.data, wine.test.data, wine.train.class, 5)
r.classification.10 = knn(wine.train.data, wine.test.data, wine.train.class, 10)


confusion.3 <- table("Real" = wine.test.class, "Prediction" = r.classification.3)
confusion.5 <- table("Real" = wine.test.class, "Prediction" = r.classification.5)
confusion.10 <- table("Real" = wine.test.class, "Prediction" = r.classification.10)

recall.3 <- sapply(1:ncol(confusion.3), function(i){ confusion.3[i,i]/sum(confusion.3[i,])})
precision.3 <- sapply(1:ncol(confusion.3), function(i){ confusion.3[i,i]/sum(confusion.3[,i])})

recall.5 <- sapply(1:ncol(confusion.5), function(i){ confusion.5[i,i]/sum(confusion.5[i,])})
precision.5 <- sapply(1:ncol(confusion.5), function(i){ confusion.5[i,i]/sum(confusion.5[,i])})

recall.10 <- sapply(1:ncol(confusion.10), function(i){ confusion.10[i,i]/sum(confusion.10[i,])})
precision.10 <- sapply(1:ncol(confusion.10), function(i){ confusion.10[i,i]/sum(confusion.10[,i])})

recall.3
precision.3

recall.5
precision.5

recall.10
precision.10


# Load the asbio library and bats data. Show on a scatter plot how bat's 
# arm length depends on its' age.
# Then create a linear model that explains this relationship.
# Think about what types of transformations to use.
# Conduct model diagnostics and provide the final formula for 
# the average arm length as a function of bat age.

library(asbio)
data("bats")
head(bats)
plot(bats)
naive.model = lm(forearm.length~days, data=bats)
summary(naive.model)
plot(bats$days, bats$forearm.length)

days.transf = log(bats$days, base = 10)
forearm.modification.1 = log(bats$forearm.length, base = 2)
forearm.modification.2 = bats$forearm.length^2

plot(bats$days, forearm.modification.1)
plot(days.transf, forearm.modification.1)

plot(bats$days, forearm.modification.2)
plot(days.transf, forearm.modification.2)

plot(days.transf, bats$forearm.length)


plot(naive.model)
improved.model = lm(bats$forearm.length~days.transf)
plot(improved.model)
summary(improved.model)

days.seq <- seq(0.001,0.2, by = 0.001)
forearm.seq <- seq(0.1,2, by = 0.01)

bat.matrix <- matrix(rep(0, length(days.seq)*length(forearm.seq)), nrow = length(forearm.seq))
colnames(bat.matrix) <- days.seq
rownames(bat.matrix) <-forearm.seq

max.square <-0
column.max <- 0
row.max <- 0

for(i in 1:ncol(bat.matrix)){
  for(j in 1:nrow(bat.matrix)){
    bat.y <- (bats$days + 5)^(days.seq[i])
    bat.x <- bats$forearm.length^(forearm.seq[j])
    bat.model <- lm(bat.y ~ bat.x)
    z <- summary(bat.model)
    bat.matrix[j,i] <- z$r.squared
    if(max.square < z$r.squared){
      max.square <- z$r.squared
      column.max <- days.seq[i]
      row.max <- forearm.seq[j]
    }
  }
  
}

max.square
column.max
row.max
max(bat.matrix)

#After experiments on the matrix, the best model due to R.squared is obtained when:
days.final.modification = (bats$days + 0)^0.001
forearm.final.modification =  bats$forearm.length^1.49
plot(days.final.modification, forearm.final.modification)

better.model = lm(forearm.final.modification~days.final.modification)
summary(better.model)
plot(better.model)


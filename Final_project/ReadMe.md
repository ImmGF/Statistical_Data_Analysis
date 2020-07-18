## Use data from two files:	
- protein.RData
- cancer.RData

### Read the data in R with the command load (“protein.RData“).

- The data.train training set contains a Y column for the explained variable - the level of a certain protein in a group of patients. others columns are explanatory variables. On it you need to train the appropriate one model.
- The data.test test file has no Y variable. Should be used on it learned model.

### Read the data in R with the command load (“cancer.RData“).

- The data.train training set contains, in addition to columns for predictors, a Y column indicating drug effect on tumor cell lines. The remaining columns are explanatory variables (gene expression in lines).
- The data.test test set does not have a Y variable. Use the learned model on it.

## Objective: Choose the best data model and choose the best subset of explanatory variables so that the test error is as small as possible.


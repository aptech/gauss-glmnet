# glmnetFit

## Purpose

Fits a generalized linear model (GLM) with LASSO or elasticnet penalty.

## Format

rslt = glmnetFit(y, X);
rslt = glmnetFit(y, X, family);
rslt = glmnetFit(y, X, family, ctl);

## Input
||
|:---- |:-----|
|y | Nx1 matrix, the dependent variable.
|X | NxP matrix, the dependent variables.
|family| Optional input, string, the distribution of the dependent variable. Currently "normal" is the only supported option.
|ctl| Optional input, an instance of a `glmnetControl` structure.

## Output
||
|:----- |:----
|rslt| An instance of a `glmnetResults` structure, containing the results of the estimation.

## Examples

### Basic example

```
// Get full path to dataset
fname = getGAUSSHome() $+ "pgks/glmnet/examples/pcancer.csv";

// Load all variables
X = loadd(fname);

lpsa = X[.,9];
X = X[.,1:8];

// To hold estimation results
struct glmnetResult rslt;

// Estimate model
rslt = glmnetFit(lpsa, X);

// Compute predictions for each lambda value
y_hat = X * rslt.betas + rslt.intercept;

// Compute MSE for each lambda value
mse = meanc((lpsa - y_hat).^2);
```

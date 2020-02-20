# glmnetFit

## Purpose

Fits a generalized linear model (GLM) with LASSO or elasticnet penalty.

## Format

rslt = glmnetFit(y, X);  
rslt = glmnetFit(y, X, family);  
rslt = glmnetFit(y, X, family, ctl);

## Input
|Name|Description|
|:---- |:-----|
|y | Nx1 matrix, the dependent variable.|
|X | NxP matrix, the dependent variables.|
|family| Optional input, string, the distribution of the dependent variable. Currently "normal" is the only supported option.|
|ctl| Optional input, an instance of a `glmnetControl` structure.|
| ctl.alpha|         The elasticnet mixing parameter. 0 <= alpha <= 1.  (1 - alpha) * L2 + alpha * L1. alpha = 1 for LASSO. alpha = 0 for Ridge. Default = 1.|
| ctl.nlam |         Maximum number of lambda values. Default = 100. Note this will be ignored if specific `lambdas` are supplied. |
| ctl.weights|       Observation weights. |
| ctl.penalties|     1xP vector, relative penalty for each predictor. |
| ctl.standardize|   1 to standardize predictors before estimation, 0 to use unstandardized predictors. Parameter estimates are always returned in terms of unstandardized parameters. Default = 1. 11/27/2018 -- currently unused. Set to 1 in the CPP code. |
| ctl.lambdas|       User supplied lambda values. If set, `nlam` will be ignored. Not recommended for use. |
| ctl.threshold|     Convergence threshold for each lambda solution. Default = 1e-5. Iterations stop when the maximum reduction in the criterion value as a result of each parameter update over a single pass is less than `threshold` times the null criterion value. |
| ctl.largest  |     Maximum number of variables allowed to enter the model. Default is all predictors. |
| ctl.flmin    |     if flmin < 1.0: Minimum lamda = flmin*(largest lamda value) if flmin >= 1.0: Use supplied lambda values (see `lambdas` member above). Default = 0.001. |
| ctl.exclude  |  Variables, columns of `X`, to exclude from consideration for the model. | 



## Output
|Name|Description|
|:----- |:----|
|rslt| An instance of a `glmnetResults` structure, containing the results of the estimation.|
|rslt.n_lambdas |    Scalar, the number of lambda values for which the model was estimated.|
|rslt.intercept |    1 x n_lambdas vector, containing the intercept estimate for each lambda.|
|rslt.betas     |    P x n_lambdas matrix, containing the parameter estimates for each lambda.|
|rslt.r_squared |    n_lambdas x 1 vector, containing the r_squared for each estimated model.|
|rslt.lambdas   |    n_lambdas x 1 vector, containing the values of the lambda path.|
|rslt.mse       |    mean squared error, currently not used.|
|rslt.df        |    n_lambdas x 1 vector, the degrees of freedom for each estimated model.|


## Examples

### Basic example with default settings

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

### Example 2: Constrain max number of possible variables in the model 

```
// Get full path to dataset
fname = getGAUSSHome() $+ "pgks/glmnet/examples/pcancer.csv";

// Load all variables
X = loadd(fname);

lpsa = X[.,9];
X = X[.,1:8];

// Declare control structure and fill with default settings
struct glmnetControl ctl;
ctl = glmnetControlcreate();

// Do not allow more than 3 predictors in the model
ctl.largest = 3;

// To hold estimation results
struct glmnetResult rslt;

// Estimate model
rslt = glmnetFit(lpsa, X, "normal", ctl);

// Compute predictions for each lambda value
y_hat = X * rslt.betas + rslt.intercept;
```

new;

// Make glmnet package functions available
library glmnet;

// Get file name with full path to
// the location of this file
fname = __FILE_DIR $+ "pcancer.csv";

// Load all variables
X = loadd(fname);

lpsa = X[.,9];
X = X[.,1:8];

// To hold output results;
struct glmnetResult rslt;

// Estimate model
rslt = glmnetFit(lpsa, X);

// Compute MSE from training set
rslt.mse = meanc((lpsa - (X * rslt.betas + rslt.intercept)).^2);

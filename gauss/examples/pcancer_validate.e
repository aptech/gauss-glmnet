new;

// Make glmnet package functions available
library glmnet;

// Get file name with full path to
// the location of this file
fname = __FILE_DIR $+ "pcancer.csv";

// Load all variables
X = loadd(fname);

// For repeatable sampling
rndseed 12849;

// Find indices for approximately 30% of the data
// (without replacement) to use as test set
test_idx = sampleData(seqa(1, 1, rows(X)), ceil(rows(X) * 0.3));

// Create training and set
X_test = X[test_idx,.];
X_train = delrows(X, test_idx);

// Declare 'ctl' to be a glmnetControl structure
// and fill with default values
struct glmnetControl ctl;
ctl = glmnetControlCreate();

// Use the ridge penalty
ctl.alpha = 0;

// Exclude variable number 3 from the model
ctl.exclude = 3;

// Split training set into dependent and indepedent vars
lpsa_train = X_train[.,9];
X_train = X_train[.,1:8];

// To hold output results;
struct glmnetResult rslt;

// Estimate model
rslt = glmnetFit(lpsa_train, X_train, "normal", ctl);

/*
** Compute MSE from test set
*/

lpsa_test = X_test[.,9];
X_test = X_test[.,1:8];

rslt.mse = meanc((lpsa_test - (X_test * rslt.betas + rslt.intercept)).^2);

/*
** Plot 1/lambda vs MSE
*/

struct plotControl myPlot;
myPlot = plotGetDefaults("xy");

font_name = "arial";

plotSetTextInterpreter(&myPlot, "Latex", "axes");
plotSetXLabel(&myPlot, "{1} / {\\lambda}", font_name, 14);
plotSetYLabel(&myPlot, "\\text{Test set MSE}");

plotSetTitle(&myPlot, "Prostate cancer lpsa prediction with ridge regression", font_name, 18);

plotCanvasSize("px", 800 | 600);

plotXY(myPlot, 1./rslt.lambdas, rslt.mse);

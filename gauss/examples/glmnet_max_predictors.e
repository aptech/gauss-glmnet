library glmnet;

X = loadd(__FILE_DIR $+ "x_glmnet.csv");
y = loadd(__FILE_DIR $+ "y_glmnet.csv");

// Declare control structure
// and will with default settings
struct glmnetControl ctl;
ctl = glmnetControlCreate();

// Only allow 10 predictors in the model
ctl.largest = 10;

struct glmnetResult rslt;
rslt = glmnetFit(y, X, "normal", ctl);

/*
** Plot Coefficient values vs L1 Norm
*/

b = rslt.betas';
l1 = sumr(abs(b));

// Find which parameters are non-zero
idx = selif(seqa(1, 1, cols(b)), sumc(b));

struct plotControl myPlot;
myPlot = plotGetDefaults("xy");

// Get 1 color for each non-zero parameter 
// from across the HSLuv color space
clrs = getHSLuvPalette(rows(idx));

plotSetLineColor(&myPlot, clrs);

plotSetXLabel(&myPlot, "L1 Norm", "arial", 14);
plotSetYLabel(&myPlot, "Coefficients");

// Specify predictor index in legend
plotSetTextInterpreter(&myPlot, "Latex", "legend");
plotSetLegend(&myPlot, "X_{"$+ntos(idx)$+"}", "top right inside", 1);

// Comment out line below to run with GAUSS 18
plotSetLegendBkd(&myPlot, 0);

// Make room on right side for legend
plotSetXRange(&myPlot, 0, 8);

plotSetTitle(&myPlot, "LASSO regression, limit 10 predictors", "arial", 18);

// Plot the path for each coefficient
// included in the final model
plotXY(myPlot, l1, b[.,idx]);


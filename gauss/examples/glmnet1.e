library glmnet;

// Load data
X = loadd(__FILE_DIR $+ "x_glmnet.csv");
y = loadd(__FILE_DIR $+ "y_glmnet.csv");

// Estimate model
struct glmnetResult rslt;
rslt = glmnetFit(y, X);

/*
** Plot Coefficient values vs L1 Norm
*/

// Compute L1 norm
b = rslt.betas';
l1 = sumr(abs(b));

struct plotControl myPlot;
myPlot = plotGetDefaults("xy");

// Get 20 colors from across the HSLuv color space
clrs = getHSLuvPalette(cols(b));

plotSetLineColor(&myPlot, clrs);

plotSetXLabel(&myPlot, "L1 Norm", "arial", 14);
plotSetYLabel(&myPlot, "Coefficients");

plotXY(myPlot, l1, b);


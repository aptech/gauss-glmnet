#ifndef GLMNET_H
#define GLMNET_H

//#ifdef __cplusplus__
//extern "C" {
//#endif

extern "C" void solns_(int *, int*, int*, double*, int*, int*, double *b);
extern "C" void elnet_(int *ka, double *parm, int *no, int *ni, double *x, double *y, double *w, int *jd, double *vp, double *cl, int *ne, int *nx, int *nlam, double *flmin, double *ulam, double *thr, int *isd, int *intr, int *maxit, int *lmu, double *a0, double *ca, int *ia, int *nin, double *rsq, double *alm, int *nlp, int *jerr);
extern "C" int gml_elnet(double *parm, double *no_ptr, double *ni_ptr, double *x, double *y, double *w, double *jd_ptr, double *jd_len_ptr, double *vp, double *cl, double *nx_ptr, double *nlam_ptr, double *flmin, double *ulam, double *thr, double *lmu_ptr, double *a0, double *ca, double *ia_ptr, double *nin_ptr, double *rsq, double *alm, double *nlp_ptr, double *jerr_ptr);

//#ifdef __cplusplus__
//}
//#endif

#endif

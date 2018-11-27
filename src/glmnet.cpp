#include <stdlib.h>
#include <memory.h>
#include "glmnet.h"

#ifndef max
#define max(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) ((a < b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a < b) ? (a) : (b))
#endif

int gml_elnet(double *parm, double *no_ptr, double *ni_ptr, double *x, double *y, double *w, double *jd_ptr, double *jd_len_ptr, double *vp, double *cl, double *nx_ptr, double *nlam_ptr, double *flmin, double *ulam, double *thr, double *lmu_ptr, double *a0, double *ca, double *ia_ptr, double *nin_ptr, double *rsq, double *alm, double *nlp_ptr, double *jerr_ptr)
{
    int ka = 1; // Covariance updating algorithm =1
    int isd = 1; // Standardize predictors?
    int no = static_cast<int>(*no_ptr); // First dim -> X_Dims == w_Dims == y_Dims
    int ni = static_cast<int>(*ni_ptr);
    int nx = static_cast<int>(*nx_ptr);     // PASS --- memlimit arg
    int ne = min(ni, nx); // min(x cols, memlimit)
    int nlam = static_cast<int>(*nlam_ptr); // PASS

    // named kw = nlam
    // ignored kw = ka, ne, isd
    //
    int jd_len = static_cast<int>(*jd_len_ptr);
    int *jd = new int[jd_len];
    for (int i = 0; i < jd_len; ++i)
        jd[i] = jd_ptr[i];

    int lmu = 0;
    int *ia = new int[nx];
    memset(ia, 0, nx * sizeof(int));
    int *nin = new int[nlam];
    memset(nin, 0, nlam * sizeof(int));
    int nlp = 0;
    int jerr = 0;
    int intr = 1;
    int maxit = 100000;

//    call elnet(ka,parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,
//               intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
    elnet_(&ka,parm,&no,&ni,x,y,w,jd,vp,cl,&ne,&nx,&nlam,flmin,ulam,thr,&isd,&intr,&maxit,
           &lmu,a0,ca,ia,nin,rsq,alm,&nlp,&jerr);

    // ret = lmu, *a0, *ca, *ia, *nin, *rsq, *alm, nlp, jerr
    // ignoring ia, nin, nlp for real use-case
    *lmu_ptr = static_cast<double>(lmu);
    *nlp_ptr = static_cast<double>(nlp);
    *jerr_ptr = static_cast<double>(jerr);

    double *b = new double[ni * nlam];
    memset(b, 0, ni * nlam * sizeof(double));
    solns_(&ni, &nx, &lmu, ca, ia, nin, b);
    memcpy(ca, b, ni * nlam * sizeof(double));
    delete[] b;

    for (int i = 0; i < nx; ++i)
        ia_ptr[i] = static_cast<double>(ia[i]);
    delete[] ia;

    for (int i = 0; i < nlam; ++i)
        nin_ptr[i] = static_cast<double>(nin[i]);
    delete[] nin;

    return 0;
}


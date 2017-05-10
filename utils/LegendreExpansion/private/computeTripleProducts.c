/******************************************************************************
 * computeTripleProducts.c
 *
 * Compute triple products of Legendre polynomials.
 * Could be improved further.
 ******************************************************************************/

#include <mex.h>
#include <math.h>
#include <matrix.h>

/******************************************************************************/
double Ilmn(int l, int m, int n)
{
    /* variables */
    int    i;
    double p1=1.0, p2=1.0, p3=1.0, p4=1.0;
    double Ilmn;
    
    /*compute*/
    for (i=1; i<=(l+m-n-1); i+=2) { p1 = ( (double)i/(double)(i+1) )*p1; }
    for (i=1; i<=(l+n-m-1); i+=2) { p2 = ( (double)i/(double)(i+1) )*p2; }
    for (i=1; i<=(m+n-l-1); i+=2) { p3 = ( (double)i/(double)(i+1) )*p3; }
    for (i=1; i<=(m+n+l-1); i+=2) { p4 = ( (double)(i+1)/(double)i )*p4; }
    
    return Ilmn = p1*p2*p3*p4;
    
    
}

/******************************************************************************/
void computeX(double *X, int rows, int cols, 
        double *nVals, double *mVals, double *lpr)
        
{
    
    /*variables*/
    int m, n, l;
    double mdb, ndb,ldb;
    int i, j;
    int r, started;
    double f1, f2, f3, f4, f5;
    double Itemp;
    
    
    /* Compute */
    
    l = lpr[0];
    ldb = lpr[0];
    
    for (j=0; j<cols; j++) {
        
        m = mVals[j];
        mdb = mVals[j];
        
        i = 0;
        started = 0;
        
        /* Loop over values of n*/
        while (started==0 & i<rows) {
            
            n   = nVals[i];
            ndb = nVals[i];
            r   = (l+m+n) % 2;
            
            if ( r==0 & (l+n-m)>=0 & (l+m-n)>=0 & (n+m-l)>=0 ) {
                
                /* we have started */
                started = 1;
                
                /* Set the first entry */
                Itemp = Ilmn(l,m,n);
                X[rows*j+i] = 2.0*Itemp/(1+ldb+mdb+ndb);
                
                /* Then increase n in steps of 2 to compute following entries */
                while (i+2 < rows) {
                    
                    n   = n+2;
                    ndb = ndb + 2.;
                    i   = i+2;
                    
                    r = (l+m+n)%2; 
                    if ( r==0 & l+n-m>=0 & l+m-n>=0 & n+m-l>=0 ) {

                        /* set factors for update */
                        f1 = (ldb+mdb+ndb-1)/(ldb+mdb+ndb+1);
                        if (l+m-n+1 > 0){f2 = (ldb+mdb-ndb+2)/(ldb+mdb-ndb+1); }
                        else {f2 = 1.0;}
                        
                        if (l+n-m > 0 ) {f3 = (ldb+ndb-mdb-1)/(ldb+ndb-mdb);}
                        else {f3 = 1.0;}
                        
                        if (m+n-l > 0) {f4 = (mdb+ndb-ldb-1)/(mdb+ndb-ldb);}
                        else {f4=1.0;}
                        
                        if (l+m+n-1 > 0) {f5 = (ldb+mdb+ndb)/(ldb+mdb+ndb-1);}
                        else {f5 = 1.0;}
                        
                        /* update */
                        X[rows*j+i] = f1*f2*f3*f4*f5*X[rows*j+i-2];
                        
                    }
                }/* end inner while */
            }/* end if */
            
            /* Increase i */
            i = i+1;
            
        }/* end while ~started */
    }/* end for over columns */
            
}




/******************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
        
{
    int rows, cols;
    double *nVals, *mVals, *lpr;
    double *X;
    
    
    /* check inputs */
    if (nrhs != 3){
        mexErrMsgTxt("smat: requires 3 input argument."); }
    else if (nlhs>1){
        mexErrMsgTxt("smat: requires 1 output argument."); }
    
    
    /***** assign pointers *****/
    nVals  = mxGetPr(prhs[0]);
    mVals  = mxGetPr(prhs[1]);
    lpr    = mxGetPr(prhs[2]);
    
    /***** get size *****/
    rows = mxGetN(prhs[0]);
    cols = mxGetN(prhs[1]);
    
    /***** create return argument *****/
    plhs[0] = mxCreateDoubleMatrix(rows,cols,mxREAL);
    X = mxGetPr(plhs[0]);
    
    /***** Do the computations in a subroutine *****/
    computeX(X,rows,cols,nVals,mVals,lpr);
    return;
}
/******************************************************************************/
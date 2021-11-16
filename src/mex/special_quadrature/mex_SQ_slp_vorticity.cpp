#include "SQ_common.h"

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
    
    double *xtar,*ytar;
    double *xsrc, *ysrc, *xpsrc, *ypsrc, *quad_weights, *wazp;
    double *xsrc32, *ysrc32, *xpsrc32, *ypsrc32, *quad_weights32, *wazp32;
    double *panel_breaks_x, *panel_breaks_y, *q1, *q2, *omegatar, *meanlen;
    int Ntar, Nsrc;
    double *out_omega, *nmodifs;
    double *pan2bndry, *bnds;
    bool periodic;
    
    if (nrhs != 19)
        mexErrMsgTxt("mex_do_vel_quad: incorrect number of input arguments.\n");
    
    xtar = mxGetPr(prhs[0]);
    ytar = mxGetPi(prhs[0]);
    Ntar = mxGetM(prhs[0]);
    xsrc = mxGetPr(prhs[1]);
    ysrc = mxGetPi(prhs[1]);
    Nsrc = mxGetM(prhs[1]);
    xpsrc = mxGetPr(prhs[2]);
    ypsrc = mxGetPi(prhs[2]);
    quad_weights = mxGetPr(prhs[3]);
    panel_breaks_x = mxGetPr(prhs[4]);
    panel_breaks_y = mxGetPi(prhs[4]);
    wazp = mxGetPr(prhs[5]);
    xsrc32 = mxGetPr(prhs[6]);
    ysrc32 = mxGetPi(prhs[6]);
    xpsrc32 = mxGetPr(prhs[7]);
    ypsrc32 = mxGetPi(prhs[7]);
    quad_weights32 = mxGetPr(prhs[8]);
    wazp32 = mxGetPr(prhs[9]);
    q1 = mxGetPr(prhs[10]);
    q2 = mxGetPi(prhs[10]);
    omegatar = mxGetPr(prhs[11]);
    meanlen = mxGetPr(prhs[12]);
    
    if (ysrc == NULL) {
        mwSize cs = Nsrc;
        ysrc = (double*) mxCalloc(cs,sizeof(double));
        ysrc32 = (double*) mxCalloc(cs*2,sizeof(double));
        
    }
    if (ypsrc == NULL) {
        mwSize cs = Nsrc;
        ypsrc = (double*) mxCalloc(cs,sizeof(double));
        ypsrc32 = (double*) mxCalloc(cs*2,sizeof(double));
    }
    if (panel_breaks_y == NULL) {
        mwSize cs = mxGetM(prhs[4]);
        panel_breaks_y = (double*) mxCalloc(cs,sizeof(double));
    }
    
    double * gridSolidmat, *gridSolidNx, *gridSolidNy;
    gridSolidmat = mxGetPr(prhs[13]);
    gridSolidNx = mxGetPr(prhs[14]);
    gridSolidNy = mxGetPr(prhs[15]);
    
    const mwSize *dims; //dims[0] is the number of maxNpans
    dims = mxGetDimensions(prhs[13]);
    
    pan2bndry = mxGetPr(prhs[16]);
    
    bnds = mxGetPr(prhs[17]);
    periodic = static_cast<int>(mxGetScalar(prhs[18]));
    
    plhs[0] = mxCreateDoubleMatrix(Ntar,1,mxREAL);
    out_omega = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    nmodifs = mxGetPr(plhs[1]);

    //We will only modify the velocities in this function. We could hack
    //our way past this, and write directly to us_re and us_im, but since
    //this memcpy is completely negligible time-wise we might as well do
    //things properly.
    memcpy(out_omega,omegatar,Ntar*sizeof(double));
    
    double xmin = bnds[0];
    double ymin = bnds[2];
    double xmax = bnds[1];
    double ymax = bnds[3];
    
    if (Nsrc == 0) //if no solids, return early
        return;
        
//#pragma omp parallel for
    for(int j = 0;j<Ntar;j++) {
        
        Complex nzpan[16], tz[16], tzp[16], tf[16], tn[16];
        Complex tz32[32], tzp32[32], tf32[32];
        Complex nzpan32[32], p32[33], q32[33], tn32[32];
        double twazp[16], twazp32[32];
        double tmpT[16], tmpb[16], tW32[32], tW[16];
        
        // for non-periodic domains only work on points inside the box
        bool in_box = (xtar[j] > xmin && xtar[j] < xmax && ytar[j] > ymin && ytar[j] < ymax);
        
        if (periodic || in_box)
        {
            // The point in the loop
            Complex z = Complex(xtar[j],ytar[j]);
            
            // Check how close the point z is to panels using precomp. boxes
            Complex zrel = Complex((xtar[j]-xmin)/(*meanlen),(ytar[j]-ymin)/(*meanlen));
            
            int midx = static_cast<int>(floor(real(zrel))) % static_cast<int>(*gridSolidNx);
            int midy = static_cast<int>(floor(imag(zrel))) % static_cast<int>(*gridSolidNy);
            
            int nind = static_cast<int>(midy*(*gridSolidNx) + midx);
            int solidind = nind*dims[0];
            int ind = 0; int pk;
            for (int npj = 0; npj<dims[0]; npj++) {
                // Go through all panels in the vector stored at gridSolidmat
                ind = solidind+npj;
                pk = static_cast<int>(gridSolidmat[ind])-1; // -1 here since C++ zero-based
                                
                if (pk > -1) 
                { // Only panels not equal to -1 will be considered
                    int b1 = static_cast<int>(pan2bndry[pk]);
                    Complex mid = Complex(0.5*(panel_breaks_x[pk+b1+1]+panel_breaks_x[pk+b1]),
                            0.5*(panel_breaks_y[pk+b1+1]+panel_breaks_y[pk+b1]));
                    Complex len = Complex(panel_breaks_x[pk+b1+1]-panel_breaks_x[pk+b1],
                            panel_breaks_y[pk+b1+1]-panel_breaks_y[pk+b1]);
                                        
                    z = Complex(xtar[j],ytar[j]);
                    bool check_sq = find_target_pt(z, mid, len, bnds);
                    
                    if (check_sq) {
                        Complex nz = 2*(z-mid)/len; // rescale z
                        
                        for (int k = 0; k<16; k++) {
                            tz[k] = Complex(xsrc[pk*16+k],ysrc[pk*16+k]);
                            nzpan[k] = 2*(tz[k]-mid)/len;
                        }
                        
                        Complex lg1 = log(1-nz);
                        Complex lg2 = log(-1-nz);
                        
                        // Is the point between the panel and the real axis?
                        if (imag(nz) > 0 && real(nz) > -1 && real(nz) < 1) {
                            int furthercheck = 0;
                            for (int k = 0; k<16; k++) {
                                if (imag(nzpan[k]) > imag(nz)) {
                                    furthercheck = 1;
                                    break;
                                }
                            }
                            if (furthercheck) {
                                for (int k=0; k<16; k++) {
                                    tmpT[k] = real(nzpan[k]);
                                    tmpb[k] = imag(nzpan[k]);
                                }
                                vandernewtonT(tmpT,tmpb,16);
                                double test = tmpb[15];
                                for (int k=14; k>=0; k--) {
                                    test = test*real(nz) + tmpb[k];
                                }
                                
                                if (test > imag(nz)) {
                                    lg1 -= pi*_i;
                                    lg2 += pi*_i;
                                }
                            }
                        }
                        
                        p32[0] = lg1-lg2;
                        
                        bool accurate = sq_necessary(lg1-lg2, 16, pk, z, xsrc,
                                ysrc, xpsrc, ypsrc, q1, q2, quad_weights, wazp, tz,
                                tzp, tW, tn, tf);
                        
                        // Does standard Q suffice? In that case, don't do anything!
                        if (!accurate) {
                            //No! First: attempt 32-point quadrature
                            Complex Ic16 = 0;
                            Complex rc;
                            double sum16;
                            Complex Ic1 = 0;
                            Complex Ic2 = 0;
                            Complex Ih1 = 0;
                            Complex Ih2 = 0;
                            
                            double r1, r2, dens1, dens2, rsq, rdotf, bdotf, rdotb;
                            
                            for (int k = 0; k<16; k++) {
                                // vorticity is real(int_C(q/(n*(z-tau))))/(2*pi)
                                rc = tz[k] - z;
                                Ic16 += tf[k]/(tn[k]*rc)*tW[k]*tzp[k];
                            }

                            sum16 = real(Ic16)/(2*pi);
                            
                            // upsample density
                            IPmultR(tf,tf32);
                            
                            accurate = sq_necessary(lg1-lg2, 32, pk, z, xsrc32,
                                    ysrc32, xpsrc32, ypsrc32, q1, q2,
                                    quad_weights32, wazp32, tz32, tzp32, tW32, tn32, tf);
                            
                            if (accurate) {
                                // 32 quad suffices!
                                Complex Ic32 = 0;
                                double sum32;
                                
                                for (int k=0; k<32; k++) {
                                    // vorticity is real(int_C(q/(n*(z-tau))))/(2*pi)
                                    rc = tz32[k] - z;
                                    Ic32 += tf32[k]/(tn32[k]*rc)*tW32[k]*tzp32[k];
                                }

                                sum32 = real(Ic32)/(2*pi);
                                
                                // add 32 point quadrature, take off existing 16 point quadrature
                                double modif = sum32 - sum16;
                                
                                out_omega[j] += modif;
                                
                            } else {
                                // Need to use SQ
                                for (int k=0; k<32; k++) {
                                    nzpan32[k] = 2*(tz32[k]-mid)/len;
                                }
                                
                                int sign = -1;
                                
                                // p[k] is the expansion of t^(k-1)/(t - z)
                                for (int k = 1; k<33; k++) {
                                    p32[k] = nz*p32[k-1] + (1.0-sign)/k;
                                    sign = -sign;
                                }
                                
                                //Solve the vandermonde systems to get the
                                //quadrature weights.
                                vandernewton(nzpan32,p32,32);
                                
                                Complex Ic_helsing = 0;
                                
                                for (int k = 0; k<32; k++) {
                                    // vorticity is -real(Ic[q/n](z))/(2*pi)
                                    Ic_helsing += p32[k]*tf32[k]/tn32[k];
                                }

                                double sq_prod = -real(Ic_helsing)/(2*pi);
                                
                                double modif = sq_prod - sum16;

                                out_omega[j] += modif;
                            }

                            nmodifs[0] += 1;
                        }
                    }
                }
            }
        }
    }
}
#include "SQ_common.h"

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
    
    double *xtar,*ytar;
    double *xsrc, *ysrc, *xpsrc, *ypsrc, *quad_weights, *wazp;
    double *bx, *by;
    double *xsrc32, *ysrc32, *xpsrc32, *ypsrc32, *quad_weights32, *wazp32;
    double *panel_breaks_x, *panel_breaks_y, *q1, *q2, *sigma1tar, *sigma2tar, *meanlen;
    int Ntar, Nsrc;
    double *out_sigma1, *out_sigma2, *nmodifs;
    double *pan2bndry, *bnds;
    bool periodic;
    
    if (nrhs != 20)
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
    bx = mxGetPr(prhs[11]);
    by = mxGetPi(prhs[11]);
    sigma1tar = mxGetPr(prhs[12]);
    sigma2tar = mxGetPi(prhs[12]);
    meanlen = mxGetPr(prhs[13]);
        
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
    gridSolidmat = mxGetPr(prhs[14]);
    gridSolidNx = mxGetPr(prhs[15]);
    gridSolidNy = mxGetPr(prhs[16]);
    
    const mwSize *dims; //dims[0] is the number of maxNpans
    dims = mxGetDimensions(prhs[14]);
    
    pan2bndry = mxGetPr(prhs[17]);
    
    bnds = mxGetPr(prhs[18]);
    periodic = static_cast<int>(mxGetScalar(prhs[19]));
    
    plhs[0] = mxCreateDoubleMatrix(Ntar,1,mxCOMPLEX);
    out_sigma1 = mxGetPr(plhs[0]);
    out_sigma2 = mxGetPi(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    nmodifs = mxGetPr(plhs[1]);
    
    //We will only modify the velocities in this function. We could hack
    //our way past this, and write directly to us_re and us_im, but since
    //this memcpy is completely negligible time-wise we might as well do
    //things properly.
    memcpy(out_sigma1,sigma1tar,Ntar*sizeof(double));
    memcpy(out_sigma2,sigma2tar,Ntar*sizeof(double));
    
    double xmin = bnds[0];
    double xmax = bnds[1];
    double ymin = bnds[2];
    double ymax = bnds[3];
    
    if (Nsrc == 0) //if no solids, return early
        return;
        
#pragma omp parallel for
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
            Complex btar = Complex(bx[j], by[j]);
            
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
                
                if (true)
                //if (pk > -1) 
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
                        
                        bool add_limit = false;
                        for (int k = 0; k<16; k++) {
                            if (real(tz[k]) == real(z) && imag(tz[k]) == imag(z)) {
                                if (imag(nz) > 0) {
                                    add_limit = true;
                                    p32[0] = p32[0] - 2*_i*pi;
                                }
                            }
                        }
                        
                        bool accurate = sq_necessary(lg1-lg2, 16, pk, z, xsrc,
                                ysrc, xpsrc, ypsrc, q1, q2, quad_weights, wazp, tz,
                                tzp, tW, tn, tf);
                        
                        // Does standard Q suffice? In that case, don't do anything!
                        if (add_limit || !accurate) {
                            //No! First: attempt 32-point quadrature
                            Complex Ic16 = 0;
                            Complex rc;
                            Complex sum16;
                            Complex Ic1 = 0;
                            Complex Ic2 = 0;
                            Complex Ih1 = 0;
                            Complex Ih2 = 0;
                            
                            for (int k = 0; k<16; k++) {
                                // stress is ...
                                if (real(tz[k]) == real(z) && imag(tz[k]) == imag(z)) {
                                    //skip self-interaction
                                } else {
                                    rc = z - tz[k];

                                    Ic1 += (tf[k]/(-tn[k]*rc))*tW[k]*tzp[k];
                                    Ic2 += conj(tf[k])/(-tn[k]*rc)*tW[k]*tzp[k];
                                    Ih1 += (conj(tz[k]*tf[k])/(tn[k]*rc*rc))*tW[k]*tzp[k];
                                    Ih2 += (conj(tf[k])/(tn[k]*rc*rc))*tW[k]*tzp[k];
                                }
                            }
                            
                            Ic16 = 2*btar*imag(Ic1) + _i*conj(btar*(Ic2+Ih1-conj(z)*Ih2));
                            sum16 = Ic16/(4*pi);
                            
                            // upsample density
                            IPmultR(tf,tf32);
                            
                            accurate = sq_necessary(lg1-lg2, 32, pk, z, xsrc32,
                                    ysrc32, xpsrc32, ypsrc32, q1, q2,
                                    quad_weights32, wazp32, tz32, tzp32, tW32, tn32, tf);
                            
                            if (!add_limit && accurate) {
                                // 32 quad suffices!
                                Complex Ic32 = 0;
                                Complex sum32;
                                Ic1 = 0.0;
                                Ic2 = 0.0;
                                Ih1 = 0.0;
                                Ih2 = 0.0;
                                
                                for (int k=0; k<32; k++) {
                                    // stress is ...
                                    if (real(tz32[k]) == real(z) && imag(tz[k]) == imag(z)) {
                                        //skip
                                    } else {
                                        rc = z - tz32[k];

                                        Ic1 += (tf32[k]/(-tn32[k]*rc))*tW32[k]*tzp32[k];
                                        Ic2 += conj(tf32[k])/(-tn32[k]*rc)*tW32[k]*tzp32[k];
                                        Ih1 += (conj(tz32[k]*tf32[k])/(tn32[k]*rc*rc))*tW32[k]*tzp32[k];
                                        Ih2 += (conj(tf32[k])/(tn32[k]*rc*rc))*tW32[k]*tzp32[k];
                                    }
                                }
                                
                                Ic32 = 2*btar*imag(Ic1) + _i*conj(btar*(Ic2+Ih1-conj(z)*Ih2));
                                sum32 = Ic32/(4*pi);
                                
                                // add 32 point quadrature, take off existing 16 point quadrature
                                Complex modif = sum32 - sum16;

                                out_sigma1[j] += real(modif);
                                out_sigma2[j] += imag(modif);
                                
                            } else {
                                // Need to use SQ
                                for (int k=0; k<32; k++) {
                                    nzpan32[k] = 2*(tz32[k]-mid)/len;
                                }
                                
                                int sign = -1;
                                
                                // p[k] is the expansion of t^(k-1)/(t - z)
                                // q[k] is the expansion of t^(k-1)/(t - z)^2
                                q32[0] = -1/(1+nz)-1/(1-nz);
                                for (int k = 1; k<33; k++) {
                                    p32[k] = nz*p32[k-1] + (1.0-sign)/k;
                                    q32[k] = nz*q32[k-1] + p32[k-1];
                                    sign = -sign;
                                }
                                
                                //Solve the vandermonde systems to get the
                                //quadrature weights.
                                vandernewton(nzpan32,p32,32);
                                vandernewton(nzpan32,q32,32);
                                
                                Ic1 = 0;
                                Ic2 = 0;
                                Ih1 = 0;
                                Ih2 = 0;
                                
                                for (int k = 0; k<32; k++) {
                                    Ic1 += p32[k]*tf32[k]/tn32[k];
                                    Ic2 += p32[k]*conj(tf32[k])/tn32[k];
                                    Ih1 += 2*q32[k]*conj(tz32[k]*tf32[k])/tn32[k]/len;
                                    Ih2 += 2*q32[k]*conj(tf32[k])/tn32[k]/len;
                                }

                                Complex sq_prod = (2*btar*imag(Ic1) + _i*conj(btar*(Ic2+Ih1-conj(z)*Ih2)))/(4*pi);
                                
                                Complex modif = sq_prod - sum16;
                                
                                out_sigma1[j] += real(modif);
                                out_sigma2[j] += imag(modif);
                            }

                            nmodifs[0] += 1;
                        }
                    }
                }
            }
        }
    }
}
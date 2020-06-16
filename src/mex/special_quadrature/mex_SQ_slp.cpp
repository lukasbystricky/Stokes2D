#include "SQ_common.h"

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
    
    double *xtar,*ytar;
    double *xsrc, *ysrc, *xpsrc, *ypsrc, *quad_weights, *wazp;
    double *xsrc32, *ysrc32, *xpsrc32, *ypsrc32, *quad_weights32, *wazp32;
    double *panel_breaks_x, *panel_breaks_y, *q1, *q2, *u1tar, *u2tar, *meanlen;
    int Ntar, Nsrc;
    double *out_u1, *out_u2, *nmodifs;
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
    u1tar = mxGetPr(prhs[11]);
    u2tar = mxGetPi(prhs[11]);
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
    
    plhs[0] = mxCreateDoubleMatrix(Ntar,1,mxCOMPLEX);
    out_u1 = mxGetPr(plhs[0]);
    out_u2 = mxGetPi(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    nmodifs = mxGetPr(plhs[1]);
    
    //We will only modify the velocities in this function. We could hack
    //our way past this, and write directly to us_re and us_im, but since
    //this memcpy is completely negligible time-wise we might as well do
    //things properly.
    memcpy(out_u1,u1tar,Ntar*sizeof(double));
    memcpy(out_u2,u2tar,Ntar*sizeof(double));
    
//    double xmax = pi; double xmin = -pi; double ymax = pi; double ymin = -pi;
    double xmin = bnds[0];
    double xmax = bnds[1];
    double ymin = bnds[2];
    double ymax = bnds[3];
    
    if (Nsrc == 0) //if no solids, return early
        return;
    
#pragma omp parallel for
    for(int j = 0;j<Ntar;j++) {
        
        Complex nzpan[16], tz[16], tzp[16], tf[16];
        Complex tz32[32], tzp32[32], tf32[32];
        Complex nzpan32[32], p32[33], r32[32], n32[32];
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
                
                if (pk > -1) { // Only panels not equal to -1 will be considered
                    
                    int b1 = static_cast<int>(pan2bndry[pk]);
                    Complex mid = Complex(0.5*(panel_breaks_x[pk+b1+1]+panel_breaks_x[pk+b1]),
                            0.5*(panel_breaks_y[pk+b1+1]+panel_breaks_y[pk+b1]));
                    Complex len = Complex(panel_breaks_x[pk+b1+1]-panel_breaks_x[pk+b1],
                            panel_breaks_y[pk+b1+1]-panel_breaks_y[pk+b1]);
                    
                    bool check_sq = find_target_pt(z, mid, len, bnds);
                    
                    if (check_sq) {
                        
                        Complex nz = 2*(z-mid)/len; // rescale z
                        
                        Complex lg1 = log(1-nz);
                        Complex lg2 = log(-1-nz);
                        
                        // SLP can be written as sum as sum of M3 and M4
                        // Expressions can be found in Pålsson and Tornberg, 2019
                        // https://arxiv.org/pdf/1909.12581.pdf
                        
                        Complex M3_old = 0, M4_old = 0, sum16;
                        Complex testsum = 0;
                        
                        for (int k = 0; k<16; k++) {
                            tz[k] = Complex(xsrc[pk*16+k],ysrc[pk*16+k]);
                            tzp[k] = Complex(xpsrc[pk*16+k],ypsrc[pk*16+k]);
                            tW[k] = quad_weights[pk*16+k];
                            nzpan[k] = 2*(tz[k]-mid)/len;
                            twazp[k] = wazp[pk*16+k];
                            tf[k] = Complex(q1[pk*16+k],q2[pk*16+k]);
                            testsum += tW[k]*tzp[k]/(tz[k]-z);
                            
                            M3_old += twazp[k]*(tz[k]-z)/conj(tz[k]-z)*conj(tf[k]);
                            M4_old += twazp[k]*log(abs(tz[k]-z))*tf[k];
                        }
                        
                        sum16 = -(-M3_old/(8*pi) + M4_old/(4*pi));
                        
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
                                    lg1 -= pi*_i; //OBS check signs here
                                    lg2 += pi*_i;
                                }
                            }
                        }
                        
                        p32[0] = lg1-lg2;
                        
                        // Does standard Q suffice? In that case, don't do anything!
                        if (abs(p32[0]-testsum) > 1e-14) {
                            //No! First: attempt 32-point quadrature
                            
                            Complex orig32[32];
                            IPmultR(tf,tf32);
                            
                            Complex o32sum = 0;
                            
                            for (int k=0; k<32; k++) {
                                tz32[k] = Complex(xsrc32[pk*32+k],ysrc32[pk*32+k]);
                                tzp32[k] = Complex(xpsrc32[pk*32+k],ypsrc32[pk*32+k]);
                                tW32[k] = quad_weights32[pk*32+k];
                                twazp32[k] = wazp32[pk*32+k];
                                orig32[k] = tW32[k]/(tz32[k]-z);
                                o32sum += tzp32[k]*orig32[k];
                            }
                            
                            if (abs(o32sum-p32[0]) < 1e-14) {
                                // 32 quad suffices!
                                
                                //Complex newsum_slp = 0, newsum_dlp = 0;
                                Complex M3_new = 0, M4_new = 0, sum32;
                                
                                for (int k=0; k<32; k++) {
                                    
                                    M3_new += twazp32[k]*(tz32[k]-z)/conj(tz32[k]-z)*conj(tf32[k]);
                                    M4_new += twazp32[k]*log(abs(tz32[k]-z))*tf32[k];
                                    
                                }
                                
                                sum32 = -(-M3_new/(8*pi) + M4_new/(4*pi));
                                
                                // add 32 point quadrature, take off existing 16 point quadrature
                                Complex modif = sum32 - sum16;
                                out_u1[j] += real(modif);
                                out_u2[j] += imag(modif);
                                
                            } else {
                                
                                // Need to use SQ
                                for (int k=0; k<32; k++) {
                                    n32[k] = -_i*tzp32[k]/abs(tzp32[k]);
                                    nzpan32[k] = 2*(tz32[k]-mid)/len;
                                }
                                
                                Complex gamma = 0.5*len;
                                double sign = -1;
                                
                                // p is the expansion of 1/z
                                // r is the expansion of log(|z|)
                                for (int k = 1; k<33; k++) {
                                    p32[k] = nz*p32[k-1] + (1.0-sign)/k;
                                    r32[k-1] = (lg1-sign*lg2-p32[k])/k + log(gamma)*(1.0-sign)/k;
                                    sign = -sign;
                                }
                                
                                //Solve the vandermonde systems to get the
                                //quadrature weights.
                                vandernewton(nzpan32,p32,32);
                                vandernewton(nzpan32,r32,32);
                                
                                //Complex newsum1 = 0, newsum2 = 0, newsum3 = 0, newsum4 = 0;
                                Complex M3_helsing = 0, M4_helsing = 0;
                                for (int k = 0; k<32; k++) {
                                    M3_helsing += conj(p32[k])*n32[k]*(tz32[k]-z)*conj(tf32[k]);
                                    M4_helsing += imag(gamma*conj(n32[k])*r32[k])*tf32[k];
                                }
                                
                                Complex modif = -(-_i*M3_helsing/(8*pi) + M4_helsing/(4*pi)) - sum16;
                                out_u1[j] += real(modif);
                                out_u2[j] += imag(modif);
                            }
                            
                            nmodifs[0] += 1;
                        }
                    }
                }
            }
        }
    }
}
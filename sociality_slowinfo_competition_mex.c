/***********************************************************************************************************
 * [t,SA,SB,IA,IB,EQFLAG] = sociality_slowinfo_competition_mex(t_max,a,b,c,d,q,alpha,beta,gamma,sigma,tau,eqtol,init_pop,strain_total)
 ***********************************************************************************************************/

#include <mex.h>
#include <math.h>

/***********************************
 * Constant parameter values
 ***********************************/
#define MAXSTEPS 1e6 /* Maximum number of steps for ODE solver */
#define INTERVAL 1e2 /* Check if the system is close to equilibrium */
#define EPS 1e-6 /* ODE solver tolerance */
#define TINY 1e-30 /* Constant value for solver */
#define TINY2 1e-30 /* Constant value for solver to avoid tolerance issues */
/* RK solver parameters */
#define b21 0.2
#define b31 3.0/40.0
#define b32 9.0/40.0
#define b41 0.3
#define b42 -0.9
#define b43 1.2
#define b51 -11.0/54.0
#define b52 2.5
#define b53 -70.0/27.0
#define b54 35.0/27.0
#define b61 1631.0/55296
#define b62 175.0/512.0
#define b63 575.0/13824.0
#define b64 44275.0/110592
#define b65 253.0/4096.0
#define c1 37.0/378.0
#define c3 250.0/621.0
#define c4 125.0/594.0
#define c6 512.0/1771.0
#define dc5 -277.00/14336

/*************************************
 * Define structure for model parameters
 *************************************/
struct PARAM{
    double t_max;
    double a;
    double b;
    double d;
    double q; 
    double alpha;
    double beta;
    double gamma;
    double sigma;
    double tau;
    double eqtol;
    int strain_total;
};

/*************************************
 * Function prototypes
 *************************************/
int my_rungkut (double *T, double *SA_out, double *SB_out, double *IA_out, double *IB_out, double *EQFLAG, double *init_pop, double *c, struct PARAM *p);
void rkqs(double *SA, double *SB, double *IA, double *IB,  double *DSADT, double *DSBDT, double *DIADT, double *DIBDT, double *h, double *hnext, double *SA_SCALE, double *SB_SCALE, double *IA_SCALE, double *IB_SCALE, double *c, struct PARAM *p);
void rkck(double *SA, double *SB, double *IA, double *IB,  double *DSADT, double *DSBDT, double *DIADT, double *DIBDT, double *SAout, double *SBout, double *IAout, double *IBout, double *SAerr, double *SBerr, double *IAerr, double *IBerr, double h, double *c, struct PARAM *p);
void dynamic(double *SA, double *SB, double *IA, double *IB,  double *DSADT, double *DSBDT, double *DIADT, double *DIBDT, double *c, struct PARAM *p);
double FMAX(double, double);
double FMIN(double, double);

/*************************************
 * Main function
 *************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *T, *SA, *SB, *IA, *IB, *EQFLAG, *init_pop, *c, *parameter;
    double *t_temp, *SA_temp, *SB_temp, *IA_temp, *IB_temp;
    int i, j, k, colLen, maxsteps;
    struct PARAM p;
    
    /* Allocate inputs */
    if(nrhs!=14){
        mexErrMsgTxt("Incorrect number of input arguments!\n");
    }
    else{
        parameter= mxGetPr(prhs[0]);
        p.t_max= *parameter;
        parameter= mxGetPr(prhs[1]);
        p.a= *parameter;
        parameter= mxGetPr(prhs[2]);
        p.b= *parameter;
        c= mxGetPr(prhs[3]);
        parameter= mxGetPr(prhs[4]);
        p.d= *parameter;   
        parameter= mxGetPr(prhs[5]);        
        p.q= *parameter;
        parameter= mxGetPr(prhs[6]);
        p.alpha= *parameter;
        parameter= mxGetPr(prhs[7]);
        p.beta= *parameter;
        parameter= mxGetPr(prhs[8]);
        p.gamma= *parameter;
        parameter= mxGetPr(prhs[9]);
        p.sigma= *parameter;
        parameter= mxGetPr(prhs[10]);
        p.tau= *parameter;
        parameter= mxGetPr(prhs[11]);
        p.eqtol= *parameter;
        init_pop= mxGetPr(prhs[12]);    
        parameter= mxGetPr(prhs[13]);
        p.strain_total= (int)*parameter;
    }
    maxsteps = (int)MAXSTEPS;
    
    /* Allocate memory */
    t_temp = malloc(maxsteps*sizeof(double));
    SA_temp = malloc(maxsteps*(p.strain_total)*sizeof(double));
    SB_temp = malloc(maxsteps*(p.strain_total)*sizeof(double));
    IA_temp = malloc(maxsteps*(p.strain_total)*sizeof(double));
    IB_temp = malloc(maxsteps*(p.strain_total)*sizeof(double));
    
    /* Initialise this output */           
    plhs[5] = mxCreateDoubleMatrix(1, 1, mxREAL);
    EQFLAG = mxGetPr(plhs[5]);
    EQFLAG[0] = 0;
    
    /* Call ODE solver */
    colLen = my_rungkut(t_temp, SA_temp, SB_temp, IA_temp, IB_temp, EQFLAG, init_pop, c, &p);
    
    /* Create outputs */
    plhs[0] = mxCreateDoubleMatrix(colLen, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(colLen, p.strain_total, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(colLen, p.strain_total, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(colLen, p.strain_total, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(colLen, p.strain_total, mxREAL);
    T = mxGetPr(plhs[0]);
    SA = mxGetPr(plhs[1]);
    SB = mxGetPr(plhs[2]);
    IA = mxGetPr(plhs[3]);
    IB = mxGetPr(plhs[4]);
    
    /* Copy data to outputs */
    for (i=0;i<colLen;i++){
        T[i] = t_temp[i];
        for (j=0;j<p.strain_total;j++) {
            SA[i + j*colLen] = SA_temp[i + j*maxsteps];
            SB[i + j*colLen] = SB_temp[i + j*maxsteps];
            IA[i + j*colLen] = IA_temp[i + j*maxsteps];
            IB[i + j*colLen] = IB_temp[i + j*maxsteps];
        }
    }
    
    /* Free memory */
    free(t_temp);
    free(SA_temp);
    free(SB_temp);
    free(IA_temp);
    free(IB_temp);
    
    return;
}

/*****************************************
 * ODE solver
 ****************************************/
int my_rungkut (double *T, double *SA_out, double *SB_out, double *IA_out, double *IB_out, double *EQFLAG, double *init_pop, double *c, struct PARAM *p){
    
    double *SA, *SB, *IA, *IB, *DSADT, *DSBDT, *DIADT, *DIBDT, *SA_SCALE, *SB_SCALE, *IA_SCALE, *IB_SCALE;
    double *SAMIN, *SAMAX, *SBMIN, *SBMAX, *IAMIN, *IAMAX, *IBMIN, *IBMAX, hnext[1], h[1];
    double t, nextcheck;
    int i, j, k, exitflag, count, maxsteps;
    
    /* Allocate memory */
    SA = malloc(p->strain_total*sizeof(double));
    SB = malloc(p->strain_total*sizeof(double));
    IA = malloc(p->strain_total*sizeof(double));
    IB = malloc(p->strain_total*sizeof(double));
    DSADT = malloc(p->strain_total*sizeof(double));
    DSBDT = malloc(p->strain_total*sizeof(double));
    DIADT = malloc(p->strain_total*sizeof(double));
    DIBDT = malloc(p->strain_total*sizeof(double));
    SA_SCALE = malloc(p->strain_total*sizeof(double));
    SB_SCALE = malloc(p->strain_total*sizeof(double));
    IA_SCALE = malloc(p->strain_total*sizeof(double));
    IB_SCALE = malloc(p->strain_total*sizeof(double));
    SAMIN = malloc(p->strain_total*sizeof(double));
    SBMIN = malloc(p->strain_total*sizeof(double));
    IAMIN = malloc(p->strain_total*sizeof(double));
    IBMIN = malloc(p->strain_total*sizeof(double));    
    SAMAX = malloc(p->strain_total*sizeof(double));
    SBMAX = malloc(p->strain_total*sizeof(double));
    IAMAX = malloc(p->strain_total*sizeof(double));
    IBMAX = malloc(p->strain_total*sizeof(double));
    
    /* Other parameters */
    exitflag = 1;
    count=0;
    k=1;
    h[0] = 1e-3;
    hnext[0] = 1e-3;
    t=0;
    nextcheck = INTERVAL;
    maxsteps = (int)MAXSTEPS;
    
    /* Initialise populations */
    for(i=0;i<p->strain_total;i++){
        SA[i] = init_pop[i*4];
        SB[i] = init_pop[i*4+1];
        IA[i] = init_pop[i*4+2];
        IB[i] = init_pop[i*4+3];
    }
    
    /* Initialise equilibrium arrays */
    for(i=0;i<p->strain_total;i++){
        SAMIN[i] = SA[i];
        SAMAX[i] = SA[i];
        SBMIN[i] = SB[i];
        SBMAX[i] = SB[i];
        IAMIN[i] = IA[i];
        IAMAX[i] = IA[i];
        IBMIN[i] = IB[i];
        IBMAX[i] = IB[i];
    }
    
    /* Update output */
    T[0]=t;
    for (i=0; i<p->strain_total; i++) {
        SA_out[i*maxsteps] = SA[i];
        SB_out[i*maxsteps] = SB[i];
        IA_out[i*maxsteps] = IA[i];
        IB_out[i*maxsteps] = IB[i];
    }
    
    /* Main loop: */
    do{
        /* This ensures the final step lands us on the final time point */
        if(1.1*hnext[0]>(p->t_max-t)){
            hnext[0] = p->t_max-t;
            h[0] = p->t_max-t;
            t=p->t_max;
            exitflag=0;
        }
        else{
            h[0] = hnext[0];
            t+=h[0];
        }
        if(t>=p->t_max) {
            t=p->t_max;
            exitflag=0;
        }
        /* This is where the equations are first solved */
        dynamic(SA, SB, IA, IB, DSADT, DSBDT, DIADT, DIBDT, c, p);

        /* Adjust the step size to maintain accuracy */
        for (i=0; i<p->strain_total; i++){
            SA[i] = FMAX(SA[i],0);
            SB[i] = FMAX(SB[i],0);
            IA[i] = FMAX(IA[i],0);
            IB[i] = FMAX(IB[i],0);
            SA_SCALE[i]=fabs(SA[i])+fabs(DSADT[i]*(*h))+TINY;
            SB_SCALE[i]=fabs(SB[i])+fabs(DSBDT[i]*(*h))+TINY;
            IA_SCALE[i]=fabs(IA[i])+fabs(DIADT[i]*(*h))+TINY;
            IB_SCALE[i]=fabs(IB[i])+fabs(DIBDT[i]*(*h))+TINY;
        }
        
        /* RK solver & adaptive step-size */
        rkqs(SA, SB, IA, IB, DSADT, DSBDT, DIADT, DIBDT, h, hnext, SA_SCALE, SB_SCALE, IA_SCALE, IB_SCALE, c, p);
        
        /* Make sure nothin has gone negative */
        for (i=0; i<p->strain_total; i++){
            SA[i] = FMAX(SA[i],0);
            SB[i] = FMAX(SB[i],0);
            IA[i] = FMAX(IA[i],0);
            IB[i] = FMAX(IB[i],0);
        }
        
        /* Update output */
        count++;
        T[count] = t;
        for (i=0; i<p->strain_total; i++) {
            SA_out[count + i*maxsteps] = SA[i];
            SB_out[count + i*maxsteps] = SB[i];
            IA_out[count + i*maxsteps] = IA[i];
            IB_out[count + i*maxsteps] = IB[i];
        }
    
        /* For equilibrium check */
        for (i=0; i<p->strain_total; i++){
            SAMIN[i] = FMIN(SAMIN[i],SA[i]);
            SBMIN[i] = FMIN(SBMIN[i],SB[i]);
            SAMAX[i] = FMAX(SAMAX[i],SA[i]);
            SBMAX[i] = FMAX(SBMAX[i],SB[i]);
            IAMIN[i] = FMIN(IAMIN[i],IA[i]);
            IBMIN[i] = FMIN(IBMIN[i],IB[i]);
            IAMAX[i] = FMAX(IAMAX[i],IA[i]);
            IBMAX[i] = FMAX(IBMAX[i],IB[i]);
        }
        
        /* Check if we're close to equilibrium */
        if(t>nextcheck){
            exitflag = 0;
            for (i=0; i<p->strain_total; i++){
                if(fabs(SAMAX[i]-SAMIN[i])>p->eqtol || fabs(SBMAX[i]-SBMIN[i])>p->eqtol || fabs(IAMAX[i]-IAMIN[i])>p->eqtol || fabs(IBMAX[i]-IBMIN[i])>p->eqtol){
                    exitflag = 1;
                    break;
                }
            }
            /* If close to equilibrium, then break */
            if(exitflag==0){
                t=p->t_max;
                T[count] = t;
                EQFLAG[0] = 1;
                break;
            }
            
            /* If not, then reset min/max values for each class */
            nextcheck+=INTERVAL;
            for (i=0; i<p->strain_total; i++){
                SAMIN[i] = SA[i];
                SBMIN[i] = SB[i];
                SAMAX[i] = SA[i];
                SBMAX[i] = SB[i];
                IAMIN[i] = IA[i];
                IBMIN[i] = IB[i];
                IAMAX[i] = IA[i];
                IBMAX[i] = IB[i];
            }
        }
    }while(count<(maxsteps-1) && t<=p->t_max && exitflag);
    count++;
    
    /* Free memory */
    free(SA);
    free(SB);
    free(IA);
    free(IB);
    free(DSADT);
    free(DSBDT);
    free(DIADT);
    free(DIBDT);
    free(SA_SCALE);
    free(SB_SCALE);
    free(IA_SCALE);
    free(IB_SCALE);
    free(SAMIN);
    free(SBMIN);
    free(IAMIN);
    free(IBMIN);
    free(SAMAX);
    free(SBMAX);
    free(IAMAX);
    free(IBMAX);
    
    return count;
}

/***************************************
 * This generates the adaptive step-size
 **************************************/
void rkqs(double *SA, double *SB, double *IA, double *IB,  double *DSADT, double *DSBDT, double *DIADT, double *DIBDT, double *h, double *hnext, double *SA_SCALE, double *SB_SCALE, double *IA_SCALE, double *IB_SCALE, double *c, struct PARAM *p)
{
    double *SA_temp, *SB_temp, *IA_temp, *IB_temp, *SA_err, *SB_err, *IA_err, *IB_err;
    double htemp, errmax;
    int i, j, count;
    
    /* Allocate memory */
    SA_temp = malloc(p->strain_total*sizeof(double));
    SB_temp = malloc(p->strain_total*sizeof(double));
    IA_temp = malloc(p->strain_total*sizeof(double));
    IB_temp = malloc(p->strain_total*sizeof(double));
    SA_err = malloc(p->strain_total*sizeof(double));
    SB_err = malloc(p->strain_total*sizeof(double));
    IA_err = malloc(p->strain_total*sizeof(double));
    IB_err = malloc(p->strain_total*sizeof(double));
    
    count = 0;
    for(;;)
    {
        rkck(SA, SB, IA, IB, DSADT, DSBDT, DIADT, DIBDT, SA_temp, SB_temp, IA_temp, IB_temp, SA_err, SB_err, IA_err, IB_err, *h, c, p);
        
        errmax= 0.0;
        for(i=0;i<p->strain_total;i++){
            errmax= FMAX(errmax, fabs(SA_err[i]/(SA_SCALE[i])));
            errmax= FMAX(errmax, fabs(SB_err[i]/(SB_SCALE[i]))); 
            errmax= FMAX(errmax, fabs(IA_err[i]/(IA_SCALE[i])));
            errmax= FMAX(errmax, fabs(IB_err[i]/(IB_SCALE[i])));   
        }
        errmax/= EPS;
        if(errmax<=1.0) break;
        htemp= 0.9*(*h)*pow(errmax, -0.25);
        *h= (*h>=0.0 ? FMAX(htemp, 0.1*(*h)) : FMIN(htemp, 0.1*(*h)));
        count++;
            
        if(count>1e4){
            printf("%f\n",errmax);
            mexErrMsgTxt("stuck in loop!\n");
            break;
        }
    }    
    if(errmax > 1.89E-4) {
        *hnext= 0.9*(*h)*pow(errmax, -0.2);
    }
    else {
        *hnext= 5.0*(*h);
    }    
    *hnext = FMAX(*hnext, p->t_max/MAXSTEPS);
    
    for(i=0;i<p->strain_total;i++){
        SA[i] = SA_temp[i];
        SB[i] = SB_temp[i];
        IA[i] = IA_temp[i];
        IB[i] = IB_temp[i];
    }
    
    /* Free memory */
    free(SA_temp);
    free(SB_temp);
    free(IA_temp);
    free(IB_temp);
    free(SA_err);
    free(SB_err);
    free(IA_err);
    free(IB_err);
}

/**************************************
 * Standard RK solver
 **************************************/
void rkck(double *SA, double *SB, double *IA, double *IB,  double *DSADT, double *DSBDT, double *DIADT, double *DIBDT, double *SAout, double *SBout, double *IAout, double *IBout, double *SAerr, double *SBerr, double *IAerr, double *IBerr, double h, double *c, struct PARAM *p){
    int i, j;
    double *SAk1, *SAk2, *SAk3, *SAk4, *SAk5, *SAk6, *SAtemp;
    double *SBk1, *SBk2, *SBk3, *SBk4, *SBk5, *SBk6, *SBtemp;
    double *IAk1, *IAk2, *IAk3, *IAk4, *IAk5, *IAk6, *IAtemp;
    double *IBk1, *IBk2, *IBk3, *IBk4, *IBk5, *IBk6, *IBtemp;
    double dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0, dc4=c4-13525.0/55296.0, 
            dc6=c6-0.25;
    
    /* Allocate memory */    
    SAk1 = malloc(p->strain_total*sizeof(double));
    SAk2 = malloc(p->strain_total*sizeof(double));
    SAk3 = malloc(p->strain_total*sizeof(double));
    SAk4 = malloc(p->strain_total*sizeof(double));
    SAk5 = malloc(p->strain_total*sizeof(double));
    SAk6 = malloc(p->strain_total*sizeof(double));
    SAtemp = malloc(p->strain_total*sizeof(double));
    
    SBk1 = malloc(p->strain_total*sizeof(double));
    SBk2 = malloc(p->strain_total*sizeof(double));
    SBk3 = malloc(p->strain_total*sizeof(double));
    SBk4 = malloc(p->strain_total*sizeof(double));
    SBk5 = malloc(p->strain_total*sizeof(double));
    SBk6 = malloc(p->strain_total*sizeof(double));
    SBtemp = malloc(p->strain_total*sizeof(double));
    
    IAk1 = malloc(p->strain_total*sizeof(double));
    IAk2 = malloc(p->strain_total*sizeof(double));
    IAk3 = malloc(p->strain_total*sizeof(double));
    IAk4 = malloc(p->strain_total*sizeof(double));
    IAk5 = malloc(p->strain_total*sizeof(double));
    IAk6 = malloc(p->strain_total*sizeof(double));
    IAtemp = malloc(p->strain_total*sizeof(double));
    
    IBk1 = malloc(p->strain_total*sizeof(double));
    IBk2 = malloc(p->strain_total*sizeof(double));
    IBk3 = malloc(p->strain_total*sizeof(double));
    IBk4 = malloc(p->strain_total*sizeof(double));
    IBk5 = malloc(p->strain_total*sizeof(double));
    IBk6 = malloc(p->strain_total*sizeof(double));  
    IBtemp = malloc(p->strain_total*sizeof(double));  
    
    for(i=0;i<p->strain_total;i++){
        SAtemp[i] = SA[i] + b21*h*DSADT[i];
        SBtemp[i] = SB[i] + b21*h*DSBDT[i];
        IAtemp[i] = IA[i] + b21*h*DIADT[i];
        IBtemp[i] = IB[i] + b21*h*DIBDT[i];
    }
    dynamic(SAtemp, SBtemp, IAtemp, IBtemp, SAk2, SBk2, IAk2, IBk2, c, p);
    
    for(i=0;i<p->strain_total;i++){
        SAtemp[i] = SA[i]+h*(b31*DSADT[i]+b32*SAk2[i]);
        SBtemp[i] = SB[i]+h*(b31*DSBDT[i]+b32*SBk2[i]);
        IAtemp[i] = IA[i]+h*(b31*DIADT[i]+b32*IAk2[i]);
        IBtemp[i] = IB[i]+h*(b31*DIBDT[i]+b32*IBk2[i]);
    }    
    dynamic(SAtemp, SBtemp, IAtemp, IBtemp, SAk3, SBk3, IAk3, IBk3, c, p);
    
    for(i=0;i<p->strain_total;i++){
        SAtemp[i] = SA[i]+h*(b41*DSADT[i]+b42*SAk2[i]+b43*SAk3[i]);
        SBtemp[i] = SB[i]+h*(b41*DSBDT[i]+b42*SBk2[i]+b43*SBk3[i]);
        IAtemp[i] = IA[i]+h*(b41*DIADT[i]+b42*IAk2[i]+b43*IAk3[i]);
        IBtemp[i] = IB[i]+h*(b41*DIBDT[i]+b42*IBk2[i]+b43*IBk3[i]);
    }
    dynamic(SAtemp, SBtemp, IAtemp, IBtemp, SAk4, SBk4, IAk4, IBk4, c, p);
    
    for(i=0;i<p->strain_total;i++){
        SAtemp[i] = SA[i]+h*(b51*DSADT[i]+b52*SAk2[i]+b53*SAk3[i]+b54*SAk4[i]);
        SBtemp[i] = SB[i]+h*(b51*DSBDT[i]+b52*SBk2[i]+b53*SBk3[i]+b54*SBk4[i]);
        IAtemp[i] = IA[i]+h*(b51*DIADT[i]+b52*IAk2[i]+b53*IAk3[i]+b54*IAk4[i]);
        IBtemp[i] = IB[i]+h*(b51*DIBDT[i]+b52*IBk2[i]+b53*IBk3[i]+b54*IBk4[i]);
    }
    dynamic(SAtemp, SBtemp, IAtemp, IBtemp, SAk5, SBk5, IAk5, IBk5, c, p);
    
    for(i=0;i<p->strain_total;i++){
        SAtemp[i] = SA[i]+h*(b61*DSADT[i]+b62*SAk2[i]+b63*SAk3[i]+b64*SAk4[i]+b65*SAk5[i]);
        SBtemp[i] = SB[i]+h*(b61*DSBDT[i]+b62*SBk2[i]+b63*SBk3[i]+b64*SBk4[i]+b65*SBk5[i]);        
        IAtemp[i] = IA[i]+h*(b61*DIADT[i]+b62*IAk2[i]+b63*IAk3[i]+b64*IAk4[i]+b65*IAk5[i]);
        IBtemp[i] = IB[i]+h*(b61*DIBDT[i]+b62*IBk2[i]+b63*IBk3[i]+b64*IBk4[i]+b65*IBk5[i]);
    }
    dynamic(SAtemp, SBtemp, IAtemp, IBtemp, SAk6, SBk6, IAk6, IBk6, c, p);
    
    for(i=0;i<p->strain_total;i++){
        SAout[i]= SA[i]+h*(c1*DSADT[i]+c3*SAk3[i]+c4*SAk4[i]+c6*SAk6[i]);
        SAerr[i]= h*(dc1*DSADT[i]+dc3*SAk3[i]+dc4*SAk4[i]+dc5*SAk5[i]+dc6*SAk6[i]);
        SBout[i]= SB[i]+h*(c1*DSBDT[i]+c3*SBk3[i]+c4*SBk4[i]+c6*SBk6[i]);
        SBerr[i]= h*(dc1*DSBDT[i]+dc3*SBk3[i]+dc4*SBk4[i]+dc5*SBk5[i]+dc6*SBk6[i]);
        IAout[i]= IA[i]+h*(c1*DIADT[i]+c3*IAk3[i]+c4*IAk4[i]+c6*IAk6[i]);
        IAerr[i]= h*(dc1*DIADT[i]+dc3*IAk3[i]+dc4*IAk4[i]+dc5*IAk5[i]+dc6*IAk6[i]);
        IBout[i]= IB[i]+h*(c1*DIBDT[i]+c3*IBk3[i]+c4*IBk4[i]+c6*IBk6[i]);
        IBerr[i]= h*(dc1*DIBDT[i]+dc3*IBk3[i]+dc4*IBk4[i]+dc5*IBk5[i]+dc6*IBk6[i]);
    }

    /* Free memory */
    free(SAk1);
    free(SAk2);
    free(SAk3);
    free(SAk4);
    free(SAk5);
    free(SAk6);
    free(SBk1);
    free(SBk2);
    free(SBk3);
    free(SBk4);
    free(SBk5);
    free(SBk6);
    free(SAtemp);
    free(SBtemp);
    free(IAk1);
    free(IAk2);
    free(IAk3);
    free(IAk4);
    free(IAk5);
    free(IAk6);
    free(IBk1);
    free(IBk2);
    free(IBk3);
    free(IBk4);
    free(IBk5);
    free(IBk6);
    free(IAtemp);
    free(IBtemp);
}

/**************************************
 * Population and evolutionary dynamics
 **************************************/
void dynamic(double *SA, double *SB, double *IA, double *IB,  double *DSADT, double *DSBDT, double *DIADT, double *DIBDT, double *c, struct PARAM *p){
    
    int i, j;
    double *Ni, N, infection_mult, information_mult;
    
    /* Allocate memory */
    Ni = malloc(p->strain_total*sizeof(double));
    
    /* Population sums */
    N = 0;
    for(i=0;i<p->strain_total;i++){
        Ni[i] = SA[i] + SB[i] + IA[i] + IB[i];
        N += Ni[i];
    }
    
    /* Infections and information */
    infection_mult = 0;
    information_mult = 0;
    for(i=0;i<p->strain_total;i++){
        infection_mult += c[i]*(IA[i]+IB[i]);
        information_mult += c[i]*(SB[i]+IB[i]);
    }
    infection_mult = p->beta*FMAX(infection_mult,TINY2)/FMAX(N,TINY2);
    information_mult = p->tau*FMAX(information_mult,TINY2)/FMAX(N,TINY2);
    
    /* ODEs */
    for(i=0;i<p->strain_total;i++){
        DSADT[i] = (p->b - p->q*N)*Ni[i] - c[i]*SA[i]*(infection_mult + information_mult) - p->d*SA[i] + p->gamma*IA[i] + p->sigma*SB[i];
        DSBDT[i] = c[i]*SA[i]*information_mult - c[i]*SB[i]*infection_mult - p->a*p->d*SB[i] + p->gamma*IB[i] - p->sigma*SB[i];
        DIADT[i] = c[i]*SA[i]*infection_mult - c[i]*IA[i]*information_mult - (p->d + p->alpha + p->gamma)*IA[i] + p->sigma*IB[i];
        DIBDT[i] = c[i]*SB[i]*infection_mult + c[i]*IA[i]*information_mult - (p->a*p->d + p->alpha + p->gamma + p->sigma)*IB[i];
    }
    
    /* Free memory */
    free(Ni);
}

/***************************************
 * Return maximum of two inputs
 ***************************************/
double FMAX(double l, double r)
{
    if(l>r)return l;
    else   return r;
}

/***************************************
 * Return minimum of two inputs
 ***************************************/
double FMIN(double l, double r)
{
    if(l<r)return l;
    else   return r;
}

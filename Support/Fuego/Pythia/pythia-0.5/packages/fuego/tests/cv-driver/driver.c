/* driver.c */

#include <stdio.h>

void burncv_( double *dt, double *phi, double *phihist,
               double *thist, double *rho, double *ene, 
               int *nsp, int * itmax, int * nsteps);

void fgindx_(int * iwrk, double *rwrk, int * mm, int * kk, 
             int * ii, int * nfit );
void fgxty_(double * x, int * iwrk, double *rwrk, double * y);
void fephity_(double * phi, int * iwrk, double *rwrk, double * y);
void feytphi_(double * y, int * iwrk, double *rwrk, double * phi);
void feeytt_(double *ene, double * y, int * iwrk, double *rwrk, 
             double * tout);

int main()
{

    double phi[100], phihist[200000], thist[2000];
    double x[100],y[100], tfinal, temperature[2000];
    double rckwrk[1];
    int    ickwrk[1];
    FILE*  output;
    
    double dt = 1000;  // integrate to equil
    int nsteps;
    double rho   = 4.58e-4;
    double ene   = 1.28e10;
    int    itmax = 2000;
    
    int id; /* Loop Counter */ 
    
    /* get number of species */
    int idum1, ns, idum2, idum3;
    fgindx_(ickwrk, rckwrk, &idum1, &ns, &idum2, &idum3);

    /* set initial conditions 
       (only valid for the ordering in HydrogenOxygen.ck2) */
    x[0] = 2; /* H2 */
    x[2] = 1; /* O2 */
    x[9] = 3.72; /* N2 */

    fgxty_( x, ickwrk, rckwrk, y);
    feytphi_( y, ickwrk, rckwrk, phi );

    /* integrate */
    burncv_( &dt, phi, phihist, thist, &rho, &ene, &ns, 
              &itmax, &nsteps);

    /* Postprocessing starts here */
    printf(" NSteps: %d \n ", nsteps);
    fephity_( phi, ickwrk, rckwrk, y); 
    feeytt_( &ene, y, ickwrk, rckwrk, &tfinal);
    printf("Equilibrium Temperature: %g\n", tfinal); 

    /* Output temperature history */
    printf("Computing temperature history\n");
    printf("and output to temperature.out\n");
    output = fopen("temperature.out","w");
    for (id = 0; id<nsteps; id++) {
        fephity_( &phihist[id*ns], ickwrk, rckwrk, y);
        feeytt_( &ene, y, ickwrk, rckwrk, &temperature[id]);
        fprintf(output,"%g %g\n",thist[id], temperature[id]);
    }

    return 0;
}

/* End of file */

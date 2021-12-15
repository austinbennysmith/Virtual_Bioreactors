#include "grid/cartesian1D.h"
#include "run.h"

scalar T[];

// T[right] = neumann(0.0);
// T[left] = neumann(0.0);

FILE * fp_params;

double mydt = 0.015;
float U = 2.0;

// See: http://basilisk.fr/sandbox/M1EMN/BASIC/heat.c
int main() {
  periodic(right);
  L0 = 10.;
  X0 = -L0/2;
  // N = 100;
  // CFL = -4.0;
  // DT = (L0/N)*(L0/N)/2 ;
// #define EPS 0.1 
  
  {
  char params[200];
  sprintf(params, "params.txt");
  fp_params=fopen(params, "w");
  }

  fprintf(fp_params, "Advection Stability Condition: |U*dt/dx|<=1, where dx=L0/N \n");
  fprintf(fp_params, "|a*dt/dx|= %g \n", fabs(U*mydt/(L0/N)));
  fprintf(fp_params, "Diffusion Stability Condition: b*mu<=0.5, where b=diffusion coefficient (=1 here) and mu=dt/dx^2, where dx=L0/N \n");
  fprintf(fp_params, "b*mu=mu= %g \n", mydt/(sq(L0/N)));
  fclose(fp_params);

  run();
}

event init (t = 0) {
  periodic(right);
  foreach()
    // T[] =  1./EPS*(fabs(x)<EPS)/2;
  	// T[] = sin((pi/5)*x)+5;
  	T[] = exp(-x*x);
  	// T[] = sin((pi/5.0)*x);
  boundary ({T});
}

event printdata (t += 1.0; t < 1000.0) {
  foreach()
    fprintf (stdout, "%g %g %g\n", x, T[], t);
  fprintf (stdout, "\n\n");
}

event integration (i++) {
  // double dt = DT;
  scalar dT[],q[];
  dt = mydt;
  foreach() {
  	// Helpful sentence from a paper (at https://www.researchgate.net/profile/Wim-Desmet/publication/251548727_Finite_Volume_Convective_Flux_Reconstruction_Using_High-Order_Skew-Symmetric-Like_Schemes/links/02e7e52f1183f86850000000/Finite-Volume-Convective-Flux-Reconstruction-Using-High-Order-Skew-Symmetric-Like-Schemes.pdf):
  	// "in the absence of sources and sinks, the rate of change in time of the sum of the conserved variable phi in a volume V equals the sum of the flux on the surface S of the volume"
  	// To get the flux q[], I want to estimate the rate of change of T (aka, get T'). To do that, I need to do a 2nd-order centered estimate for T'.
  	// Info on approximations of derivatives: https://www.dam.brown.edu/people/alcyew/handouts/numdiff.pdf

  // My ACTUAL advection-diffusion scheme (FINITE VOLUMES, forward time, centered space for diffusion, Lax-Wendroff for advection):
    q[] = (T[-1,0]-T[0,0])/Delta - U*(T[0,0]+T[-1,0])/2.0 - ((sq(U)*dt)/(2.0*Delta))*(T[0,0]-T[-1,0]);

    // Lax-Wendroff Advection:
    // q[] = - U*(T[0,0]+T[-1,0])/2.0 - ((sq(U)*dt)/(2.0*Delta))*(T[0,0]-T[-1,0]);

 	// Centered space diffusion:
 	// q[] = (T[-1,0]-T[0,0])/Delta;
  }
  boundary ({q});
  foreach() {
    dT[] = ( q[0,0]  - q[1,0] )/Delta;
    // Alternative method for advection-diffusion (FINITE DIFFERENCES, forward time, centered space for diffusion, Lax-Wendroff for advection):
    // dT[] = (T[1,0]-2*T[0,0]+T[-1,0])/(sq(Delta)) + (U/(2*Delta))*(T[1,0]-T[-1,0]) + ((sq(U)*dt)/(2*sq(Delta)))*(T[1,0]-2*T[0,0]+T[-1,0]);

    // Lax-Wendroff Advection:
    // dT[] = (U/(2*Delta))*(T[1,0]-T[-1,0]) + ((sq(U)*dt)/(2*sq(Delta)))*(T[1,0]-2*T[0,0]+T[-1,0]);

    // Centered space diffusion:
    // dT[] = (T[1,0]-2*T[0,0]+T[-1,0])/(sq(Delta));
  }
  foreach()
    T[] += dt*dT[];
  boundary ({T});
}
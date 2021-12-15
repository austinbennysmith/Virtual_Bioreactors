#include "run.h"
#include "view.h"

scalar T[];

// T[right] = 0.0;
// T[left] = 0.0;
// T[top] = 0.0;
// T[bottom] = 0.0;

// T[right] = neumann(0.0);
// T[left] = neumann(0.0);
// T[top] = neumann(0.0);
// T[bottom] = neumann(0.0);

FILE * fp_params;

const double tmax = 100.0;

double dt;
float U = 2.0;

// See: http://basilisk.fr/sandbox/M1EMN/BASIC/heat.c
int main() {
  // init_grid(1<<LEVEL);
  periodic(right);
  periodic(top);
  L0 = 100.;
  X0 = Y0 = -L0/2;
  N=128;
  // N = 100;
  // CFL = -4.0;
  // DT = (L0/N)*(L0/N)/2 ;
// #define EPS 0.1 
  
  {
  char params[200];
  sprintf(params, "params.txt");
  fp_params=fopen(params, "w");
  }

  fprintf(fp_params, "Delta: %g \n", L0/N);
  fprintf(fp_params, "dt: %g \n", dt);
  fprintf(fp_params, "CFL: %g \n", CFL);
  fclose(fp_params);

  run();
}

event init (t = 0) {
  // mask (y > 5.0 ? top : none);
  // mask (y < -5.0 ? bottom : none);
  // periodic(right);
  foreach()
    // T[] =  1./EPS*(fabs(x)<EPS)/2;
  	// T[] = sin((pi/5)*x)+5;
  	T[] = exp(-(0.1*x*x+0.1*y*y));
  	// T[] = sin((pi/5.0)*x);
  boundary ({T});
}

event printdata (t += 1.0; t < tmax) {
  foreach()
    fprintf (stdout, "%g %g %g %g\n", x, y, T[], t);
  fprintf (stdout, "\n\n");
}

scalar dT[],qh[],qv[];
event integration (i++) {
  // double dt = DT;
  dt = 0.05;
  foreach() {
  // // 	// Helpful sentence from a paper (at https://www.researchgate.net/profile/Wim-Desmet/publication/251548727_Finite_Volume_Convective_Flux_Reconstruction_Using_High-Order_Skew-Symmetric-Like_Schemes/links/02e7e52f1183f86850000000/Finite-Volume-Convective-Flux-Reconstruction-Using-High-Order-Skew-Symmetric-Like-Schemes.pdf):
  // // 	// "in the absence of sources and sinks, the rate of change in time of the sum of the conserved variable phi in a volume V equals the sum of the flux on the surface S of the volume"
  // // 	// To get the flux q[], I want to estimate the rate of change of T (aka, get T'). To do that, I need to do a 2nd-order centered estimate for T'.
  // // 	// Info on approximations of derivatives: https://www.dam.brown.edu/people/alcyew/handouts/numdiff.pdf
  // //   // q[]=(T[1,0]-T[-1,0])/(2*Delta); // heat equation alone
  // //   // q[]=(T[0,0] - T[-1,0])/(Delta)- U*(T[1, 0]+T[-1, 0])/2.;
  // //   // q[]=-(T[0,0] - T[-1,0])/Delta;

    // My ACTUAL advection-diffusion scheme:
    qh[] = (T[-1,0]-T[0,0])/Delta - U*(T[0,0]+T[-1,0])/2.0 - ((sq(U)*dt)/(2.0*Delta))*(T[0,0]-T[-1,0]);
    qv[] = (T[0,-1]-T[0,0])/Delta;

 	// // // My OLD advection-diffusion scheme:
  	// qh[] = -(T[0,0]-T[-1,0])/Delta + U*((T[0,0]+T[-1,0])/2.0 - ((T[0,0]-T[-1,0])/2));
   //  qv[] = -(T[0,1]-T[0,-1])/Delta;

  // // 	// // Advection alone:
  	// qh[] = U*((T[0,0]+T[-1,0])/2.0 - ((T[0,0]-T[-1,0])/2));
   //  qv[] = U*((T[0,0]+T[0,-1])/2.0 - ((T[0,0]-T[0,-1])/2));

  // // 	//Diffusion alone:
  // 	// qh[] = -(T[0,0]-T[-1,0])/Delta;
  //  //  qv[] = -(T[0,0]-T[0,-1])/Delta;

  // // // q[]=U*(T[0, 0]+T[-1, 0])/2.;
  }
  boundary ({qh});
  boundary({qv});
  foreach() {
    dT[] = ( qh[0,0]  - qh[1,0] )/Delta + ( qv[0,0]  - qv[0,1] )/Delta;

    //   // Alternative method for advection-diffusion (FINITE DIFFERENCES, forward time, centered space for diffusion, Lax-Wendroff for advection):
    // dT[] = (T[1,0]-2*T[0,0]+T[-1,0])/(sq(Delta)) + (U/(2*Delta))*(T[1,0]-T[-1,0]) + ((sq(U)*dt)/(2*sq(Delta)))*(T[1,0]-2*T[0,0]+T[-1,0]) + (T[0,1]-2*T[0,0]+T[0,-1])/(sq(Delta)) + (U/(2*Delta))*(T[0,1]-T[0,-1]) + ((sq(U)*dt)/(2*sq(Delta)))*(T[0,1]-2*T[0,0]+T[0,-1]);
  }

  foreach()
    T[] = T[] + dt*dT[];
  boundary ({T});
}

event profiles (t = 0; t+=1.0; t<=tmax)
{
  FILE * fp = fopen("xprof", "a");
  for (double y = -L0/2; y <= L0/2; y += 0.01)
    fprintf (fp, "%g %g %g\n", t, y, interpolate (T, 0, y));
  fclose (fp);
  
  fp = fopen("yprof", "a");
  for (double x = -L0/2; x <= L0/2; x += 0.01)
    fprintf (fp, "%g %g %g\n", t, x, interpolate (T, x, 0));
  fclose (fp);
}

event Tmovie (t+=10.0, t<tmax)
{
 clear();
 // cells(lc={0.5,0.5,0.5}, lw=0.5);
 // draw_vof ("cs", "fs", filled=-1, fc = {1.0,1.0,1.0}, lw=2);
 squares("T", min=0.0, max=0.5, linear=true, map=cool_warm);
 // cells();
 save ("Tmovie.mp4");
}

// event ymovie (t+=0.1)
// {
//  clear();
//  // cells(lc={0.5,0.5,0.5}, lw=0.5);
//  // draw_vof ("cs", "fs", filled=-1, fc = {1.0,1.0,1.0}, lw=2);
//  squares("u.y", spread=-1, linear=true, map=cool_warm);
//  // cells();
//  save ("ymovie.mp4");
// }

event gfsview (t += 1.0, t<=tmax) { // RC
    char name_gfs[200];
    sprintf(name_gfs,"Slice-%0.1f.gfs",t);

    FILE* fp_gfs = fopen (name_gfs, "w");
    output_gfs(fp_gfs);
    fclose(fp_gfs);
}
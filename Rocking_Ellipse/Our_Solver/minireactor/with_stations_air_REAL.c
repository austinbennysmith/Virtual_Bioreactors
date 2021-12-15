#define FILTERED //RC
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2)) //RC

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
#include "fractions.h"

scalar TT[];
// scalar * tracers = {T};
double b = 0.00002; // Diffusion coefficient

double thetaNOW;
double omegaNOW;

// The following doubles will be used in the acceleration event to avoid really long lines of code:
double gravityX;
double gravityY;
double coriolisX;
double coriolisY;
double centripetalX;
double centripetalY;

FILE *fp1 ;

#define MAXLEVEL 9 // RC was 4, needs to be bigger to capture the setup

// Defining Froude separately from the other nondimensional quantities because it will be used to set up the reference Velocity:
#define Fr 1.0 // Froude number

// DIMENSIONAL QUANTITIES:
#define rhoWater 1000.0 // kg/m^3
#define rhoAir 1.225 // kg/m^3
#define muWater 0.001 // approximatley the viscosity of water
// #define muAir 1.81e-5 // approximateley the viscosity of air
#define muAir 1.81e-4 // visc. of air *10
// #define muAir 1.81e-3 // visc. of air *10
#define sig 0.0728  //surface tension of water
const double semiminor = 1.0; // semiminor axis (cm)
const double semimajor = semiminor*3.0; // semimajor axis (cm)
const double maxDegrees = 7.0; // degrees through which the reactor rotates
double maxRads = maxDegrees*(3.14159265/180.0);
#define dimensional_period 1.0 // rocking period in actual seconds (will convert to nondimensional time units below for the simulation)
#define refLength 0.01  // semiminor axis length (m). Could define in terms of semiminor? nah
#define refVelocity (Fr*sqrt(9.8*refLength))  // Reference length, defined in terms of the Froude number
#define refTime (refLength/refVelocity)

// DIMENSIONLESS QUANTITIES:
const double period = (dimensional_period/refTime); // how many NONDIMENSIONAL TIME UNITS it takes to go through a complete rocking cycle. DEFINE IN TERMS OF dimensional_period
double BB = (2.0*3.14159265)/period; // constant used in rocking motion equations below
#define Re (rhoWater*refVelocity*refLength/muWater)  // Reynolds number
#define We (rhoWater*pow(refVelocity,2)*refLength/sig)
#define rho_ratio (rhoAir/rhoWater)
#define mu_ratio (muAir/muWater)
// Should include Peclet numbers eventually

const double tmax = (5.0*period); // Runs for 5 oscillation periods (defined in dimensionless time units)
const double tplus = (period/100.0); // events that run at set t invervals will have those intervals defined by this quantity

FILE *fp_params;

// RC
FILE * fp_stats;

FILE * fp_interface;

u.t[top] = dirichlet(0.0);
u.n[top] = dirichlet(0.0);

u.t[bottom] = dirichlet(0.0);
u.n[bottom] = dirichlet(0.0);

u.t[left] = dirichlet(0.0);
u.n[left] = dirichlet(0.0);

u.t[right] = dirichlet(0.0);
u.n[right] = dirichlet(0.0);

// I use the following variables to determine the cell length/width
const double ymax = semiminor+0.5; // ymax of domain (used in masking & profiling, etc.)
const double ymin = -semiminor-0.5; // ymin of domain (used in masking & profiling, etc.)
const double xmax = semimajor+3.0; // xmax of domain (used in profiling)
const double xmin = -semimajor-3.0; // xmin of domain (used in profiling)
// long int sizeNow; // Total # of cells at t=0.0 (as long int)
// double sizeNowDouble; // Total # of cells at t=0.0 (as double)
// double numCellsX; // # of cells in x (long) direction
// double numCellsY; // # of cells in y (short) direction
double DeltaX; // Length of smallest cell in x direction
// double DeltaY; // Length of smallest cell in y direction
// double ylength; // Domain length in y direction
// double xlength; // Domain length in x direction
double Diffusion_Stability; // test for stability of the FTCS diffusion scheme
double Advection_Stability; // test for stability of the Lax-Wendroff advection scheme

// scalar s1[];

// void draw_frame(char * fname)
// {
// 	cells();
// 	squares("s1");
// 	save("fname.png");
// }

scalar circle[];

int main() {
  L0 = 2.0*semimajor+2.0;
  origin(-L0/2., -semiminor-0.5);// Change -0.5 to -L0/16??
  // periodic(right);
  init_grid (1 << MAXLEVEL);

  // RC
  {
    char name[200];
    sprintf(name, "logstats.dat");
    fp_stats = fopen(name, "w");
  }

  // I rescale the system so that water has unit density, and the air properties are just defined in terms of water
  rho1 = 1.0; // water density
  rho2 = rho_ratio; // air density
  mu1 = 1/Re; // water dynamic viscosity
  mu2 = mu_ratio*mu1; // air dynamic viscosity
  f.sigma=1/We; // change this later

  DT = 1.0e-3; // RC
  NITERMIN = 1; // default 1
  NITERMAX = 500; // default 100
  TOLERANCE = 1e-5; // default 1e-3

  char params[200];
  sprintf(params, "params.txt");
  fp_params=fopen(params, "w");
  {
    char params[200];
	sprintf(params, "params.txt");
	fp_params=fopen(params, "w");
  }

  fprintf(fp_params, "rho1: %g \n", rho1);
  fprintf(fp_params, "rho2: %g \n", rho2);
  fprintf(fp_params, "mu1: %g \n", mu1);
  fprintf(fp_params, "mu2: %g \n", mu2);
  fprintf(fp_params, "sigma: %g \n", f.sigma);
  fprintf(fp_params, "Dimensional period: %g \n", dimensional_period);
  fprintf(fp_params, "reference Time scale: %g \n", refTime);
  fprintf(fp_params, "Dimensionless period: %g \n", period);
  fprintf(fp_params, "tplus: %g \n", tplus);
  fprintf(fp_params, "refLength: %g \n", refLength);
  fprintf(fp_params, "refVelocity: %g \n", refVelocity);
  fprintf(fp_params, "refTime: %g \n", (refLength/refVelocity));
  fprintf(fp_params, "Reynolds Number: %g \n", Re);
  fprintf(fp_params, "Weber Number: %g \n", We);
  fprintf(fp_params, "Froude Number: %g \n", Fr);
  fprintf(fp_params, "Dimensionless gravity: %g \n", 1/sq(Fr));
  fprintf(fp_params, "L0: %g \n", L0);
  fclose(fp_params);

  {
    char intName[200];
    sprintf(intName, "loginterface.dat");
    fp_interface = fopen(intName, "w");
  }
  run();

  fclose(fp_interface);
  fclose(fp_stats);
}

event acceleration (i++)
{
  thetaNOW = maxRads*sin(BB*t); // Derivation in notebook. Should write this up at some point.
  omegaNOW = BB*maxRads*cos(BB*t); // Derivative of omegaNOW w/respect to t
  face vector av = a;
  // Which of the following 2 ways of doing it is better??
  foreach_face(x) {
    gravityX = (1/sqrt(2))*sq(1/Fr)*sin(thetaNOW);
    coriolisX = 2.0*((u.y[]+u.y[-1])/2.0)*BB*maxRads*cos(BB*t);
    centripetalX = x*sq(BB)*sq(maxRads)*sq(cos(BB*t));
    av.x[] -= gravityX - coriolisX + centripetalX;
    // av.x[] -= (1/sqrt(2))*sq(1/Fr)*sin(thetaNOW) + 2.0*((u.y[]+u.y[-1])/2.0)*BB*maxRads*cos(BB*t) + x*sq(BB)*sq(maxRads)*sq(cos(BB*t));
  }
  foreach_face(y) {
    gravityY = (1/sqrt(2))*sq(1/Fr)*cos(thetaNOW);
    coriolisY = -2.0*((u.x[]+u.x[-1])/2.0)*BB*maxRads*cos(BB*t);
    centripetalY = sq(BB)*sq(maxRads)*sq(cos(BB*t))*(y+semiminor);
    av.y[] -= gravityY - coriolisY + centripetalY;
    // av.y[] -= (1/sqrt(2))*sq(1/Fr)*cos(thetaNOW) - 2.0*((u.x[]+u.x[-1])/2.0)*BB*maxRads*cos(BB*t) + sq(BB)*sq(maxRads)*sq(cos(BB*t))*(y+semiminor);
    // av.x[] -= (1/sqrt(2))*sq(1/Fr)*sin(thetaNOW);
  }
}

event init(t = 0) {
 
  // initially velocity is 0 everywhere
  foreach () {
    u.x[] = 0.;
    u.y[] = 0.;
  }

  // RC cut the top and bottom of domain
  // Could try changing this to -0.5 and 0.5 just to see
  mask (y > ymax ? top : none);
  mask (y < ymin ? bottom : none);
  
  // mask (circle < 1.0 ? top : none);

  refine(level<MAXLEVEL && (1.0/sq(semimajor))*sq(x) + (1.0/sq(semiminor))*sq(y) < 1.1);
  unrefine(level<MAXLEVEL && (1.0/sq(semimajor))*sq(x) + (1.0/sq(semiminor))*sq(y) > sq(1.1));
  // adapt_wavelet((scalar *){f, u.x, u.y}, (double[]){1e-6, 1e-2, 1e-2}, MAXLEVEL, 4);

  fraction (f, -y); // Could automate fill height here by putting a dimensional quantity at the beginning if I wanted

  // foreach() {
  //   TT[] = exp(-(10*x*x+10*y*y));
  // }
  // fraction (TT, -(0.05*sq(x) + sq(y+0.7) - sq(0.2)));
  fraction (TT, intersection((y), (-((1.0/sq(semimajor))*sq(x) + (1.0/sq(semiminor))*sq(y) - 1.0))));
  // fraction (TT, -(sq(x) + sq(y) - sq(0.5)));
  boundary({TT});

  // boundary conditions
  // u.t[top] = dirichlet(1.);
  // u.t[bottom] = dirichlet(0.);

}

event circle_flow (i++) {
  fraction (circle, ((1.0/sq(semimajor))*sq(x) + (1.0/sq(semiminor))*sq(y) - 1.0));
  foreach() {
    foreach_dimension() {
      u.x[] = (1. - circle[])*u.x[];
      u.y[] = (1. - circle[])*u.y[];
    }
  }
  boundary ((scalar *){u});
}

// event adapt (i++)
// {
//   adapt_wavelet({u, f}, (double[]){1e-2, 1e-2, 1e-2}, LEVEL);
// }
// event adapt(i++)
// {
// 	// Order matters for the following (in which the regions called by refine() and unrefine() overlap. If they don't overlap, you can at least switch refine() and unrefine(), idk about other rearrangements.)
// 	// unrefine(y<0.1);
// 	// refine(y<0.1&&level<9);
// 	// adapt_wavelet((scalar *){f, u.x, u.y}, (double[]){1e-6, 1e-2, 1e-2}, (LEVEL+2), (LEVEL-3));
//   // refine(circle[]>0.0&&level<9);
//   // Going fast to get video editing working:
//   adapt_wavelet((scalar *){f, u.x, u.y}, (double[]){1e-6, 1e-2, 1e-2}, MAXLEVEL, 4);
  
// 	// // refine(y>0.4&&x>1.0&&level<=8);
// 	// unrefine(y<0.1);
// 	// refine(y<0.1&&level<9);
// 	// adapt_wavelet((scalar *){f, u.x, u.y}, (double[]){1e-6, 1e-2, 1e-2}, (LEVEL+2), (LEVEL-3));

// 	// adapt_wavelet((scalar *){u}, (double[]){3e-2, 3e-2}, 9, 4);
// 	// unrefine(x>2.0);

// 	// Order doesn't matter for the following:
// 	// adapt_wavelet((scalar *){f, u.x, u.y}, (double[]){1e-6, 1e-2, 1e-2}, 8, 4);
// 	// refine(y<-0.45 && level<8);

// 	// Order doesn't matter for the following:
// 	// adapt_wavelet((scalar *){u}, (double[]){3e-2, 3e-2}, 9, 4);
// 	// unrefine(y<0.1);

// 	// Order doesn't matter for the following:
// 	// The following two lines manually set the refinement level to be higher near the interface (y=0.3):
// 	// refine(y<0.35 && y>0.25 && level<9);
// 	// unrefine(y<0.25 && y>0.35);

// 	// The following line is how I probably want to actually do the AMR for this code:
// 	// adapt_wavelet((scalar *){f, u.x, u.y}, (double[]){1e-6, 1e-2, 1e-2}, 8, 4);
	
// 	// The following line does just u. If you don't put (scalar *) in front of it, an error is thrown since u is a vector but Basilisk expects a scalar
// 	// adapt_wavelet((scalar *){u}, (double[]){3e-2, 3e-2}, 9, 4);

// 	// Test Radu suggested:
// 	// refine(y>0.45&&level<=8);
// 	// adapt_wavelet((scalar *){f}, (double[]){1e-6}, 8, 4);
// 	// unrefine(y>0.45&&x>3.0&&level<=8);
// }

scalar dTT[],qh[],qv[];
event integration (i++) {
  // Setting up scalar field to represent the stability condition for diffusion:
  // foreach() {
  //  diffStability[] = b*
  // }

  foreach() {
    // advection-diffusion:
      qh[] = b*(TT[-1,0]-TT[0,0])/Delta + u.x[]*(TT[0,0]+TT[-1,0])/2.0 - ((sq(u.x[])*dt)/(2.0*Delta))*(TT[0,0]-TT[-1,0]);
      qv[] = b*(TT[0,-1]-TT[0,0])/Delta + u.y[]*(TT[0,0]+TT[0,-1])/2.0 - ((sq(u.y[])*dt)/(2.0*Delta))*(TT[0,0]-TT[-1,0]);;

  // double dt = DT;
  // scalar dT[];
  // dt = dtnext (dt);

    // diffusion scheme:
    // qh[] = b*(TT[-1,0]-TT[0,0])/Delta;
    // qv[] = b*(TT[0,-1]-TT[0,0])/Delta;

    // advection scheme:
    // qh[] = u.x[]*(TT[0,0]+TT[-1,0])/2.0 - ((sq(u.x[])*dt)/(2.0*Delta))*(TT[0,0]-TT[-1,0]);
    // qv[] = u.y[]*(TT[0,0]+TT[0,-1])/2.0 - ((sq(u.y[])*dt)/(2.0*Delta))*(TT[0,0]-TT[-1,0]);
  }
  boundary ({qh});
  boundary({qv});
  foreach() {
    dTT[] = ( qh[0,0]  - qh[1,0] )/Delta + ( qv[0,0]  - qv[0,1] )/Delta;
    }
  //   // Alternative method for advection-diffusion (FINITE DIFFERENCES, forward time, centered space for diffusion, Lax-Wendroff for advection):
  // dT[] = (T[1,0]-2*T[0,0]+T[-1,0])/(sq(Delta)) + (u.x[]/(2*Delta))*(T[1,0]-T[-1,0]) + ((sq(u.x[])*dt)/(2*sq(Delta)))*(T[1,0]-2*T[0,0]+T[-1,0]) + (T[0,1]-2*T[0,0]+T[0,-1])/(sq(Delta)) + (u.y[]/(2*Delta))*(T[0,1]-T[0,-1]) + ((sq(u.y[])*dt)/(2*sq(Delta)))*(T[0,1]-2*T[0,0]+T[0,-1]);
  boundary ({qh});
  boundary({qv});
  foreach() {
    dTT[] = ( qh[0,0]  - qh[1,0] )/Delta + ( qv[0,0]  - qv[0,1] )/Delta;

  //   // Alternative method for advection-diffusion (FINITE DIFFERENCES, forward time, centered space for diffusion, Lax-Wendroff for advection):
  // dT[] = (T[1,0]-2*T[0,0]+T[-1,0])/(sq(Delta)) + (u.x[]/(2*Delta))*(T[1,0]-T[-1,0]) + ((sq(u.x[])*dt)/(2*sq(Delta)))*(T[1,0]-2*T[0,0]+T[-1,0]) + (T[0,1]-2*T[0,0]+T[0,-1])/(sq(Delta)) + (u.y[]/(2*Delta))*(T[0,1]-T[0,-1]) + ((sq(u.y[])*dt)/(2*sq(Delta)))*(T[0,1]-2*T[0,0]+T[0,-1]);
  }

  // THESE ARE THE TWO LINES THAT CAUSE IT TO INITIALIZE WEIRDLY
  foreach()
    TT[] += DT*dTT[];
  boundary ({TT});
}

event end (t = tmax) { // RC restricted to 400
  printf ("i = %d t = %g\n", i, t);
}

// RC added this for profiling
event logstats (t += tplus) {

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // i, timestep, no of cells, real time elapsed, cpu time
    fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
    fflush(fp_stats);

    if (t<1.0) { // Could change this condition in the future but won't for now

      {
      char params[200];
      sprintf(params, "params.txt");
      fp_params=fopen(params, "a");
      }

      // sizeNow = grid->n;
      // xlength = L0;
      // ylength = ymax-ymin;
      // sizeNowDouble = (double) sizeNow;
      // numCellsY = sqrt(sizeNow/xlength);
      // numCellsX = (sizeNow/numCellsY);
      // DeltaX = (xlength/numCellsX);
      // DeltaY = (ylength/numCellsY);
      DeltaX = (L0/pow(2, MAXLEVEL));

      // Here I perform a stability test for the FTCS diffusion discretization.
      // I must have Diffusion_Stability<=0.5
      // Source: http://math.tifrbng.res.in/~praveen/notes/cm2013/heat_2d.pdf
      Diffusion_Stability = (2*b*DT)/sq(DeltaX); // DeltaX = DeltaY so it doesn't matter which one I put here

      // Here I should perform a stability test for the Lax-Wendroff advection discretization.
      // I haven't actually set it up yet, but it should involve a scalar field for the entire domain.
      // I must have Advection_Stability<=1. Stability condition is the same as in 1D.
      // Source: http://pages.erau.edu/~snivelyj/ep711sp12/EP711_7.pdf
      
      // Putting grid setup parameters in the params.txt file. I'm taking the MAXIMUM # of cells in order to get the MINIMUM cell size.
      // I'm doing this because I want the scheme to be stable for every grid cell, so I test the stability for the smallest grid cell.
      // If it's stable for the smallest grid cell, it's stable for all of them.
      // fprintf(fp_params, "(Max) Total #Cells: %g \n", sizeNowDouble);
      // fprintf(fp_params, "(Max) #Cells along width: %g \n", numCellsY);
      // fprintf(fp_params, "(Max) #Cells along length: %g \n", numCellsX);
      fprintf(fp_params, "(Min) DeltaX: %g \n", DeltaX);
      // fprintf(fp_params, "(Min) DeltaY: %g \n", DeltaY);
      fprintf(fp_params, "Diffusion Stability condition (must be <=0.5, see code for source): %g \n", Diffusion_Stability);
      // fprintf(fp_params, "Advection Stability condition (must be <=1, see code for source): %g \n", Advection_Stability);

      fclose(fp_params);
    }
}

double uxmax;
double uymax;
double UMAX;
event advection_test (t+=tplus) {
  uxmax = statsf(u.x).max;
  uymax = statsf(u.y).max;
  UMAX = sqrt(sq(uxmax)+sq(uymax));
  Advection_Stability = (UMAX*DT)/DeltaX;

  FILE * fp_uxmax = fopen("ux_max", "a");
  FILE * fp_uymax = fopen("uy_max", "a");
  FILE * fp_UMAX = fopen("UMAX", "a");
  FILE * fp_adv = fopen("advection_stability", "a");

  fprintf(fp_uxmax, "%g %g \n", t, uxmax);
  fprintf(fp_uymax, "%g %g \n", t, uymax);
  fprintf(fp_UMAX, "%g %g \n", t, UMAX);
  fprintf(fp_adv, "%g %g \n", t, Advection_Stability);

  fclose(fp_uxmax);
  fclose(fp_uymax);
  fclose(fp_UMAX);
  fclose(fp_adv);
}

event gfsview (t += tplus*10.0) { // RC
    char name_gfs[200];
    sprintf(name_gfs,"Slice-%0.1f.gfs",t);

    FILE* fp_gfs = fopen (name_gfs, "w");
    output_gfs(fp_gfs);
    fclose(fp_gfs);
}

event xmovie (t+=tplus)
{
 view (fov=9, width=800, height=350);
 clear();
 // cells(lc={0.5,0.5,0.5}, lw=0.5);
 squares("u.x", min=-0.2, max=0.2, linear=true, map=cool_warm);
 draw_vof ("f", lc = {0.0,0.0,0.0}, lw=2);
 draw_vof("circle", lc = {0.0,0.0,0.0}, lw=2);
 // cells();
 save ("xmovie.mp4");
}

event ymovie (t+=tplus)
{
 view (fov=9, width=800, height=350);
 clear();
 squares("u.y", min=-0.5, max=0.5, linear=true, map=cool_warm);
 draw_vof ("f", lc = {0.0,0.0,0.0}, lw=2);
 draw_vof("circle", lc = {0.0,0.0,0.0}, lw=2);
 // draw_vof("circle", lc = {0.0,0.0,0.0}, lw=2);
 // cells(); // Movie is black when this line is included
 save ("ymovie.mp4");
}

event tracemovie (t+=tplus)
{
 view (fov=9, width=800, height=350);
 clear();
 // cells(lc={0.5,0.5,0.5}, lw=0.5);
 squares("TT", spread=-1, linear=true, map=cool_warm);
 draw_vof ("f", lc = {0.0,0.0,0.0}, lw=2);
 draw_vof("circle", lc = {0.0,0.0,0.0}, lw=2);
 // cells();
 save ("tracemovie.mp4");
}

event shearmovie (t+=tplus)
{
  scalar shear[];
  foreach()
  	shear[] = (mu(f[]))*(u.x[0, 1]-u.x[0, -1])/(2.*Delta);
  // foreach()
  //   if (y<0.3)
  //   {
  //     shear[] = (mu1)*(u.x[0, 1]-u.x[0, -1])/(2.*Delta);
  //   }
  // foreach()
  //   if (y>=0.3)
  //   {
  //     shear[] = (mu2)*(u.x[0, 1]-u.x[0, -1])/(2.*Delta);
  //   }
  // foreach()
  // shear[] = (u.x[0, 1]-u.x[0, -1])/(2.*Delta);
  boundary ({shear});
  view (fov=9, width=800, height=350);
  clear();
  squares("shear", spread=1, linear=true, map=cool_warm);
  draw_vof("f", lc = {0.0,0.0,0.0}, lw=2);
  draw_vof("circle", lc = {0.0,0.0,0.0}, lw=2);
  save("shearMovie.mp4");
  FILE * fp_shear = fopen("shear", "w");
  output_field ({shear}, fp_shear, linear=true);
  fclose(fp_shear);
}

event interface (t+=tplus*10.0)
{
  char name_interface[100];
  sprintf(name_interface, "interface_%g.dat", t);

  FILE * fp2 = fopen(name_interface, "w");
  output_facets (f, fp2);
}

event loginterface (t += tplus) {    
    scalar posX[],posY[];
    position (f, posX, {1, 0});
    position (f, posY, {0,1});
    fprintf(fp_interface, "%i %g %1.4f %1.4f %1.4f %1.4f %1.4f\n", i, t, statsf(f).sum, statsf(posX).min, statsf(posX).max, statsf(posY).min, statsf(posY).max);
    fflush(fp_interface);
}

event profiles (t = 0; t+=tplus; t<=tmax) // RC restricted the output a little, don't overdo it at first!
{
  FILE * fp = fopen("xprof", "a");
  for (double y = ymin; y <= ymax; y += 0.01)
    fprintf (fp, "%g %g %g\n", t, y, interpolate (u.x, 0, y));
  fclose (fp);
  
  fp = fopen("yprof", "a");
  for (double x = xmin; x <= xmax; x += 0.01)
    fprintf (fp, "%g %g %g\n", t, x, interpolate (u.y, x, 0));
  fclose (fp);
  
  // scalar shear[];
  // foreach()
  //   if (y<0.3)
  //   {
  //     shear[] = mu1.*(u.x[0, 1]-u.x[0, -1])/(2.*Delta);
  //   }
  // foreach()
  //   if (y>=0.3)
  //   {
  //     shear[] = mu2.*(u.x[0, 1]-u.x[0, -1])/(2.*Delta);
  //   }
  // // fp=fopen("shearprof", "a");
  // // for (double y = -0.5; y <= 0.5; y += 0.01)
  // //   for (double x = -4; x <= 4; x += 0.01)
  // //     fprintf(fp, "%g ", shear);
  // //   fprintf(fp, "\n");
  // fp=fopen("Xshearprof", "a");
  // for (double y = -0.5; y <= 0.5; y += 0.01)
  //   fprintf (fp, "%g %g %g\n", t, y, shear);
  // fclose (fp);

  // fp=fopen("Yshearprof", "a");
  // for (double x = -4; x <= 4; x += 0.01)
  //   fprintf (fp, "%g %g %g\n", t, x, shear);
  // fclose (fp);
}

// ybelow and yabove will be used to mark points on the ellipse for each x value
double ybelow;
double yabove;
event stations (t = 0; t+=tplus; t<=tmax) // RC restricted the output a little, don't overdo it at first!
{
  FILE * fp = fopen("stations", "a");
  for (double x = -semimajor; x <= semimajor+0.1; x+=0.1) {
    // I do the rounding stuff below because I need a consistent grid for Python postprocessing
    yabove = sqrt(1-sq(x/3));
    yabove = (double)floor(yabove*10.0)/10.0;
    ybelow = -yabove;
    for (double y = ybelow; y <= yabove; y+=0.1) {
      fprintf (fp, "%g %g %g %g\n", t, x, y, interpolate(TT, x, y));
    }
  }
  fclose(fp);
}

event vortmovie (t+=tplus, t<=10.0*tplus)
{
  scalar omega[];
  vorticity(u, omega);

  view (fov=9, width=800, height=350);
  clear();
  squares("omega", spread=-1.0, linear=true, map=cool_warm);
  draw_vof ("f", lc = {0.0,0.0,0.0}, lw=2);
  draw_vof("circle", lc = {0.0,0.0,0.0}, lw=2);
  // draw_vof("fs",fc = {0.0,0.0,0.0}, lw=1); 
  // draw_vof("f", lc = {0.0,0.0,0.0}, lw=1); // For some reason Basilisk throws an error if you don't put spaces between 'lc='
  save("Vorticity.mp4");
}
// See "Bioreactor_Equations_2.pdf" for derivations of the equations and nondimensionalizations used in this code.

// The following 2 lines are helpful for interface stuff:
#define FILTERED
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))

// Header files to be included (you can read more about them at the Basilisk C website)
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
#include "fractions.h"
#include "tracer.h"
#include "diffusion.h"
#include "tag.h"

scalar T[]; // This will be the Oxygen tracer
scalar * tracers = {T};

// Defining scalar fields to keep track of how much of the oxygen tracer is in the water, air, and in and out of the ellipse:
scalar T_water[]; // In the water
scalar T_air[]; // In the air
scalar T_ellipseIN[]; // In the ellipse
scalar T_ellipseOUT[]; // Outside of the ellipse

scalar WATER[]; // Scalar to keep track of the mass of the water

scalar UV[]; // Scalar field for the magnitude of velocity
// double b = 0.00002; // Diffusion coefficient

double thetaNOW; // Angle of rotation
double omegaNOW; // Angular velocity of the ellipse

// The following doubles will be used in the acceleration event to avoid really long lines of code:
double gravityX; // x-component of gravitational acceleration
double gravityY; // y-component of gravitational acceleration
double coriolisX; // x-component of coriolis acceleration
double coriolisY; // y-component of coriolis acceleration
double centripetalX; // x-component of centripetal acceleration
double centripetalY; // y-component of centripetal acceleration

FILE *fp1 ; // pointer file to be used later

#define MAXLEVEL 9 // Maximum level of refinement. I will define grid refinement in the init function

// Defining Froude number separately from the other nondimensional quantities because it will be used to set up the reference Velocity:
#define Fr 1.0 // Froude number

// DIMENSIONAL QUANTITIES:
#define rhoWater 1000.0 // water density, kg/m^3
#define rhoAir 1.225 // air density, kg/m^3
#define muWater 0.001 // water viscosity, Pa*s
#define muAir 1.81e-5 // air viscosity, Pa*s

// If you want to make the setup computationally gentler, you make the air viscosity one or two magnitudes greater:
// #define muAir 1.81e-4 // visc. of air *10
// #define muAir 1.81e-3 // visc. of air *100

#define sig 0.0728  //surface tension of water, N/m
const double semiminor = 1.0; // semiminor axis (dm)
const double semimajor = semiminor*3.0; // semimajor axis (dm)
const double maxDegrees = 7.0; // degrees through which the reactor rotates
double maxRads = maxDegrees*(3.14159265/180.0); // radians through which the reactor rotates
#define dimensional_period 1.0 // rocking period in actual seconds (will convert to nondimensional time units below for the simulation)
#define refLength 0.1  // semiminor axis length (m)
#define refVelocity (Fr*sqrt(9.8*refLength))  // Reference length, defined in terms of the Froude number
#define refTime (refLength/refVelocity) // Reference time scale (not actually used for the simulation, but useful to have)

// DIMENSIONLESS QUANTITIES:
const double period = (dimensional_period/refTime); // how many NONDIMENSIONAL TIME UNITS it takes to go through a complete rocking cycle
double BB = (2.0*3.14159265)/period; // constant used in rocking motion equations below
#define Re (rhoWater*refVelocity*refLength/muWater)  // Reynolds number
#define We (rhoWater*pow(refVelocity,2)*refLength/sig) // Weber number
#define rho_ratio (rhoAir/rhoWater) // Density ratio
#define mu_ratio (muAir/muWater) // Viscosity ratio 
// Should include Peclet numbers eventually

const double tmax = (50.0*period); // You probably want to change this to a 5, not 50 (50 oscillations is a LONG simulation)
const double tplus = (period/100.0); // events that run at set t invervals will have those intervals defined by this quantity

// Pointer files for outputting data:
FILE *fp_params;
FILE * fp_stats;
FILE * fp_interface;
FILE * fp_mass;
FILE * fp_mass_water;
FILE * fp_mass_air;
FILE * fp_mass_IN;
FILE * fp_mass_OUT;
FILE * fp_howmuch_water;

// Boundary conditions:
u.t[top] = dirichlet(0.0);
u.n[top] = dirichlet(0.0);

u.t[bottom] = dirichlet(0.0);
u.n[bottom] = dirichlet(0.0);

u.t[left] = dirichlet(0.0);
u.n[left] = dirichlet(0.0);

u.t[right] = dirichlet(0.0);
u.n[right] = dirichlet(0.0);

const double ymax = semiminor+0.5; // ymax of domain (used in masking & profiling, etc.)
const double ymin = -semiminor-0.5; // ymin of domain (used in masking & profiling, etc.)
const double xmax = semimajor+3.0; // xmax of domain (used in profiling)
const double xmin = -semimajor-3.0; // xmin of domain (used in profiling)
// The following is from an old attempt to do some stability tests. It's probably not relevant anymore but I'm leaving it here in case anyone wants to build on it.
double DeltaX; // Length of smallest cell in x direction
double Diffusion_Stability; // test for stability of the FTCS diffusion scheme
double Advection_Stability; // test for stability of the Lax-Wendroff advection scheme

// Defining a scalar field for keeping track of inside/outside the ellipse:
scalar circle[];

int main() {
  L0 = 2.0*semimajor+2.0; // Length of the simulation domain
  origin(-L0/2., -semiminor-0.5); // Defining the origin of the setup
  init_grid (1 << MAXLEVEL); // Setting up the grid

  // Setting up some data files:
  {
    char name[200];
    sprintf(name, "logstats.dat");
    fp_stats = fopen(name, "w");
  }

  {
    char name[200];
    sprintf(name, "mass.dat");
    fp_mass = fopen(name, "w");
  }
  {
    char name[200];
    sprintf(name, "mass_water.dat");
    fp_mass_water = fopen(name, "w");
  }
  {
    char name[200];
    sprintf(name, "mass_air.dat");
    fp_mass_air = fopen(name, "w");
  }
  {
    char name[200];
    sprintf(name, "mass_IN.dat");
    fp_mass_IN = fopen(name, "w");
  }
  {
    char name[200];
    sprintf(name, "mass_OUT.dat");
    fp_mass_OUT = fopen(name, "w");
  }
  {
    char name[200];
    sprintf(name, "howmuch_water.dat");
    fp_howmuch_water = fopen(name, "w");
  }

  // I rescale the system so that water has unit density, and the air properties are just defined in terms of water
  rho1 = 1.0; // water density
  rho2 = rho_ratio; // air density
  mu1 = 1/Re; // water dynamic viscosity
  mu2 = mu_ratio*mu1; // air dynamic viscosity
  f.sigma=1/We; // change this later

  // Some computational parameters:
  DT = 1.0e-3;
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
  fclose(fp_mass);
  fclose(fp_mass_water);
  fclose(fp_mass_air);
  fclose(fp_mass_IN);
  fclose(fp_mass_OUT);
  // fclose(fp_howmuch_water);
}

// Event to define gravitational, coriolis and centripetal accelerations:
// See "Rotating_Reference_Frame.pdf" for derivations of the accelerations written out here
event acceleration (i++)
{
  thetaNOW = maxRads*sin(BB*t); // Angle defining how much the ellipse has rotated relative to horizontal
  omegaNOW = BB*maxRads*cos(BB*t); // Derivative of omegaNOW w/respect to t
  face vector av = a;
  foreach_face(x) {
    gravityX = (1/sqrt(2))*sq(1/Fr)*sin(thetaNOW); // x component of gravity
    coriolisX = 2.0*((u.y[]+u.y[-1])/2.0)*BB*maxRads*cos(BB*t); // x component of Coriolis
    centripetalX = x*sq(BB)*sq(maxRads)*sq(cos(BB*t)); // x component of Centripetal
    av.x[] -= gravityX - coriolisX + centripetalX; // Combining the accelerations
  }
  foreach_face(y) {
    gravityY = (1/sqrt(2))*sq(1/Fr)*cos(thetaNOW); // y component of gravity
    coriolisY = -2.0*((u.x[]+u.x[-1])/2.0)*BB*maxRads*cos(BB*t); // y component of Coriolis
    centripetalY = sq(BB)*sq(maxRads)*sq(cos(BB*t))*(y+semiminor); // y component of Centripetal
    av.y[] -= gravityY - coriolisY + centripetalY; // Combining the accelerations
  }
}

event init(t = 0) {
 
  // initially velocity is 0 everywhere
  foreach () {
    u.x[] = 0.;
    u.y[] = 0.;
  }

  // Could try changing this to -0.5 and 0.5 just to see
  mask (y > ymax ? top : none);
  mask (y < ymin ? bottom : none);

  refine(level<MAXLEVEL && (1.0/sq(semimajor))*sq(x) + (1.0/sq(semiminor))*sq(y) < 1.0); // Refining to the maximum level inside the ellipse
  unrefine(level<MAXLEVEL && (1.0/sq(semimajor))*sq(x) + (1.0/sq(semiminor))*sq(y) > sq(1.0)); // Unrefining outside of the ellipse

  fraction (f, -y); // Could automate fill height here by putting a dimensional quantity at the beginning if I wanted

  // Different initializations for the tracer position:
  fraction (T, -(0.05*sq(x) + sq(y+0.7) - sq(0.2))); // Blob in the water
  // fraction (T, intersection((y), (-((1.0/sq(semimajor))*sq(x) + (1.0/sq(semiminor))*sq(y) - 1.0)))); // Headspace (air region) filled with tracer
  // fraction (T, -(sq(x) + sq(y) - sq(0.5))); // Circle of tracer at interface
}

// Event for diffusion of the Oxygen tracer:
// example: http://basilisk.fr/src/test/kh-ns.c
event tracer_diffusion(i++) { // not sure if it should be i++ or +=dt here
  diffusion (T, dt, mu);
}

// The following removes small droplets. Without using this, the simulation blows up (fun to watch)
event small_droplet_removal (i++) {
/* Removes any small droplets that have formed, that are smaller than a specific
    size */
    remove_droplets(f, 8);       // Removes droplets of diameter 8 cells or less
    remove_droplets(f, 8, true); // Removes bubbles of diameter 8 cells or less
}

// This event sets up the ellipse & boundary conditions
event circle_flow (i++) {
  fraction (circle, ((1.0/sq(semimajor))*sq(x) + (1.0/sq(semiminor))*sq(y) - 1.0)); // Defining the elliptical reactor outline
  foreach() {
    foreach_dimension() { // Boundary conditions (velocities are 0 outside the ellipse):
      u.x[] = (1. - circle[])*u.x[];
      u.y[] = (1. - circle[])*u.y[];
    }
  }
  boundary ((scalar *){u});
}

// Defining values in the scalar field for the magnitude of water velocity:
event vel_Mag (t+=tplus) {
  foreach() {
    if (f[]==1) { // This is just for water. If you change the conditionas you can get fluid velocity for air, or both air & water.
      UV[]=sqrt(sq(u.x[])+sq(u.y[]));
    }
    else {
      UV[]=0;
    }
  }
}

// Defining values for the scalar fields representing the Oxygen tracer in various domains of the simuation:
event Tracer_Fields (t+=tplus) {
  // Tracer field in water inside ellipse:
  foreach() {
    if (f[]==1 && circle[]==0) {
      T_water[]=T[];
    }
    else {
      T_water[]=0;
    }
  }
  // Tracer field in air inside ellipse:
  foreach() {
    if (f[]==0 && circle[]==0) {
      T_air[]=T[];
    }
    else {
      T_air[]=0;
    }
  }
  // Tracer field inside the ellipse:
  foreach() {
    if (circle[]==0) {
      T_ellipseIN[]=T[];
    }
    else {
      T_ellipseIN[]=0;
    }
  }
  // Tracer field outside the ellipse:
  foreach() {
    if (circle[]==1) {
      T_ellipseOUT[]=T[];
    }
    else {
      T_ellipseOUT[]=0;
    }
  }
}

// Measuring the mass of the water:
event howmuch_water_isthere (t+=tplus) {
  foreach() {
    if (f[]==1 && circle[]==0) {
      WATER[]=f[];
    }
  }
}

// Simulation ends at t=tmax, prints out a few numbers.
event end (t = tmax) {
  printf ("i = %d t = %g\n", i, t);
}

// Printing some profiles:
event logstats (t += tplus) {

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
    fflush(fp_stats);
}

// Measuring the tracer mass in air, water, in the ellipse, outside of the ellipse (outputing to data files):
event tracer_mass (t+=tplus) {
  stats s1 = statsf(T);
  fprintf(fp_mass, "t: %g Total mass: %g \n", t, s1.sum);
  fflush(fp_mass);

  stats s2 = statsf(T_water);
  fprintf(fp_mass_water, "t: %g Total mass: %g \n", t, s2.sum);
  fflush(fp_mass_water);

  stats s3 = statsf(T_air);
  fprintf(fp_mass_air, "t: %g Total mass: %g \n", t, s3.sum);
  fflush(fp_mass_air);

  stats s4 = statsf(T_ellipseIN);
  fprintf(fp_mass_IN, "t: %g Total mass: %g \n", t, s4.sum);
  fflush(fp_mass_IN);

  stats s5 = statsf(T_ellipseOUT);
  fprintf(fp_mass_OUT, "t: %g Total mass: %g \n", t, s5.sum);
  fflush(fp_mass_OUT);
}

// Measuring the mass of the water:
event measure_water (t+=tplus) {
  stats swater = statsf(WATER);
  fprintf(fp_howmuch_water, "t: %g Water mass: %g \n", t, swater.sum);
  fflush(fp_howmuch_water);
}

// Advection stability test (this is pretty old, I'm not sure if it still works, but feel free to play with it)
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

// Outputing gfsview files. Gfsview is a useful functionality for looking at simulation data.
// You can look at any of the scalar fields, grid cells, etc. Just type "gfsview2D FILENAME" in the terminal
event gfsview (t += tplus*10.0) { // RC
    char name_gfs[200];
    sprintf(name_gfs,"Slice-%0.1f.gfs",t);

    FILE* fp_gfs = fopen (name_gfs, "w");
    output_gfs(fp_gfs);
    fclose(fp_gfs);
}

// Making movies:
event xmovie (t+=tplus)
{
 view (fov=9, width=800, height=350); // Setting up field of view & size
 clear();
 squares("u.x", spread=-1, linear=true, map=cool_warm); // Graphing x component of velocity
 draw_vof ("f", lc = {0.0,0.0,0.0}, lw=2); // Drawing a curve to represent the air-water interface
 draw_vof("circle", lc = {0.0,0.0,0.0}, lw=2); // Curve to represent boundary of the ellipse
 save ("xmovie.mp4");
}

event ymovie (t+=tplus)
{
 view (fov=9, width=800, height=350);
 clear();
 squares("u.y", spread=-1, linear=true, map=cool_warm);
 draw_vof ("f", lc = {0.0,0.0,0.0}, lw=2);
 draw_vof("circle", lc = {0.0,0.0,0.0}, lw=2);
 save ("ymovie.mp4");
}

event MagMovie (t+=tplus)
{
 view (fov=9, width=800, height=350);
 clear();
 squares("UV", min=0.0, max=0.5, linear=true, map=cool_warm);
 draw_vof ("f", lc = {0.0,0.0,0.0}, lw=2);
 draw_vof("circle", lc = {0.0,0.0,0.0}, lw=2);
 // draw_vof("circle", lc = {0.0,0.0,0.0}, lw=2);
 // cells(); // Movie is black when this line is included
 save ("VelMagmovie.mp4");
}

// Tracer movies
// This movie sets the color scale based on minimum & maximum values of the scalar field.
event TmovieMinMax (t+=tplus)
{
 view (fov=9, width=800, height=350);
 clear();
 squares("T", min=0.0, max=0.5, linear=true, map=cool_warm);
 draw_vof ("f", lc = {0.0,0.0,0.0}, lw=2);
 draw_vof("circle", lc = {0.0,0.0,0.0}, lw=2);
 save ("TmovieMinMax.mp4");
}

// This movie sets the color scale based on the variance of the scalar field.
event TmovieSpread (t+=tplus)
{
 view (fov=9, width=800, height=350);
 clear();
 squares("T", spread=-1, linear=true, map=cool_warm);
 draw_vof ("f", lc = {0.0,0.0,0.0}, lw=2);
 draw_vof("circle", lc = {0.0,0.0,0.0}, lw=2);
 save ("TmovieSpread.mp4");
}

// Movie to track shear stress:
// event shearmovie (t+=tplus)
// {
//   scalar shear[];
//   foreach()
//   	shear[] = (mu(f[]))*(u.x[0, 1]-u.x[0, -1])/(2.*Delta);
//   // foreach()
//   //   if (y<0.3)
//   //   {
//   //     shear[] = (mu1)*(u.x[0, 1]-u.x[0, -1])/(2.*Delta);
//   //   }
//   // foreach()
//   //   if (y>=0.3)
//   //   {
//   //     shear[] = (mu2)*(u.x[0, 1]-u.x[0, -1])/(2.*Delta);
//   //   }
//   // foreach()
//   // shear[] = (u.x[0, 1]-u.x[0, -1])/(2.*Delta);
//   boundary ({shear});
//   view (fov=9, width=800, height=350);
//   clear();
//   squares("shear", spread=1, linear=true, map=cool_warm);
//   draw_vof("f", lc = {0.0,0.0,0.0}, lw=2);
//   draw_vof("circle", lc = {0.0,0.0,0.0}, lw=2);
//   save("shearMovie.mp4");
//   FILE * fp_shear = fopen("shear", "w");
//   output_field ({shear}, fp_shear, linear=true);
//   fclose(fp_shear);
// }

// Movie of vorticity:
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

// Outputting data files with air-water interface location that can be analyzed in Python, Matlab, etc.
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

// Profiling the x and y velocities along particular lines in the ellipes:
event profiles (t = 0; t+=tplus; t<=tmax)
{
  FILE * fp = fopen("xprof", "a");
  for (double y = ymin; y <= ymax; y += 0.01)
    fprintf (fp, "%g %g %g\n", t, y, interpolate (u.x, 0, y));
  fclose (fp);
  
  fp = fopen("yprof", "a");
  for (double x = xmin; x <= xmax; x += 0.01)
    fprintf (fp, "%g %g %g\n", t, x, interpolate (u.y, x, 0));
  fclose (fp);
}

double ybelow;
double yabove;
event stations (t = 0; t+=tplus; t<=tmax)
{
  FILE * fp = fopen("stations", "a");
  for (double x = -semimajor; x <= semimajor+0.1; x+=0.1) {
    yabove = sqrt(1-sq(x/3));
    yabove = (double)floor(yabove*10.0)/10.0;
    ybelow = -yabove;
  }
  fclose(fp);
}
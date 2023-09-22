//  ================================================================================
//  ||                                                                            ||
//  ||              esPIC                                                         ||
//  ||              ------------------------------------------------------        ||
//  ||              E L E C T R O S T I C   P A R T I C L E - I N - C E L L       ||
//  ||                                                                            ||
//  ||              D E M O N S T R A T I O N   C O D E                           ||
//  ||              ------------------------------------------------------        ||
//  ||                                                                            ||
//  ||       Developed by: Scott R. Runnels, Ph.D.                                ||
//  ||                     University of Colorado Boulder                         ||
//  ||                                                                            ||
//  ||                For: CU Boulder CSCI 4576/5576 and associated labs          ||
//  ||                                                                            ||
//  ||           Copyright 2020 Scott Runnels                                     ||
//  ||                                                                            ||
//  ||                     Not for distribution or use outside of the             ||
//  ||                     this course.                                           ||
//  ||                                                                            ||
//  ================================================================================


class particles
{

 public:
  int n;
  VD x;
  VD y;
  VD vx;
  VD vy;
  VI active;
  VD xf, yf;
  double fluxTimer;
  double Qp;
  double previousInjectionTime;

  // ==
  // ||
  // ||  Constructor
  // ||
  // ==

  particles(int _n)
    {
      n  = _n;
      x.resize(n+1);
      y.resize(n+1);
      vx.resize(n+1);
      vy.resize(n+1);
      xf.resize(n+1);
      yf.resize(n+1);
      active.resize(n+1);

      for ( int i = 1 ; i < active.size() ; ++i ) active[i] = 0;
      fluxTimer = 0.;

      Qp = 1.; // Extensive physical quantity represented by each particle
      previousInjectionTime = 0.;
    }

  // ==
  // ||
  // ||  Adds new particles to the list of particles
  // ||
  // ==

  void add(VD &xAdd , VD &yAdd, VD &vxAdd , VD &vyAdd)
  {

    // Add new particles in spaces left inactive (active[j] != 1)
    // to the existing arrays.

    int count = 0;
    for ( int j = 1 ; j <= x.size()-1 ; ++j )
      {
	if ( active[j] == 0 )
	  {
	    ++count;
	    
	    if ( count > xAdd.size()-1 )
	      {
		return;
	      }
	    
	    x [j] = xAdd [count];
	    y [j] = yAdd [count];
	    vx[j] = vxAdd[count];
	    vy[j] = vyAdd[count];
	    active[j] = 1;

	  }
      }
    
    // If we are out of space in the existing arrays, resize
    // them and try again.

    FatalError("More room needed for particles. Adjust your -nPtcl value on the command line.");

  }

  
  // ==
  // ||
  // ||  Moves particles the (provided) dt
  // ||
  // ==

  void move(double dt)
  {
    for ( int i = 1 ; i < x.size() ; ++i )
      {
      if ( active[i] == 1 )
	{

	  // Apply force to update velocity with acceleration

	  vx[i] += dt * xf[i];
	  vy[i] += dt * yf[i];

	  // Apply velocity to update position
	  
	  x [i] += dt * vx[i];
	  y [i] += dt * vy[i];
	}
      }
  }

  // ==
  // ||
  // ||  ------------------------
  // ||  Function addFlux
  // ||  ------------------------
  // ||
  // ||  Regarding particles:
  // ||
  // ||  (o) Each particle represents an extensive physical quantity, e.g., mass in [kg], or
  // ||      charge in [Coulombs].  
  // ||  (o) We will call that physical quantity, Q, and its units will be [Q], i.e., the
  // ||      extensive quantity.  On a per-particle basis, it is Qp = [Q/particle]
  // ||  (o) Np is the number of particles.  So the density of particles in a volume, V, is
  // ||      rho [EQ/m^3] = Np * Qp / V
  // ||  (o) The flux of our extensive quantity [EQ] is the velocity of the particles,
  // ||      times their density, or
  // ||
  // ||             flux, f [EQ/sec*m^2] = v [m/sec] * rho [Np*Qp/m^3 ]
  // ||
  // ||  (o) For a given area, A [m^2], the number of particles crossing the area is the
  // ||      particle flux, or p in units of [particles/sec].
  // ||
  // ||      p[particles/sec] = f [Q/sec*m^2] * A [m^2] = v [m/sec] * rho [Np*Qp/m^3] * A [m^2]
  // ||
  // ||      As part of the setup, we specify v and rho, so that for an area provided to this
  // ||      routine, we can calculate p:
  // ||
  // ||                    p [particles/sec] = v * rho * A
  // ||
  // ||  (o) It is helpful to consider the inverse, called the particle interval, or I, in
  // ||      units of [seconds/particle]:
  // ||
  // ||                    I [sec/particle] = 1/p = 1/(v rho A)
  // ||
  // ||      The value I gives us the time between each set of particles crossing area, A.
  // ||
  // ||  (o) In the 2D geometry considered here, A = dx * dy = (1) * dy = y1 - y0.
  // ||
  // ||
  // ==

  void addFlux(double time, double y0 , double y1 , double density, double vx_bdy, double vy_bdy)
  {

    double Area     = y1 - y0;                      // Flux area
    double Interval = 1. / (vx_bdy*density*Area);   // Particle interval

    // See if it is time to inject new particles.  Assume for the moment that it is
    // time, in fact, it is past time to inject new particles (time steps might not be
    // just right):
    
    double howLateWeAre = (time - previousInjectionTime) - Interval;

    // If it's not time yet, then were early (i.e., negative late-ness):      

    if ( howLateWeAre < 0. ) return;

    // So we are late!  Calculate how far into the mesh the particles have traveled
    // since entering.

    double distanceIn = howLateWeAre * vx_bdy;

    // Calculate how many particles have crossed the cross section

    double Np = density*Area/Qp ;

    // Calculate the vertical spacing of the particles

    double dy = (y1-y0) / Np;  
                               
    // Allocate memory to store the new particles
    
    VD xAdd ;  xAdd.resize(Np+1);    VD yAdd ;  yAdd.resize(Np+1);
    VD vxAdd; vxAdd.resize(Np+1);    VD vyAdd; vyAdd.resize(Np+1);

    // Compute their initial position

    for ( int i = 1; i <=  Np ; ++i )
      {
	xAdd[i] = 0.;
	yAdd[i] = y0 + dy/2. + dy*(i-1);
      }

    // Store their initial velocity

    for ( int i = 1 ; i <= Np ; ++i )
      {
	vxAdd[i] = vx_bdy;
	vyAdd[i] = vy_bdy;
      }
    
    // Record the time that this set of particles entered the mesh.  Note that
    // we are recording the actual time, which is the current time minus how
    // late we are.

    previousInjectionTime = time - howLateWeAre;

    // Add the particles

    add(xAdd , yAdd, vxAdd , vyAdd);
    
  }

  


#include "particle_plotter.h"


};



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

//  ==
//  ||
//  ||   M E S H  +  F I E L D   V A R I A B L E   N U M B E R I N G   S Y S T E M
//  ||
//  ==
//
//  The mesh for a 3x3 mesh is shown below. 
//
//  Some key variables
//  ------------------------------------
//
//  nRealx = number of **physical** nodes in the x-direction
//  nRealy =    "   "       "         "    "  "  y-    "
//
//  In the figure below,
//
//  nRealx = 3   nRealy = 3
//
//  Field variable/natural node numbering
//  -------------------------------------
//
//  The ID numbers of the field variables, i.e., variable we're solving for, is shown
//  in () and [] in the figure.  Note:
//
//  [...] are Physical nodes --> those that lie in the physical doman for this processor
//  (...) are Ghost    nodes --> those that lie outside the physical doman for this processor
//
//  i-j node numbering
//  -------------------------------------
//
//  In addition to natural node numbering, used in the linear system solver and in the
//  field variable, nodes are numbered using an i-j logical numbering scheme, also as
//  shown in the figure below.  
//
//  Physical nodes --> range from 1,1 to nRealx,nRealy
//  Ghost    nodes --> are the candy-coating, with lower left at 0,0 and upper right at nRealx+1,nRealy+1
//
//    (21)----(22)----(23)----(24)-----(25)  <--- nRealx+1 , nRealy+1
//     |       |       |       |        |   
//     |       |       |       |        |  
//     |       |       |       |        |  
//    (16)----[17]----[18]----[19]-----(20)
//     |    1,3|    2,3|    3,3|        |  
//     |       |       |       |        |  
//     |       |       |       |        |  
//    (11)----[12]----[13]----[14]-----(15)
//     |    1,2|    2,2|    3,2|        |  
//     |       |       |       |        |  
//     |       |       |       |        |  
//    (6)-----[7]-----[8]-----[9]------(10) 
//     |    1,1|    2,1|    3,1|        |  
//     |       |       |       |        |  
//     |       |       |       |        |  
//    (1)-----(2)-----(3)-----(4)------(5)
//  0,0                              4,0
//


//  ==
//  ||
//  ||      C L A S S:   L A P L A C I A N O N G R I D 
//  ||
//  == 

class LaplacianOnGrid
{

public:

  double x0, x1, y0, y1;
  VD x,y;
  int nRealx    , nRealy   , nField;
  double dx, dy;
  VDD Acoef;
  VII Jcoef ;
  VD  phi ;  VD  b ;
  int bandwidth;
  int myPE;
  VD Qval;

  //  ==
  //  ||
  //  ||  Constructor:  Initialize values
  //  ||
  //  ==

  LaplacianOnGrid(double _x0   , double _x1,  double _y0, double _y1 , int ncell_x , int ncell_y, mpiInfo &myMPI )
  {

    x0 = _x0;    x1 = _x1;
    y0 = _y0;    y1 = _y1;
    
    nRealx     = ncell_x + 1; 
    nRealy     = ncell_y + 1; 
    nField  = (nRealx+2)*(nRealy+2);  
    dx     = (x1-x0)/ncell_x;
    dy     = (y1-y0)/ncell_y;

    // Allocate memory -- Note that the node numbers and field variable numbers must
    // match.  So even though we only will be caring about real nodes, their node
    // numbers are naturally ordered, so must be of size nField.

    x.resize(nField+1); y.resize(nField+1); Qval.resize(nField+1);
		       
    for ( int i = 1 ; i <= nRealx ; ++i )
      for ( int j = 1 ; j <= nRealy ; ++j )
	{
	  int p = pid(i,j);
	  x[p] = x0 + (i-1)*dx;
	  y[p] = y0 + (j-1)*dy;
	}

    bandwidth = 5;
    
    Acoef.resize(nField+1 ); rLOOP Acoef[r].resize(bandwidth+1);
    Jcoef.resize(nField+1 ); rLOOP Jcoef[r].resize(bandwidth+1);
    b.resize(nField+1 );
    phi.resize(nField+1);

    rLOOP phi[r] = 0.;
    
    myPE = myMPI.myPE;

  }

  //  ==
  //  ||
  //  ||  Form Linear System Ax = b
  //  ||
  //  ==

  void FormLS(mpiInfo &myMPI)
  {
    
    rLOOP cLOOP Acoef[r][c] = 0.;  // Initialize linear system
    rLOOP cLOOP Jcoef[r][c] = 0.;  //
    rLOOP       Jcoef[r][1] = r;   //
    rLOOP b[r] = 0.;
    
    double dx2 = dx*dx;           // Form matrix entries for the interior grid points
    double dy2 = dy*dy;           // Form matrix entries for the interior grid points
    
    iLOOP      // A*phi = b has extra rows and columns to accommodate nodes outside the physical domain.
      jLOOP    // That is why we can loop all the from 1 to nRealx and 1 to nRealy.
      {
	int p = pid(i,j);
	Acoef[ p ][ 1 ] = -2./dx2 - 2./dy2;
	Acoef[ p ][ 2 ] =  1./dx2;
	Acoef[ p ][ 3 ] =  1./dx2;
	Acoef[ p ][ 4 ] =  1./dy2;
	Acoef[ p ][ 5 ] =  1./dy2;

	Jcoef[ p ][ 1 ] =  pid( i  , j  );
	Jcoef[ p ][ 2 ] =  pid( i+1, j  );
	Jcoef[ p ][ 3 ] =  pid( i-1, j  );
	Jcoef[ p ][ 4 ] =  pid( i  , j+1);
	Jcoef[ p ][ 5 ] =  pid( i  , j-1);
	
      }

    // Populate RHS

    double c1 = 10.;
    double c2 = 10.;

    // Q IS INTERPRETED AS NEGATIVE CHARGE HERE
    // rLOOP b[r] = -( -c1*Qval[r] + c2*phi[r] );

    // Q IS INTERPRETED AS POSITIVE CHARGE HERE
    rLOOP b[r] = -( c1*Qval[r] + c2*phi[r] );

    // Populate arrays, one for each side, containing the Dirichlet BCs

    VD phiL, phiR, phiT, phiB; phiL.resize(nRealy+2);  phiR.resize(nRealy+2); phiB.resize(nRealx+2); phiT.resize(nRealx+2);

    for ( int i = 0 ; i <= nRealx + 1 ; ++i ) { phiB[i] = 10.; phiT[i] = -10.; }
    for ( int j = 0 ; j <= nRealy + 1 ; ++j ) { phiL[j] = 1.; phiR[j] = -1.; }

    // Put a diagonal on each ghost cell, regardless of PE.  The actual phi values used here
    // will not be used in serial cases and will get overwritten during parallel communication.

    ApplyBCs(   0     , nRealx+1 ,        0 ,        0 , phiB );
    ApplyBCs(   0     , nRealx+1 , nRealy+1 , nRealy+1 , phiT );
    ApplyBCs(   0     ,        0 ,        0 , nRealy+1 , phiL );
    ApplyBCs(nRealx+1 , nRealx+1 ,        0 , nRealy+1 , phiR );
    
    // Now, apply actual BCs -- Dirichlet on top and bottom and left
    
    if ( myMPI.jPE  == 0             ) ApplyBCs(  0      ,   nRealx+1 ,   1      , 1          , phiB);
    if ( myMPI.jPE  == myMPI.nPEy-1  ) ApplyBCs(  0      ,   nRealx+1 ,   nRealy , nRealy     , phiT);
    if ( myMPI.iPE  == 0             ) ApplyBCs(  1      ,   1        ,   0      , nRealy+1   , phiL);

      
    // Now, apply actual BCs -- Neumann on right
    
    if ( myMPI.iPE  == myMPI.nPEx-1  ) ApplyNeumannBCs(  nRealx ,   nRealx   ,   0      , nRealy+1   , -1 ,   0 );

    // Other BCs, not being used currently

    // if ( myMPI.iPE  == myMPI.nPEx-1  ) ApplyBCs(  nRealx ,   nRealx   ,   0      , nRealy+1   , phiR);
    // if ( myMPI.jPE  == 0             ) ApplyNeumannBCs(  0      ,   nRealx+1 ,   1      , 1          ,  0 ,   1 );
    // if ( myMPI.jPE  == myMPI.nPEy-1  ) ApplyNeumannBCs(  0      ,   nRealx+1 ,   nRealy , nRealy     ,  0 ,  -1 );
    
      
  }

  void ApplyBCs(int iMin, int iMax , int jMin, int jMax, VD &phiValues)
  {
    int count = 0;
    for ( int i = iMin ; i <= iMax ; ++i )
      for ( int j = jMin ; j <= jMax ; ++j )
	{
	  int p = pid(i,j);
	  cLOOP Acoef[p][c] = 0.;
	  Acoef[ p ] [ 1 ] = 1. ;  b[ p ] =  phiValues[count];
	  ++count;
	}

    return;
  }

  void ApplyNeumannBCs(int iMin, int iMax , int jMin, int jMax , int dirx , int diry)
  {
    for ( int i = iMin ; i <= iMax ; ++i )
      for ( int j = jMin ; j <= jMax ; ++j )
	{
	  int p  = pid(i,j);
	  int p2 = pid(i+dirx,j+diry);
	  
	  cLOOP Acoef[p][c] = 0.;
	  
	  Acoef[ p ] [ 1 ] = 1. ;

	  for ( int k = 1 ; k <= bandwidth ; ++k )
	    if ( Jcoef[p][k] == pid(i+dirx,j+diry) ) Acoef[ p ] [ k ] = -1.;

	  b[ p ] = 0.;
	}

    return;
  }

  //  ==
  //  ||
  //  ||  ParticlesOnMesh
  //  ||
  //  ||  Quantify the particles' interaction with the mesh
  //  ||
  //  ||  References:  [1] https://www.particleincell.com/2010/es-pic-method/
  //  ||
  //  ==

  void ParticlesOnMesh(particles &PTCL, mpiInfo &myMPI)
  {
    double hx, hy;
    double w[5];
    int    p[5];
    int    iL, iR, jB, jT;
    int    iPEnew, jPEnew;            // These store the i-j indicies of the processor receiving a particle,
                                      //    if that particle is leaving the mesh.
    VI ptcl_send_list, ptcl_send_PE;  // These collect information about particles that
                                      //    have left this processor and are heading onto
                                      //    another processor.

    // -
    // |
    // | Determine which particles are still on this mesh and which have left
    // |
    // -

    for ( int k = 1 ; k <= PTCL.n ; ++k )
      {
	// First, check to be sure the particle is still in the mesh.  If it is not, set its
	// "active" flag to zero, and note the processor to which it is going for MPI exchange.
	
	if ( PTCL.active[k] == 1 )
	  {
	    iPEnew = myMPI.iPE;
	    jPEnew = myMPI.jPE;

	    
	    if ( PTCL.x[k] < x0     ) { PTCL.active[k] = -1;  iPEnew = myMPI.iPE - 1 ; } 
	    if ( PTCL.x[k] > x1     ) { PTCL.active[k] = -1;  iPEnew = myMPI.iPE + 1 ; }
	    if ( PTCL.y[k] < y0     ) { PTCL.active[k] = -1;  jPEnew = myMPI.jPE - 1 ; }
	    if ( PTCL.y[k] > y1     ) { PTCL.active[k] = -1;  jPEnew = myMPI.jPE + 1 ; }
	    
	  }

	// The particle is not in the mesh.  Collect this particle into a holding array that will 
	// be sent to the neighboring processor.

	if ( PTCL.active[k] == -1)
	  {
	    if ( iPEnew >= 0 && iPEnew < myMPI.nPEx )
	    if ( jPEnew >= 0 && jPEnew < myMPI.nPEy )
	      {
		ptcl_send_list.push_back(k);
		ptcl_send_PE  .push_back(iPEnew + jPEnew * myMPI.nPEx); 
	      }

	    PTCL.active[k] = 0;  // Remove it from the list of active particles
	  }
      }

    // -
    // |
    // | Give and receive particles to/with other processors
    // |
    // -

    myMPI.ParticleExchange( ptcl_send_list , ptcl_send_PE , PTCL);
    
    // -
    // |
    // | Accumulate particles to the nodes
    // |
    // -

    for ( int k = 1 ; k <= nField ; ++k ) Qval[k] = 0.;

    for ( int k = 1 ; k <= PTCL.n ; ++k )
      {
	if ( PTCL.active[k] == 1 )
	  {
	    iL = int ( (PTCL.x[k]-x0) / dx ) + 1;    // point to the left
	    jB = int ( (PTCL.y[k]-y0) / dy ) + 1;    // point below 
	    iR = iL + 1;                             // point to the right
	    jT = jB + 1;                             // point above
	    
	    // 1-2-3-4 are temporary numbers that refer to the lower left, lower right,
	    // upper right, upper left, in that order.
	    
     	    p[1] = pid(iL,jB);
     	    p[2] = pid(iR,jB);
     	    p[3] = pid(iR,jT);
     	    p[4] = pid(iL,jT);
	    
     	    // Compute weights for spreading particle k to the four surrounding nodes.
     	    // Here, hx is the fractional x-distance from particle k to node 1.  Similar for hy.
     	    // See [1].
	    
     	    hx = (PTCL.x[k] - x[ p[1] ])/dx;
     	    hy = (PTCL.y[k] - y[ p[1] ])/dy;
	    
     	    w[1] = ( 1. - hx ) * ( 1. - hy );
     	    w[2] =        hx   * ( 1. - hy );
     	    w[3] =        hx   *        hy  ;
     	    w[4] = ( 1. - hx ) *        hy  ;
	    
     	    // Spread particle k to the four surrounding points
	
	    for ( int i = 1 ; i <= 4 ; ++i ) Qval[ p[i] ] +=  w[i] * PTCL.Qp;
     	  }
       }

    // -
    // |
    // | Sum Qval on processor boundaries
    // |
    // -

    // For debugging -- set value to myPE
    // iLOOP jLOOP { int p = pid(i,j) ; Qval[ p ] = myPE + 1; }
    
    myMPI.PEsum(Qval); 
    
    // -
    // |
    // | Compute forces on particles
    // |
    // -
    
    for ( int k = 1 ; k <= PTCL.n ; ++k ) PTCL.xf[k] = PTCL.yf[k] = 0.;

    for ( int k = 1 ; k <= PTCL.n ; ++k )
      {
	if ( PTCL.active[k] == 1 )
	  {
	    iL = int ( (PTCL.x[k]-x0) / dx ) + 1;    // point to the left
	    jB = int ( (PTCL.y[k]-y0) / dy ) + 1;    // point below 
	    iR = iL + 1;                             // point to the right
	    jT = jB + 1;                             // point above

	    // Compute a phi value for the top, bottom right, and left sides
	    // of the cell in which this particle lies by averaging the
	    // endpoints that define that side.
	    //
	    //     -->   phiT  <--   (average of top two node)
	    //     +-------------+
	    //     |             |
	    //     |        *    |
	    //     |    particle | phiR  (is the average of the two right-hand side nodes)
	    //     |             |
	    //     |             |
	    //     +-------------+
	    //     -->   phiB  <--  (average of bottom two nodes)
	    //
	    //

	    double phiT = ( phi[ pid(iL,jT) ] + phi[ pid(iR,jT) ] )/2.;
	    double phiB = ( phi[ pid(iL,jB) ] + phi[ pid(iR,jB) ] )/2.;
	    double phiL = ( phi[ pid(iL,jT) ] + phi[ pid(iL,jB) ] )/2.;
	    double phiR = ( phi[ pid(iR,jT) ] + phi[ pid(iR,jB) ] )/2.;

	    // Compute the gradient of the electric field using the average phi
	    // values on each of the four sides
	    
	    PTCL.xf[ k ] = -0.01*( phiR - phiL ) / dx;
	    PTCL.yf[ k ] = -0.01*( phiT - phiB ) / dy;

      // Bug fix?
	    // PTCL.xf[ k ] *= -1;
	    // PTCL.yf[ k ] *= -1;
     	  }
       }


  }
  

  //  ==
  //  ||
  //  ||  Utility routines
  //  ||
  //  ==

  int pid(int i,int j) { return (i+1) + (j)*(nRealx+2); }  // Given i-j, return point ID.  Here i-j is the physical grid.
                                                           // The row/col numbers must include the candy-coating.
  #include "plotter.h"
  #include "gauss_seidel.h"
};



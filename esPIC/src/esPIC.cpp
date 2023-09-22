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

#include "mpi.h"
#include "esPIC.h"
#include "particles.h"
#include "mpiInfo.h"
#include "LaplacianOnGrid.h"


//  ==
//  ||
//  ||
//  ||  Main Program
//  ||
//  ||
//  ==

int main(int argc, char *argv[])
{
  
   mpiInfo myMPI;
   MPI_Init     (&argc         , &argv       );
   MPI_Comm_size(MPI_COMM_WORLD, &myMPI.numPE);
   MPI_Comm_rank(MPI_COMM_WORLD,&myMPI.myPE  );


   int nPEx, nPEy, nCellx, nCelly, nPtcl ;
   double tEnd, dt, tPlot;
   double flux;


    if ( myMPI.myPE == 0 )
      {
       cout << "\n";
       cout << "------------------------------------------------------\n";
       cout << "\n";
       cout << " E L E C T R O S T I C   P A R T I C L E - I N - C E L L\n";
       cout << " D E M O   C O D E                           \n";
       cout << "\n";
       cout << " Running on " << myMPI.numPE << " processors \n";
       cout << "\n";
       cout << "------------------------------------------------------\n";
       cout << "\n";
      }
    

     
   for (int count =  0 ; count < argc; ++count)
     {
       if ( !strcmp(argv[count],"-nPEx"    ) ) nPEx   = atoi(argv[count+1]);
       if ( !strcmp(argv[count],"-nPEy"    ) ) nPEy   = atoi(argv[count+1]);
       if ( !strcmp(argv[count],"-nCellx"  ) ) nCellx = atoi(argv[count+1]);
       if ( !strcmp(argv[count],"-nCelly"  ) ) nCelly = atoi(argv[count+1]);
       if ( !strcmp(argv[count],"-nPtcl"   ) ) nPtcl  = atoi(argv[count+1]);
       if ( !strcmp(argv[count],"-flux"    ) ) flux   = atof(argv[count+1]);
       if ( !strcmp(argv[count],"-tEnd"    ) ) tEnd   = atof(argv[count+1]);
       if ( !strcmp(argv[count],"-dt"      ) ) dt     = atof(argv[count+1]);
       if ( !strcmp(argv[count],"-tPlot"   ) ) tPlot  = atof(argv[count+1]);
     }

   if ( myMPI.myPE == 0 )
     {
       cout << endl;
       cout << "Input Summary: " << endl;
       cout << "--------------------------- " << endl;
       cout << "No. PE   in  x-direction: " << nPEx   << endl;
       cout << "             y-direction: " << nPEy   << endl;
       cout << "No. Cells in x-direction: " << nCellx << endl;
       cout << "             y-direction: " << nCelly << endl;
       cout << "Flux density            : " << flux   << endl;
       cout << "End Time                : " << tEnd   << endl;
       cout << "Time Step               : " << dt     << endl;
       cout << "Plot Interval           : " << tPlot  << endl;
       cout << endl;
     }

   myMPI.GridDecomposition(nPEx,nPEy,nCellx,nCelly);

   // -
   // |
   // | Parallel Grid Generation and Laplace Solver
   // |
   // -
   
   double totalLength = 1.;
   double eachPElength_x = totalLength / nPEx;
   double eachPElength_y = totalLength / nPEy;

   double x0 = eachPElength_x * myMPI.iPE;   double x1 = x0 + eachPElength_x;
   double y0 = eachPElength_y * myMPI.jPE;   double y1 = y0 + eachPElength_y;
   
   LaplacianOnGrid MESH(x0,x1,y0,y1,nCellx,nCelly, myMPI );


   // -
   // |
   // | Set up Particles
   // |
   // -
   
   particles PTCL(nPtcl); 

   int    count   = 0;
   double density = flux/(MESH.y1-MESH.y0);
   double vx_bdy  = .5;
   double vy_bdy  = .0;

   

   // -
   // |
   // | Time Marching Loop
   // |
   // -

   struct timespec t0, t1;
   StartTimer(t0);

   double timeSinceLastPlot = 0.;
   int latestIterCount;
   
   for ( double t = 0. ; t <= tEnd ; t += dt )
     {
       

       // Inject particles
       
       if ( myMPI.iPE == 0 ) PTCL.addFlux( t , MESH.y0, MESH.y1, density, vx_bdy , vy_bdy );

       // Move particles
       
       PTCL.move( dt );

       // Map between particles and the mesh

       MESH.ParticlesOnMesh(PTCL,myMPI);

       // Solve Poisson Equation

       MPI_Barrier(MPI_COMM_WORLD);
       
       MESH.FormLS(myMPI);

       MPI_Barrier(MPI_COMM_WORLD);

       MESH.GS_or_Jacobi(5000, MESH.b , MESH.phi , myMPI , 2, latestIterCount);

       // Plot
 
       timeSinceLastPlot += dt;
       if ( timeSinceLastPlot >= tPlot )
	 {
	   if ( myMPI.myPE == 0 ) cout << "Plotting output at time = " << t << " Latest matrix solver iteration count = " << latestIterCount << endl;
	   PTCL.plot( "ptcl" ,            count , myMPI.myPE );
	   MESH.plot( "mesh" , MESH.Qval, count , myMPI      );
	   MESH.plot( "phi"  , MESH.phi , count , myMPI      );
	   timeSinceLastPlot = 0.;
	   ++count;
	 }
     }


   EndTimer("main", t0, t1);


   MPI_Finalize();
   return 0;

}

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

void plot(string descriptor, int timeIdx , int myPE)
  {
    // Convert my processor number to a string:
    
    std::stringstream myPE_str;  myPE_str << myPE;
    std::stringstream Idx_str ;  Idx_str  << timeIdx;
    
    // Create plot filename:
    //    Filename is the following concatenated information:
    //    descriptor + "_" + myPE + "_" + timeIdx + ".plt"

    string filename = descriptor + "_" + myPE_str.str() + "_" + Idx_str.str() + ".plt";

    // Open plot file
    
    std::fstream f;
    f.open(filename.c_str(),std::ios::out);

    // Write to plot file

    for ( int i = 1 ; i <= n ; ++i )
      {
	if ( active[i] == 1 )
	  {
	    f << x[i] << " " << y[i] << "   0.  " << endl;
	  }
      }

    f.close();
  }


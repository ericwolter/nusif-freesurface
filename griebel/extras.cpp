# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>

using namespace std;

# include "datadef.hpp"
# include "extras.hpp"

//****************************************************************************80

void output_one_real ( REAL **field, int imin, int imax, int jmin, int jmax, 
  char *output_filename )

//****************************************************************************80
//
//  Purpose:
//
//    OUTPUT_ONE_REAL outputs a single real field.
//
//  Modified:
//
//    16 March 2004
//
//  Author:
//
//    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer
//
//  Reference:
//
//    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
//    Numerical Simulation in Fluid Dynamics,
//    SIAM, 1998.
//
//  Parameters:
//
//    Input, REAL **FIELD, the field to be written out.
//
//    Input, int IMIN, IMAX, JMIN, JMAX, specify that
//    FIELD(IMIN:IMAX,JMIN:JMAX) is to be written out.
//
//    Input, char *OUTPUT_FILENAME, the name of the file
//    into which the data is to be written.
//
{
  ofstream output_file;
  int i;
  int j;

  output_file.open ( output_filename );

  if ( !output_file )
  {
    cout << "\n";
    cout << "OUTPUT_ONE_REAL - Fatal error!\n";
    cout << "  Cannot open the input file \"" << output_filename << "\".\n";
    return;
  }

  for ( j = jmin; j <= jmax; j++ )
  {
    for ( i = imin; i <= imax; i++)
    {
      output_file << setw(6)  << i           << "  "
                  << setw(6)  << j           << "  "
                  << setw(12) << field[i][j] << "\n";
    }
  }

  output_file.close ( );

  return;
}

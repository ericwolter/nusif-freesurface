# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>

using namespace std;

# include "datadef.hpp"
# include "boundary.hpp"

//****************************************************************************80

void SETBCOND ( REAL **U, REAL **V, REAL **P, REAL **TEMP, int **FLAG,
  int imax, int jmax, int wW, int wE, int wN, int wS )

//****************************************************************************80
//
//  Purpose:
//
//    SETBCOND sets the boundary conditions at the boundary strip.
//
//  Discussion:
//
//  The flags wW,wE,wN, and wS can have the values:
//  1 = slip               2 = no-slip
//  3 = outflow            4 = periodic
//  Moreover, no-slip conditions are set at internal obstacle cells
//  by default.
//  For temperature, adiabatic boundary conditions are set.
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
{
  int i;
  int j;

//
//  western and eastern boundary.
//
  for(j=0;j<=jmax+1;j++)
  {
    if(wW == 1 )
    {
      U[0][j] = 0.0;
      V[0][j] = V[1][j];
    }
    if(wW == 2 )
    {
      U[0][j] = 0.0;
      V[0][j] = (-1.0)*V[1][j];
    }
    if(wW == 3)
    {
      U[0][j] = U[1][j];
      V[0][j] = V[1][j];
    }
    if(wW == 4 )
    {
      U[0][j] = U[imax-1][j];
      V[0][j] = V[imax-1][j];
      V[1][j] = V[imax][j];
      P[1][j] = P[imax][j];
    }

    TEMP[0][j] = TEMP[1][j];

    if(wE == 1 )
    {
      U[imax][j] = 0.0;
      V[imax+1][j] = V[imax][j];
    }
    if(wE == 2 )
    {
      U[imax][j] = 0.0;
      V[imax+1][j] = (-1.0)*V[imax][j];
    }
    if(wE == 3)
    {
      U[imax][j] = U[imax-1][j];
      V[imax+1][j] = V[imax][j];
    }
    if(wE == 4 )
    {
      U[imax][j] = U[1][j];
      V[imax+1][j] = V[2][j];
    }

   TEMP[imax+1][j] = TEMP[imax][j];
  }
//
//  northern and southern boundary.
//
  for(i=0;i<=imax+1;i++)
  {
    if(wN == 1 )
    {
      V[i][jmax] = 0.0;
      U[i][jmax+1] = U[i][jmax];
    }
    if(wN == 2 )
    {
      V[i][jmax] = 0.0;
      U[i][jmax+1] = (-1.0)*U[i][jmax];
    }
    if(wN == 3)
    {
      V[i][jmax] = V[i][jmax-1];
      U[i][jmax+1] = U[i][jmax];
    }
    if(wN == 4 )
    {
      V[i][jmax] = V[i][1];
      U[i][jmax+1] = U[i][2];
    }

    TEMP[i][0] = TEMP[i][1];

    if(wS == 1 )
    {
      V[i][0] = 0.0;
      U[i][0] = U[i][1];
    }
    if(wS == 2 )
    {
      V[i][0] = 0.0;
      U[i][0] = (-1.0)*U[i][1];
    }
    if(wS == 3)
    {
      V[i][0] = V[i][1];
      U[i][0] = U[i][1];
    }
    if(wS == 4 )
    {
      V[i][0] = V[i][jmax-1];
      U[i][0] = U[i][jmax-1];
      U[i][1] = U[i][jmax];
      P[i][1] = P[i][jmax];
    }

   TEMP[i][jmax+1] = TEMP[i][jmax];
  }
//
//  setting the boundary values at inner obstacle cells
//  (only no-slip)
//
 for(i=1;i<=imax;i++)
   for(j=1;j<=jmax;j++)
     if(FLAG[i][j] & 0x000f)
//
//  The mask 0x000f filters the obstacle cells adjacent to fluid cells.
//
       switch (FLAG[i][j])
	 {
          case B_N:  {
		       V[i][j]   = 0.0;
                       U[i][j]   = -U[i][j+1];
                       U[i-1][j] = -U[i-1][j+1];
                       TEMP[i][j] = TEMP[i][j+1];
                       break;
	             }
          case B_O:  {
		       U[i][j]   = 0.0;
                       V[i][j]   = -V[i+1][j];
                       V[i][j-1] = -V[i+1][j-1];
                       TEMP[i][j] = TEMP[i+1][j];
                       break;
	             }
          case B_S:  {
		       V[i][j-1] = 0.0;
                       U[i][j]   = -U[i][j-1];
                       U[i-1][j] = -U[i-1][j-1];
                       TEMP[i][j] = TEMP[i][j-1];
                       break;
	             }
          case B_W:  {
		       U[i-1][j] = 0.0;
                       V[i][j]   = -V[i-1][j];
                       V[i][j-1] = -V[i-1][j-1];
                       TEMP[i][j] = TEMP[i-1][j];
                       break;
	             }
          case B_NO: {
		       V[i][j]   = 0.0;
                       U[i][j]   = 0.0;
                       V[i][j-1] = -V[i+1][j-1];
                       U[i-1][j] = -U[i-1][j+1];
                       TEMP[i][j] = 0.5*(TEMP[i][j+1]+TEMP[i+1][j]);
                       break;
	             }
          case B_SO: {
		       V[i][j-1] = 0.0;
                       U[i][j]   = 0.0;
                       V[i][j]   = -V[i+1][j];
                       U[i-1][j] = -U[i-1][j-1];
                       TEMP[i][j] = 0.5*(TEMP[i][j-1]+TEMP[i+1][j]);
                       break;
                     }
          case B_SW: {
		       V[i][j-1] = 0.0;
                       U[i-1][j] = 0.0;
                       V[i][j]   = -V[i-1][j];
                       U[i][j]   = -U[i][j-1];
                       TEMP[i][j] = 0.5*(TEMP[i][j-1]+TEMP[i-1][j]);
                       break;
	             }
          case B_NW: {
		       V[i][j]   = 0.0;
                       U[i-1][j] = 0.0;
                       V[i][j-1] = -V[i-1][j-1];
                       U[i][j]   = -U[i][j+1];
                       TEMP[i][j] = 0.5*(TEMP[i][j+1]+TEMP[i-1][j]);
                       break;
	             }
	  default : break;
  }
}
//****************************************************************************80

void SETSPECBCOND ( char* problem, REAL **U, REAL **V, REAL **P, REAL **TEMP,
  int imax, int jmax, REAL UI, REAL VI )

//****************************************************************************80
//
//  Purpose:
//
//    SETSPECBCOND sets specific boundary conditions, depending on "problem".
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
{
  int i;
  int j;
//
//  Flow past a backward facing step without free boundary.
//  U = 1.0 at the left boundary.
//
  if ( strcmp(problem, "backstep" ) == 0 )
  {
    for(j=(jmax/2)+1;j<=jmax;j++)
    {
      U[0][j] = 1.0;
    }
    return;
  }
//
//  Flow past an obstacle: U = 1.0 at left boundary.
//
  else if ( strcmp(problem, "circle")==0 )
  {
     V[0][0] = 2*VI-V[1][0];
     for(j=1;j<=jmax;j++)
       {
        U[0][j] = UI;
        V[0][j] = 2*VI-V[1][j];
       }
     return;
  }
//
//  natural convection; left T = 0.5 right T = -0.5
//  left wall heated, right wall cooled.
//  upper and lower wall adiabatic
//
  else if ( strcmp ( problem, "convection" ) == 0 )
  {
    for ( j = 0; j <= jmax+1; j++ )
    {
	TEMP[0][j] = 2*(0.5)-TEMP[1][j];
	TEMP[imax+1][j] = 2*(-0.5)-TEMP[imax][j];
    }
    for(i=0;i<=imax+1;i++)
    {
      TEMP[i][0] = TEMP[i][1];
      TEMP[i][jmax+1] = TEMP[i][jmax];
    }
    return;
  }
//
//  The breaking dam.
//
  else if ( strcmp(problem,"dam") == 0 )
  {
    return;
  }
//
//  Driven Cavity: U = 1.0 at the upper boundary.
//
  else if ( strcmp(problem, "dcavity")==0 )
  {
    for(i=0;i<=imax;i++)
    {
      U[i][jmax+1] = 2.0 - U[i][jmax];
    }
    return;
  }
//
//  The fluid drop.
//
  else if ( strcmp(problem,"drop")==0 )
  {
    return;
  }
//
//  fluidtrap: left T = 0.5 right T = -0.5
//  left wall heated, right wall cooled,
//  upper and lower wall adiabatic
//
  else if( strcmp(problem, "fluidtrap")==0)
  {
     for(j=0; j<=jmax+1; j++)
       {
	TEMP[0][j] = 2*(0.5)-TEMP[1][j];
	TEMP[imax+1][j] = 2*(-0.5)-TEMP[imax][j];
       }
     for(i=0;i<=imax+1;i++)
       {
	TEMP[i][0] = TEMP[i][1];
	TEMP[i][jmax+1] = TEMP[i][jmax];
       }
     return;
  }
//
//  Inflow for injection molding: U = 1.0 in the mid of left boundary.
//
  else if(strcmp(problem, "molding")==0)
  {
     for (j=(int)(0.4*jmax)+1;j<=(int)(0.6*jmax);j++)
        U[0][j] = 1.0;
     return;
  }
//
//  Flow past an obstacle: U = 1.0 at left boundary.
//
  else if( strcmp(problem, "plate")==0 )
  {
     V[0][0] = 2*VI-V[1][0];
     for(j=1;j<=jmax;j++)
       {
        U[0][j] = UI;
        V[0][j] = 2*VI-V[1][j];
       }
     return;
  }
//
//  Rayleigh-Benard flow: top T = -0.5 bottom T = 0.5
//  left and right walls adiabatic,
//  lower wall heated, upper wall cooled.
//
  else if ( strcmp(problem, "rayleigh")==0 )
  {

     for(j=0; j<=jmax+1; j++)
       {
	TEMP[0][j] = TEMP[1][j];
	TEMP[imax+1][j] = TEMP[imax][j];
       }

     for(i=0;i<=imax+1;i++)
       {
	TEMP[i][0] = 2*(0.5)-TEMP[i][1];
	TEMP[i][jmax+1] = 2*(-0.5)-TEMP[i][jmax];
       }
     return;
  }
//
//  Flow past a backward facing step, with free boundary.
//  U = 1.0 at the left boundary.
//
  else if(strcmp(problem, "wave")==0)
  {
     for(j=(jmax/2)+1;j<=jmax;j++)
        U[0][j]    =   1.0;
     return;
  }

//
//  Else
//
  else
  {
    cout << "\n";
    cout << "SETSPECBCOND - Fatal error!\n";
    cout << "  Problem " << problem << " not defined!\n";
    exit ( 1 );
  }

  return;
}

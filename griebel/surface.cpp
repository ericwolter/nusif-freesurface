# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>

using namespace std;

# include "datadef.hpp"
# include "surface.hpp"
# include "visual.hpp"

//****************************************************************************80

struct particleline *INIT_PARTICLES ( int *N, int imax, int jmax,
  REAL delx, REAL dely, int ppc, char *problem, REAL **U, REAL **V )

//****************************************************************************80
//
//  Purpose:
//
//    INIT_PARTICLES initializes particles for free boundary problems.
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
 int i,j,ip,jp;
 struct particleline *Particlelines;
 REAL x,y;
//
//  Sizes for the "drop" problem.
//
 REAL height=0,rad=0,mpx=0,mpy=0,vstart=0;
//
//  Initialization of some parameters.
//
 if(strcmp(problem, "dam")==0)
    *N = 1;
//
//  Height of the basin, radius of the drop, coordinates of the center of the drop,
//  and the initial vertical velocity of the drop:
//
 if(strcmp(problem, "drop")==0)
   {
    *N = 2;
    height = 1./2.*jmax*dely;
    rad    = 0.1*jmax*dely;
    mpx    = 0.5*imax*delx;
    mpy    = 2./3.*jmax*dely;
    vstart = -2.0;
   }

 if((Particlelines=(struct particleline *)
             malloc((unsigned)(*N) * sizeof(struct particleline))) == NULL)
   {
    cout << "\n";
    cout << "INIT_PARTICLES - Fatal error!\n";
    cout << "  No memory!\n";
    exit(0);
   }
 Particlelines -= 1;
//
//  Particlelines from 1 to N.
//
 for (i=1;i<=*N;i++)
   {
    Particlelines[i].length = 0;
    Particlelines[i].Particles = PARTALLOC(-1.,-1.);
   }

//
//  Set the particles.
//
 for (i=1;i<=imax;i++)
    for (j=1;j<=jmax;j++)
       for (ip=1;ip<=ppc;ip++)
         {
          x = (i-1)*delx+(ip-.5)/((REAL)ppc)*delx;
          for (jp=1;jp<=ppc;jp++)
	    {
             y = (j-1)*dely+(jp-.5)/((REAL)ppc)*dely;

             if(strcmp(problem, "dam")==0)
                if (x<0.2*imax*delx)
                   SET_PART(&Particlelines[1],x,y);

             if(strcmp(problem, "drop")==0)
	       {
                if (y<height)
		  {
                   SET_PART(&Particlelines[1],x,y);
		  }
                else if ((x-mpx)*(x-mpx)+(y-mpy)*(y-mpy) <= rad*rad)
		  {
                   SET_PART(&Particlelines[2],x,y);
                   V[i][j] = vstart;
		  }
	       }
	    }
	 }
 return (Particlelines);
}
//****************************************************************************80

void SET_PART(struct particleline *Partline, REAL x,REAL y)

//****************************************************************************80
//
//  Purpose:
//
//    SET_PART adds particle to "Partline" at (x,y).
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
 struct particle *part;
//
//  Create a particle.
//
 part =  PARTALLOC(x,y);
//
//  Add it to PARTLINE in the first position.
//
 part->next = (*Partline).Particles->next;
 (*Partline).Particles->next = part;
 (*Partline).length++;
}
//****************************************************************************80

void MARK_CELLS(int **FLAG,int imax,int jmax,REAL delx,REAL dely,
                int *ifull,int *isurf,
                int N, struct particleline *Particlelines)

//****************************************************************************80
//
//  Purpose:
//
//    MARK_CELLS marks the cells of the fluid domain.
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
 int i,j,n;
 REAL x,y;
 struct particle *part,*help;
//
//  set all cells which are not obstacle cells to empty cells.
//
 for (i=0;i<=imax+1;i++)
    for (j=0;j<=jmax+1;j++)
       if (FLAG[i][j] >= C_F)
          FLAG[i][j] = (FLAG[i][j] | C_E) & ~C_NSWO;
//
//  Mark cells containing particles as fluid cells (loop over particles) .
//
 for (n=1;n<=N;n++)
    for(part=Particlelines[n].Particles; part->next != NULL; part=part->next)
      {
       x = part->next->x;
       y = part->next->y;

       i = (int)(x/delx)+1;
       j = (int)(y/dely)+1;

       if (FLAG[i][j] < C_F)
         {
//
//  delete particles in obstacle cells.
//
          help = part->next->next;
          free(part->next);
          part->next = help;
          Particlelines[n].length--;
	 }
       else
          FLAG[i][j] = FLAG[i][j] & ~C_E;
      }
//
//  Mark surface cells.
//
 (*ifull) = 0;
 (*isurf) = 0;
 for (j=1;j<=jmax;j++)
    for (i=1;i<=imax;i++)
      {
       if( (FLAG[i][j] & C_F) && (FLAG[i][j] < C_E) )
         {
          if (FLAG[i-1][j] & C_E)
             FLAG[i][j] = FLAG[i][j] | C_W;
          if (FLAG[i+1][j] & C_E)
             FLAG[i][j] = FLAG[i][j] | C_O;
          if (FLAG[i][j-1] & C_E)
             FLAG[i][j] = FLAG[i][j] | C_S;
          if (FLAG[i][j+1] & C_E)
             FLAG[i][j] = FLAG[i][j] | C_N;

          if (FLAG[i][j] < 0x0100)
             (*ifull)++;
          else
             (*isurf)++;
         }
      }
// printf ("\nGeometry of the fluid domain:\n\n");
//  for(j=jmax+1;j>=0;j--)
//    {
//     for(i=0;i<=imax+1;i++)
//        printf("%.4x ",FLAG[i][j]);
//     printf ("\n");
//    }
//
}
//****************************************************************************80

void SET_UVP_SURFACE(REAL **U,REAL **V,REAL **P,int **FLAG,REAL GX,REAL GY,
                     int imax,int jmax,REAL Re,REAL delx,REAL dely,REAL delt)

//****************************************************************************80
//
//  Purpose:
//
//    SET_UVP_SURFACE sets the boundary values at a free surface.
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
 int i,j;
//
//  Set velocity values in empty cells to zero.
//
 for (j=1;j<=jmax;j++)
    for (i=1;i<=imax-1;i++)
       if ((FLAG[i][j] & C_E) && FLAG[i+1][j] & C_E)
          U[i][j] = 0.0;
 for (j=1;j<=jmax-1;j++)
    for (i=1;i<=imax;i++)
       if ((FLAG[i][j] & C_E) && FLAG[i][j+1] & C_E)
          V[i][j] = 0.0;

 for (j=1;j<=jmax;j++)
    for (i=1;i<=imax;i++){
//
//  treat only surface cells.
//
       if ( (!(FLAG[i][j] & C_E)) || (FLAG[i][j] < 0x0100) )
         switch (FLAG[i][j] & C_NSWO){
//
//  mask NSWO_E=0x0f00 filters surface cells.
//
          case C_N   :{ V[i][j] = V[i][j-1]-dely/delx*(U[i][j]-U[i-1][j]);
                        if (FLAG[i-1][j+1] & C_E)
                          U[i-1][j+1] =
                                U[i-1][j]-dely/delx*(V[i][j]-V[i-1][j]);
                      } break;
          case C_S   :{ V[i][j-1] = V[i][j]+dely/delx*(U[i][j]-U[i-1][j]);
                        if (FLAG[i-1][j-1] & C_E)
                           U[i-1][j-1] =
                                 U[i-1][j]+dely/delx*(V[i][j-1]-V[i-1][j-1]);
	              } break;
	  case C_O   :{ U[i][j] = U[i-1][j]-delx/dely*(V[i][j]-V[i][j-1]);
                        if (FLAG[i+1][j-1] & C_E)
                           V[i+1][j-1] =
                                 V[i][j-1]-delx/dely*(U[i][j]-U[i][j-1]);
		      } break;
	  case C_W :{ U[i-1][j] = U[i][j]+delx/dely*(V[i][j]-V[i][j-1]);
                        if (FLAG[i-1][j-1] & C_E)
                           V[i-1][j-1] =
                                 V[i][j-1]+delx/dely*(U[i-1][j]-U[i-1][j-1]);
		      } break;

          case C_NO  :{ U[i][j]     = U[i-1][j];
                        V[i][j]     = V[i][j-1];
                        if (FLAG[i-1][j+1] & C_E)
                           U[i-1][j+1] =
                                 U[i-1][j]-dely/delx*(V[i][j]-V[i-1][j]);
                        if (FLAG[i+1][j+1] & C_E)
			  {
                           U[i][j+1]   = U[i][j];
                           V[i+1][j]   = V[i][j];
			  }
                        if (FLAG[i+1][j-1] & C_E)
                           V[i+1][j-1] =
                                 V[i][j-1]-delx/dely*(U[i][j]-U[i][j-1]);
		      } break;
          case C_NW  :{ U[i-1][j]   = U[i][j];
                        V[i][j]     = V[i][j-1];
                        if (FLAG[i-1][j+1] & C_E)
			  {
                           U[i-1][j+1] = U[i-1][j];
                           V[i-1][j]   = V[i][j];
			  }
                        if (FLAG[i-1][j-1] & C_E)
                           V[i-1][j-1] =
                                 V[i][j-1]+delx/dely*(U[i-1][j]-U[i-1][j-1]);
		      } break;
          case C_SW  :{ U[i-1][j]   = U[i][j];
                        V[i][j-1]   = V[i][j];
                        if (FLAG[i-1][j-1] & C_E)
			  {
                           U[i-1][j-1] = U[i-1][j];
                           V[i-1][j-1] = V[i][j-1];
			  }
		      } break;
          case C_SO  :{ U[i][j]     = U[i-1][j];
                        V[i][j-1]   = V[i][j];
                        if (FLAG[i-1][j-1] & C_E)
                           U[i-1][j-1] =
                                 U[i-1][j]+dely/delx*(V[i][j-1]-V[i-1][j-1]);
                        if (FLAG[i+1][j-1] & C_E)
			  {
                           U[i][j-1]   = U[i][j];
                           V[i+1][j-1] = V[i][j-1];
			  }
		      } break;
          case C_WO  :{ U[i][j]     += delt*GX;
                        U[i-1][j]   += delt*GX;
                        if (FLAG[i-1][j-1] & C_E)
                           V[i-1][j-1]  =
                                 V[i][j-1]+delx/dely*(U[i-1][j]-U[i-1][j-1]);
                        if (FLAG[i+1][j-1] & C_E)
                           V[i+1][j-1]  =
                                 V[i][j-1]-delx/dely*(U[i][j]-U[i][j-1]);
		      } break;
          case C_NS  :{ V[i][j]     += delt*GY;
                        V[i][j-1]   += delt*GY;
                        if (FLAG[i-1][j+1] & C_E)
                           U[i-1][j+1]  =
                                 U[i-1][j]-dely/delx*(V[i][j]-V[i-1][j]);
                        if (FLAG[i-1][j-1] & C_E)
                           U[i-1][j-1]  =
                                 U[i-1][j]+dely/delx*(V[i][j-1]-V[i-1][j-1]);
		      } break;

          case C_NWO :{ V[i][j]      = V[i][j-1]-dely/delx*(U[i][j]-U[i-1][j]);
                        U[i][j]     += delt*GX;
                        U[i-1][j]   += delt*GX;
                        if (FLAG[i-1][j-1] & C_E)
                           V[i-1][j-1]  =
                                 V[i][j-1]+delx/dely*(U[i-1][j]-U[i-1][j-1]);
                        if (FLAG[i+1][j-1] & C_E)
                           V[i+1][j-1]  =
                                 V[i][j-1]-delx/dely*(U[i][j]-U[i][j-1]);
                        if (FLAG[i-1][j+1] & C_E)
			  {
                           V[i-1][j]    = V[i][j];
                           U[i-1][j+1]  = U[i-1][j];
			  }
                        if (FLAG[i+1][j+1] & C_E)
			  {
                           V[i+1][j]    = V[i][j];
                           U[i][j+1]    = U[i][j];
                          }
		      } break;
          case C_NSW :{ U[i-1][j]    = U[i][j]+delx/dely*(V[i][j]-V[i][j-1]);
                        V[i][j]     += delt*GY;
                        V[i][j-1]   += delt*GY;
                        if (FLAG[i-1][j-1] & C_E)
			  {
                           V[i-1][j-1]  = V[i][j-1];
                           U[i-1][j-1]  = U[i-1][j];
			  }
                        if (FLAG[i-1][j+1] & C_E)
			  {
                           V[i-1][j]    = V[i][j];
                           U[i-1][j+1]  = U[i-1][j];
			  }
		      } break;
          case C_SWO :{ V[i][j-1]    = V[i][j]+dely/delx*(U[i][j]-U[i-1][j]);
                        U[i][j]     += delt*GX;
                        U[i-1][j]   += delt*GX;
                        if (FLAG[i-1][j-1] & C_E)
			  {
                           U[i-1][j-1]  = U[i-1][j];
                           V[i-1][j-1]  = V[i][j-1];
			  }
                        if (FLAG[i+1][j-1] & C_E)
			  {
                           U[i][j-1]    = U[i][j];
                           V[i+1][j-1]  = V[i][j-1];
			  }
		      } break;
          case C_NSO :{ U[i][j]      = U[i-1][j]-delx/dely*(V[i][j]-V[i][j-1]);
                        V[i][j]     += delt*GY;
                        V[i-1][j]   += delt*GY;
                        if (FLAG[i-1][j+1] & C_E)
                           U[i-1][j+1]  =
                                 U[i-1][j]-dely/delx*(V[i][j]-V[i-1][j]);
                        if (FLAG[i-1][j-1] & C_E)
                           U[i-1][j-1]  =
                                 U[i-1][j]+dely/delx*(V[i][j-1]-V[i-1][j-1]);
                        if (FLAG[i+1][j-1] & C_E)
			  {
                           U[i][j-1]    = U[i][j];
                           V[i+1][j-1]  = V[i][j-1];
			  }
                        if (FLAG[i+1][j+1] & C_E)
			  {
                           U[i][j+1]    = U[i][j];
                           V[i+1][j]    = V[i][j];
			  }
		      } break;
          case C_NSWO:{ U[i][j]     += delt*GX;
                        U[i-1][j]   += delt*GX;
                        V[i][j]     += delt*GY;
                        V[i][j-1]   += delt*GY;
                        if (FLAG[i-1][j+1] & C_E)
			  {
                           U[i-1][j+1]  = U[i-1][j];
                           V[i-1][j]    = V[i][j];
			  }
                        if (FLAG[i+1][j+1] & C_E)
			  {
                           U[i][j+1]    = U[i][j];
                           V[i+1][j]    = V[i][j];
			  }
                        if (FLAG[i-1][j-1] & C_E)
			  {
                           U[i-1][j-1]  = U[i-1][j];
                           V[i-1][j-1]  = V[i][j-1];
			  }
                        if (FLAG[i+1][j-1] & C_E)
			  {
                           U[i][j-1]    = U[i][j];
                           V[i+1][j-1]  = V[i][j-1];
			  }
		      } break;
	  default     : break;
	 }
      }
//
//  Second loop for pressure boundary values.
//
 for (j=1;j<=jmax;j++)
    for (i=1;i<=imax;i++)
      {
       if (! ((FLAG[i][j] & C_E) || (FLAG[i][j] < 0x0100 )))
         switch (FLAG[i][j] & C_NSWO){
	  case C_N   : P[i][j] = 2./Re/dely*(V[i][j]-V[i][j-1]); break;
	  case C_S   : P[i][j] = 2./Re/dely*(V[i][j]-V[i][j-1]); break;
	  case C_O   : P[i][j] = 2./Re/delx*(U[i][j]-U[i-1][j]); break;
	  case C_W   : P[i][j] = 2./Re/delx*(U[i][j]-U[i-1][j]); break;

          case C_NO  : P[i][j] = 1./Re/2.*
                              ((U[i][j]+U[i-1][j]-U[i][j-1]-U[i-1][j-1])/dely+
                               (V[i][j]+V[i][j-1]-V[i-1][j]-V[i-1][j-1])/delx);
                       break;
          case C_NW  : P[i][j] = -1./Re/2.*
                              ((U[i][j]+U[i-1][j]-U[i][j-1]-U[i-1][j-1])/dely+
                               (V[i+1][j]+V[i+1][j-1]-V[i][j]-V[i][j-1])/delx);
                       break;
          case C_SW  :{ P[i][j] = 1./Re/2.*
                              ((U[i][j+1]+U[i-1][j+1]-U[i][j]-U[i-1][j])/dely+
                               (V[i+1][j]+V[i+1][j-1]-V[i][j]-V[i][j-1])/delx);
		      }
                       break;
          case C_SO  :{ P[i][j] = -1./Re/2.*
                              ((U[i][j+1]+U[i-1][j+1]-U[i][j]-U[i-1][j])/dely+
                               (V[i][j]+V[i][j-1]-V[i-1][j]-V[i-1][j-1])/delx);
		      }
                       break;
          case C_WO  : P[i][j] = 0.; break;
          case C_NS  : P[i][j] = 0.; break;
          case C_NWO : P[i][j] = 0.; break;
          case C_NSW : P[i][j] = 0.; break;
          case C_SWO : P[i][j] = 0.; break;
          case C_NSO : P[i][j] = 0.; break;
          case C_NSWO: P[i][j] = 0.; break;
          default    : break;
	 }
     }
}

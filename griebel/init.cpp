# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>
#include <stdio.h>
using namespace std;

# include "datadef.hpp"

//****************************************************************************80

int READ_PARAMETER ( char *Inputfile, char *problem,
		   REAL *xlength,REAL *ylength,int *imax,int *jmax,
		   REAL *delx,REAL *dely,
		   REAL *t_end,REAL *delt, REAL *tau,
		   REAL *del_trace,REAL *del_inj,REAL *del_streak,REAL *del_vec,
                   char *vecfile,char *tracefile,char *streakfile,
                   char *infile,char *outfile,
                   int *N,REAL *pos1x,REAL *pos1y,REAL *pos2x, REAL *pos2y,
		   int *itermax,REAL *eps,REAL *omg,REAL *gamma,int *p_bound,
		   REAL *Re,REAL *Pr,REAL *beta,REAL *GX,REAL *GY,
		   REAL *UI,REAL *VI,REAL *TI,
		   int *wW,int *wE,int *wN,int *wS)

//****************************************************************************80
//
//  Purpose:
//
//    READ_PARAMETER reads the input parameters from "Inputfile".
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
  char c;
  FILE *fp;

  if ((fp = fopen(Inputfile, "r")) == NULL)
  {
    cout << "\n";
    cout << "READ_PARAMETER - Fatal error!\n";
    cout << "Error while opening inputfile " << Inputfile << "\n";
    return(1);
  }
//
//  Reading the type of the problem and checking if defined or not.
//
  fscanf(fp, "%s", problem);
//
//  Seek the end of the line.
//
     for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));

  if( strcmp(problem, "convection") && strcmp(problem, "rayleigh") &&
      strcmp(problem, "fluidtrap") &&
      strcmp(problem, "dcavity") && strcmp(problem, "backstep") &&
      strcmp(problem, "plate") && strcmp(problem, "circle") &&
      strcmp(problem, "dam") && strcmp(problem, "drop") &&
      strcmp(problem, "molding") && strcmp(problem, "wave"))
  {
    cout << "Problem " << problem << " not defined!\n";
    cout << "Choose dcavity\n";
    cout << "       convection\n";
    cout << "       rayleigh\n";
    cout << "       fluidtrap\n";
    cout << "       backstep\n";
    cout << "       plate\n";
    cout << "	    circle\n";
    cout << "	    drop\n";
    cout << "       dam\n";
    cout << "       molding\n";
    cout << "       wave\n";
    return(1);
  }
//
//  reading "Inputfile" line for line.
//
  fscanf(fp, "%s", infile);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, "%s", outfile);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));

  fscanf(fp, INREAL, xlength);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, ylength);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, "%d", imax);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, "%d", jmax);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));

  (*delx) = (*xlength)/(*imax); (*dely) = (*ylength)/(*jmax);

  fscanf(fp, INREAL, t_end);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, delt);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, tau);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));

  fscanf(fp, INREAL, del_trace);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, del_inj);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, del_streak);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, del_vec);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));

  fscanf(fp, "%s", vecfile);
     for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, "%s", tracefile);
     for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, "%s", streakfile);
     for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));

  fscanf(fp, "%d", N);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, pos1x);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, pos1y);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, pos2x);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, pos2y);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));

  fscanf(fp, "%d", itermax);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, eps);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, omg);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, gamma);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, "%d", p_bound);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));

  fscanf(fp, INREAL, Re);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, Pr);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, beta);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, GX);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, GY);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, UI);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, VI);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, TI);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));

  fscanf(fp, "%d", wW);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, "%d", wE);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, "%d", wN);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, "%d", wS);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
//
//  Printing problem parameters.
//
  cout << "\n";
  cout << "Problem: " << problem << "\n";

  cout << "\n";
  cout << "xlength=      " << *xlength    << "  "
       << "ylength=      " << *ylength    << "  "
       << "imax=         " << *imax       << "  "
       << "jmax=         " << *jmax       << "\n";

  cout << "delx=         " << *delx       << "  "
       << "dely=         " << *dely       << "\n";

  cout << "t_end=        " << *t_end      << "  "
       << "delt=         " << *delt       << "  "
       << "tau=          " << *tau        << "\n";

  cout << "del_trace=    " << *del_trace  << "  "
       << "del_inj=      " << *del_inj    << "  "
       << "del_streak=   " << *del_streak << "  "
       << "del_vec=      " << *del_vec    << "\n";

  cout << "vecfile=      " << vecfile     << "  "
       << "tracefile=    " << tracefile   << "  "
       << "streakfile=   " << streakfile  << "\n";

  cout << "infile=       " << infile      << "  "
       << "outfile=      " << outfile     << "\n";

  cout << "N=            " << *N          << "  "
       << "pos1x=        " << *pos1x      << "  "
       << "pos1y=        " << *pos1y      << "  "
       << "pos2x=        " << *pos2x      << "  "
       << "pos2y=        " << *pos2y      << "\n";

  cout << "itermax=      " << *itermax    << "  "
       << "eps=          " << *eps        << "  "
       << "omg=          " << *omg        << "  "
       << "gamma=        " << *gamma      << "  "
       << "p_bound=      " << *p_bound    << "\n";

  cout << "Re=           " << *Re         << "  "
       << "Pr=           " << *Pr         << "  "
       << "beta=         " << *beta       << "  "
       << "GX=           " << *GX         << "  "
       << "GY=           " << *GY         << "\n";

  cout << "UI=           " << *UI         << "  "
       << "VI=           " << *VI         << "  "
       << "TI=           " << *TI         << "\n";

  cout << "wW=           " << *wW         << "  "
       << "wE=           " << *wE         << "  "
       << "wN=           " << *wN         << "  "
       << "wS=           " << *wS         << "\n";
//
//  Closing "Inputfile".
//
  fclose(fp);

  if ((*p_bound !=1) && (*p_bound != 2)){
    cout << "p_bound must be 1 or 2!\n";
    return(1);
  }
  if ((*wW > 4)||(*wW < 1)){
    cout << "wW must be 1,2,3, or 4\n";
    return(1);
  }
  if ((*wE > 4)||(*wE < 1)){
    cout << "wE must be 1,2,3, or 4\n";
    return(1);
  }
  if ((*wN > 4)||(*wN < 1)){
    cout << "wN must be 1,2,3, or 4\n";
    return(1);
  }
  if ((*wS > 4)||(*wS < 1)){
    cout << "wS must be 1,2,3, or 4\n";
    return(1);
  }
  if (((*wW==4) && (*wE!=4)) || (*wW!=4) && (*wE==4)){
    cout << "Periodic boundary conditions need wW=wE=4\n";
    return(1);
  }
  if (((*wS==4) && (*wN!=4)) || (*wS!=4) && (*wN==4)){
    cout << "Periodic boundary conditions need wS=wN=4\n";
    return(1);
  }

  return(0);
}
//****************************************************************************80

void INIT_UVP(char *problem,
              REAL **U,REAL **V,REAL **P,REAL **TEMP,int imax,int jmax,
              REAL UI,REAL VI,REAL TI)

//****************************************************************************80
//
//  Purpose:
//
//    INIT_UVP sets the initial values for U,V,P, and TEMP.
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

  for(i=0;i<=imax+1;i++)
    for(j=0;j<=jmax+1;j++)
      {
	U[i][j] = UI;
	V[i][j] = VI;
	P[i][j] = 0.;
	TEMP[i][j] = TI;
      }
//
//  Set U=0.0 in the lower half for the flow past a backward facing step.
//
  if(strcmp(problem, "backstep")==0)
     for(i=0;i<=imax+1;i++)
        for(j=0;j<=jmax/2;j++)
           U[i][j] = 0.0;
}


//****************************************************************************80

void INIT_FLAG(char *problem,int **FLAG,int imax,int jmax,REAL delx,REAL dely,
               int *ibound)

//****************************************************************************80
//
//  Purpose:
//
//    INIT_FLAG initializes the integer array FLAG, dependent of the problem type.
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
  int low,up;
  REAL mx,my,x,y,rad1;
//
//  boundary strip to C_B
//
  for(i=0;i<=imax+1;i++)
    {
     FLAG[i][0]      = C_B;
     FLAG[i][jmax+1] = C_B;
    }
  for(j=1;j<=jmax;j++)
    {
     FLAG[0][j]      = C_B;
     FLAG[imax+1][j] = C_B;
    }
//
//  all inner cells fluid cells.
//
  for(i=1;i<=imax;i++)
     for(j=1;j<=jmax;j++)
        FLAG[i][j] = C_F;
//
//  problem dependent obstacle cells in the interior.
//
  if(strcmp(problem,"fluidtrap")==0)
    {
     for(i=9*imax/22+1;i<=13*imax/22;i++)
       {
	for(j=1;j<=4*jmax/11;j++)
           FLAG[i][j] = C_B;
        for(j=8*jmax/11+1;j<=jmax;j++)
           FLAG[i][j] = C_B;
       }
    }

  if ( strcmp(problem,"plate") == 0 )
  {
//
//  flow past an inclined plate.
//
     low = 2*jmax/5;
//
//  lower and upper bound of the plate
//
     up  = 3*jmax/5;
     FLAG[low][low]     = C_B;
     FLAG[low][low+1]   = C_B;
     FLAG[up][up-1]     = C_B;
     FLAG[up][up]       = C_B;
     for (i=low+1;i<=up-1;i++)
        for (j=i-1;j<=i+1;j++)
           FLAG[i][j] = C_B;
    }
//
//  flow past a backward facing step.
//
  if ( strcmp(problem,"backstep") == 0 )
    {

     for (i=1;i<=jmax;i++)
        for (j=1;j<=jmax/2;j++)
           FLAG[i][j] = C_B;
    }
//
//  flow past a cylinder/circle.
//
  if(strcmp(problem,"circle")==0)
    {

     mx = 20.0/41.0*jmax*dely;
     my = mx;
     rad1 = 5.0/41.0*jmax*dely;
     for (i=1;i<=imax;i++)
        for (j=1;j<=jmax;j++)
          {
           x = (i-0.5)*delx;
           y = (j-0.5)*dely;
           if ((x-mx)*(x-mx)+(y-my)*(y-my) <= rad1*rad1)
              FLAG[i][j] = C_B;
	  }
    }

  if(strcmp(problem,"molding")==0)
    {
     mx = jmax*dely/2;
     my = jmax*dely/2;
     rad1 = jmax*dely/6;
     for (i=1;i<=imax;i++)
        for (j=1;j<=jmax;j++)
          {
           x = (i-0.5)*delx;
           y = (j-0.5)*dely;
           if ((x-mx)*(x-mx)+(y-my)*(y-my) <= rad1*rad1)
              FLAG[i][j] = C_B;
	  }
    }
//
//  flow past a backward facing step.
//
  if ( strcmp(problem,"wave") == 0 )
    {

     for (i=1;i<=jmax;i++)
        for (j=1;j<=jmax/2;j++)
           FLAG[i][j] = C_B;
    }
//
//   Printing the geometry of the fluid domain.
//
  cout << "\n";
  cout << "Geometry of the fluid domain:\n";
  cout << "\n";

  for(j=jmax+1;j>=0;j--)
    {
     for(i=0;i<=imax+1;i++)
        if (!(FLAG[i][j] & C_F))
           cout << "**";
        else
           cout << "  ";
     cout << "\n";
    }
  cout << "\n";
  cout << "\n";
//
//  FLAGs for boundary cells.
//
  (*ibound) = 0;
  for(i=1;i<=imax;i++)
     for(j=1;j<=jmax;j++){
        if (!(FLAG[i][j] & C_F))
           (*ibound)++;
        FLAG[i][j] += ((FLAG[i-1][j] & C_F)*B_W + (FLAG[i+1][j] & C_F)*B_O +
                      (FLAG[i][j-1] & C_F)*B_S + (FLAG[i][j+1] & C_F)*B_N)/C_F;
        switch (FLAG[i][j]){
           case 0x0003:
           case 0x0007:
           case 0x000b:
           case 0x000c:
           case 0x000d:
           case 0x000e:
           case 0x000f:
           {
             cout << "Illegal obstacle cell [" << i << "][" << j << "]\n";
             exit(0);
           }
	 }
      }
}
//****************************************************************************80

void WRITE_bin ( REAL **U, REAL **V, REAL **P, REAL **TEMP, int **FLAG,
		      int imax, int jmax, char* file )

//****************************************************************************80
//
//  Purpose:
//
//    WRITE_BIN writes U,V,P,TEMP,FLAG into a file for subsequent calculations.
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
 FILE *fp;

 fp = fopen(file, "wb");

 fwrite(&imax, sizeof(int), 1, fp);
 fwrite(&jmax, sizeof(int), 1, fp);

 for(i=0;i<=imax+1;i+=1)
   fwrite(U[i], sizeof(REAL), jmax+2, fp);
 for(i=0;i<=imax+1;i+=1)
   fwrite(V[i], sizeof(REAL), jmax+2, fp);
 for(i=0;i<=imax+1;i+=1)
   fwrite(P[i], sizeof(REAL), jmax+2, fp);
 for(i=0;i<=imax+1;i+=1)
   fwrite(TEMP[i], sizeof(REAL), jmax+2, fp);
 for(i=0;i<=imax+1;i+=1)
   fwrite(FLAG[i], sizeof(int), jmax+2, fp);
 fclose(fp);
}
//****************************************************************************80

int READ_bin(REAL **U,REAL **V,REAL **P,REAL **TEMP,int **FLAG,
		 int imax,int jmax,char* file)

//****************************************************************************80
//
//  Purpose:
//
//    READ_BIN reads initial values from a file.
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
  FILE *fp;

  if(strcmp(file, "none") == 0)
  {
    return(-1);
  }

  if( (fp = fopen(file,"rb")) == NULL)
  {
    cout << "\n";
    cout << "READ_BIN - Fatal error!\n";
    cout << "  Cannot open file " << file << "\n";
   return(1);
  }

  fread(&i, sizeof(int), 1, fp);
  fread(&j, sizeof(int), 1, fp);

  if(i!=imax || j!=jmax)
  {
    cout << "\n";
    cout << "READ_BIN - Fatal error!\n";
    cout << "  IMAX or JMAX has the wrong values in " << file << "\n";
    cout << "  IMAX = " << i << "\n";
    cout << "  JMAX = " << j << "\n";
    return(2);
  }

 for(i=0;i<=imax+1;i+=1)
   fread(U[i], sizeof(REAL), jmax+2, fp);
 for(i=0;i<=imax+1;i+=1)
   fread(V[i], sizeof(REAL), jmax+2, fp);
 for(i=0;i<=imax+1;i+=1)
   fread(P[i], sizeof(REAL), jmax+2, fp);
 for(i=0;i<=imax+1;i+=1)
   fread(TEMP[i], sizeof(REAL), jmax+2, fp);
 for(i=0;i<=imax+1;i+=1)
   fread(FLAG[i], sizeof(int), jmax+2, fp);
 fclose(fp);

 return(0);
}

//****************************************************************************80

REAL **RMATRIX(int nrl,int nrh,int ncl,int nch)

//****************************************************************************80
//
//  Purpose:
//
//    RMATRIX allocates memory for a [nrl,nrh]x[ncl,nch]-array of REAL-type.
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
  REAL **m;
  if((m = (REAL**) malloc((unsigned) (nrh-nrl+1)*sizeof(double*))) == NULL)
     {
        cout << "\n";
        cout << "RMATRIX - Fatal error!\n";
        cout << "  No memory.\n";
      exit(0);
     }
  m -= nrl;
  for(i=nrl;i<=nrh;i++)
    {
     if((m[i] = (REAL*) malloc((unsigned) (nch-ncl+1)*sizeof(double)))==NULL)
       {
        cout << "\n";
        cout << "RMATRIX - Fatal error!\n";
        cout << "  No memory.\n";
        exit(0);
       }
     m[i] -= ncl;
    }
  return m;
}
//****************************************************************************80

void FREE_RMATRIX(REAL** m,int nrl,int nrh,int ncl,int nch)

//****************************************************************************80
//
//  Purpose:
//
//    FREE_RMATRIX frees the memory of an array allocated with RMATRIX.
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
  for (i=nrh;i>=nrl;i--) free((void*) (m[i]+ncl));
  free((char*) (m+nrl));
}
//****************************************************************************80

int **IMATRIX(int nrl,int nrh,int ncl,int nch)

//****************************************************************************80
//
//  Purpose:
//
//    IMATRIX allocates memory for a [nrl,nrh]x[ncl,nch]-array of integer-type.
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
  int **m;
  if((m = (int**) malloc((unsigned) (nrh-nrl+1)*sizeof(int*))) == NULL)
     {
     cout << "\n";
     cout << "IMATRIX - Fatal error!\n";
     cout << "  No memory.\n";
      exit(0);
     }
  m -= nrl;
  for(i=nrl;i<=nrh;i++)
    {
     if((m[i] = (int*) malloc((unsigned) (nch-ncl+1)*sizeof(int)))==NULL)
       {
        cout << "\n";
        cout << "IMATRIX - Fatal error!\n";
        cout << "  No memory.\n";
        exit(0);
       }
     m[i] -= ncl;
    }
  return m;
}
//****************************************************************************80

void FREE_IMATRIX ( int** m, int nrl, int nrh, int ncl, int nch )

//****************************************************************************80
//
//  Purpose:
//
//    FREE_IMATRIX frees the memory of an array allocated with IMATRIX.
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

  for ( i=nrh; i>=nrl; i-- )
  {
    free((char*) (m[i]+ncl));
  }

  free((char*) (m+nrl));
}

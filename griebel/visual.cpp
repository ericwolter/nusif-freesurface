# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "datadef.hpp"
# include "visual.hpp"

//****************************************************************************80

void OUTPUTVEC_bin ( REAL **U, REAL **V, REAL **P, REAL **TEMP,
  REAL **PSI, REAL **ZETA, REAL **HEAT, int **FLAG,
  REAL xlength, REAL ylength, int imax, int jmax, char* vecfile )

//****************************************************************************80
//
//  Purpose:
//
//    OUTPUTVEC_bin writes U, V, P, PSI, and ZETA into "vecfile" for visualization.
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
  FILE *fp;
  int i;
  int j;
  float temp;

  fp = fopen ( vecfile, "wb" );

  temp = ( float ) xlength;
  fwrite ( &temp, sizeof(float), 1, fp );
  temp = ( float ) ylength;
  fwrite ( &temp, sizeof(float), 1, fp );
  temp = ( float ) imax;
  fwrite ( &temp, sizeof(float), 1, fp );
  temp = ( float ) jmax;
  fwrite ( &temp, sizeof(float), 1, fp );

  for(j=1;j<=jmax;j+=1)
  {
    for(i=1;i<=imax;i+=1)
    {
      if ( (FLAG[i][j] & C_F) && (FLAG[i][j] < C_E) ) 
      {
        temp = (U[i][j]+U[i-1][j])/2.0;
      }
      else
      {
	temp = 0.0;
      }
      fwrite(&temp, sizeof(float), 1, fp);
    }
  }

  for(j=1;j<=jmax;j+=1)
  {
    for(i=1;i<=imax;i+=1)
    {
      if ( (FLAG[i][j] & C_F) && (FLAG[i][j] < C_E) )
      {
        temp = (V[i][j]+V[i][j-1])/2.0;
      }
      else
      {
        temp = 0.0;
      }
      fwrite(&temp, sizeof(float), 1, fp);
    }
  }

  for(j=1;j<=jmax;j+=1)
  {
    for(i=1;i<=imax;i+=1)
    {
      if( (FLAG[i][j] & C_F) && (FLAG[i][j] < C_E) ) 
      {
        temp = P[i][j];
      }
      else
      {
        temp = 0.0;
      }
      fwrite(&temp, sizeof(float), 1, fp);
    }
  }

  for(j=1;j<=jmax;j+=1)
  {
    for(i=1;i<=imax;i+=1)
    {
      if( (FLAG[i][j] & C_F) && (FLAG[i][j] < C_E) ) 
      {
        temp = TEMP[i][j];
      }
      else
      {
        temp = -0.5;
      }
      temp = TEMP[i][j];
      fwrite(&temp, sizeof(float), 1, fp);
    }
  }

  for(j=1;j<=jmax-1;j+=1)
  {
    for(i=1;i<=imax-1;i+=1)
    {
      temp = ZETA[i][j];
      fwrite(&temp, sizeof(float), 1, fp);
    }
  }

  for(j=0;j<=jmax;j+=1)
  {
    for(i=0;i<=imax;i+=1)
    {
      temp = PSI[i][j];
      fwrite(&temp, sizeof(float), 1, fp);
    }
  }

  for(j=0;j<=jmax;j+=1)
  {
    for(i=0;i<=imax;i+=1)
    {
      temp = HEAT[i][j];
      fwrite(&temp, sizeof(float), 1, fp);
    }
  }

  fclose ( fp );
}
//****************************************************************************80

void COMPPSIZETA(REAL **U,REAL **V,REAL **PSI,REAL **ZETA,int **FLAG,
  int imax,int jmax,REAL delx,REAL dely)

//****************************************************************************80
//
//  Purpose:
//
//    COMPPSIZETA computes the stream function and vorticity
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
//  Computation of the vorticity zeta at the upper right corner
//  of cell (i,j) (only if the corner is surrounded by fluid cells).
//  
 for (i=1;i<=imax-1;i++)
    for (j=1;j<=jmax-1;j++)
        if( ((FLAG[i][j] & C_F) && (FLAG[i][j] < C_E))  &&
            ((FLAG[i+1][j] & C_F) && (FLAG[i+1][j] < C_E))  &&
            ((FLAG[i][j+1] & C_F) && (FLAG[i][j+1] < C_E))  &&
            ((FLAG[i+1][j+1] & C_F) && (FLAG[i+1][j+1] < C_E)) )
          ZETA[i][j] = (U[i][j+1]-U[i][j])/dely - (V[i+1][j]-V[i][j])/delx;
       else
          ZETA[i][j] = 0.0;
//  
//  Computation of the stream function at the upper right corner 
//  of cell (i,j) (only if bother lower cells are fluid cells).
//  
 for (i=0;i<=imax;i++)
   {
    PSI[i][0] = 0.0;
    for(j=1;j<=jmax;j++)
        if( ((FLAG[i][j] & C_F) && (FLAG[i][j] < C_E))  ||
            ((FLAG[i+1][j] & C_F) && (FLAG[i+1][j] < C_E)) )
          PSI[i][j] = PSI[i][j-1] + U[i][j]*dely;
       else
          PSI[i][j] = PSI[i][j-1];
  } 
}
//****************************************************************************80

void COMP_HEAT(REAL **U,REAL **V,REAL **TEMP,REAL **HEAT,int **FLAG,
               REAL Re,REAL Pr,int imax,int jmax,REAL delx,REAL dely)

//****************************************************************************80
//
//  Purpose:
//
//    COMP_HEAT computes the heat function. 
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
//  Computation at the upper right corner of cell (i,j).
//  
 for (i=0;i<=imax;i++)
   {
    HEAT[i][0] = 0.0;
    for(j=1;j<=jmax;j++)
       if( ((FLAG[i][j] & C_F) && (FLAG[i][j] < C_E))  ||
           ((FLAG[i+1][j] & C_F) && (FLAG[i+1][j] < C_E)) )
          HEAT[i][j] = HEAT[i][j-1] + 
                    dely*(U[i][j]*0.5*(1.0+TEMP[i+1][j]+TEMP[i][j])*Re*Pr-
                          (TEMP[i+1][j]-TEMP[i][j])/delx );
   }
}
//****************************************************************************80

struct particle *PARTALLOC ( REAL x, REAL y )

//****************************************************************************80
//
//  Purpose:
//
//    PARTALLOC allocate memory for a particle and set the coordinates to (x,y).
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

  if((part=(struct particle *)malloc(sizeof(struct particle))) == NULL)
  {
    cout << "\n";
    cout << "PARTALLOC - Fatal error!\n";
    cout << "  No memory!\n";	
    exit(0);
  }

  part->x = x; part->y = y;
  part->next = NULL;
  return( part );
}
//****************************************************************************80

struct particleline *SET_PARTICLES(int N,REAL pos1x,REAL pos1y,
				   	 REAL pos2x,REAL pos2y)

//****************************************************************************80
//
//  Purpose:
//
//    SET_PARTICLES allocates memory for a particleline where particles are injected.
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
  REAL hx,hy;
  struct particleline *Partlines;

  if((Partlines=(struct particleline *)
                 malloc((unsigned)(N) * sizeof(struct particleline))) == NULL)
  {	
    cout << "\n";
    cout << "SET_PARTICLES - Fatal error!\n";
    cout << "  No memory!\n";	
    exit(0);
  }

  Partlines -= 1;  

   if (N>=2)
     {
      hx  = (pos2x-pos1x)/(N-1);
      hy  = (pos2y-pos1y)/(N-1);
      for(i=1; i<=N; i++){
         Partlines[i].length = 0;
         Partlines[i].Particles = 
             PARTALLOC(pos1x+hx*(i-1),pos1y+hy*(i-1));
         Partlines[i].Particles->next = 
             PARTALLOC(pos1x+hx*(i-1),pos1y+hy*(i-1));
         Partlines[i].length++;
       }
     }
   return(Partlines);
}

//****************************************************************************80

void ADVANCE_PARTICLES(int imax,int jmax,REAL delx,REAL dely,REAL delt,
                       REAL **U,REAL **V,int **FLAG,
                       int N,struct particleline *Partlines)

//****************************************************************************80
//
//  Purpose:
//
//    ADVANCE_PARTICLES updates the position of particles.
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
  int i, j, k;
  REAL x, y, x1, y1, x2, y2, u, v; 
  struct particle *part,*help;

  for(k=1;k<=N;k++){
     for(part=Partlines[k].Particles; part->next != NULL; part=part->next)
    {
//  
//  first element is only a dummy element.
//  
        x = part->next->x; y = part->next->y;
//  
//  Computation of new x-coordinates by discretizing dx/dt=u.
//  
        i = (int)(x/delx)+1;  j = (int)((y+0.5*dely)/dely)+1;

        x1 = (i-1)*delx;    y1 = ((j-1)-0.5)*dely;
        x2 = i*delx;        y2 = (j-0.5)*dely;
//  
//  bilinear interpolation.
//  
        u = ((x2-x)*(y2-y)*U[i-1][j-1] +
  	     (x-x1)*(y2-y)*U[i][j-1]   +
	     (x2-x)*(y-y1)*U[i-1][j]   +
	     (x-x1)*(y-y1)*U[i][j])/delx/dely;
//  
//  Computation of new y-coordinates by discretizing dy/dt=v.
//  
        i = (int)((x+0.5*delx)/delx)+1; j = (int)(y/dely)+1;

        x1 = ((i-1)-0.5)*delx;    y1 = (j-1)*dely;
        x2 = (i-0.5)*delx;        y2 = j*dely;
//
//  Bilinear interpolation.
//
        v = ((x2-x)*(y2-y)*V[i-1][j-1] +
	     (x-x1)*(y2-y)*V[i][j-1]   +
	     (x2-x)*(y-y1)*V[i-1][j]   +
	     (x-x1)*(y-y1)*V[i][j])/delx/dely;

        x += delt*u;   y += delt*v; 
//
//  Determine new cell for the particle.
//
        i = (int)(x/delx)+1;   j = (int)(y/dely)+1;
//
//  if particle left the fluid domain, delete it.
//
        if (x>=imax*delx || y>=jmax*dely || x<=0 || y<=0){
          help = part->next->next;
          free(part->next);
          part->next = help;
          Partlines[k].length--;
	  if (help == NULL)
	    break; 
        }
        else{
//  
//  special treatment if particle would be in an inner obstacle cell.
//  
          if (FLAG[i][j] < C_F)
             ADVANCE_AT_BOUND(i,j,&x,&y,u,v,U,V,FLAG,delx,dely,delt);
            
          part->next->x = x; part->next->y = y;
         }
      }
   }
}

//****************************************************************************80

void ADVANCE_AT_BOUND(int i,int j,REAL *x,REAL *y,REAL u,REAL v,
                      REAL **U,REAL **V,int **FLAG,
                      REAL delx,REAL dely,REAL delt)

//****************************************************************************80
//
//  Purpose:
//
//    ADVANCE_AT_BOUND updates particle positions near bound.
//
//  Discussion:
//
//    Computation of new particle location of a particle near a no-slip wall,
//    guaranteeing, that the new position is not in the obstacle cell.
//    Here a modified interpolation algorithm is applied, using the fact that
//    at no-skip walls, the velocity is not only given at the midpoint of the
//    edge but on the whole edge.
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
 int ialt,jalt;
 REAL xalt,yalt;
 REAL ul,ur,vo,vu;     
 REAL x1,x2,y1,y2;

//
//  get old particle position.
//
 xalt = (*x)-delt*u;          yalt = (*y)-delt*v;
 ialt = (int)(xalt/delx)+1;   jalt = (int)(yalt/dely)+1; 

  if (i != ialt)
  {
   if (FLAG[ialt+1][jalt] < C_F)
      ur = 0.0;      else{
      if (yalt>= (jalt-0.5)*dely)
         if (FLAG[ialt+1][jalt+1] < C_F){
            y2 = jalt *dely;
            ur = U[ialt][jalt]*(y2-yalt)*2.0/dely;
           }
         else{  
            y1 = (jalt-0.5)*dely;
            y2 = (jalt+0.5)*dely;
            ur = (U[ialt][jalt]*(y2-yalt)+U[ialt][jalt+1]*(yalt-y1))/dely;
	   }
      else      
         if (FLAG[ialt+1][jalt-1] < C_F){
            y1 = (jalt-1.0)*dely;
            ur = U[ialt][jalt]*(yalt-y1)*2.0/dely;
	   }
         else{ 
            y1 = (jalt-1.5)*dely;
            y2 = (jalt-0.5)*dely;
            ur = (U[ialt][jalt-1]*(y2-yalt)+U[ialt][jalt]*(yalt-y1))/dely;  
	   }
     }
   if (FLAG[ialt-1][jalt] < C_F)
      ul = 0.0;  
   else{
      if (yalt>= (jalt-0.5)*dely) 
         if (FLAG[ialt-1][jalt+1] < C_F){
            y2 = jalt *dely;
            ul = U[ialt-1][jalt]*(y2-yalt)*2.0/dely;
           }
         else{   
            y1 = (jalt-0.5)*dely;
            y2 = (jalt+0.5)*dely;
            ul = (U[ialt-1][jalt]*(y2-yalt)+U[ialt-1][jalt+1]*(yalt-y1))/dely;
	   }
      else       
         if (FLAG[ialt-1][jalt-1] < C_F){
            y1 = (jalt-1.0)*dely;
            ul = U[ialt-1][jalt]*(yalt-y1)*2.0/dely;
	   }
         else{ 
            y1 = (jalt-1.5)*dely;
            y2 = (jalt-0.5)*dely;
            ul = (U[ialt-1][jalt-1]*(y2-yalt)+U[ialt-1][jalt]*(yalt-y1))/dely;
	   }
     }
   u = (ul*(ialt*delx-xalt)+ur*(xalt-(ialt-1)*delx))/delx;
   (*x) = xalt+u*delt;
  }

  if (j != jalt)
  {
   if (FLAG[ialt][jalt+1] < C_F)
      vo = 0.0;   
   else{
      if (xalt>= (ialt-0.5)*delx) 
         if (FLAG[ialt+1][jalt+1] < C_F){
            x2 = ialt*delx;
            vo = V[ialt][jalt]*(x2-xalt)*2.0/delx;
           }
         else{  
            x1 = (ialt-0.5)*delx;
            x2 = (ialt+0.5)*delx;
            vo = (V[ialt][jalt]*(x2-xalt)+V[ialt+1][jalt]*(xalt-x1))/delx;
	   }
      else      
         if (FLAG[ialt-1][jalt+1] < C_F){
            x1 = (ialt-1.0)*delx;
            vo = V[ialt][jalt]*(xalt-x1)*2.0/delx;
	   }
         else{   
            x1 = (ialt-1.5)*delx;
            x2 = (ialt-0.5)*delx;
            vo = (V[ialt-1][jalt]*(x2-xalt)+V[ialt][jalt]*(xalt-x1))/delx;
	   }
     }
   if (FLAG[ialt][jalt-1] < C_F)
      vu = 0.0;  
   else{
      if (xalt>= (ialt-0.5)*delx) 
         if (FLAG[ialt+1][jalt-1] < C_F){
            x2 = ialt*delx;
            vu = V[ialt][jalt-1]*(x2-xalt)*2.0/delx;
           }
         else{   
            x1 = (ialt-0.5)*delx;
            x2 = (ialt+0.5)*delx;
            vu = (V[ialt][jalt-1]*(x2-xalt)+V[ialt+1][jalt-1]*(xalt-x1))/delx;
	   }
      else       
         if (FLAG[ialt-1][jalt-1] < C_F){
            x1 = (ialt-1.0)*delx;
            vu = V[ialt][jalt-1]*(xalt-x1)*2.0/delx;
	   }
         else{ 
            x1 = (ialt-1.5)*delx;
            x2 = (ialt-0.5)*delx;
            vu = (V[ialt-1][jalt-1]*(x2-xalt)+V[ialt][jalt-1]*(xalt-x1))/delx;
	   }
     }
   v = (vu*(jalt*dely-yalt)+vo*(yalt-(jalt-1)*dely))/dely;
   (*y) = yalt+v*delt;
  }
}
//****************************************************************************80

void INJECT_PARTICLES(int N, struct particleline *Partlines)

//****************************************************************************80
//
//  Purpose:
//
//    INJECT_PARTICLES injects new particles for streaklines.
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
  struct particle *part;

  for(i=1; i<=N; i++){
    part = PARTALLOC(Partlines[i].Particles->x,Partlines[i].Particles->y);
    part->next = Partlines[i].Particles->next;
    Partlines[i].Particles->next = part;
    Partlines[i].length++;       
  }
}
//****************************************************************************80

void WRITE_PARTICLES(char *partfile, int N, struct particleline *Partlines)

//****************************************************************************80
//
//  Purpose:
//
//    WRITE_PARTICLES appends particles to an ASCII file.
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
  struct particle *part;

  fp = fopen(partfile,"ab");

  for(i=1; i<=N; i++){
    fprintf(fp,"%d\n",Partlines[i].length);
    for(part=Partlines[i].Particles; part->next != NULL; part=part->next)
      fprintf(fp,"%3.3f %3.3f\n", part->next->x, part->next->y);
  }

  fclose(fp);
}
//****************************************************************************80

void WRITE_PARTICLES_bin(char *partfile, int N, struct particleline *Partlines)

//****************************************************************************80
//
//  Purpose:
//
//    WRITE_PARTICLES appends particles to a file in binary format.
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
 float temp, temp2[2];
 struct particle *part;

 fp = fopen(partfile, "ab");
 for(i=1; i<=N; i++){
    temp=Partlines[i].length;
    fwrite(&temp, sizeof(float), 1, fp);
    part=Partlines[i].Particles;
    for(; part->next != NULL; part=part->next){
      temp2[0]=part->next->x;
      temp2[1]=part->next->y;
      fwrite(temp2, sizeof(float), 2, fp);
    }
 }
 fclose(fp);
}
//****************************************************************************80

void PARTICLE_TRACING(char* outputfile,REAL t,int imax,int jmax,
		      REAL delx,REAL dely,REAL delt,
                      REAL **U,REAL **V,int **FLAG,
                      int N, struct particleline *Partlines, int write)

//****************************************************************************80
//
//  Purpose:
//
//    PARTICLE_TRACING moves particles and appends them to a file if wanted.
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
  FILE *fp;

  if(t == 0)
  {
    fp = fopen(outputfile, "wb");
    fprintf(fp,"%d\n%d\n%f\n%f\n%d\n", imax, jmax, delx, dely, N);
    fclose(fp);
    WRITE_PARTICLES(outputfile,N,Partlines);
  }

  ADVANCE_PARTICLES(imax,jmax,delx,dely,delt,U,V,FLAG,N,Partlines);

  if(write & 1)
  {
    WRITE_PARTICLES(outputfile,N,Partlines);
  }    
}
//****************************************************************************80

void STREAKLINES ( char* streakfile,int write,
                 int imax, int jmax, REAL delx, REAL dely,REAL delt,REAL t, 
                 REAL **U,REAL **V,int **FLAG,
                 int N, struct particleline *Partlines)

//****************************************************************************80
//
//  Purpose:
//
//    STREAKLINES manages streakline particles.
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
  FILE *fp;

  if ( t == 0 )
  {
    fp = fopen(streakfile, "wb");
    fprintf(fp, "%d\n", imax);
    fprintf(fp, "%d\n", jmax);
    fprintf(fp, "%g\n", delx);
    fprintf(fp, "%g\n", dely);
    fprintf(fp, "%d\n", N);
    fclose(fp);
    WRITE_PARTICLES_bin(streakfile,N,Partlines);
   }

  ADVANCE_PARTICLES(imax,jmax,delx,dely,delt,U,V,FLAG,N,Partlines);

  if(write & 2)  INJECT_PARTICLES(N,Partlines);
  if(write & 4)  WRITE_PARTICLES_bin(streakfile,N,Partlines);

}

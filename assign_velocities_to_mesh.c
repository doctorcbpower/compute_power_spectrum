#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "header.h"
#include "power.h"

// Inputs:
// NumPart - number of particles on task
// (x,y,z) - positions
// BoxSize - length of periodic box
// Nx - dimension of FFT mesh on task
// Nmesh - dimension of full FFT mesh
//
// Outputs:
// local_mesh - contains values of mesh cells

void assign_velocities_to_mesh(long long NumPart, long long NumPartTot, float *x, float *y, float *z, float *vx, float *vy, float *vz, double BoxSize, int Nx, int Nmesh, fftw_real *local_mesh)
{
    int j,k,l,jj,kk,ll;
    long long i,indx;
	float dx,dy,dz;

    double *vxmesh,*vymesh,*vzmesh;
    double *vx2mesh,*vy2mesh,*vz2mesh;

    double sigma_x,sigma_y,sigma_z;

    float hx0,hy0,hz0;
    
#if defined(CIC) || defined(TSC)
	float hx1,hy1,hz1;
#endif
    
    float fx0,fx1;  // Weight values for particles close to
                    // task boundaries
    
    double xmin=0.0;
    double SlabSize=BoxSize/NTask;

    if(NTask>1) xmin=ThisTask*BoxSize/NTask;
	    
    vxmesh=(double*)malloc(sizeof(double)*Nx*Nmesh*Nmesh);
/*
    vymesh=(double*)malloc(sizeof(double)*Nx*Nmesh*Nmesh);
    vzmesh=(double*)malloc(sizeof(double)*Nx*Nmesh*Nmesh);
    vx2mesh=(double*)malloc(sizeof(double)*Nx*Nmesh*Nmesh);
    vy2mesh=(double*)malloc(sizeof(double)*Nx*Nmesh*Nmesh);
    vz2mesh=(double*)malloc(sizeof(double)*Nx*Nmesh*Nmesh);
 */
    for(i=0;i<NumPart;i++)
	{
        // What cell does particle lie in?
		j=floor(Nx*(x[i]-xmin)/SlabSize);
        if(j<0) j=0;
        if(j>Nx) j=Nx-1;

        k=floor(Nmesh*y[i]/BoxSize);
		
		if(k>=Nmesh) k-=Nmesh;
		if(k<0) k+=Nmesh;

        l=floor(Nmesh*z[i]/BoxSize);

		if(l>=Nmesh) l-=Nmesh;
		if(l<0) l+=Nmesh;

        indx = l+Nmesh*(Nx*k+j);
        
        vxmesh[indx]+=vx[i];
/*
        vx2mesh[indx]+=vx[i]*vx[i];

        vymesh[indx]+=vy[i];
        vy2mesh[indx]+=vy[i]*vy[i];

        vzmesh[indx]+=vz[i];
        vz2mesh[indx]+=vz[i]*vz[i];
 */
    }
     
    double sumx,sumy,sumz;
        
    sumx=0.0;
    sumy=0.0;
    sumz=0.0;

    for(int j=0;j<Nx;j++)
    {
      for(int k=0;k<Nmesh;k++)
        {
         for(int l=0;l<Nmesh;l++)
          {
            indx = l+Nmesh*(Nx*k+j);
              printf("%lld %d\n",indx,Nx*Nmesh*Nmesh-1);
              sumx+=vxmesh[indx];
/*
              sumy+=vymesh[indx];
            sumz+=vzmesh[indx];
  */
//            sigma_x=vx2mesh[indx]-vxmesh[indx]*vxmesh[indx];
 //           sigma_y=vy2mesh[indx]-vymesh[indx]*vymesh[indx];
  //          sigma_z=vz2mesh[indx]-vzmesh[indx]*vzmesh[indx];
//            local_mesh[indx]=sqrt(sigma_x+sigma_y+sigma_z);
          }
        }
    }
        
    
    fprintf(stdout,"Check - summation over mesh : %18.8g %18.8g %18.8g\n",sumx,sumy,sumz);
    fflush(stdout);
    
    free(vxmesh);
    /*
    free(vymesh);
    free(vzmesh);
    free(vx2mesh);
    free(vy2mesh);
    free(vz2mesh);
     */
    return;
        
 }


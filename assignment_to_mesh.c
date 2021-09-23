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

void assign_to_mesh(long long NumPart, long long NumPartTot, float *x, float *y, float *z, double BoxSize, int Nx, int Nmesh, fftw_real *local_mesh)
{
    int j,k,l,jj,kk,ll;
    long long i,indx;
	float dx,dy,dz;
	float hx0,hy0,hz0;
    
#if defined(CIC) || defined(TSC)
	float hx1,hy1,hz1;
#endif
    
    float fx0,fx1;  // Weight values for particles close to
                    // task boundaries
    
    double xmin=0.0;
    double SlabSize=BoxSize/NTask;
    
    if(NTask>1) xmin=ThisTask*BoxSize/NTask;
	
	for(i=0;i<NumPart;i++)
	{
	// What cell does particle lie in?
		j=floor(Nx*(x[i]-xmin)/SlabSize);

        k=floor(Nmesh*y[i]/BoxSize);
		
		if(k>=Nmesh) k-=Nmesh;
		if(k<0) k+=Nmesh;

        l=floor(Nmesh*z[i]/BoxSize);

		if(l>=Nmesh) l-=Nmesh;
		if(l<0) l+=Nmesh;
		
	// Assign NGP weights as default
		hx0=1.0;
		hy0=1.0;
		hz0=1.0;

#if defined(CIC) || defined(TSC)
		dx=x[i]-xmin;

        // Account for cases where there are particles
        // across the periodic boundaries of the box
        if(NTask>1)
        {
            // Particles is on task NTask-1 and is needed
            // by task 0
			if(dx>BoxSize/2)dx-=BoxSize;
			// Particles is on tasks 0 and is needed by task
			// NTask-1
			if(dx<-(BoxSize/2)) dx+=BoxSize;
		}
        
        j=floor(Nx*dx/SlabSize);

		dx=(float)Nx*dx/SlabSize-j-0.5;
        dy=(float)Nmesh*(y[i]/BoxSize)-k-0.5;
		dz=(float)Nmesh*(z[i]/BoxSize)-l-0.5;
		
#ifdef CIC
		hx0-=fabs(dx);
		hy0-=fabs(dy);
		hz0-=fabs(dz);

		hx1=fabs(dx);
		hy1=fabs(dy);
		hz1=fabs(dz);

        if(dx>0)
			jj=j+1;
		else
			jj=j-1;

        fx0=fx1=1.0;
		
        if(NTask==1)
        {
            if(jj>=Nx) jj-=Nx;
            if(jj<0) jj+=Nx;
        } else {
			if(j<0 || j>=Nx)
			{
				fx0=0.;
				if(j<0) j=0;
				if(j>=Nx) j=Nx-1;
			}
			if(jj<0 || jj>=Nx)
			{
				fx1=0.;
				if(jj<0) jj=0;
				if(jj>=Nx) jj=Nx-1;
			}
		}

        hx0*=fx0;
        hx1*=fx1;
        
        // Compute values in y-direction
		if(dy>0)
			kk=k+1;
		else
			kk=k-1;

        if(kk>=Nmesh) kk-=Nmesh;
        if(kk<0) kk+=Nmesh;

        // Compute values in z-direction
		if(dz>0)
			ll=l+1;
		else
			ll=l-1;
 
		if(ll>=Nmesh) ll-=Nmesh;
		if(ll<0) ll+=Nmesh;
#endif
#endif
        indx = l + (2*(Nmesh/2+1))*(k+Nmesh*j);
        local_mesh[indx]+=(fftw_real)(hx0*hy0*hz0);
#ifdef CIC
        // Add to mesh points
        indx = ll + (2*(Nmesh/2+1))*(kk+Nmesh*j);
        local_mesh[indx]+=(fftw_real)(hx0*hy1*hz1);
        
        indx = l + (2*(Nmesh/2+1))*(kk+Nmesh*j);
        local_mesh[indx]+=(fftw_real)(hx0*hy1*hz0);
        
        indx = ll + (2*(Nmesh/2+1))*(k+Nmesh*j);
        local_mesh[indx]+=(fftw_real)(hx0*hy0*hz1);
            
        indx = l + (2*(Nmesh/2+1))*(k+Nmesh*jj);
        local_mesh[indx]+=(fftw_real)(hx1*hy0*hz0);
        
        indx = ll + (2*(Nmesh/2+1))*(k+Nmesh*jj);
        local_mesh[indx]+=(fftw_real)(hx1*hy0*hz1);
        
        indx = l + (2*(Nmesh/2+1))*(kk+Nmesh*jj);
        local_mesh[indx]+=(fftw_real)(hx1*hy1*hz0);
        
        indx = ll + (2*(Nmesh/2+1))*(kk+Nmesh*jj);
        local_mesh[indx]+=(fftw_real)(hx1*hy1*hz1);
#endif
		 
	}
    
    double sum=0.0;
    
    for(int j=0;j<Nx;j++)
    {
        for(int k=0;k<Nmesh;k++)
        {
            for(int l=0;l<Nmesh;l++)
            {
                indx=l+(2*(Nmesh/2+1))*(k+Nmesh*j);
                local_mesh[indx]*=(pow(Nmesh,3.)/(double)NumPartTot);
                local_mesh[indx]-=1.0;
                sum+=local_mesh[indx];

            }
        }
    }
        
    double global_sum;

#ifdef ENABLE_MPI
    MPI_Reduce(&sum,&global_sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
    global_sum=sum;
#endif
    
    if(ThisTask==0)
    {
        fprintf(stdout,"Check - summation over mesh : %18.8g\n",global_sum);
        fflush(stdout);
    }
}


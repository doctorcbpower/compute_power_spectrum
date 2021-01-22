#include <math.h>

#include "power.h"
#include "header.h"

void calculate_spectrum(int Nmesh, int nslab_y, double BoxSize, fftw_complex *fft_of_mesh, char *filename)
{
    double *local_avpower, *local_avks;
    long long *local_nmodes;

    double *global_avpower, *global_avks;
    long long *global_nmodes;

    double dx,dy;
    double kx,ky,kz,k2;
    double krad;
    double fcorrect;
    int i, nx, ny, nz;
    int p;
    int ks;
    FILE *outfile;

    //
    local_avpower = (double*)malloc(sizeof(double) * 3*Nmesh);
    local_avks = (double*)malloc(sizeof(double) * 3*Nmesh);
    local_nmodes = (long long*)malloc(sizeof(long long) *3*Nmesh);
    
    for(i=0;i<3*Nmesh;i++)
    {
        local_nmodes[i]=0;
        local_avpower[i]=local_avks[i]=0.0;
    }
    
    for(ny=0;ny<nslab_y;++ny)
      for(nx=0;nx<Nmesh;++nx)
        for(nz=0;nz<(Nmesh/2+1);++nz)
        {
            i=nz+(Nmesh/2+1)*(nx+Nmesh*ny);
    
            if(nx>Nmesh/2)
                kx=nx-Nmesh;
            else
                kx=nx;

            if(ny+ThisTask*nslab_y>Nmesh/2)
                ky=ny+ThisTask*nslab_y-Nmesh;
            else
                ky=ny+ThisTask*nslab_y;

            if(nz>Nmesh/2)
                kz=nz-Nmesh;
            else
                kz=nz;

            dx=fft_of_mesh[i].re/pow(Nmesh,3.);
            dy=fft_of_mesh[i].im/pow(Nmesh,3.);
            
            krad=sqrt(kx*kx+ky*ky+kz*kz);
            ks=(int)krad;
            
                // Here we need to correct for the mass assignment scheme; need
                //          k = PI * k_(x/y/z)/2 * k_Ny
                // where
                //          k_(x/y/z) = 2*PI*i_(x/y/z)/(N*dx)
                // and
                //          k_Ny = PI/dx
                // which means that
                //          k = PI * 2*PI*i/(N*dx) * dx/PI * 1/2
                // which gives
                //          k = PI * i/N
            
            kx*=M_PI/Nmesh;
            ky*=M_PI/Nmesh;
            kz*=M_PI/Nmesh;
            
                // Correct for interpolation scheme...
            fcorrect = 1.0;
            
#ifdef CORRECT_SMOOTHING
            if(kx!=0.) fcorrect *= sin(kx)/kx;
            if(ky!=0.) fcorrect *= sin(ky)/ky;
            if(kz!=0.) fcorrect *= sin(kz)/kz;
            
#ifndef CIC
            p=1.;
#else
            p=2.;
#endif
            fcorrect = pow(fcorrect,2.*p);
#endif
            local_nmodes[ks]++;
            local_avks[ks]+=krad;
            local_avpower[ks]+=(dx*dx+dy*dy)/fcorrect;
        }
        
    global_avpower = (double*)malloc(sizeof(double) * 3*Nmesh);
    global_avks = (double*)malloc(sizeof(double) * 3*Nmesh);
    global_nmodes = (long long*)malloc(sizeof(long long) *3*Nmesh);
    
    for(i=0;i<3*Nmesh;i++)
    {
        global_nmodes[i]=0;
        global_avpower[i]=global_avks[i]=0.0;
    }
    
    for(i=0;i<Nmesh;i++)
    {
        MPI_Reduce(&local_nmodes[i],&global_nmodes[i],1,MPI_LONG_LONG_INT,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(&local_avks[i],&global_avks[i],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(&local_avpower[i],&global_avpower[i],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    }
    
    if(ThisTask==0)
    {
        outfile=fopen(filename,"w");
        fprintf(outfile,"# BoxSize: %g\n",BoxSize);
        fprintf(outfile,"# Nmesh: %d\n",Nmesh);
        fprintf(outfile,"# %6s %18s %18s %s\n","Number","Log10 k [h/Mpc]","Log10 P(k) [h^3/Mpc^3]","Number of modes");
        for(i=0;i<3*Nmesh;i++)
#ifdef DEBUG
        if(global_nmodes[i]>1)
        {
            fprintf(stdout,"%5d %18.8g %18.8g %10lld\n",i,log10(2*M_PI/BoxSize)+log10(global_avks[i]/global_nmodes[i]),log10(global_avpower[i]/global_nmodes[i])+3*log10(BoxSize),global_nmodes[i]);
            fflush(stdout);
        }
#endif
        for(i=0;i<3*Nmesh;i++)
        if(global_nmodes[i]>1)
            fprintf(outfile,"%5d %18.8g %18.8g %10lld\n",i,log10(2*M_PI/BoxSize)+log10(global_avks[i]/global_nmodes[i]),log10(global_avpower[i]/global_nmodes[i])+3*log10(BoxSize),global_nmodes[i]);
        fclose(outfile);
    }
    
    if(ThisTask==0)
        fprintf(stdout,"\nWritten output to %s...\n",filename);
    
    free(global_nmodes);
    free(global_avks);
    free(global_avpower);
    
    free(local_nmodes);
    free(local_avks);
    free(local_avpower);
}

#include "power.h"

void create_fftplan(int Nmesh, rfftwnd_mpi_plan fftw_forward_plan, int slabstart_x, int nslab_x, int slabstart_y, int nslab_y, int fftsize)
{
    fft_forward_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, Nmesh, Nmesh, Nmesh, FFTW_FORWARD, FFTW_ESTIMATE);
                            
    rfftwnd_mpi_local_sizes(fft_forward_plan, &nslab_x, &slabstart_x, &nslab_y, &slabstart_y, &fftsize);
}

void compute_fft(fftw_real *local_mesh, fftw_complex *fftw_mesh)
{
    workspace = (fftw_real*) malloc(sizeof(fftw_real) * fftsize);
  
    for(i=0;i<fftsize;i++) workspace[i]=0.0;
    
    rfftwnd_mpi(fft_forward_plan,1,local_mesh,workspace,FFTW_TRANSPOSED_ORDER);
    
    fftw_complex *fft_of_mesh;
    
    fft_of_mesh = (fftw_complex *) local_mesh;
    
    rfftwnd_mpi_destroy_plan(fft_forward_plan);
}


#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <hdf5.h>

#include "power.h"
#include "header.h"

#define BUFSIZE 500

int main(int argc, char *argv[])
{
  char file_root[BUFSIZE],filename[BUFSIZE],output_file[BUFSIZE];
  
  int i;
  int isHDF5=0;
  int isDistributed=0;
	
  int Nslab=0, Nmesh=0;
  
  long long NumPart=0, NumPartRead=0, NumPartOnTask=0;
  struct sim_info header;
  
  float *x,*y,*z;
  int *ptype;
    
  ThisTask=0;
  NTask=1;

#ifdef ENABLE_MPI
  // First initialise MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD,&NTask);
#endif
    
  // Read information from the command line...
    if(argc<2)
    {
        if(ThisTask==0)
        {
            fprintf(stdout,"Usage: %s -input <filename> -output <filename> -meshsize NNN (<isHDF5>)\n",argv[0]);
            fflush(stdout);
        }
#ifdef ENABLE_MPI
        MPI_Finalize();
#endif
        exit(0);
    }
    
    i=1;
    
    while(i<argc)
    {
        if(strcmp(argv[i],"-input")==0) sprintf(file_root,"%s",argv[++i]);
        if(strcmp(argv[i],"-output")==0) sprintf(output_file,"%s",argv[++i]);
        if(strcmp(argv[i],"-meshsize")==0) Nmesh=atoi(argv[++i]);
        if(strcmp(argv[i],"-isHDF5")==0) isHDF5=1;
        
        i++;
    }
    
    // Assume slab decomposition along the x-dimension
    Nslab=1;
#ifdef ENABLE_MPI
    Nslab*=NTask;
#endif
        
    // Get filename and figure out if it's distributed across multiple files
    check_input_filenames(filename,file_root,isHDF5,&isDistributed);
    
    if(ThisTask==0)
    {
        fprintf(stdout,"Reading header...\n");
        fflush(stdout);
    }
    
    if(isHDF5==1)
        read_hdf5_header(filename, &header, &NumPart);
    else
        read_gadget_binary_header(filename, &header, &NumPart);
    
    if(ThisTask==0)
    {
        fprintf(stdout,"Number of files: %d\n",header.NumFiles);
        fprintf(stdout,"Number of particles: %lld\n",NumPart);
        fprintf(stdout,"Time/Expansion Factor: %g\n",header.time);
        fprintf(stdout,"BoxSize: %g\n",header.BoxSize);
        fflush(stdout);
    }
    
        // Allocate memory....
    
    NumPartOnTask=NumPart;
    
    if(NTask>1)
    {
        NumPartOnTask/=NTask;
        NumPartOnTask*=2;
    }
    
    x = (float*)malloc(sizeof(float)*NumPartOnTask);
    y = (float*)malloc(sizeof(float)*NumPartOnTask);
    z = (float*)malloc(sizeof(float)*NumPartOnTask);
    ptype = (int*)malloc(sizeof(int)*NumPartOnTask);
    
    if(ThisTask==0)
    {
        fprintf(stdout,"Assigning %lld particles per processor...\n",NumPartOnTask);
        fflush(stdout);
    }
    
    if(ThisTask==0)
    {
        fprintf(stdout,"Reading particles...\n");
        fflush(stdout);
    }
    
    if(isHDF5==1)
        read_particles_from_hdf5(file_root,x,y,z,ptype,header.NumFiles,&NumPartRead);
    else
        read_particles_from_gadget_binary(file_root,x,y,z,ptype,header.NumFiles,&NumPartRead);
    
    if(NTask>1)
    {
        if(ThisTask==0)
        {
            fprintf(stdout,"Splitting particles across tasks...\n");
            fflush(stdout);
        }
        
        split_across_processors_by_slabs(&NumPartOnTask,NumPartRead,x,y,z,header.BoxSize, Nmesh, Nslab);
    }
    
        // Start setting up the FFT
    
    fft_forward_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, Nmesh, Nmesh, Nmesh, FFTW_FORWARD, FFTW_ESTIMATE);
    
    rfftwnd_mpi_local_sizes(fft_forward_plan, &nslab_x, &slabstart_x, &nslab_y, &slabstart_y, &fftsize);
    
    local_mesh = (fftw_real *) malloc(sizeof(fftw_real) * fftsize);
    
    for(i=0;i<fftsize;i++) local_mesh[i]=0.0;
    
        // Assign particles to mesh
    
    if(ThisTask==0)
    {
        fprintf(stdout,"\nAssigning particles to mesh...\n");
        fflush(stdout);
    }
    
    assign_to_mesh(NumPartOnTask,NumPart,x,y,z,header.BoxSize,nslab_x,Nmesh,local_mesh);
    
    if(ThisTask==0)
    {
        fprintf(stdout,"\nCalculating FFT...\n");
        fflush(stdout);
    }
        
    workspace = (fftw_real*) malloc(sizeof(fftw_real) * fftsize);
    
    for(i=0;i<fftsize;i++) workspace[i]=0.0;
    
    rfftwnd_mpi(fft_forward_plan,1,local_mesh,workspace,FFTW_TRANSPOSED_ORDER);
    
    fftw_complex *fft_of_mesh;
    
    fft_of_mesh = (fftw_complex *) local_mesh;
    
    rfftwnd_mpi_destroy_plan(fft_forward_plan);
        
    calculate_spectrum(Nmesh, nslab_y, header.BoxSize, fft_of_mesh, output_file);
    
    if(ThisTask==0)
    {
        fprintf(stdout,"\nFinished...\n");
        fflush(stdout);
    }
    
    MPI_Finalize();
    
    return(0);
}



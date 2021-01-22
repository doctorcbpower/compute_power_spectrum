#ifdef DOUBLE_FFTW
#include <drfftw_mpi.h>
#else
#include <srfftw_mpi.h>
#endif

static rfftwnd_mpi_plan fft_forward_plan, fft_inverse_plan;
static int fftsize, maxfftsize;
static fftw_complex *fft_of_rhogrid;

static int slabstart_x, nslab_x, slabstart_y, nslab_y, smallest_slab;

fftw_real *workspace;
fftw_real *local_mesh;

int nx,ny,nz;

void assign_to_mesh(long long, long long, float *, float *, float *, double, int, int, fftw_real *);
void calculate_spectrum(int, int, double, fftw_complex *, char *);

/*
#define rhocrit 27.755

#define BUFSIZE 500

hsize_t dims[2],start[2],count[2];
hid_t hdf5_file, hdf5_headergrp, hdf5_grp, hdf5_attribute, hdf5_datatype, hdf5_dataset, hdf5_dataspace_in_file, hdf5_dataspace_in_memory;

struct sim_info header;

int j,k,l,m,n,buf[BUFSIZE];
MPI_Status status;
char file_root[128],filename[128],image_file[128],debug_file[128];
FILE *infile,*outfile;
char buffer[BUFSIZE];
double SliceSize,ImageSliceSize,f_Image=0.1;
double mp;

int nsend,nrec;

float *send;
float *rec;

char label[4];
char junk[256-24];
long long NumPart=0, NumPartInFile=0;
long long fileoffset;

float pos[3];
float *x,*y,*z;
int *ptype;
float *posx,*posy,*posz;
float xmin,ymin,zmin;
float xmax,ymax,zmax;
double sum,sumx,sumy,sumz;
long long NumPartRead=0;

#ifdef VELOCITIES
float vel[3];
float *vx,*vy,*vz;
float *velx,*vely,*velz;
float vxmin,vymin,vzmin;
float vxmax,vymax,vzmax;
double sumvx,sumvy,sumvz;
#endif

int PMGRID;

int NSlab;
int Nmesh;
int nx,ny,nz;
int ix,iy,iz;
long long num_x[BUFSIZE]={0},num_gather[BUFSIZE]={0};
long long num_proc=0, num_tot=0;

double dx,dy,dz;

int rank;
float *dummy;

float hx0,hy0,hz0,hx1,hy1,hz1;
*/


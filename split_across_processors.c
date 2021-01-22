#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

#include "header.h"

void split_across_processors_by_slabs(long long *NumPartOnTask,long long NumPart,float *x, float *y, float *z, double lbox, int Nmesh, int Nslab)
{
	int i,j;
	int nleft[Nslab],nright[Nslab];
	int ix,iix,left_nx,right_nx;
#ifdef CIC
	int keep_nx=0;
#elif TSC
	int keep_nx=1;
#endif
	float dx, ilbox;
  
	int **ProcessorID;
	int ProcessorCount[Nslab];
  
	ilbox=1.0/lbox;
	
  // Allocate memory for 2D array hosting particle IDs
	ProcessorID = (int**)malloc(sizeof(int*)*Nslab);
  
	for(i=0;i<Nslab;i++)
	{
      ProcessorID[i] = (int*)malloc(2*sizeof(int)*NumPart);
    }
  
  // Initialise the count arrays
	for(i=0;i<Nslab;i++)
    {
      ProcessorCount[i]=0;
      nleft[i]=0;
      nright[i]=0;
    }
  
	// Loop over particles and assign to appropriate task
	for(i=0;i<NumPart;i++)
    {
      // Figure out which slab the particle is on
      ix=(int)Nslab*(x[i]*ilbox);

#if defined(CIC) || defined(TSC)
	  if(ix<0) ix+=Nslab;
	  if(ix>=Nslab) ix-=Nslab;
#else
	  if(ix<0) ix=0;
	  if(ix>=Nslab) ix=Nslab-1;
#endif
      // Now assign particle ID to processor...
      ProcessorID[ix][ProcessorCount[ix]]=i;
      // ... and increment the number of particles on the processor...
      ProcessorCount[ix]+=1;
      
#if defined(CIC) || defined(TSC)
	  // For CIC or TSC, need to know if it is needed by a neigbouring
	  // slab - so compute offset with respect to slab boundary
      dx=x[i]*ilbox-(float)ix/Nslab;
      
	  // Where does it lie relative to left boundary...
      left_nx=Nmesh*dx;    // Cell with respect to
      // ... and right boundary
      right_nx=-Nmesh*(dx-1.0/(float)Nslab);
	  
      if(left_nx<=keep_nx) {
		iix=ix-1;
		if(iix<0) iix+=Nslab;
	// Now assign particle ID to processor...
		ProcessorID[iix][ProcessorCount[iix]]=i;
	// ... and increment the number of particles on the processor...
		ProcessorCount[iix]+=1;
		nleft[ix]++;
      }
      
      if(right_nx<=keep_nx) {
		iix=ix+1;
		if(iix>=Nslab) iix-=Nslab;
	// Now assign particle ID to processor...
		ProcessorID[iix][ProcessorCount[iix]]=i;
	// ... and increment the number of particles on the processor...
		ProcessorCount[iix]+=1;
		nright[ix]++;
      }
#endif
    }
  
#ifdef DEBUG
    for(i=0;i<Nslab;i++)
    {
  	  if(nright[i]>0 || nleft[i]>0)
	    {
	    printf("ThisTask: %d PE: %d nleft: %d nright: %d Count: %d\n",
	    ThisTask,i,nleft[i],nright[i],ProcessorCount[i]);
	    }
    }
#endif
	
  // Move particles to the appropriate task
  
  int ProcessorCount_Received[Nslab];
  
  for(i=0;i<Nslab;i++)
    {
        ProcessorCount_Received[i]=0;
        if(i!=ThisTask)
        {
            MPI_Send(&ProcessorCount[i],1,MPI_INT,i,ThisTask,MPI_COMM_WORLD);
            MPI_Recv(&ProcessorCount_Received[i],1,MPI_INT,i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    }
  
  // Determine total number of particles on the task
  int ProcessorCount_Total=0;
  
  for(i=0;i<Nslab;i++)  // Number of particles that will be moved to task
    {
      if(i!=ThisTask) ProcessorCount_Total+=ProcessorCount_Received[i];
    }
  
  ProcessorCount_Total+=ProcessorCount[ThisTask];  // Add particles already on task
    	
  for(i=0;i<Nslab;i++)
    {
      if(i!=ThisTask)
		printf("ThisTask: %d PE: %d Sending: %d Receiving: %d\n",
	       ThisTask,i,ProcessorCount[i],ProcessorCount_Received[i]);
    }
  
  fprintf(stdout,"Allocating memory for %d particles on PE %d\n\n",ProcessorCount_Total,ThisTask);
  fflush(stdout);
  
  float *posx, *posy, *posz;
  
  posx = (float*)malloc(sizeof(float)*ProcessorCount_Total);
  posy = (float*)malloc(sizeof(float)*ProcessorCount_Total);
  posz = (float*)malloc(sizeof(float)*ProcessorCount_Total);
    
  for(i=0;i<ProcessorCount[ThisTask];i++)
    {
      posx[i]=x[ProcessorID[ThisTask][i]];
      posy[i]=y[ProcessorID[ThisTask][i]];
      posz[i]=z[ProcessorID[ThisTask][i]];
    }
	
  int ProcessorCount_Offset[Nslab];
	
  //    ProcessorCount_Offset[ThisTask]=ProcessorCount[ThisTask];
  for(i=0;i<NTask;i++)
	{
	ProcessorCount_Offset[i]=ProcessorCount_Received[i];
	}
    
  int nsend, nrec;
  float *send, *rec;
  int noffset;
  
  for(i=0;i<NTask;i++)
    {
      if(i!=ThisTask)
	{
	  nsend=ProcessorCount[i];
	  nrec=ProcessorCount_Received[i];
	  
	  send = (float*)malloc(sizeof(float)*nsend);
	  rec = (float*)malloc(sizeof(float)*nrec);
	  
	  noffset=ProcessorCount[ThisTask];
	  
	  for(j=0;j<i;j++)
	    noffset+=ProcessorCount_Offset[j];
	  
	  // First do x
	  for(j=0;j<nsend;j++)
	    send[j]=x[ProcessorID[i][j]];
	  
	  MPI_Sendrecv(send,nsend,MPI_FLOAT,i,ThisTask,
		       rec,nrec,MPI_FLOAT,i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  
	  for(j=0;j<nrec;j++)
	    posx[j+noffset]=rec[j];
	  
	  // ... then do y
	  for(j=0;j<nsend;j++)
	    send[j]=y[ProcessorID[i][j]];
	  
	  MPI_Sendrecv(send,nsend,MPI_FLOAT,i,ThisTask,
		       rec,nrec,MPI_FLOAT,i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  
	  for(j=0;j<nrec;j++)
	    posy[j+noffset]=rec[j];
	  
	  // ... then do z
	  for(j=0;j<nsend;j++)
	    send[j]=z[ProcessorID[i][j]];
	  
	  MPI_Sendrecv(send,nsend,MPI_FLOAT,i,ThisTask,
		       rec,nrec,MPI_FLOAT,i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  
	  for(j=0;j<nrec;j++)
	    posz[j+noffset]=rec[j];
	  
	  free(send);
	  free(rec);
	}
    }
	  
  for(i=0;i<ProcessorCount_Total;i++)
    {
      x[i]=posx[i];
      y[i]=posy[i];
      z[i]=posz[i];
    }
  
  *NumPartOnTask=ProcessorCount_Total;

  free(posx);
  free(posy);
  free(posz);

  free(ProcessorID);

  double xc=0.0,yc=0.0,zc=0.0;
  double xmin=1.e10,xmax=0.0,ymin=1.e10,ymax=0.0,zmin=1.e10,zmax=0.0;
	
  for(i=0;i<ProcessorCount_Total;i++)
	{
		xc+=x[i];
		yc+=y[i];
        zc+=z[i];
        
        if(x[i]<xmin) xmin=x[i];
        if(x[i]>xmax) xmax=x[i];
        if(y[i]<ymin) ymin=y[i];
        if(y[i]>ymax) ymax=y[i];
        if(z[i]<zmin) zmin=z[i];
        if(z[i]>zmax) zmax=z[i];
    }

	xc/=ProcessorCount_Total;
  yc/=ProcessorCount_Total;
  zc/=ProcessorCount_Total;
    
  fprintf(stdout,"Number of particles on Task %d: %d \n",ThisTask,ProcessorCount_Total);
  fprintf(stdout,"Average position on Task %d: (%f,%f,%f) \n",ThisTask,xc,yc,zc);
  fprintf(stdout,"Boundaries on Task %d: (%g,%g) \n",ThisTask,xmin,xmax);
  fflush(stdout);        
}

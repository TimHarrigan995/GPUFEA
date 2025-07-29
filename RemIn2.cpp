//============================================================================
// Name        : RemIn2.cpp
// Author      : Tim
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
#include "stdio.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <malloc.h>
#include <errno.h>

#include <iostream>
#include <vector>
#include <array>
#include <bits/stdc++.h>
#include <numeric> // For std::accumulate (optional, for calculating total_size)
#include <fstream>
#include <string>
#include <cuda_runtime.h>
#include "cudss.h"

#define EXIT_FREE \
    do { \
        free(csr_offsets_h); \
        free(csr_columns_h); \
        free(csr_values_h); \
        free(x_values_h); \
        free(b_values_h); \
        free(diag_h);       \
        free(row_scale_h); \
        free(col_scale_h); \
        cudaFree(csr_offsets_d); \
        cudaFree(csr_columns_d); \
        cudaFree(csr_values_d); \
        cudaFree(x_values_d); \
        cudaFree(b_values_d); \
        cudaFree(diag_d); \
        cudaFree(row_scale_d); \
        cudaFree(col_scale_d); \
    } while(0);


using namespace std;

static FILE *output,*ebfile;
static FILE *input,*disdens,*loadfile;
double *matptr,*gamptr;
double modexpnt,transexpnt;
int nmatp,mtype,matlsize;
int *eldat;
double *volume;
double* loadptr;


void elementDOFconnectivity(int elgdofs[], int elsiz,vector <int>* globaldofs,int aSize);
int twodsolid(int nel,int nnodes,int itype,int nint,double matprops[],double xx[][2],double s[],FILE *output,
		double *volume,double gammar,double gexpnt);
int threedsolid(int nel,int nnodes,int nint,double matprops[3],double xx[][3],double s[],FILE *output,double *volume,double gammar,double gexpnt);

//    threedsolid(int, int, int, double*, double (*) [3], double*, _IO_FILE*, double*, double, double)'



/* **************************************************************************


     INPUT PHASE:

   get input data:
       number of nodes,
       DOF per node,
       number of load cases,
       number of element groups

   **********************************************************************   */
int main() {

//	int **elarry,*eldat;

double *stifptr, *workptr, *matptr, *nsptr;
double modexpnt,transexpnt,nulstate,dtime,gamnorm,ratenorm,dgamnorm,tol,itol;
double oblnorm;
double *gamptr,*dgamptr,*oldgamptr,*rateptr,*relvol,*ebmatrow,*s;
double x,y,z,lmag,e,nu,ux,uy,uz,thic,ddmy;
int maxmelem=0;
int netmaxdel=0;
int itermax=20;
int iter=0;
int blalcdone=0;
int c,c2,c3,ntnodes,masterdof,nlcase,nelgps;
int ndoftot,n,dofx,dofy,dofz,dmy,eldmy,eldof,netdof;
int rind,numel,numrel,eltyp,offset,offset2;
int maxnodes,remflag,maxdel,elsize,memelem;
char title[80];
std::vector<int>* elarry;
//std::vector<int> eldata;


std::vector<int>* gldofs;
std::vector<int> StrungDofs;
std::vector<double>* glStif;
std::vector<double> StrungStif;

int *remdiags,*dofptr;
double xx2d[4][2],xx3d[8][3];
double *coorptr;
double gammar,gexpnt;
int eNodDF[50];

/* open files */

auto input=fopen("/home/tim/FEINFILE", "r");
auto output=fopen("/home/tim/FEOUTFILE", "w");
auto disdens=fopen("FEDISDENS", "w");
auto ebfile=fopen("EBMAT", "w+");
//loadfile=fopen("FELOADS", "w+");

perror(" file open error status (initial):  ");

  /* read in overal control data */
c=0;
while((title[c++]=fgetc(input))!= '\n');

fscanf(input,"%d %d %d %d",&ntnodes,&masterdof,&nlcase,&nelgps);



/* memory for element groups */
//elarry is an array of pointers to element group data arrays (2 per element group)

//int* elarry[nelgps*2];

elarry = new std::vector<int>[nelgps*2];

//elarry = calloc(nelgps*2,sizeof(int *));

/* print up a nice looking header for the output file */

fprintf(output," ********************************************************\n");
fprintf(output,"\n\n      A LINEAR ELASTIC FINITE ELEMENT PROGRAM\n");
fprintf(output,"                FOR BONE REMODELLING STUDIES \n");
fprintf(output,"     AT TRUMAN MEDICAL CENTER/UMKC MED SCHOOL, KC MO\n");
fprintf(output,"     WRITTEN BY TIM HARRIGAN - FIRST REVISION 6/6/91\n\n");

/* *************************************************************************
   get memory for nodal point coordinates
              and degrees of freedom
   **********************************************************************    */

coorptr = new double[ntnodes*3];

dofptr = new int[ntnodes*3];

/* **********************************************************************
   read and echo in nodal points and initial dofs
   **********************************************************************    */

fprintf(output," Nodal Point Degrees of Freedom and Coordinates\n\n");
fprintf(output,
   "  n  xdof  ydof  zdof         x             y            z\n\n");

ndoftot=1;
for(c=0;c<ntnodes;c++){  //each node
  fscanf(input, " %d %d %d %d %lf %lf %lf",&n,&dofx,&dofy,&dofz,&x,&y,&z);
    dmy=3*(n-1);
  if(n<=ntnodes){
    coorptr[dmy]=x;
    coorptr[dmy+1]=y;
    coorptr[dmy+2]=z;
    dofptr[dmy]=dofx;
    dofptr[dmy+1]=dofy;
    dofptr[dmy+2]=dofz;
    if(dofx>0)ndoftot++;
    if(dofy>0)ndoftot++;
    if(dofz>0)ndoftot++;
  }  //if(n<=ntnodes

fprintf(output," %5d %5d %5d %5d %12.9lf %12.9lf,%12.9lf\n",
   n,dofptr[dmy],dofptr[dmy+1],dofptr[dmy+2],x,y,z);
}  //for (c=0;
ndoftot=ndoftot-1;
loadptr= new double[ndoftot*nlcase];
//allocate vectors for stiffness matrices

gldofs = new std::vector<int>[ndoftot+1];
glStif = new std::vector<double>[ndoftot+1];

int arrSize=ndoftot;
int diagindex[ndoftot+1];

/*
   LOOP THROUGH ELEMENT GROUPS

   read in element information
      to get the matrix profile
      and remodelling element indices
      then store the information on disk

   ***************************************************************    */

numrel = 0;  //number of remodeling elements
rind=0; //remodeling element index

// loop through element groups
for(c=0;c<nelgps;c++){  //each element group

  fscanf(input," %d %d %d %d", &numel,&eltyp,&maxnodes,&remflag);

  if(remflag == 1)numrel+= numel; //all elements in a group remodel

  elarry[c*2].push_back(numel);
  elarry[c*2].push_back(eltyp);
  elarry[c*2].push_back(maxnodes);
  elarry[c*2].push_back(remflag);

  printf("%d %d %d %d\n",elarry[0][0],elarry[0][1],elarry[0][2],elarry[0][3]);
/* define netmaxdel = maximum element degrees of freedom
   for later use in remodeling iterations */

  if(eltyp/10 == 2){  //planar elements
      maxdel=maxnodes*2;}
  if(eltyp == 3){ // volume elements
      maxdel=maxnodes*3;}

  netmaxdel = (maxdel >netmaxdel) ? maxdel : netmaxdel;

/* *****************************************************************
   get memory for element info
   ******************************************************************   */

  if(eltyp == 3)elsize= 4 + maxnodes*4;  //3d elements
  if(eltyp/10 == 2)elsize= 4 + maxnodes*3;  //2d elements

  memelem = numel*elsize;

  //allocating memory for element data

  elarry[c*2+1].reserve(memelem);

  eldat = (int *)calloc(memelem, sizeof(int));

//fill up eldat and then copy it over to elarray
/* ********************************************************************
     read in elements,
       send info to column height routine,
       store elements for later
    *******************************************************************  */
  for(c2=0;c2<numel;c2++){  //for each element in the group

  offset=c2*elsize;
  //offset is an offset within eldat - eldat points to an array of individual element data
  /* ********************************************************************
     for each element
     eldat[offset]=element number
     eldat[offset+1]=material property set number
     eldat[offset+2]=number of nodes
     eldat[offset+3]=element density index (-1 if no remodelling)
     ********************************************************************   */

  fscanf(input," %d %d %d", &eldat[offset], &eldat[offset+1],
          &eldat[offset+2]);

  eldat[offset+3] = (remflag == 1)? rind++ : -1; // stores remodeling data for an individual element
  //set up to echo element data
  fprintf(output,"\nElement Number %d degrees of freedom\n", c2+1);
  fprintf(output,"local node, global node, xdof, ydof, zdof\n\n");

  offset2=offset+maxnodes+4; //offset2 is the offset to local dofs in element data

  //assigning global degrees of freedom to element dof arrays

  eldmy = 0;
  netdof=0;
  for(c3=0;c3<eldat[offset+2];c3++){  //for each node in the element

      fscanf(input," %d", &eldat[offset+c3+4]);// get the next element node
      dmy=3*(eldat[offset+c3+4]-1); /* index to element nodal point dof in node array*/

      if(eltyp/10 != 2)eldat[offset2+eldmy++]=dofptr[dmy];
      eldat[offset2+eldmy++]=dofptr[dmy+1];
      eldat[offset2+eldmy++]=dofptr[dmy+2];

    fprintf(output," %5d %5d %5d %5d %5d\n",c3, eldat[offset+c3+4],dofptr[dmy],dofptr[dmy+1],dofptr[dmy+2]);

    }  // each element node

  int elemdfnum=0;
  if(eltyp/10 ==2)elemdfnum=eldat[offset+2]*2;
  if(eltyp == 3) elemdfnum=eldat[offset+2]*3;

  //set up the column dof vector arrrays

  elementDOFconnectivity(&eldat[offset2],elemdfnum,gldofs,8);

  } /* each element in  group */

  // add eldat into the array in elarry
  for (int i=0;i<memelem;i++){
  elarry[c*2+1].push_back(eldat[i]);}

} /* loop through element groups */

/* *******************************************
   Read in number of material property sets.
      If the element is to be remodelled, the elastic modulus
      in its property set must be the modulus for volumetric
      density equal to one.
   *******************************************    */
fscanf(input," %d", &nmatp);  // number of material property sets


fprintf(output,"\n\n Material Property Sets:\n\n");

matlsize=4*nmatp;

matptr=new double[matlsize]; //array for overall material property data

/* **************************************
   read in and echo material property information
   **************************************    */
int pn=0;

for(c=0;c<nmatp;c++){
fscanf(input," %d", &mtype);
if(mtype == 2)
	{fscanf(input," %d %lf %lf %lf",&pn, &e, &nu, &thic);
	fprintf(output,"\n Set number= %d (%dD) E = %lf Nu =%lf Thickness=%lf\n",pn,mtype,e,nu,thic);
	}
if(mtype == 3)
	{fscanf(input," %d %lf %lf",&pn, &e, &nu);
	fprintf(output,"\n Set number= %d (%dD) E = %lf Nu =%lf\n",pn,mtype,e,nu);
	}

dmy=4*c;
matptr[dmy]=pn;
matptr[dmy+1]=e;
matptr[dmy+2]=nu;
if(mtype == 2)matptr[dmy+3]=thic;
} //for each material

/* ***************************************
   get memory for volumetric density and rate info
   *************************************** */
    gamptr=new double[numrel];
//   gamptr=calloc(numrel, sizeof(double)); /* actual density in iterations */
//   if(gamptr == NULL)printf("density allocation failed (main)\n");
//rateptr=new double[numrel];
//   rateptr=calloc(numrel, sizeof(double)); /* rate of change */
//   if(rateptr == NULL)printf("rate allocation failed (main)\n");
//dgamptr=new double[numrel];
//   dgamptr=calloc(numrel, sizeof(double));  /* increment in density */
//   if(dgamptr == NULL)printf("density increment allocation failed (main)\n");
//oldgamptr=new double[numrel];
//   oldgamptr=calloc(numrel, sizeof(double));
		/* converged density from previous iteration */
//   if(oldgamptr == NULL)printf("old density allocation failed (main)\n");
relvol=new double[numrel];
//   relvol=calloc(numrel,sizeof(double)); /* element volumes */
//   if(relvol == NULL)printf("element volume allocation failed (main)\n");
remdiags=new int[numrel+1]; /* addresses of remodelling matrix diagonals */
//   remdiags=calloc(numrel+1,sizeof(int));
//   if(remdiags == NULL)printf("remodelling diagonal index allocation failed (main)\n");
volume = new double[100];  //check - total elements?
/* set up remodelling matrix as full */

for(c=0;c<=numrel;c++)remdiags[c]=c*(c+1)/2;

/* print out number of remodeling elements */

fprintf(output,"number of remodeling elements = %d\n",numrel);

/* ************************************
   read in initial density distribution
   ************************************  */

  for(c=0;c<numrel;c++)
        fscanf(input," %lf",&gamptr[c]);

//   for(c=0;c<numrel;c++)oldgamptr[c]=gamptr[c];

/* *****************************************
   read in loading information:
       for each case,
       create the load vector and store it.
   ***************************************** */


fprintf(output,"\n\nLoad Cases:\n");
for(c=0;c<nlcase;c++){

int ncase,nloads;


fscanf(input," %d %d", &ncase, &nloads);
fprintf(output,"\n load case %d has %d nodal loads\n", ncase, nloads);
fprintf(output,"\n node    direction    value\n");

int lnode,ldof;
double lmag;
for(c2=0;c2<nloads;c2++){
    fscanf(input," %d %d %lf", &lnode,&ldof,&lmag);
    netdof=dofptr[3*(lnode-1)+(ldof-1)];
    loadptr[ndoftot*c+netdof-1]= lmag; //subtract 1 from netdof to get to row/column in equations
    printf(" %5d %5d %5d %lf %5d %5d\n", lnode,ldof,netdof,lmag, ndoftot,c);
    fprintf(output," %5d %5d %lf\n", lnode,ldof,lmag);
  }
}

/* ************************************
   read in remodelling simulation data
   ************************************ */
fscanf(input," %lf %lf %lf %lf %lf %lf"
   ,&modexpnt,&transexpnt,&nulstate,&dtime,&tol,&itol);

/* ****************************************************
   calculate addresses of elements in finite element matrix
   ************************************************  */


// concatenate row data to get the overall matrix index array for gpu routine

size_t total_size = 0;
for (int c=0;c<(ndoftot+1);c++) { //each dof
	diagindex[c]=total_size;
    total_size += gldofs[c].size();
printf("Index array %d\n",c);
    //print out dofs to check/debug
std::cout << "Vector elements: ";
for (int num : gldofs[c]) {
    std::cout << num << " ";
	fprintf(output,"%d ",num);
}
std::cout << std::endl;
fprintf(output,"\n");

} //each row


StrungDofs.reserve(total_size); // Pre-allocate memory for efficiency
for (int c=0;c<(ndoftot+1);c++) {
    StrungDofs.insert(StrungDofs.end(), gldofs[c].begin(), gldofs[c].end());}
std::cout << "Vector elements: ";
for (int num : StrungDofs) {
    std::cout << num << " ";
}
std::cout << std::endl;



// generate double-precision vectors with the same length as the DOF vectors
// so we can add in stiffness components
int dofsize;
for (int c=0;c<ndoftot;c++) { //each dof
	glStif[c].resize(gldofs[c].size(),0.0);
	fprintf(output,"%d \n",glStif[c].size());}
//*************************************************************
// end of set-up phase
//  next section of the code will
//  a) generate element stiffness matrices and load them in the global matrix
//  b) define teh load vector
//  c) set up and solve or displacements
//**************************************************



// loop through element groups
// element groups all have the same type of element
// (plane strain, plane stress, axisymmetric, 3d)
// all elements in an element group are either remodeling or not remodeling

for(c=0;c<nelgps;c++){  //each element group

	  //  eldat[0]=numel - number of elements in the group
	  //  eldat[1]=eltyp - element type
	  //  eldat[2]=maxnodes - maximum nodes in any element in the group
	  //  eldat[3]=remflag - flag for whether to let the elements change density over time);

    numel=elarry[c*2][0];
    eltyp=elarry[c*2][1];
    maxnodes=elarry[c*2][2];
    remflag=elarry[c*2][3];


    if(eltyp == 3)elsize= 4 + maxnodes*4;  //3d elements
    if(eltyp/10 == 2)elsize= 4 + maxnodes*3;  //2d elements

	  for(c2=0;c2<numel;c2++){  //each element

	  offset=c2*elsize;
	  offset2=offset+maxnodes+4; //offset2 is the offset to local dofs in element data
	  //offset is an offset within eldat - eldat points to an array of individual element data
	  /* ********************************************************************
	     for each element, setting up the data for 2D solids
	     eldat[offset]=element number
	     eldat[offset+1]=material property set number
	     eldat[offset+2]=number of nodes
	     eldat[offset+3]=element density index (-1 if the element is not remodeling)
	     ********************************************************************   */
	  int elnumber=elarry[c*2+1][offset];
	  int matnumber=elarry[c*2+1][offset+1];
	  int enodes=elarry[c*2+1][offset+2];
	  int denind=elarry[c*2+1][offset+3];
	  int off2nodes=offset+4;
			  //for each node, get coordinates
	  int elnp;
	  int dfind=0;
	  for(c3=0;c3<enodes;c3++){
		  elnp=elarry[c*2+1][off2nodes+c3];
//		  printf("%d %d\n",c3,elnp);

		  if(eltyp/10 == 2){ //eltyp=20, 21, 22 are 2d (defined in the y-z plane)
		    xx2d[c3][0]=coorptr[(elnp-1)*3+1];
            xx2d[c3][1]=coorptr[(elnp-1)*3+2];  //2d analysis occurrs in the y-z plane
		  eNodDF[dfind++]=dofptr[(elnp-1)*3+1];
		  eNodDF[dfind++]=dofptr[(elnp-1)*3+2];}
		  if(eltyp == 3){ // 3d elements
				xx3d[c3][0]=coorptr[(elnp-1)*3+0];
				xx3d[c3][1]=coorptr[(elnp-1)*3+1];
				xx3d[c3][2]=coorptr[(elnp-1)*3+2];
		        eNodDF[dfind++]=dofptr[(elnp-1)*3+0];
		        eNodDF[dfind++]=dofptr[(elnp-1)*3+1];
		        eNodDF[dfind++]=dofptr[(elnp-1)*3+2];}

	  } //each node

	  double matpass[4];
	  int elemdof=0;
	  if (eltyp/10==2){
      elemdof=2*enodes;
	  s=new double[elemdof*(elemdof+1)/2];
	  matpass[0]=matptr[(matnumber-1)*4];
	  matpass[1]=matptr[(matnumber-1)*4+1];
	  matpass[2]=matptr[(matnumber-1)*4+2];
	  matpass[3]=matptr[(matnumber-1)*4+3];
	  }
	  if(eltyp==3){
	   elemdof=3*enodes;
	   s=new double[elemdof*(elemdof+1)/2];
	   matpass[0]=matptr[(matnumber-1)*4];
	   matpass[1]=matptr[(matnumber-1)*4+1];
	   matpass[2]=matptr[(matnumber-1)*4+2];
	  }

	  gammar= denind>=0 ? gamptr[denind] : -1;
	  gexpnt=3.1; //dummy value

	  if(eltyp/10 == 2){
	  int Elret=twodsolid(elnumber,enodes,eltyp,3,matpass,xx2d,s,output,volume,gammar,gexpnt);
	  printf("elem %d\n",elnumber);
	  }
	  if(eltyp == 3){
	  int Elret=threedsolid(elnumber,enodes,3,matpass,xx3d,s,output,volume,gammar,gexpnt);
      printf("elem %d\n, elnumber");
	  }
//*************************************
// adding the local stiffness matrix components to the global stiffness matrix
//*************************************
	  int sindex=0;
	  int AddOffset=0;
	  int growindex=0;
	  int gcolindex=0;
      sindex=0;
	  for(int rr=0;rr<elemdof;rr++){fprintf(output,"%d ",eNodDF[rr]);}
	  fprintf(output,"\n\n");

	  for(int rr=0; rr<elemdof;rr++){
		  for (int cc=0;cc<=rr;cc++){ //rr is rows in element matrix

			sindex=rr*(rr+3)/2 - (rr-cc);

				  growindex=eNodDF[rr]-1;
				  gcolindex=eNodDF[cc]-1;
	//			  fprintf(output, " add %d %d %f \n",growindex,gcolindex,s[sindex]);

				  if(gcolindex>=growindex){ //column greater than row, index as usual
					  	  	  	  	  	  //find the dof in the index array, add to the corresponding spot in the value array
				  if(growindex>=0){	  //(growindex will be -1 for fixed DOF)
				  auto it= std::find(gldofs[growindex].begin(),gldofs[growindex].end(),gcolindex);
			      AddOffset=std::distance(gldofs[growindex].begin(),it);
			      glStif[growindex][AddOffset] += s[sindex];} //if growindex>=0
				  } //if gcolindex>=growindex
				  //this takes care of negative gcolindex because gcolindex >= growindex and growindex >= 0

				  else{ //column less than the row = we add the stiffness from the transposed position
				  if(gcolindex>=0){ //that dof not fixed
				  auto it= std::find(gldofs[gcolindex].begin(),gldofs[gcolindex].end(),growindex);
				  AddOffset=std::distance(gldofs[gcolindex].begin(),it);
				  glStif[gcolindex][AddOffset] += s[sindex];  } //if gcolindex>=0
				  } //else

		  } //cc
		  fprintf(output,"\n\n");
	  } //rr

	  //print out dofs to check/debug
	  for (int c7=0;c7<ndoftot;c7++) { //each dof
//	  std::cout << "Vector elements: \n";
		  for (int c8=0;c8<glStif[c7].size();c8++) {
			  fprintf(output,"%f ",glStif[c7][c8]);
		  }
	  fprintf(output,"\n");
	  }


	  } //each element


     } //each element group

//print out dofs to check/debug
for (int c=0;c<ndoftot;c++) { //each dof
std::cout << "Vector entries: \n";
for (c2=0;c2<glStif[c].size();c2++) {
fprintf(output,"%f ",glStif[c][c2]);
}
fprintf(output,"\n");
}

// at this point in the code, the global stiffness matrix exists as a set of vectors.
// now they are concatenated to put together a stiffness matrix for the solver
// this is similar to the assembly of the column matrix.

StrungStif.reserve(total_size); // Pre-allocate memory for efficiency
for (int c=0;c<ndoftot;c++) {
    StrungStif.insert(StrungStif.end(), glStif[c].begin(), glStif[c].end());}
std::cout << "Stiffness matrix entries: ";
for (int c=0;c<total_size;c++) {
    printf("%f ",StrungStif[c]);
    fprintf(output,"%f ",StrungStif[c]);
}
std::cout << std::endl;

//****************************
//  Solver phase - sets up the GPU,
//                 sends over the matrices
//                 renumbers as available
//                 Scales if appropriate
//                 solves
//                 sends back data
//                 frees up memory
// ****************************
// this code largely follows the Nvidia example
int major, minor, patch;
cudssGetProperty(MAJOR_VERSION, &major);
cudssGetProperty(MINOR_VERSION, &minor);
cudssGetProperty(PATCH_LEVEL,   &patch);
printf("CUDSS Version (Major,Minor,PatchLevel): %d.%d.%d\n", major, minor, patch);

cudaError_t cuda_error = cudaSuccess;
cudssStatus_t status = CUDSS_STATUS_SUCCESS;

n = ndoftot; // degrees of freedom
int nnz = total_size; //number of non-zero entries
int nrhs = 1; //number or right hand sides

int *csr_offsets_h = NULL;  //offsets in the matrix file to diagonal entries
int *csr_columns_h = NULL; // map to the matrix entries
double *csr_values_h = NULL; // values of the matrix entries
double *x_values_h = NULL;  //values of the rhs entries
double *b_values_h = NULL; // values of the solution entries
double *diag_h = NULL;
double *row_scale_h = NULL;
double *col_scale_h = NULL;
//these are pointers to the data structures on the GPU device
//  _d refers to "device"
int *csr_offsets_d = NULL;
int *csr_columns_d = NULL;
double *csr_values_d = NULL;
double *x_values_d = NULL;
double *b_values_d = NULL;
double *diag_d = NULL;
double *row_scale_d = NULL;
double *col_scale_d = NULL;

/* Allocate host memory for the sparse input matrix A,
   right-hand side x and solution b*/

csr_offsets_h = (int*)malloc((n + 1) * sizeof(int));
csr_columns_h = (int*)malloc(nnz * sizeof(int));
csr_values_h = (double*)malloc(nnz * sizeof(double));
x_values_h = (double*)malloc(nrhs * n * sizeof(double));
b_values_h = (double*)malloc(nrhs * n * sizeof(double));

diag_h = (double*)malloc(n * sizeof(double));
row_scale_h = (double*)malloc(n * sizeof(double));
col_scale_h = (double*)malloc(n * sizeof(double));

if (!csr_offsets_h || ! csr_columns_h || !csr_values_h ||
    !x_values_h || !b_values_h || !diag_h || !row_scale_h || !col_scale_h) {
    printf("Error: host memory allocation failed\n");
    return -1;
}

//get memory allocated for matrices and on the GPU
cudaMalloc(&csr_offsets_d, (n + 1) * sizeof(int));
cudaMalloc(&csr_columns_d, nnz * sizeof(int));
cudaMalloc(&csr_values_d, nnz * sizeof(double));
cudaMalloc(&b_values_d, nrhs * n * sizeof(double));
cudaMalloc(&x_values_d, nrhs * n * sizeof(double));
cudaMalloc(&diag_d, n * sizeof(double));
cudaMalloc(&row_scale_d, n * sizeof(double));
cudaMalloc(&col_scale_d, n * sizeof(double));

for(int ii=0;ii<n;ii++){
	csr_offsets_h[ii]=diagindex[ii];
	b_values_h[ii]=loadptr[ii];  // only one load case for the first try
}
csr_offsets_h[n]=nnz;

//echo b for checking
    for (int i = 0; i < n; i++)
    fprintf(output,"bvalues[%d] = %f\n", i, b_values_h[i]);


//fprintf(output,"overall matrix\n");
for(int ii=0;ii<nnz;ii++){
	csr_columns_h[ii]=StrungDofs[ii];
	csr_values_h[ii]=StrungStif[ii];
	//fprintf(output,"%d %f\n",StrungDofs[ii],StrungStif[ii]);
}

// TO DO: set up values for arrays

cudaMemcpy(csr_offsets_d, csr_offsets_h, (n + 1) * sizeof(int),cudaMemcpyHostToDevice);
cudaMemcpy(csr_columns_d, csr_columns_h, nnz * sizeof(int),cudaMemcpyHostToDevice);
cudaMemcpy(csr_values_d, csr_values_h, nnz * sizeof(double),cudaMemcpyHostToDevice);
cudaMemcpy(b_values_d, b_values_h, nrhs * n * sizeof(double),cudaMemcpyHostToDevice);

//Setting up the cudss variables and data structures

/* Create a cuDSS library handle */
cudssHandle_t handle;

cudssCreate(&handle);
/* Set up a custom stream for the library handle */
cudaStream_t stream = NULL;
cudssSetStream(handle, stream);

/* Create a cuDSS solver configuration and initialize data objects */
cudssConfig_t solverConfig;
cudssConfigCreate(&solverConfig);

cudssData_t solverData;
cudssDataCreate(handle, &solverData);

/* Setting algorithmic parameters */
cudssAlgType_t reorder_alg = CUDSS_ALG_DEFAULT;
cudssConfigSet(solverConfig, CUDSS_CONFIG_REORDERING_ALG, &reorder_alg, sizeof(cudssAlgType_t));
int ione = 1;
cudssConfigSet(solverConfig, CUDSS_CONFIG_USE_MATCHING, &ione, sizeof(int));

cudssAlgType_t matching_alg = CUDSS_ALG_DEFAULT; // matching with scaling, same as CUDSS_ALG_5
cudssConfigSet(solverConfig, CUDSS_CONFIG_MATCHING_ALG, &matching_alg, sizeof(cudssAlgType_t));

/* Create matrix objects for the right-hand side b and solution x (as dense matrices). */
cudssMatrix_t xc, b;

int64_t nrows = n, ncols = n;
int ldb = ncols, ldx = nrows;
cudssMatrixCreateDn(&b, ncols, nrhs, ldb, b_values_d, CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR);
cudssMatrixCreateDn(&xc, nrows, nrhs, ldx, x_values_d, CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR);

/* Create a matrix object for the sparse input matrix. */
// uses device data
cudssMatrix_t A;
cudssMatrixType_t mtype     = CUDSS_MTYPE_SYMMETRIC;
cudssMatrixViewType_t mview = CUDSS_MVIEW_UPPER;
cudssIndexBase_t base       = CUDSS_BASE_ZERO;

cudssMatrixCreateCsr(&A, nrows, ncols, nnz, csr_offsets_d, NULL, csr_columns_d, csr_values_d, CUDA_R_32I, CUDA_R_64F, mtype, mview, base);

/* Symbolic factorization */
cudssExecute(handle, CUDSS_PHASE_ANALYSIS, solverConfig, solverData, A, xc, b);
size_t sizeWritten;
int perm[ndoftot];
cudssDataGet(handle, solverData, CUDSS_DATA_PERM_REORDER_ROW, &perm, sizeof(perm), &sizeWritten);
for (int i = 0; i < n; i++)
printf("reorder row perm[%d] = %d\n", i, perm[i]);

cudssDataGet(handle, solverData, CUDSS_DATA_PERM_REORDER_COL, &perm, sizeof(perm), &sizeWritten);
for (int i = 0; i < n; i++)
    printf("reorder col perm[%d] = %d\n", i, perm[i]);

int used_matching = 0;
cudssConfigGet(solverConfig, CUDSS_CONFIG_USE_MATCHING, &used_matching, sizeof(int), &sizeWritten);

if (used_matching) {

cudssDataGet(handle, solverData, CUDSS_DATA_PERM_MATCHING, &perm, sizeof(perm), &sizeWritten);
}

/*  Note: currently these features are only supported for CUDSS_ALG_1 and CUDSS_ALG_2
    reordering algorithms. */

if (reorder_alg == CUDSS_ALG_1 || reorder_alg == CUDSS_ALG_2) {

cudssDataGet(handle, solverData, CUDSS_DATA_PERM_ROW, &perm, sizeof(perm), &sizeWritten);

for (int i = 0; i < ndoftot; i++)
    printf("final row perm[%d] = %d\n", i, perm[i]);

cudssDataGet(handle, solverData, CUDSS_DATA_PERM_COL, &perm, sizeof(perm), &sizeWritten);
for (int i = 0; i < ndoftot; i++)
    printf("final col perm[%d] = %d\n", i, perm[i]);
}

int64_t memory_estimates[6] = {0};
cudssDataGet(handle, solverData, CUDSS_DATA_MEMORY_ESTIMATES, &memory_estimates, sizeof(memory_estimates), &sizeWritten);

printf("memory estimates: device: %ld (stable) %ld (peak)\n",
        memory_estimates[0], memory_estimates[1]);
printf("memory estimates: host: %ld (stable) %ld (peak)\n",
        memory_estimates[2], memory_estimates[3]);
printf("memory estimates: hybrid peak: %ld = (GPU) %ld (CPU)\n",
        memory_estimates[4], memory_estimates[5]);
fflush(0);

/* Factorization */
cudssExecute(handle, CUDSS_PHASE_FACTORIZATION, solverConfig, solverData, A, xc, b);
/* (optional) Recommended: getting runtime errors from device side
    Note: cudssDataGet is a synchronous API.
    Note: per cuDSS documentation, CUDSS_DATA_INFO is always returned
    as a pointer to int.
*/
int info;

cudssDataGet(handle, solverData, CUDSS_DATA_INFO, &info, sizeof(info), &sizeWritten);

printf("cuDSS info = %d\n", info);
int npivots;
cudssDataGet(handle, solverData, CUDSS_DATA_NPIVOTS, &npivots, sizeof(npivots), &sizeWritten);
printf("cuDSS npivots = %d\n", npivots);

int inertia[2];
cudssDataGet(handle, solverData, CUDSS_DATA_INERTIA, &inertia, sizeof(inertia), &sizeWritten);

printf("cuDSS inertia = %d %d\n", inertia[0], inertia[1]);
/* (optional) Recommended: getting data back when the user does not know the size
    Note: cudssDataGet is a synchronous API.
*/
cudssDataGet(handle, solverData, CUDSS_DATA_LU_NNZ, NULL, 0, &sizeWritten);

printf("cuDSS requests %zu bytes for returning the #nnz in the factors\n", sizeWritten);
/* Note: per cuDSS documentation, CUDSS_DATA_LU_NNZ is always returned as a
   pointer to int64_t.
   In the general case, the user can allocate the returned #bytes and pass a
   pointer to the allocated piece of memory into the cudssDataGet API and
   reinterpret the results with correct types (defined in the documentation).
*/
int64_t lu_nnz;
cudssDataGet(handle, solverData, CUDSS_DATA_LU_NNZ, &lu_nnz, sizeof(int64_t), NULL);

printf("cuDSS #nnz in LU = %ld\n", lu_nnz);

/* (optional) Getting data back when the user knows the size
*/
cudssDataGet(handle, solverData, CUDSS_DATA_DIAG, diag_d, n*sizeof(double), &sizeWritten);
cudaMemcpy(diag_h, diag_d, n * sizeof(double),cudaMemcpyDeviceToHost);

for (int i = 0; i < n; i++)
    printf("diag[%d] = %f\n", i, diag_h[i]);
/* Note: scales are always real (never complex) and positive */

if (matching_alg == CUDSS_ALG_DEFAULT || matching_alg == CUDSS_ALG_5) {

	cudssDataGet(handle, solverData, CUDSS_DATA_SCALE_ROW, row_scale_d, n*sizeof(double), &sizeWritten);
	cudaMemcpy(row_scale_h, row_scale_d, n * sizeof(double), cudaMemcpyDeviceToHost);

	for (int i = 0; i < n; i++)
        printf("row scale[%d] = %f\n", i, row_scale_h[i]);

	cudssDataGet(handle, solverData, CUDSS_DATA_SCALE_COL, col_scale_d, n*sizeof(double), &sizeWritten);
	cudaMemcpy(col_scale_h, col_scale_d, n * sizeof(double), cudaMemcpyDeviceToHost);

    for (int i = 0; i < n; i++)
        printf("col scale[%d] = %f\n", i, col_scale_h[i]);
}

int iter_ref_nsteps = 2;
cudssConfigSet(solverConfig, CUDSS_CONFIG_IR_N_STEPS, &iter_ref_nsteps, sizeof(iter_ref_nsteps));

/* Solving */
cudssExecute(handle, CUDSS_PHASE_SOLVE, solverConfig, solverData, A, xc, b);
/* Destroying opaque objects, matrix wrappers and the cuDSS library handle */
cudssMatrixDestroy(A);
cudssMatrixDestroy(b);
cudssMatrixDestroy(xc);
cudssDataDestroy(handle, solverData);
cudssConfigDestroy(solverConfig);
cudssDestroy(handle);
cudaStreamSynchronize(stream);
/* Print the solution and compare against the exact solution */
cudaMemcpy(x_values_h, x_values_d, nrhs * n * sizeof(double), cudaMemcpyDeviceToHost);

int passed = 1;
for (int i = 0; i < n; i++) {
    printf("x[%d] = %1.4f\n", i, x_values_h[i]);
    fprintf(output,"x[%d] = %1.4f \n", i, x_values_h[i]);

}
//Free up allocated memory structures
EXIT_FREE;


}





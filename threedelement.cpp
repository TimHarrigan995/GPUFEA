/* c subroutine to do stiffness calculations for 3D solid finite elements */

#include <stdio.h>
#include <math.h>
double stdm3d(int nnodes,double xx[],double b[][6],double r,double s,double t,int nel);
int matstf3d(double d[6][6],double matprops[4],double gammar,double gexpnt);

//void matstf3d(double [6][6],int,int);
double pow(double,double);

static double xg[4][4]={{0.,0.,0.,0.},
                        {-.5773502691896,0.5773502691896,0.,0.},
                        {-.7745966692415,0.,0.7745966692415,0.},
                        {-.8611363115941,-.3399810435849,
                               0.3399810435849,0.8611363115941}};
static double wgt[4][4]={{2.,0.,0.,0.},
                         {1.,1.,0.,0.},
                         {.5555555555556,.888888888889,.5555555555556,0.},
                         {.3478548451375,.6521451548625,
                                         .6521451548625,.3478548451375}};

static double *elstf,*gblstf,*coords;
static int dofpernode,*nodofs,*eldofs;
static double ym;
static double pr;

int threedsolid(int nel,int nnodes,int nint,double matprops[3],double xps[][3],double s[],FILE *output,double *volume,double gammar,double gexpnt){

   double d[6][6],db[6],wt,xbar,f,g,h,det,a;
   double ri,si,ti,stiff;
   int c,c2,c3,c4,j,k,i,ist,l,ni,dof,ind,lstif;
   double b[24][6];
   *volume=0.;

   ni=nint -1;
   dof=3*nnodes;

   double xx[3*nnodes];

      for(int c=0;c<nnodes;c++){
   	  for(int cc=0;cc<3;cc++)xx[c*3+cc]=xps[c][cc];
   		}


   lstif=dof*(dof+1)/2;

/* get stress-strain law */
//int matstf3d(double d[][6],double matprops[4],double gammar,double gexpnt);


   matstf3d(d,matprops,gammar,gexpnt);
   ist=6;
/* calculate element stiffness matrix */

   for(c2=0;c2<lstif;c2++){
      s[c2]=0.; }


   for(c=0;c<nint;c++){
     ri = xg[ni][c];
     for(c2=0;c2<nint;c2++){
       si = xg[ni][c2];
       for(c3=0;c3<nint;c3++){
         ti = xg[ni][c3];

   /* evaluate strain-displacement matrix in stdm3d */
   
//double stdm3d(int nnodes,double xx[][3],double b[][6],double r,double s,double t,int nel);


         if((det=stdm3d(nnodes,xx,b,ri,si,ti,nel)) < 1.0e-8){
/*    fprintf(output,"\n***negative jacobian determinant in el # %d\n",nel);
    fprintf(output,"***determinant = %15.12lf\n",det);
*/
    return(nel);}

         wt= wgt[ni][c] * wgt[ni][c2] * wgt[ni][c3] * det;

         *volume += wt;

         for(j=0;j<dof;j++){

           for(k=0;k<ist;k++){
             db[k]=0.;
             for(l=0;l<ist;l++)db[k] += d[k][l]*b[j][l];
           } /* for(k... */

           for(i=0;i<j+1;i++){
           stiff=0.;
             for(l=0;l<ist;l++)stiff += b[i][l]*db[l];
             ind=j*(j+1)/2 +i;
             s[ind] += stiff*wt;
           } /* for (i... */

         } /* for(j... */

        } /* for(c3... */
      } /* for(c2... */
    } /* for(c... */
return(0);
} /* end of routine */

double stdm3d(int nnodes,double xx[],double b[][6],double r,double s,double t,int nel)
{
double h[27],p[27][3],xj[3][3],xji[3][3],det,rp,sp,tp,rm,sm,tm,dum;
int i,j,k;
rp=1.0+r;
sp=1.0+s;
tp=1.0+t;
rm=1.0-r;
sm=1.0-s;
tm=1.0-t;


/* coordinate derivatives of interpolation functions */

/* with respect to r */
p[0][0]=  0.125*sp*tp;
p[1][0]= -p[0][0];
p[2][0]= -0.125*sm*tp;
p[3][0]= -p[2][0];
p[4][0]=  0.125*sp*tm;
p[5][0]= -p[4][0];
p[6][0]= -0.125*sm*tm;
p[7][0]= -p[6][0];

if(nnodes >8){
p[8][0] = -0.5*sp*r*tp;
p[0][0] -= 0.5*p[8][0];
p[1][0] -= 0.5*p[8][0];}
if(nnodes >9){
p[9][0] = -0.25*(1.0-s*s)*tp;
p[1][0] -= 0.5*p[9][0];
p[2][0] -= 0.5*p[9][0];}
if(nnodes >10){
p[10][0] = -0.5*sm*r*tp;
p[2][0] -= 0.5*p[10][0];
p[3][0] -= 0.5*p[10][0];}
if(nnodes >11){
p[11][0] = 0.25*(1.0-s*s)*tp;
p[3][0] -= 0.5*p[11][0];
p[0][0] -= 0.5*p[11][0];}

if(nnodes >12){
p[12][0] = -0.5*sp*r*tm;
p[4][0] -= 0.5*p[12][0];
p[5][0] -= 0.5*p[12][0];}
if(nnodes >13){
p[13][0] = -0.25*(1.0-s*s)*tm;
p[5][0] -= 0.5*p[13][0];
p[6][0] -= 0.5*p[13][0];}
if(nnodes >14){
p[14][0] = -0.5*sm*r*tm;
p[6][0] -= 0.5*p[14][0];
p[7][0] -= 0.5*p[14][0];}
if(nnodes >15){
p[15][0] = 0.25*(1.0-s*s)*tm;
p[7][0] -= 0.5*p[15][0];
p[4][0] -= 0.5*p[15][0];}

if(nnodes >16){
p[16][0] = 0.25*(1.0-t*t)*sp;
p[0][0] -= 0.5*p[16][0];
p[4][0] -= 0.5*p[16][0];}
if(nnodes >17){
p[17][0] = -0.25*(1.0-t*t)*sp;
p[1][0] -= 0.5*p[17][0];
p[5][0] -= 0.5*p[17][0];}
if(nnodes >18){
p[18][0] = -0.25*(1.0-t*t)*sm;
p[2][0] -= 0.5*p[18][0];
p[6][0] -= 0.5*p[18][0];}
if(nnodes >19){
p[19][0] = 0.25*(1.0-t*t)*sm;
p[3][0] -= 0.5*p[19][0];
p[7][0] -= 0.5*p[19][0];}

/* add nodes 21-27 here  */

/* with respect to s */
p[0][1]= 0.125*rp*tp;
p[1][1]= 0.125*rm*tp;
p[2][1]= -p[1][1];
p[3][1]= -p[0][1];
p[4][1]= 0.125*rp*tm;
p[5][1]= 0.125*rm*tm;
p[6][1]= -p[5][1];
p[7][1]= -p[4][1];

if(nnodes >8){
p[8][1]= 0.25*(1.0-r*r)*tp;
p[0][1] -= 0.5*p[8][1];
p[1][1] -= 0.5*p[8][1];}
if(nnodes >9){
p[9][1]= -0.5*s*rm*tp;
p[1][1] -= 0.5*p[9][1];
p[2][1] -= 0.5*p[9][1];}
if(nnodes >10){
p[10][1]= -0.25*(1.0-r*r)*tp;
p[2][1] -= 0.5*p[10][1];
p[3][1] -= 0.5*p[10][1];}
if(nnodes >11){
p[11][1] = -0.5*s*rp*tp;
p[3][1] -= 0.5*p[11][1];
p[0][1] -= 0.5*p[11][1];}

if(nnodes >12){
p[12][1]= 0.25*(1.0-r*r)*tm;
p[4][1] -= 0.5*p[12][1];
p[5][1] -= 0.5*p[12][1];}
if(nnodes >13){
p[13][1]= -0.5*s*rm*tm;
p[5][1] -= 0.5*p[13][1];
p[6][1] -= 0.5*p[13][1];}
if(nnodes >14){
p[14][1]= -0.25*(1.0-r*r)*tm;
p[6][1] -= 0.5*p[14][1];
p[7][1] -= 0.5*p[14][1];}
if(nnodes >15){
p[15][1] = -0.5*s*rp*tm;
p[7][1] -= 0.5*p[15][1];
p[4][1] -= 0.5*p[15][1];}

if(nnodes >16){
p[16][1] = 0.25*(1.0-t*t)*rp;
p[0][1] -= 0.5*p[16][1];
p[4][1] -= 0.5*p[16][1];}
if(nnodes >17){
p[17][1] = 0.25*(1.0-t*t)*rm;
p[1][1] -= 0.5*p[17][1];
p[5][1] -= 0.5*p[17][1];}
if(nnodes >18){
p[18][1] = -0.25*(1.0-t*t)*rm;
p[2][1] -= 0.5*p[18][1];
p[6][1] -= 0.5*p[18][1];}
if(nnodes >19){
p[19][1]= -0.25*(1.0-t*t)*rp;
p[3][1] -= 0.5*p[19][1];
p[7][1] -= 0.5*p[19][1];}

/* derivatives with respect to t */

p[0][2]= 0.125*rp*sp;
p[1][2]= 0.125*rm*sp;
p[2][2]= 0.125*rm*sm;
p[3][2]= 0.125*rp*sm;
p[4][2]= -p[0][2];
p[5][2]= -p[1][2];
p[6][2]= -p[2][2];
p[7][2]= -p[3][2];

if(nnodes >8){
p[8][2]= 0.25*(1-r*r)*sp;
p[0][2] -= 0.5*p[8][2];
p[1][2] -= 0.5*p[8][2];}
if(nnodes >9){
p[9][2]= 0.25*(1-s*s)*rm;
p[1][2] -= 0.5*p[9][2];
p[2][2] -= 0.5*p[9][2];}
if(nnodes >10){
p[10][2]= 0.25*(1.0-r*r)*sm;
p[2][2] -= 0.5*p[10][2];
p[3][2] -= 0.5*p[10][2];}
if(nnodes >11){
p[11][2] = 0.25*(1-s*s)*rp;
p[3][2] -= 0.5*p[11][2];
p[0][2] -= 0.5*p[11][2];}

if(nnodes >12){
p[12][2]= -p[8][2];
p[4][2] -= 0.5*p[12][2];
p[5][2] -= 0.5*p[12][2];}
if(nnodes >13){
p[13][2]= -p[9][2];
p[5][2] -= 0.5*p[13][2];
p[6][2] -= 0.5*p[13][2];}
if(nnodes >14){
p[14][2]= -p[10][2];
p[6][2] -= 0.5*p[14][2];
p[7][2] -= 0.5*p[14][2];}
if(nnodes >15){
p[15][2] = -p[11][2];
p[7][2] -= 0.5*p[15][2];
p[4][2] -= 0.5*p[15][2];}

if(nnodes >16){
p[16][2] = -0.5*t*rp*sp;
p[0][2] -= 0.5*p[16][2];
p[4][2] -= 0.5*p[16][2];}
if(nnodes >17){
p[17][2] = -0.5*t*rm*sp;
p[1][2] -= 0.5*p[17][2];
p[5][2] -= 0.5*p[17][2];}
if(nnodes >18){
p[18][2] = -0.5*t*rm*sm;
p[2][2] -= 0.5*p[18][2];
p[6][2] -= 0.5*p[18][2];}
if(nnodes >19){
p[19][2] = -0.5*t*rp*sm;
p[3][2] -= 0.5*p[19][2];
p[7][2] -= 0.5*p[19][2];}


/* add code for nodes 21-27 here */

/* jacobian at r,s,t */
   for(i=0;i<3;i++){
      for(j=0;j<3;j++){
         dum=0.;
         for(k=0;k<nnodes;k++)dum += p[k][i]*xx[k*3+j];
         xj[i][j]=dum;
      } /* for(j... */
   } /* for(i... */

/* determinant of jacobian */

   det= xj[0][0]*(xj[1][1]*xj[2][2]-xj[2][1]*xj[1][2])
      - xj[0][1]*(xj[1][0]*xj[2][2]-xj[2][0]*xj[1][2])
      + xj[0][2]*(xj[1][0]*xj[2][1]-xj[2][0]*xj[1][1]);



if(det < 1.0e-8)return det;

/* inverse of the jacobian */

   dum=1.0/ det;
   xji[0][0]=  (xj[1][1]*xj[2][2] - xj[2][1]*xj[1][2])*dum;
   xji[0][1]= -(xj[0][1]*xj[2][2] - xj[2][1]*xj[0][2])*dum;
   xji[0][2]=  (xj[0][1]*xj[1][2] - xj[1][1]*xj[0][2])*dum;

   xji[1][0]= -(xj[1][0]*xj[2][2] - xj[2][0]*xj[1][2])*dum;
   xji[1][1]=  (xj[0][0]*xj[2][2] - xj[2][0]*xj[0][2])*dum;
   xji[1][2]= -(xj[0][0]*xj[1][2] - xj[1][0]*xj[0][2])*dum;

   xji[2][0]=  (xj[1][0]*xj[2][1] - xj[2][0]*xj[1][1])*dum;
   xji[2][1]= -(xj[0][0]*xj[2][1] - xj[2][0]*xj[0][1])*dum;
   xji[2][2]=  (xj[0][0]*xj[1][1] - xj[1][0]*xj[0][1])*dum;

/* evaluate B */

   j= -1;
   for(k=0;k<nnodes;k++){
      j+=3;

      b[j-2][0]=0.;
      b[j-1][0]=0.;
      b[j][0]=0.;

      b[j-2][1]=0.;
      b[j-1][1]=0.;
      b[j][1]=0.;

      b[j-2][2]=0.;
      b[j-1][2]=0;
      b[j][2]=0.;


     for(i=0;i<3;i++){
         b[j-2][0] += xji[0][i]*p[k][i];
         b[j-1][1] += xji[1][i]*p[k][i];
         b[j][2]   += xji[2][i]*p[k][i];
      } /* for (i... ) */

     b[j][5]=b[j-1][3]=b[j-2][0];
     b[j][4]=b[j-2][3]=b[j-1][1];
     b[j-2][5]=b[j-1][4]=b[j][2];
     b[j-1][5]=b[j-2][4]=b[j][3]=0.;
   } /* for (k...) */

return det;

} /* end of routine */

//   matstf3d(d,matprops,gammar,gexpnt);
int matstf3d(double d[][6],double matprops[3],double gammar,double gexpnt)
//int matstf3d(d,nmat,phindex)
//double d[6][6];
//int nmat,phindex;
{
   int i,j,k;
//   extern double *matptr,*gamptr;
//   extern double modexpnt,transexpnt;
   double a,b,c,frac;
   double ym,pr;

//frac=modexpnt/transexpnt;

//   k=(nmat-1)*3;
//   ym=matptr[k];
//   pr=matptr[k+1];

	ym=matprops[1];
	pr=matprops[2];

	if(gammar>=0)ym *= pow(gammar,gexpnt);

   for(i=0;i<6;i++){
      for(j=0;j<6;j++)d[i][j]=0.;
      } /* for(i... */

      a=ym*(1.-pr)/((1.+pr)*(1.-2.*pr));
      b=pr/(1.-pr);
      c=(1-2.*pr)/(2-2*pr);
      d[0][0]=d[1][1]=d[2][2]=a;
      d[1][0]=d[0][1]=d[2][0]=d[0][2]=d[2][1]=d[1][2]=a*b;
      d[3][3]=d[4][4]=d[5][5]=a*c;
      return(0);
} /* end of routine */


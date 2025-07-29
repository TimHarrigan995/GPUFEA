/* c program to do stiffness calculations for 4-sided finite elements */
/* itype tells what type of element

  itype=0: axisymmetric
  itype=1: plane strain
  itype=2: plane stress

int matstf2d(int itype,double d[][4],double *thic,int matprp,int phindex)

*/

#include <stdio.h>
#include <math.h>
//int twodsolid(int,int,int,int,int,double [][2],double [][4],double *,
//                 FILE *,double *,int);
double stdm(int,double [][2],double [][4],double,double,double *,int,int);

int matstf2d(int itype,double d[4][4],double matprops[],double gammar,double gexpnt);


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
extern double modexpnt,transexpnt;

int twodsolid(int nel,int nnodes,int itype,int nint,double matprops[4],double xx[][2],double s[],FILE *output,
		double *volume,double gammar,double gexpnt)
{
double d[4][4],db[4],wt,xbar,f,g,h,det,a,b[16][4];
double ri,si,stiff,ym,pr,thic;
int c,c2,j,k,i,ist,l,ni,dof,ind,stifsiz;

for (int i=0;i<nnodes;i++)fprintf (output,"%f %f\n",xx[i][0],xx[i][1]);

*volume=0.;
ni=nint -1;
dof=2*nnodes;
stifsiz=dof*(dof+1)/2;
for(c=0;c<4;c++){
   for(c2=0;c2<4;c2++)d[c][c2]=0.;}


/* get stress-strain law */
matstf2d(itype,d,matprops,gammar,gexpnt);

thic=matprops[3];

/* calculate element stiffness matrix */

for(c2=0;c2<stifsiz;c2++){
s[c2]=0.; }

ist= (itype == 0) ? 4 : 3;

for(c=0;c<nint;c++){
  ri = xg[ni][c];
  for(c2=0;c2<nint;c2++){
    si = xg[ni][c2];

    /* evaluate strain-displacement matrix in stdm */
    if((det=stdm(nnodes,xx,b,ri,si,&xbar,nel,itype)) < 1.0e-8){

    return(nel);}  //returns for invalid elements

    if(itype > 0)xbar=thic;

    wt= wgt[ni][c] * wgt[ni][c2] * xbar * det;

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

  } /* for(c2... */
} /* for(c... */
return(0);
} /* end of routine */

double stdm(int nnodes,double xx[4][2],double b[][4],double r,double s,double *xbar,int nel,int itype)
{
double h[9],p[9][2],xj[2][2],xji[2][2],det,rp,sp,rm,sm,dum;
int i,j,k;
rp=1.0+r;
sp=1.0+s;
rm=1.0-r;
sm=1.0-s;
/* interpolation functions */
h[0]=0.25*rp*sp;
h[1]=0.25*rm*sp;
h[2]=0.25*rm*sm;
h[3]=0.25*rp*sm;
if(nnodes>4) {
h[4]= 0.5*(1.0-r*r)*sp;
h[0]-= 0.5*h[4];
h[1]-= 0.5*h[4];}
if(nnodes >5){
h[5]= 0.5*(1.0-s*s)*rm;
h[1] -= 0.5*h[5];
h[2] -= 0.5*h[5];}
if(nnodes >6){
h[6] = 0.5*(1.0-r*r)*sm;
h[2] -= 0.5*h[6];
h[3] -= 0.5*h[6];}
if(nnodes >7){
h[7] = 0.5*(1.0-s*s)*rp;
h[3] -= 0.5*h[7];
h[0] -= 0.5*h[7];}
if(nnodes >8){
h[8] = (1.0-r*r)*(1.0-s*s);
h[0] += 0.25*h[8];
h[1] += 0.25*h[8];
h[2] += 0.25*h[8];
h[3] += 0.25*h[8];
h[4] -= 0.5*h[8];
h[5] -= 0.5*h[8];
h[6] -= 0.5*h[8];
h[7] -= 0.5*h[8];}

/* coordinate derivatives of interpolation functions */
/* with respect to r */
p[0][0]=0.25*sp;
p[1][0]= -p[0][0];
p[2][0]= -0.25*sm;
p[3][0]= -p[2][0];
if(nnodes >4){
p[4][0] = -sp*r;
p[0][0] -= 0.5*p[4][0];
p[1][0] -= 0.5*p[4][0];}
if(nnodes >5){
p[5][0] = -0.5*(1.0-s*s);
p[1][0] -= 0.5*p[5][0];
p[2][0] -= 0.5*p[5][0];}
if(nnodes >6){
p[6][0] = -sm*r;
p[2][0] -= 0.5*p[6][0];
p[3][0] -= 0.5*p[6][0];}
if(nnodes >7){
p[7][0] = 0.5*(1.0-s*s);
p[3][0] -= 0.5*p[7][0];
p[0][0] -= 0.5*p[7][0];}
if(nnodes >8){
p[8][0] = -2*r*(1.0-s*s);
p[0][0] += 0.25*p[8][0];
p[1][0] += 0.25*p[8][0];
p[2][0] += 0.25*p[8][0];
p[3][0] += 0.25*p[8][0];
p[4][0] -= 0.5*p[8][0];
p[5][0] -= 0.5*p[8][0];
p[6][0] -= 0.5*p[8][0];
p[7][0] -= 0.5*p[8][0];}

/* with respect to s */
p[0][1]= 0.25*rp;
p[1][1]= 0.25*rm;
p[2][1]= -p[1][1];
p[3][1]= -p[0][1];
if(nnodes >4){
p[4][1]= 0.5*(1.0-r*r);
p[0][1] -= 0.5*p[4][1];
p[1][1] -= 0.5*p[4][1];}
if(nnodes >5){
p[5][1]= -s*rm;
p[1][1] -= 0.5*p[5][1];
p[2][1] -= 0.5*p[5][1];}
if(nnodes >6){
p[6][1]= -0.5*(1.0-r*r);
p[2][1] -= 0.5*p[6][1];
p[3][1] -= 0.5*p[6][1];}
 if(nnodes >7){
p[7][1] = -s*rp;
p[3][1] -= 0.5*p[7][1];
p[0][1] -= 0.5*p[7][1];}
if(nnodes >8){
p[8][1] = -2*s*(1.0-r*r);
p[0][1] += 0.25*p[8][1];
p[1][1] += 0.25*p[8][1];
p[2][1] += 0.25*p[8][1];
p[3][1] += 0.25*p[8][1];
p[4][1] -= 0.5*p[8][1];
p[5][1] -= 0.5*p[8][1];
p[6][1] -= 0.5*p[8][1];
p[7][1] -= 0.5*p[8][1];}

/* jacobian at r,s */
for(i=0;i<2;i++){
  for(j=0;j<2;j++){
    dum=0.;
    for(k=0;k<nnodes;k++)dum += p[k][i]*xx[k][j];
    xj[i][j]=dum;
  } /* for(j... */
} /* for(i... */

/* determinant of jacobian */
det= xj[0][0]*xj[1][1] - xj[0][1]*xj[1][0];
if(det < 1.0e-8)return det;

/* inverse of the jacobian */

dum=1.0/ det;
xji[0][0]= xj[1][1]*dum;
xji[0][1]= -xj[0][1]*dum;
xji[1][0]= -xj[1][0]*dum;
xji[1][1]= xj[0][0]*dum;

/* evaluate B */

j= -1;
for(k=0;k<nnodes;k++){
  j+=2;
  b[j-1][0]=0.;
  b[j][0]=0.;
  b[j-1][1]=0.;
  b[j][1]=0.;
  for(i=0;i<2;i++){
    b[j-1][0] += xji[0][i]*p[k][i];
    b[j][1] += xji[1][i]*p[k][i];
  } /* for (i... ) */
  b[j][2]=b[j-1][0];
  b[j-1][2]= b[j][1];
} /* for (k...) */

/* if plane stress or strain, return */

if (itype >0)return det;

/* compute the radius at r,s for axisymmetric elements */

*xbar=0.0;

for(k=0;k<nnodes;k++)*xbar += h[k]*xx[k][0];


/* for zero radius, equate hoop and radial strains */
if (*xbar < 1.0e-8){
  for(k=0;k<(nnodes*2);k++)b[k][3]=b[k][0];
  return det;
} /* if(xbar... ) */

/* for nonzero radius, compute hoop strain stiffness */

dum =1./ *xbar;
j= -1;

for(k=0;k<nnodes;k++){
  j +=2;
  b[j][3]=0;
  b[j-1][3]=h[k]*dum;
} /* for(k...)*/
return det;
} /* end of routine */

/********************************************

  routine to get 2D material stiffness matrices

****************************************** */


int matstf2d(int itype,double d[][4],double matprops[4],double gammar,double gexpnt)
{
int c,c2,k;
double ym,pr,f,g,h,a,frac,thic;

ym=matprops[1];
pr=matprops[2];
thic=matprops[3];

if(gammar>=0)ym *= pow(gammar,gexpnt);

f=ym/(1.+pr);
g=f*pr/(1.-2.*pr);
h=f+g;
d[0][0]=h;
d[1][0]=d[0][1]=g;
d[1][1]=h;
d[2][0]=d[0][2]=d[2][1]=d[1][2]=0.;
d[2][2]=f/2.;

if(itype == 1)thic=1.0;

if(itype == 0){  /* axisymmetric */

d[3][0]=d[0][3]=g;
d[3][1]=d[1][3]=g;
d[3][2]=d[2][3]=0.;
d[3][3]=h;

} /* if(itype... */

if(itype == 2){

/* in plane stress condense the stress-strain matrix */
  for(c=0;c<3;c++){
    a=d[c][3]/d[3][3];
    for(c2=0;c2<3;c2++){
      d[c][c2] -= d[3][c2] * a;
      d[c2][c]=d[c][c2];
    }/* for(c2... ) */
  }/* for(c... ) */
} /* if(itype */

return (0);
}

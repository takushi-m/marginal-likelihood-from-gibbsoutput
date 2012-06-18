#include"./dSFMT/dSFMT.h"
#include<cmath>
const double PI = 3.14159265358979324;
dsfmt_t dsfmt;

double loggamma(double x){
  const int NN = 8;
  //const double B0 = 1.0;
  //const double B1 = -1.0 / 2.0;
  const double B2 = 1.0 / 6.0;
  const double B4 = -1.0 / 30.0;
  const double B6 = 1.0 / 42.0;
  const double B8 = -1.0 / 30.0;
  const double B10 = 5.0 / 66.0;
  const double B12 = -691.0 / 2730.0;
  const double B14 = 7.0 / 6.0;
  const double B16 = -3617.0 / 510.0;

  double v,w;

  v = 1.0;
  while(x < NN){
    v = v*x;
    x = x+1.0;
  }

  w = 1.0/(x*x);

  return ((((((((B16/(16.0*15.0))*w
                +(B14/(14.0*13.0)))*w
               +(B12/(12.0*11.0)))*w
              +(B10/(10.0*9.0)))*w
             +(B8/(8.0*7.0)))*w
            +(B6/(6.0*5.0)))*w
           +(B4/(4.0*3.0)))*w
          +(B2/(2.0*1.0)))/x
    +0.5*log(2.0*PI)-log(v)-x+(x-0.5)*log(x);
}

double rand_double_01(){
  return dsfmt_genrand_close_open(&dsfmt);//[0,1)
}

unsigned int rand_uint(){
  return dsfmt_genrand_uint32(&dsfmt);
}


double logsumexp2(double x,double y){
  /*
    x   :log(A)
    y   :log(B)
    ret :log(A+B)
   */

  if(x>y){
    return x + log(1.0+exp(y-x));
  }
  else{
    return y + log(exp(x-y)+1.0);
  }
}

double logsumexpn(double *l,int n){
  double ret = l[0];

  for(int i=1;i<n;i++){
    ret = logsumexp2(ret,l[i]);
  }

  return ret;
}

double vec_inner(double *x,double *y,int size){
  double sum = 0.0;
  int i;
  for(i = 0;i<size;i++){
    sum += x[i]*y[i];
  }

  return sum;
}

void vec_add(double *a,double *b,double *ret,int size){
  for(int i=0;i<size;i++){
    ret[i] = a[i] + b[i];
  }
}

void vec_sub(double *a,double *b,double *ret,int size){
  for(int i=0;i<size;i++){
    ret[i] = a[i] - b[i];
  }
}

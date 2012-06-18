#include"./mylib.h"
#include<iostream>
#include<ctime>
using namespace std;



const int M = 1;
const int K = 2;
const int n = 20;

double **DATAX;

const double sigma = 1.0;
const double phi = 1.0;
const double beta = 0.01;

typedef unsigned char uchar;

double logZ(uchar *Yn);
double exact_F();
void MCMC_Yn(uchar **ret,const int NUM,const int BURN_IN,const int INTERVAL);
double gibbs_output(uchar **Yn,double **mu,double *pi,const int num);

int main(){
  //preprocess
  dsfmt_init_gen_rand(&dsfmt,(unsigned long)time(NULL));
  
  DATAX = new double*[n];
  for(int i=0;i<n;i++){
    DATAX[i] = new double[M];
  }

  for(int i=0;i<n;i++){
    for(int j=0;j<M;j++){
      cin>>DATAX[i][j];
    }
  }
  ////////////

  double exactf = exact_F();


  //cout<<"MCMC\n";
  int num = 5000;
  int burn_in = 5000;
  int interval = 100;
  uchar **gibbs_out = new uchar*[num];
  for(int i=0;i<num;i++){
    gibbs_out[i] = new uchar[n];
  }
  MCMC_Yn(gibbs_out,num,burn_in,interval);



  double **mu = new double*[K];
  for(int k=0;k<K;k++){
    mu[k] = new double[M];
  }
  double *pi = new double[K];

  mu[0][0] = 0;

  double f;
  for(double r = -3.0;r<6.0;r=r+0.1){
    for(double p=0.025;p<1.0;p=p+0.025){
      mu[1][0] = r;
      pi[0] = p;
      pi[1] = 1-p;

      f = gibbs_output(gibbs_out,mu,pi,num);

      cout<<r<<" "<<p<<" "<<(f-exactf)/exactf<<endl;
    }
    cout<<"\n";
  }


  //postprocess
  for(int i=0;i<n;i++){
    delete [] DATAX[i];
  }
  delete [] DATAX;
  for(int i=0;i<num;i++){
    delete [] gibbs_out[i];
  }
  delete [] gibbs_out;
  for(int k=0;k<K;k++){
    delete [] mu[k];
  }
  delete [] mu;
  delete [] pi;
  

  return 0;
}

double logZ(uchar *Yn){
  //pre process
  double *Y = new double[K];
  double *S = new double[K];
  double **T = new double*[K];
  double *U = new double[K];

  for(int k=0;k<K;k++){
    Y[k] = 0.0;
    S[k] = beta;
    U[k] = 0.0;

    T[k] = new double[M];
    for(int i=0;i<M;i++){
      T[k][i] = 0.0;
    }
  }
  /////////////

  for(int i=0;i<n;i++){
    int idx = (int)Yn[i];

    Y[idx] += 1.0; 
    S[idx] += 1.0/(sigma*sigma);
    vec_add(DATAX[i],T[idx],T[idx],M);
    U[idx] += vec_inner(DATAX[i],DATAX[i],M)/(2.0*sigma*sigma);
  }


  double ret = 0;

  ret += -M*n/2.0*log(2.0*PI*sigma*sigma);
  ret += loggamma(K*phi) - K*loggamma(phi) - loggamma(n+K*phi);
  for(int k=0;k<K;k++){
    ret += loggamma(Y[k]+phi);
    ret += M/2.0*log(beta/S[k]);
    ret += vec_inner(T[k],T[k],M)/(2*S[k]) - U[k];
  }
  

  //post process
  delete [] Y;
  delete [] S;
  delete [] U;
  for(int k=0;k<K;k++){
    delete [] T[k];
  }
  delete [] T;

  return ret;
}


int next_Yn(uchar *Yn){
  /*
    (K-1)(K-1)(K-1)...(K-1)から000...0になるタイミングで0を返す
   */

  for(int i=0;i<n;i++){
    if(Yn[i] < K-1){
      Yn[i] += 1;
      return 1;
    }
    else{
      Yn[i] = 0;
    }
  }

  return 0;
}

double exact_F(){
  uchar *Yn = new uchar[n];
  for(int i=0;i<n;i++){
    Yn[i] = 0;
  }

  double ret = logZ(Yn);

  while(next_Yn(Yn) == 1){
    ret = logsumexp2(ret,logZ(Yn));
  }

  delete [] Yn;
  return -ret;
}


void next_by_gibbs(uchar *Yn){
  double *p = new double[K];


  for(int i=0;i<n;i++){
    for(int k=0;k<K;k++){
      Yn[i] = k;
      p[k] = logZ(Yn);
    }

    double logsum = logsumexpn(p,K);
    for(int k=0;k<K;k++){
      p[k] = exp(p[k] - logsum);
    }

    double r = rand_double_01();
    double pp = 0;
    for(int k=0;k<K;k++){
      pp += p[k];
      if(r < pp){
	Yn[i] = k;
	break;
      }
    }
  }

  delete [] p;
}
void MCMC_Yn(uchar **ret,const int NUM,const int BURN_IN,const int INTERVAL){
  uchar *Yn = new uchar[n];
  for(int i=0;i<n;i++){
    Yn[i] = rand_uint()%K;
  }

  for(int i=0;i<BURN_IN;i++){
    for(int j=0;j<INTERVAL;j++){
      next_by_gibbs(Yn);
    }
  }

  for(int i=0;i<NUM;i++){
    for(int j=0;j<INTERVAL;j++){
      next_by_gibbs(Yn);
    }

    for(int j=0;j<n;j++){
      ret[i][j] = Yn[j];
    }
  }
  
  delete [] Yn;
}


double log_posterior_comp(double **mu,double *pi,uchar *Yn){
  double *Y = new double[K];
  double *U = new double[K];
  double *tmp = new double[M];


  for(int i=0;i<n;i++){
    int idx = (int)Yn[i];

    Y[idx] += 1.0;
    vec_sub(DATAX[i],mu[idx],tmp,M);
    U[idx] += vec_inner(tmp,tmp,M)/(2.0*sigma*sigma);
  }

  double ret = loggamma(K*phi) - K*loggamma(phi);

  for(int k=0;k<K;k++){
    ret += -M*Y[k]/2.0*log(2.0*PI*sigma*sigma);
    ret += -U[k];

    ret += 1.0/2.0*log(beta) - M/2.0*log(2.0*PI);
    ret += -beta/2.0*vec_inner(mu[k],mu[k],M);

    ret += (Y[k]+phi-1.0)*log(pi[k]);
  }

  delete [] Y;
  delete [] U;
  delete [] tmp;

  return ret - logZ(Yn);
}


double log_likelihood(double **mu,double *pi){
  double ret = 0;
  double *tmp_p = new double[K];
  double *tmp_vec = new double[M];

  for(int i=0;i<n;i++){
    for(int k=0;k<K;k++){
      vec_sub(DATAX[i],mu[k],tmp_vec,M);

      tmp_p[k] = log(pi[k]) - M/2.0*log(2.0*PI*sigma*sigma) - vec_inner(tmp_vec,tmp_vec,M)/(2.0*sigma*sigma);
    }
    ret += logsumexpn(tmp_p,K);
  }

  delete [] tmp_p;
  delete [] tmp_vec;

  return ret;
}

double log_prior(double **mu,double *pi){
  double ret = loggamma(K*phi) - K*loggamma(phi);

  for(int k=0;k<K;k++){
    ret += 1.0/2.0*log(beta) - M/2.0*log(2*PI) - beta/2.0*vec_inner(mu[k],mu[k],M) + (phi-1.0)*log(pi[k]);
  }

  return ret;
}

double gibbs_output(uchar **Yn,double **mu,double *pi,const int num){
  double post = log_posterior_comp(mu,pi,Yn[0]);
  for(int i=1;i<num;i++){
    post = logsumexp2(post,log_posterior_comp(mu,pi,Yn[i]));
  }
  return -(log_prior(mu,pi) + log_likelihood(mu,pi) - (post - log(num)));
}


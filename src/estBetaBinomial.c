/*
 * add correction for missing zeros ie  1 - p(0)
 *
 *
 *
 */

#include<R.h>
#include<Rinternals.h>
#include<Rmath.h>
#include<math.h>
#include<stdio.h>

#define epsilon (long double)1e-10
#define MAT_RF(i,j,nx) i+nx*j

SEXP estBetaBinomial(SEXP R_n, SEXP R_x, SEXP R_alphaI, SEXP R_alphaStar, SEXP R_b1, SEXP R_b2);

SEXP estBetaBinomialVariance(SEXP R_n, SEXP R_x, SEXP R_alphaI, SEXP R_alphaStar, SEXP R_b1, SEXP R_b2);

SEXP estBetaBinomialQuantile(SEXP R_n, SEXP R_x, SEXP R_alphaI, SEXP R_alphaStar, SEXP R_b1, SEXP R_b2);

double SumDbetabinom( int j,  long double n, int x, long double alphaI, long double alphaStar, long double b1, long double b2);

long double calc_phi(long double j, long double b1, long double b2);

long double ld_beta_bin(long double n, long double x, long double alphaI, long double alphaStar); // long double beta binomial log likelihood

long double ld_bin(long double n, long double x, long double alphaI, long double alphaStar);// long double binomial log likelihood

long double T_ld_beta_bin(long double n, long double x, long double alphaI, long double alphaStar);// long double beta binomial log likelihood truncated 0

long double T_ld_bin(long double n, long double x, long double alphaI, long double alphaStar);// long double binomial likelihood log truncated 0

SEXP estllTrunc(SEXP R_Y, SEXP R_od, SEXP R_tau, SEXP R_X, SEXP R_offset, SEXP R_N,  SEXP R_dist);


/*
 *
*/

SEXP estBetaBinomial(SEXP R_n, SEXP R_x, SEXP R_alphaI, SEXP R_alphaStar, SEXP R_b1, SEXP R_b2){

  double *t_n,*t_x,*t_alphaI,*t_alphaStar, *t_ans, *b1, *b2;
  long double n, x, alphaI, alphaStar, phi;
  long double sum=0,t1,t2;
  int j,s,trunc;
  SEXP ans;
  t_n = REAL(R_n);
  t_x = REAL(R_x);
  t_alphaI = REAL(R_alphaI);
  t_alphaStar = REAL(R_alphaStar);
  b1 = REAL(R_b1);
  b2 = REAL(R_b2);
  s = LENGTH(R_n);

  // convert all double vector to long double vectors


  ans=allocVector(REALSXP, s);
  t_ans=REAL(ans);
  //*t_ans= (double) sum;
  /*
   * as phi becomes large the overdispersion tends to 1
   * yeilding the standard binomial distribution
   * var(y) = np(1-p)(1+(n-1)/(alphaI+alphaStar+1))
   * alphaI + alphaStar > (n-1-epsilon)/epsilon  | goes to binomial distribution 
   * 
   * epsilon (cutoff point) is set to 1e-15 , the limit of doubles
   */
  for(j=0;j<s;j++){
    phi = calc_phi( (long double)(j+1), (long double) *b1, (long double) *b2);
    n=(long double)t_n[j];
    x=(long double)t_x[j];
    alphaI = ( (long double)t_alphaI[j] ) * phi;
    alphaStar= ( (long double)t_alphaStar[j] ) * phi;
      if(alphaI + alphaStar > (n-1-epsilon)/epsilon){
      	//t_ans[j] = (double)T_ld_bin( n, x, (long double)t_alphaI[j],(long double)t_alphaStar[j-1]);
	t_ans[j] = (double)ld_bin( n, x, (long double)t_alphaI[j],(long double)t_alphaStar[j-1]);
      }
      else{
	// Rprintf("%Lf, %Lf, %Lf\n", alphaI, alphaStar,phi);
	//t_ans[j] = (double)T_ld_beta_bin( n, x, alphaI, alphaStar);
       t_ans[j] = (double)ld_beta_bin( n, x, alphaI, alphaStar);
       
      }
  }
  //printf("\n");
  return(ans);
}

long double calc_phi(long double j, long double b1, long double b2){
  // return( 1/expl(b1 + j*b2) );
  //return( 1/expl(b1 + j*j*b2) );;
  //return( 1/expl(b1 + logl(j)*b2) );
  return(1/expl(logl(j)*b2));
   //return( expl(b1 + logl(j)*b2)/(1+expl(b1 + logl(j)*b2)) );
}


SEXP estBetaBinomialVariance(SEXP R_n, SEXP R_x, SEXP R_alphaI, SEXP R_alphaStar, SEXP R_b1, SEXP R_b2){
  double *t_n,*t_x,*t_alphaI,*t_alphaStar, *t_ans, *b1, *b2;
  long double n, x, alphaI, alphaStar, phi;
  int s,j;
  SEXP ans;

  t_n = REAL(R_n);
  t_x = REAL(R_x);
  t_alphaI = REAL(R_alphaI);
  t_alphaStar = REAL(R_alphaStar);
  b1 = REAL(R_b1);
  b2 = REAL(R_b2);
  s = LENGTH(R_x);
 
  ans=allocVector(REALSXP, s);
  t_ans=REAL(ans);
  n=(long double)t_n[0];
  // alpha star = sum alphaI - alphaI
  for(j=0;j<s;j++){
    phi = calc_phi( (long double)(j+1), (long double) *b1, (long double) *b2);

    x=(long double)t_x[j];
    alphaI = ( (long double)t_alphaI[j] ) * phi;
    alphaStar= ( (long double)t_alphaStar[j] ) * phi;
    if(alphaI + alphaStar > (n-1-epsilon)/epsilon){
      t_ans[j] =  (double) (n * (long double)t_alphaI[j]  / (long double)t_alphaStar[j] * (1 - (long double)t_alphaI[j]  / (long double)t_alphaStar[j] ) );
    }
    else{
      t_ans[j] = (double)(n * alphaI *(alphaStar)*(n + alphaStar + alphaI)/((alphaStar+ alphaI)*(alphaStar+ alphaI))/(1+alphaStar+ alphaI));

    }
  }
   return(ans);
}


/*


 */

SEXP estBetaBinomialQuantile(SEXP R_n, SEXP R_x, SEXP R_alphaI, SEXP R_alphaStar, SEXP R_b1, SEXP R_b2){

  double *t_n,*t_x,*t_alphaI,*t_alphaStar, *t_ans, *b1, *b2;
  long double n, x, alphaI, alphaStar, phi;
  int s,j;
  SEXP ans;

  t_n = REAL(R_n); //single
  t_x = REAL(R_x); //vector
  t_alphaI = REAL(R_alphaI); //vector
  t_alphaStar = REAL(R_alphaStar); //vector
  b1 = REAL(R_b1); // single
  b2 = REAL(R_b2); // single
  s = LENGTH(R_x);
 
  ans=allocVector(REALSXP, s);
  t_ans=REAL(ans);
  //n=(long double)t_n[j];
  n=(long double) *t_n;
  for(j=0;j<s;j++){
    x=(long double)t_x[j];
    alphaI = ( (long double)t_alphaI[j] );
    //alphaStar = sum(alphaI) - alphaI
    alphaStar= ( (long double)t_alphaStar[j] );

    //t_ans[j] = qnorm(SumDbetabinom(j, n, (int)x, alphaI, alphaStar, (long double) *b1, (long double) *b2),0,1,1,0);
    t_ans[j] = qnorm(SumDbetabinom(j+1, n, (int)x, alphaI, alphaStar, (long double) *b1, (long double) *b2),0,1,1,0);
  }

  //  Rprintf("\n");

  return(ans);
}


/*


 */
double SumDbetabinom( int j,  long double n, int x, long double alphaI, long double alphaStar, long double b1, long double b2){ // single values of x only
  long double phi, sum=0;
  int i;
  phi = calc_phi( (long double)j, b1, b2);

  if(alphaI*phi + alphaStar*phi > (n-1-epsilon)/epsilon){
    
    for(i=0;i<=x;i++){ 
      //sum += expl(T_ld_bin(n,(long double)i,alphaI,alphaStar) );
      if(i==x){sum += expl(ld_bin(n,(long double)i,alphaI,alphaStar) )/2;}
      else{sum += expl(ld_bin(n,(long double)i,alphaI,alphaStar) );}

    }

  }
  else{
    alphaI*=phi;
    alphaStar*=phi;
    // Rprintf("%Lf,%Lf,%Lf | ",phi,alphaI,alphaStar);
    for(i=0;i<=x;i++){ 
      if(i==x){sum += expl(ld_beta_bin(n,(long double)i,alphaI,alphaStar))/2;}
      else{sum += expl(ld_beta_bin(n,(long double)i,alphaI,alphaStar));}
      //sum += expl(T_ld_beta_bin(n,(long double)i,alphaI,alphaStar));

    }
  }

  return((double)sum);
}

/*
 * distribution functions
 *
 */
long double ld_beta_bin(long double n, long double x, long double alphaI, long double alphaStar){ // long double beta binomial log likelihood

  return(( (alphaI-0.5)*logl(1+x/alphaI) + (alphaStar-0.5)*logl(1+(n-x)/alphaStar) - (alphaI+alphaStar-0.5)*logl(1+n/(alphaI+alphaStar)) + x*logl(x+alphaI) + (n-x)*logl(n-x+alphaStar) - n*logl(n+alphaI+alphaStar) ) + (long double)lchoose(n,x));

}


long double ld_bin(long double n, long double x,long double alphaI, long double alphaStar){// long double binomial log likelihood

  return( (long double)lchoose(n,x) + ( x * logl( alphaI / alphaStar ) + (n-x)*logl(1 - alphaI  / alphaStar ) ) );
}

long double T_ld_beta_bin(long double n, long double x, long double alphaI, long double alphaStar){// long double beta binomial log likelihood truncated 0
  long double lbb=0,trunc0=0;
  int i;

  lbb=ld_beta_bin(n,x,alphaI,alphaStar);
  //for(i=1;i<=n;i++) trunc0+=expl(ld_beta_bin(n,i,alphaI,alphaStar));
  trunc0=expl(ld_beta_bin(n,0,alphaI,alphaStar));
  //lbb=lbb-logl(1-expl(ld_beta_bin(n,0,alphaI,alphaStar)));
  lbb=lbb-logl(1-trunc0);
  return(lbb);

}

long double T_ld_bin(long double n, long double x, long double alphaI, long double alphaStar){// long double binomial log likelihood truncated 0
  long double lb=0;

  lb=ld_bin(n,x,alphaI,alphaStar);
  
  lb=lb-logl(1-expl(ld_bin(n,0,alphaI,alphaStar)));
  
  return(lb);

}

SEXP estllTrunc(SEXP R_Y, SEXP R_od, SEXP R_tau, SEXP R_X, SEXP R_offset, SEXP R_N,  SEXP R_dist){

  double *lp,*summy, *prob, *lnb;
  double *offset, *tau, *od, *X, *N, *y;
  int x_r,x_v,i,j,n;
  int *dist;
  SEXP dimensions;
  SEXP R_ans;
  double *ans;

  dimensions=getAttrib(R_X, R_DimSymbol);
  x_r=INTEGER(dimensions)[0]; // number of observations 
  x_v=INTEGER(dimensions)[1]; // number of covariates

  od = REAL(R_od);
  X = REAL(R_X);
  tau= REAL(R_tau);
  offset = REAL(R_offset);  
  dist= INTEGER(R_dist); //set dist 1 == NB 2 == poisson;
  N= REAL(R_N);
  y = REAL(R_Y);

  lp = (double *)R_alloc(x_r,sizeof(double));
  prob = (double *)R_alloc(x_r,sizeof(double));
  summy = (double *)R_alloc(x_r,sizeof(double));
  lnb = (double *)R_alloc(x_r,sizeof(double));

  R_ans = allocVector(REALSXP, 1);
  ans = REAL(R_ans);

  for(i=0;i<x_r;i++){
    lp[i]=0;
    for(j=0;j<x_v;j++){
      lp[i]+= tau[j]*X[MAT_RF(i,j,x_r)];
     } 
    lp[i] = exp(lp[i]+offset[i]);
    prob[i] = *od/(*od+lp[i]); //define prob from dnbinom; lp is mean
  }  
    
  if ( *dist == 1) {
    for(i=0;i<x_r;i++){
      lnb[i] = dnbinom( y[i], *od, prob[i], 1);
       
      summy[i] = 0;
      for(n=0;n<=N[i];n++){ summy[i] += dnbinom(n, *od,prob[i],0);}
      // Rprintf("%f,",summy[i]);
    }
    //Rprintf("\n");
  }

  if ( *dist == 2) {
    //cut out poisson for now
    /*    lnb <- dpois( y, lam, log=TRUE)
    for( ii in 1:nrow( X)){
      summy[ii] <- sum( dpois( 0:t[ii], lam[ii], log=FALSE))
      }*/
  }
  *ans = 0;
  for(i=0;i<x_r;i++) *ans += lnb[i]-log(summy[i]);

  return(R_ans);
}

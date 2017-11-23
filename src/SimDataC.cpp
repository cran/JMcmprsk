#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>

#include <string>
#include <iostream>
#include <fstream>
using namespace Rcpp;
static void MulM(const gsl_matrix *XX, const gsl_vector *X, gsl_vector *beta) /******** ####### MulM function ############# ********/
{
    int p = XX->size1;
    int q = XX->size2;

    int i,j;
    double temp;

    for(i=0;i<p;i++)
    {
        temp=0;
        for(j=0;j<q;j++)  temp+=gsl_matrix_get(XX,i,j)*gsl_vector_get(X,j);
        gsl_vector_set(beta,i,temp);
    }

}
   
static int GetN(double t)   /******** ####### GetN function ############# ********/
{

    return (int)(t/0.5+1);

}

static double Min(const double t1, const double t2)  /******** ####### Min function ############# ********/
{
    if(t1<t2) return t1;
    else return t2;
}


// [[Rcpp::export]]
void  SimDataC(SEXP k_val,SEXP p1_val,SEXP p2_val, SEXP g_val, SEXP yfn, SEXP cfn, SEXP mfn)
{	

  size_t k=as<int> (k_val);
  size_t p1=as<int> (p1_val);
  size_t p2=as<int> (p2_val);
  size_t g=as<int> (g_val);
  size_t p1a=1;
  size_t n=0,n1,array_size=k*20;
  size_t i,j;
  
  std::string yfile=as<std::string>(yfn);
  std::string cfile=as<std::string>(cfn);
  std::string mfile=as<std::string>(mfn);
  
  //  int k=500,n,n1,p1=4,p2=2,g=2,p1a=1,,i,j;  /** ######## p1a = 1, array_size=k*20 ######  ****/

    /* allocate space for data */
    gsl_matrix *C = gsl_matrix_alloc(k,p2+2);
    gsl_matrix *FY= gsl_matrix_alloc(array_size, p1+1);          /* initially assign a size for Y, call it fake Y ****/
    gsl_vector *M1= gsl_vector_alloc(k);
    gsl_matrix_set_zero(FY);


    gsl_matrix *VC= gsl_matrix_alloc(p1a+1,p1a+1);
    gsl_vector *RA= gsl_vector_alloc(p1a+1);
    gsl_vector *RI= gsl_vector_alloc(p1a+1);
    gsl_vector *S=gsl_vector_alloc(p1a+1);
    gsl_matrix *V=gsl_matrix_alloc(p1a+1,p1a+1);
    gsl_vector *W=gsl_vector_alloc(p1a+1);

    /* allocate space for true parameters */
    gsl_vector *tbeta = gsl_vector_alloc(p1),
               *tvee = gsl_vector_alloc(g-1);
    gsl_matrix *tgamma = gsl_matrix_alloc(g,p2);



    /* assign the true value to parameters */

    double tsigma=5, tsigmab0=0.5, tsigmau=0.5, rho=0.9, sigmax=1, prob=0.5;

    double tsigmab0u=sqrt(tsigmab0*tsigmau)*rho;

    gsl_vector_set(tbeta,0,10);
    gsl_vector_set(tbeta,1,1);
    gsl_vector_set(tbeta,2,-1.5);
    gsl_vector_set(tbeta,3,-0.6);

    gsl_matrix_set(tgamma,0,0,0.8);
    gsl_matrix_set(tgamma,0,1,-1);
    gsl_matrix_set(tgamma,1,0,0.5);
    gsl_matrix_set(tgamma,1,1,-1.5);

    gsl_vector_set(tvee,0,0.5);

    gsl_matrix_set(VC,0,0,tsigmab0);
    gsl_matrix_set(VC,1,1,tsigmau);
    gsl_matrix_set(VC,1,0,tsigmab0u); 

    for(i=0;i<p1a+1;i++)
    {
        for(j=i+1;j<p1a+1;j++)    gsl_matrix_set(VC,i,j,gsl_matrix_get(VC,j,i));

    }


    gsl_linalg_SV_decomp(VC,V,S,W);

    gsl_matrix_set_zero(V);
    for(i=0;i<p1a+1;i++) gsl_matrix_set(V,i,i,sqrt(gsl_vector_get(S,i)));


    double lamda01=0.1, lamda02=0.15;                          /* baseline hazard; constant over time */
    double temp,x1,x2,t1,t2,censor,u,b0,error;



    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    int point_loc;
    double crate,rate1,rate2,max_censor=5;


    crate=0;
    rate1=0;
    rate2=0;
    point_loc=0;

    for(j=0;j<k;j++)
    {

        for(i=0;i<p1a+1;i++)    gsl_vector_set(RI,i,gsl_ran_gaussian(r,1));
        MulM(V,RI,RA);
        MulM(VC,RA,RI);
  
        b0=gsl_vector_get(RI,0);
        u=gsl_vector_get(RI,1);      

        x1=gsl_ran_gaussian(r,sqrt(sigmax))+2;
        x2=gsl_ran_bernoulli(r,prob);
        censor=gsl_ran_exponential(r, 25);

        gsl_matrix_set(C,j,2,x1);
        gsl_matrix_set(C,j,3,x2);
         
        temp=gsl_matrix_get(tgamma,0,0)*x1+gsl_matrix_get(tgamma,0,1)*x2;
        temp=exp(temp+u);
        temp=lamda01*temp;
        t1=gsl_ran_exponential(r, 1/temp);

        temp=gsl_matrix_get(tgamma,1,0)*x1+gsl_matrix_get(tgamma,1,1)*x2;
        temp=exp(temp+gsl_vector_get(tvee,0)*u);
        temp=lamda02*temp;
        t2=gsl_ran_exponential(r, 1/temp);


        if(t1<=t2 && t1<=censor && t1<=max_censor)
        {
            rate1+=1;
            gsl_matrix_set(C,j,0,t1);
            gsl_matrix_set(C,j,1,1);
            n=GetN(t1);               
        }

        if(t2<=t1 && t2<=censor && t2<=max_censor)
        {
            rate2+=1;
            gsl_matrix_set(C,j,0,t2);
            gsl_matrix_set(C,j,1,2);
            n=GetN(t2);               
        }

        if(Min(censor, max_censor)<=t1 && Min(censor, max_censor)<=t2)
        {
            crate+=1;
            gsl_matrix_set(C,j,0,Min(max_censor,censor));
            gsl_matrix_set(C,j,1,0);
            n=GetN(Min(max_censor,censor));
        }  


        for(i=point_loc;i<point_loc+n;i++)
        {
            gsl_matrix_set(FY,i,1,1.0);
            gsl_matrix_set(FY,i,2,(double)(i-point_loc)*0.5);
            gsl_matrix_set(FY,i,3,x2);
            gsl_matrix_set(FY,i,4,(double)(i-point_loc)*0.5*x2);


            error=gsl_ran_gaussian(r,sqrt(tsigma));
            
            temp=gsl_vector_get(tbeta,0)+gsl_vector_get(tbeta,1)*gsl_matrix_get(FY,i,2)+gsl_vector_get(tbeta,2)*x2
                 +gsl_vector_get(tbeta,3)*x2*gsl_matrix_get(FY,i,2);
            temp=temp+error+b0*gsl_matrix_get(FY,i,2);
            gsl_matrix_set(FY,i,0,temp);
        }

        point_loc+=n;
        gsl_vector_set(M1,j,(double)n);

    }

    n1 = point_loc;
    gsl_matrix *Y = gsl_matrix_alloc(n1,p1+p1a+1);      /* #### now Y's dimension is consistent with the jmc function input format ####### **/
   
    for(i=0;i<n1;i++)                                   /* #### this part is newly added to make Y in the correct format ####### **/
    {
        gsl_matrix_set(Y,i,0,gsl_matrix_get(FY,i,0));
        gsl_matrix_set(Y,i,1,gsl_matrix_get(FY,i,2));
        for(j=0;j<p1;j++)  gsl_matrix_set(Y,i,j+2,gsl_matrix_get(FY,i,j+1));
    }





    /***** ###### output data ############## *******/

  //write Y file
  
  FILE *output_F; 
  
  output_F=fopen(yfile.c_str(),"w");
  if (output_F == NULL)
  {printf("Can't write y file\n");
  }
  
    for(i=0;i<Y->size1;i++)
    {
        for(j=0;j<Y->size2;j++)   fprintf(output_F,"%f   ", gsl_matrix_get(Y,i,j)); 
        if(i!=Y->size1-1)  fprintf(output_F,"\n");                      
    }
 //close writing files  
  fclose(output_F);
  
  //write C file
  
  output_F=fopen(cfile.c_str(),"w");
  if (output_F == NULL)
  {printf("Can't write C file\n");
    //return R_NilValue;
  }
  
    for(i=0;i<C->size1;i++)
    {
        for(j=0;j<C->size2;j++)   fprintf(output_F,"%f   ", gsl_matrix_get(C,i,j)); 
        if(i!=C->size1-1)  fprintf(output_F,"\n");                      
    }
  //close writing files  
  fclose(output_F);
  
  //write M file
  
  output_F=fopen(mfile.c_str(),"w");
  if (output_F == NULL)
  {printf("Can't write M file\n");
    // return R_NilValue;
  }
  
 for(i=0;i<M1->size;i++)
    {
        fprintf(output_F,"%f", gsl_vector_get(M1,i)); 
        if(i!=M1->size-1)  fprintf(output_F,"\n");                      
    }
	
   fclose(output_F);


    gsl_matrix_free(C);
    gsl_matrix_free(FY);          
    gsl_vector_free(M1);

    gsl_matrix_free(VC);
    gsl_vector_free(RA);
    gsl_vector_free(RI);
    gsl_vector_free(S);
    gsl_matrix_free(V);
    gsl_vector_free(W);

    gsl_vector_free(tbeta);
    gsl_vector_free(tvee);
    gsl_matrix_free(tgamma);
	
	
    printf("censoring rate=%f\n",crate/(double)k);
    printf("rate1=%f\n",rate1/(double)k);
    printf("rate2=%f\n",rate2/(double)k);



}


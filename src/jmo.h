/*** joint analysis of ordinal repeated measures and Competing risks (Prentice's cause-specific hazard functions) ****/
/*** use proportional odds model ***/
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
#include <gsl/gsl_blas.h>



#include <Rcpp.h>
using namespace Rcpp;

namespace jmospace {

static const double xs[] = {
    1.94840741569399326713e-01,    5.84978765435932448449e-01,
    9.76500463589682838499e-01,    1.37037641095287183817e+00,
    1.76765410946320160465e+00,    2.16949918360611217335e+00,
    2.57724953773231745414e+00,    2.99249082500237420621e+00,
    3.41716749281857073593e+00,    3.85375548547144464390e+00,
    4.30554795335119844506e+00,    4.77716450350259639289e+00,
    5.27555098651588012760e+00,    5.81222594951591383294e+00,
    6.40949814926966041214e+00,    7.12581390983072757292e+00
};

static const double ws[] = {
    3.75238352592802392864e-01,    2.77458142302529898131e-01,
    1.51269734076642482578e-01,    6.04581309559126141860e-02,
    1.75534288315734303030e-02,    3.65489032665442807915e-03,
    5.36268365527972045989e-04,    5.41658406181998255789e-05,
    3.65058512956237605727e-06,    1.57416779254559402923e-07,
    4.09883216477089661816e-09,    5.93329146339663861478e-11,
    4.21501021132644757306e-13,    1.19734401709284866582e-15,
    9.23173653651829223381e-19,    7.31067642738416239302e-23
};


/*

 static const double xs[] = {
    2.73481046138152452172e-01,    8.22951449144655892596e-01,
    1.38025853919888079639e+00,    1.95178799091625397740e+00,
    2.54620215784748136221e+00,    3.17699916197995602682e+00,
    3.86944790486012269869e+00,    4.68873893930581836465e+00
};

 static const double ws[] = {
    5.07929479016613741923e-01,    2.80647458528533675357e-01,
    8.38100413989858294132e-02,    1.28803115355099736832e-02,
    9.32284008624180529895e-04,    2.71186009253788151199e-05,
    2.32098084486521065344e-07,    2.65480747401118224476e-10
};


*/

/*

 static const double xs[] = {
    2.45340708300901249903e-01,    7.37473728545394358719e-01,
    1.23407621539532300786e+00,    1.73853771211658620678e+00,
    2.25497400208927552311e+00,    2.78880605842813048055e+00,
    3.34785456738321632688e+00,    3.94476404011562521040e+00,
    4.60368244955074427298e+00,    5.38748089001123286199e+00
};

 static const double ws[] = {
    4.62243669600610089640e-01,    2.86675505362834129720e-01,
    1.09017206020023320014e-01,    2.48105208874636108814e-02,
    3.24377334223786183217e-03,    2.28338636016353967260e-04,
    7.80255647853206369398e-06,    1.08606937076928169398e-07,
    4.39934099227318055366e-10,    2.22939364553415129254e-13
};



static  const double xs[] = {
    2.08067382690736869143e-01,    6.24836719505209227388e-01,
    1.04353527375420827481e+00,    1.46553726345740918961e+00,
    1.89236049683768534600e+00,    2.32574984265644102172e+00,
    2.76779535291359376332e+00,    3.22111207656145554770e+00,
    3.68913423846167949074e+00,    4.17663674212926835927e+00,
    4.69075652394311792522e+00,    5.24328537320293598147e+00,
    5.85701464138285064504e+00,    6.59160544236774254391e+00
};

 static const double ws[] = {
    3.98604717826451445875e-01,    2.82561391259388735122e-01,
    1.41394609786954930701e-01,    4.95148892898981398798e-02,
    1.19684232143548187864e-02,    1.95733129440898961560e-03,
    2.10618100024032330920e-04,    1.43455042297144153131e-05,
    5.85771972099298173396e-07,    1.32568250154170781520e-08,
    1.47585316827768861419e-10,    6.63943671490966363565e-13,
    8.31593795120682988192e-16,    1.14013934790367616383e-19
};



 static const double xs[] = {
    1.94840741569399326713e-01,    5.84978765435932448449e-01,
    9.76500463589682838499e-01,    1.37037641095287183817e+00,
    1.76765410946320160465e+00,    2.16949918360611217335e+00,
    2.57724953773231745414e+00,    2.99249082500237420621e+00,
    3.41716749281857073593e+00,    3.85375548547144464390e+00,
    4.30554795335119844506e+00,    4.77716450350259639289e+00,
    5.27555098651588012760e+00,    5.81222594951591383294e+00,
    6.40949814926966041214e+00,    7.12581390983072757292e+00
};

 static const double ws[] = {
    3.75238352592802392864e-01,    2.77458142302529898131e-01,
    1.51269734076642482578e-01,    6.04581309559126141860e-02,
    1.75534288315734303030e-02,    3.65489032665442807915e-03,
    5.36268365527972045989e-04,    5.41658406181998255789e-05,
    3.65058512956237605727e-06,    1.57416779254559402923e-07,
    4.09883216477089661816e-09,    5.93329146339663861478e-11,
    4.21501021132644757306e-13,    1.19734401709284866582e-15,
    9.23173653651829223381e-19,    7.31067642738416239302e-23
};


*/



 double HAZ(const gsl_matrix *H, const double t);

 double CH(const gsl_matrix *H, const double t);

 double MulVV(const gsl_vector *Z,const gsl_vector *beta);

 void MulV(const gsl_vector *Z,gsl_matrix *ZZ);

 void MulM(const gsl_matrix *XX, const gsl_vector *X, gsl_vector *beta);

 void MulMM(const gsl_matrix *A,const gsl_matrix *B,gsl_matrix *AB);

 int inv_matrix(gsl_matrix *x_square);

 double Abs(const double a, const double b);

 int DiffM1(const gsl_matrix *matrixa, const gsl_matrix *matrixb);

 double DiffM(const gsl_matrix *matrixa, const gsl_matrix *matrixb);

 double DiffV(const gsl_vector *veca, const gsl_vector *vecb);

 double Min(const double t1, const double t2);

 void STAT(gsl_matrix *store,int i,double *mean,double *sd);

 int GetN(double t);

 int DiffM2(const gsl_matrix *preH1,const gsl_matrix *H1,const gsl_matrix *preH2,const gsl_matrix *H2);

 void TransM(const gsl_matrix *A, gsl_matrix *B);

 int Sbeta(gsl_vector *beta, double *sigma, const gsl_matrix *Y, const int p1a);



 int EM(
       gsl_vector *beta,
       gsl_matrix *beta2,
       gsl_vector *theta,
       gsl_matrix *gamma,
       gsl_vector *vee,
       gsl_matrix *H01,
       gsl_matrix *H02,
       gsl_matrix *sig,
       const gsl_matrix *Y,
       const gsl_matrix *C,
       const gsl_vector *M1,
       const int p1a,
       const int bq,
	   const int K_num,
	   const int j_max,
	   const int point
       );


 int GetCov(
           gsl_matrix *Cov,
           const gsl_vector *beta,
           const gsl_matrix *beta2,
           const gsl_vector *theta,
           const gsl_matrix *gamma,
           const gsl_vector *vee,
           const gsl_matrix *H01,
           const gsl_matrix *H02,
           const gsl_matrix *sig,
           const gsl_matrix *Y,
           const gsl_matrix *C,
           const gsl_vector *M1,
           const int p1a,
           const int bq,
		   const int K_num,
		   const int j_max,
		   const int point
           );





 int GetE(
          gsl_vector *FUNU,
          gsl_vector *FUNUS,
          gsl_matrix *FUNB,
          gsl_matrix *FUNBS,
          gsl_matrix *FUNBU,
          gsl_matrix *FUNE,
          gsl_matrix *FUNUSE,
          gsl_matrix *FUNUE,
          gsl_matrix *FUNP,
          gsl_matrix *FUNPS,
          gsl_matrix *FUNPR1,
          gsl_matrix *FUNPR2,
          const gsl_vector *beta,
          const gsl_matrix *beta2,
          const gsl_vector *theta,
          const gsl_matrix *gamma,
          const gsl_vector *vee,
          const gsl_matrix *H01,
          const gsl_matrix *H02,
          const gsl_matrix *sig,
          const gsl_matrix *Y,
          const gsl_matrix *C,
          const gsl_vector *M1,
          const int p1a,
          const int bq,
		   const int K_num,
		   const int j_max,
		   const int point
          );



 double Getloglik(
          const gsl_vector *beta,
          const gsl_matrix *beta2,
          const gsl_vector *theta,
          const gsl_matrix *gamma,
          const gsl_vector *vee,
          const gsl_matrix *H01,
          const gsl_matrix *H02,
          const gsl_matrix *sig,
          const gsl_matrix *Y,
          const gsl_matrix *C,
          const gsl_vector *M1,
          const int p1a,
          const int bq,
		  const int K_num,
		   const int point
          );



 int Diff(
         const gsl_vector *prebeta,
         const gsl_vector *beta,
         const gsl_matrix *prebeta2,
         const gsl_matrix *beta2,
         const gsl_vector *pretheta,
         const gsl_vector *theta,
         const gsl_matrix *pregamma,
         const gsl_matrix *gamma,
         const gsl_vector *prevee,
         const gsl_vector *vee,
         const gsl_matrix *preH01,
         const gsl_matrix *H01,
         const gsl_matrix *preH02,
         const gsl_matrix *H02,
         const gsl_matrix *presig,
         const gsl_matrix *sig
         );

#define array_size 100000
 
		 
//declare before use it
Rcpp::List jmo_cmain(int k, int n1,int p1,int p2, int p1a, int bq,int K_num, int j_max,int point, std::vector<double>beta_val, std::vector<double>theta_val, int maxiterations,std::string yfile, std::string cfile, std::string mfile, int trace);
}

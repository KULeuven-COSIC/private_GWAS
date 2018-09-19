//set 100 centers with input values up to 1000
#include <cstddef>
#include <iostream>
#include "../nfllib/include/nfl.hpp"
#include <thread>
#include <vector>

int ti=0;
/// include the FV homomorphic encryption library
namespace FV {
namespace params {
  using poly_t = nfl::poly_from_modulus<uint64_t, 4096, 186>;
  template <typename T>
  struct plaintextModulus;
  template <>
  struct plaintextModulus<mpz_class> {
    static mpz_class value() { return mpz_class(std::to_string(ti)); } 
    };
  using gauss_struct = nfl::gaussian<uint16_t, uint64_t, 2>;
  using gauss_t = nfl::FastGaussianNoise<uint16_t, uint64_t, 2>;
  gauss_t fg_prng_sk(265.0, 80, 1 << 12);
  gauss_t fg_prng_evk(265.0, 80, 1 << 12);
  gauss_t fg_prng_pk(265.0, 80, 1 << 12);
  gauss_t fg_prng_enc(265.0, 80, 1 << 12);
  }
} //namespace FV


#include <FV.hpp>
#include "ConvertBalancedBBase.hpp"
#include "CorrelationProgram.hpp"
#include "CRT.hpp"

FV::params::poly_p ConvertTowNIBNAFPol(long double input, double base, mpz_t mod, double precision)
{
  FV::params::poly_p output;
  const size_t degree=4096;
  std::array<mpz_t,degree> coeff;
  mpz_t rest,sum,sigmampz;
  mpz_inits(rest,sum,sigmampz,nullptr);
  for (size_t i=0; i<degree;i++)
    {
      mpz_inits(coeff[i],nullptr);
    }
  while(abs(input)>precision)
    {
      int sigma = sgn(input);
      long double t = sigma*input;
      int r = ceil(log(t)/log(base));
      if ((pow(base,r)-t)>(t-pow(base,r-1)))
       {
         r=r-1;
       }
      if(r>-1)
      	{ 
          mpz_set_si(sigmampz,sigma);
          mpz_add(sum,sigmampz,mod);
          mpz_mod(rest,sum,mod);
          mpz_set(coeff[r],rest);
	} 
      else 
	{
          mpz_set_si(sigmampz,-sigma);
          mpz_add(sum,sigmampz,mod);
          mpz_mod(rest,sum,mod);
          mpz_set(coeff[degree+r],rest);
	}
      input=input-sigma*pow(base,r);
    }   
  output.mpz2poly(coeff);
// Clean
     for (size_t i = 0; i < degree; i++) 
     {
       mpz_clears(coeff[i],nullptr);
     }
mpz_clears(rest,sum,sigmampz,nullptr);
  return output;
  }

template<class I>
long double ConvertToNumber(long double base, mpz_t mod,int split, I &input)
{
  long double result=0;
  const size_t dim = sizeof(input)/sizeof(input[0]);
  const size_t neg_dim=split;
  mpz_t test,div,two,min;
  mpz_inits(test,div,two,min,nullptr);
  mpz_set_ui(two,2);
  mpz_cdiv_q(div,mod,two);
  for (unsigned long int i=0; i<dim; i++)
    {
      mpz_mod(test,input[i],mod);
	if (mpz_cmp(test,div)<0)
	{
	  mpz_set(input[i],test);
	}
      else
	{
	  mpz_sub(min,test,mod);
	  mpz_set(input[i],min);
	}
      long double coeffi=mpz_get_d(input[i]);
      if (coeffi !=0)
      {
	if(i<neg_dim)
	  {
	    result=result+coeffi*pow(base,i);
	  }
	else
	  {
	    int neg_power=i-dim;
	    result=result-coeffi*pow(base,neg_power);
	  }
      }
    }
mpz_clears(test,div,two,min,nullptr);
  return result;
}

namespace FV{
ciphertext_t poly_to_ciphertext(pk_t &pk, params::poly_p const &poly)
  {
    ciphertext_t ct;
    ct.pk = &pk;
    ct.c0 = poly;
    ct.c0.ntt_pow_phi();
    ct.c0 = nfl::shoup(ct.c0 * ct.pk->delta, ct.pk->delta_shoup);
    ct.isnull = false;

    return ct;
  } 
}

long double randTobound(long double M)
{
    return (random() / ( RAND_MAX / M ) ) ;  
}

void diff(timespec &start, timespec &stop, timespec &difference)
{
  if ((stop.tv_nsec-start.tv_nsec)<0)
    {
      difference.tv_sec = stop.tv_sec - start.tv_sec -1;
      difference.tv_nsec = 1000000000+stop.tv_nsec-start.tv_nsec;
    }
  else
    {
      difference.tv_sec=stop.tv_sec-start.tv_sec;
      difference.tv_nsec = stop.tv_nsec-start.tv_nsec;
    }
}

void divandadd(timespec &start, timespec &stop,int nb_hosp)
{
  start.tv_sec=start.tv_sec+(stop.tv_sec)/nb_hosp;
  start.tv_nsec=start.tv_nsec+(stop.tv_nsec)/nb_hosp;
}

void add(timespec &start, timespec &stop)
{
  start.tv_sec=start.tv_sec+stop.tv_sec;
  start.tv_nsec=start.tv_nsec+stop.tv_nsec;
}

void significance_function(long double input[4][100], const size_t nb_centers, int primes[], size_t nb_primes, int w, long double basew, int power_rand, int power_frac, int split, timespec &time_hosp, timespec &time_comp_server, timespec &time_decr)
{
  const size_t degree=4096;
  mpz_t mpz_primes[nb_primes];
  mpz_t coeff_result[degree][nb_primes];
  for (size_t j=0; j<nb_primes;j++)
  {
    for (size_t i=0; i<degree; i++)
    {
      mpz_inits(coeff_result[i][j],nullptr);
    }  
    mpz_init_set_ui(mpz_primes[j],primes[j]);
  }
  long double M=pow(2,power_rand); 
  int bound=floor(pow(2,power_rand)/powl(basew,7*power_frac));
  long double random_double=randTobound(M);
  long double random_double2=randTobound(bound);
  for (size_t j=0;j<nb_primes;j++)
  {
    ti=primes[j];

     // KeyGen
    FV::sk_t secret_key;

    FV::evk_t evaluation_key(secret_key, 32);

    FV::pk_t public_key(secret_key, evaluation_key);

    FV::params::poly_p mess_hosp[4][nb_centers];
    FV::params::poly_p treshold;
    FV::ciphertext_t c[4][nb_centers];
    FV::ciphertext_t treshold_encr;
    mpz_t modi;
    mpz_inits(modi,nullptr);
    mpz_set_ui(modi,ti);
    timespec start_hosp, end_hosp;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&start_hosp);
    for (size_t j=0; j<4;j++)
      {
      for (size_t i=0;i<nb_centers;i++)
	  {
	    mess_hosp[j][i]=ConvertTowNIBNAFPol(input[j][i],basew,modi,0.1);
	    FV::encrypt_poly(c[j][i], public_key,mess_hosp[j][i]);
	  }
      }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&end_hosp);
    timespec difference;
    diff(start_hosp,end_hosp,difference);
    divandadd(time_hosp,difference,nb_centers);
   treshold=ConvertTowNIBNAFPol(50,basew,modi,0.1);
   treshold_encr=FV::poly_to_ciphertext(public_key, treshold);
   timespec start_comp,end_comp;
   clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&start_comp);
   FV::ciphertext_t total_cum_table[4];
   FV::ciphertext_t temp;
   for (size_t j=0; j<4;j++)
     {
       //LIBRARY CAN NOT ITERATIVELY ADD CIPHERTEXTS AT THE MOMENT
       total_cum_table[j]=c[j][0]+c[j][1]+c[j][2]+c[j][3]+c[j][4]+c[j][5]+c[j][6]+c[j][7]+c[j][8]+c[j][9]+c[j][10]+c[j][11]+c[j][12]+c[j][13]+c[j][14]+c[j][15]+c[j][16]+c[j][17]+c[j][18]+c[j][19]+c[j][20]+c[j][21]+c[j][22]+c[j][23]+c[j][24]+c[j][25]+c[j][26]+c[j][27]+c[j][28]+c[j][29]+c[j][30]+c[j][31]+c[j][32]+c[j][33]+c[j][34]+c[j][35]+c[j][36]+c[j][37]+c[j][38]+c[j][39]+c[j][40]+c[j][41]+c[j][42]+c[j][43]+c[j][44]+c[j][45]+c[j][46]+c[j][47]+c[j][48]+c[j][49]+c[j][50]+c[j][51]+c[j][52]+c[j][53]+c[j][54]+c[j][55]+c[j][56]+c[j][57]+c[j][58]+c[j][59]+c[j][60]+c[j][61]+c[j][62]+c[j][63]+c[j][64]+c[j][65]+c[j][66]+c[j][67]+c[j][68]+c[j][69]+c[j][70]+c[j][71]+c[j][72]+c[j][73]+c[j][74]+c[j][75]+c[j][76]+c[j][77]+c[j][78]+c[j][79]+c[j][80]+c[j][81]+c[j][82]+c[j][83]+c[j][84]+c[j][85]+c[j][86]+c[j][87]+c[j][88]+c[j][89]+c[j][90]+c[j][91]+c[j][92]+c[j][93]+c[j][94]+c[j][95]+c[j][96]+c[j][97]+c[j][98]+c[j][99];
	}
    
    FV::ciphertext_t alpha;
    FV::ciphertext_t beta;

    FV::ciphertext_t R0,R1,C0,C1,N;
    R0=total_cum_table[0]+total_cum_table[1];
    R1=total_cum_table[2]+total_cum_table[3];
    C0=total_cum_table[0]+total_cum_table[2];
    C1=total_cum_table[1]+total_cum_table[3];

    N = total_cum_table[0]+total_cum_table[1]+total_cum_table[2]+total_cum_table[3];
    alpha =(R1*C1)*(N*total_cum_table[0]-R0*C0)*(N*total_cum_table[0]-R0*C0)+(R1*C0)*(N*total_cum_table[1]-R0*C1)*(N*total_cum_table[1]-R0*C1)+(R0*C1)*(N*total_cum_table[2]-R1*C0)*(N*total_cum_table[2]-R1*C0)+(R0*C0)*(N*total_cum_table[3]-R1*C1)*(N*total_cum_table[3]-R1*C1);
     beta = (N*R0)*(R1*C0)*(C1*treshold_encr);
     
     
     FV::params::poly_p random1,random2;
     random1=ConvertTowNIBNAFPol(random_double,basew,modi,0.1);
     random2=ConvertTowNIBNAFPol(random_double2,basew,modi,0.1);
     FV::ciphertext_t c_random1,c_random2;
     c_random1=FV::poly_to_ciphertext(public_key,random1);
     c_random2=FV::poly_to_ciphertext(public_key,random2);

     FV::ciphertext_t result;
     
     result=(alpha-beta)*c_random1+c_random2;
     clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&end_comp);
     timespec comp_difference;
     diff(start_comp, end_comp, comp_difference);
     add(time_comp_server,comp_difference);
     std::array<mpz_t, degree> poly_result;  
     for (size_t i = 0; i <degree; i++) 
     {
       mpz_inits(poly_result[i], nullptr);
     }

     timespec start_decr, end_decr;
     clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&start_decr);
     FV::decrypt_poly(poly_result, secret_key, public_key,result);
     clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&end_decr);
     timespec decr_difference;
     diff(start_decr,end_decr,decr_difference);
     add(time_decr,decr_difference);
     for (size_t i = 0; i <degree; i++) 
     {
       mpz_set(coeff_result[i][j],poly_result[i]);
     }
     // Clean
     for (size_t i = 0; i < degree; i++) 
     {
       mpz_clears(poly_result[i],nullptr);
     }

     mpz_clears(modi,nullptr);
  }

  
  mpz_t coeffN[degree];
  for (size_t k=0;k<degree;k++)
  {
    mpz_inits(coeffN[k],nullptr);
  }
    timespec start_dec, end_dec;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&start_dec);
   for(size_t k=0; k<degree; k++)
   {
     findMinX(coeffN[k],mpz_primes,coeff_result[k],nb_primes);
   }
 

    mpz_t prod;
    mpz_inits(prod,nullptr);
    mpz_set_ui(prod,1);
    for (size_t k=0; k<nb_primes; k++)
    { 
      mpz_mul(prod,prod,mpz_primes[k]);
    }

   long double result_num=ConvertToNumber(basew,prod,split,coeffN);
   clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&end_dec);
   timespec dec_difference;
   diff(start_dec, end_dec,dec_difference);
   add(time_decr,dec_difference);
   std::cout<<"decr_time"<<time_decr.tv_nsec<<"\n";

   std::cout.precision(20);
   std::cout<<"result"<<result_num<<"\n";

  for (size_t i = 0; i<nb_primes;i++)
  {
    for(size_t k=0; k<degree;k++)
     {
        mpz_clears(coeff_result[k][i],nullptr);
     }
  }
  mpz_clears(prod,nullptr);  
  for (size_t i=0; i<nb_primes; i++)
  {
    mpz_clears(mpz_primes[i],nullptr);
  }

  for(size_t k=0; k<degree;k++)
  {
     mpz_clears(coeffN[k],nullptr);
  }
}

int main()
{
  gmp_randstate_t state;
  gmp_randinit_default(state);
  size_t nb_centers=100;
  int max_input1=250;
  long double input1[4][100];
  for (size_t j=0; j<4;j++)
      {
      for (size_t i=0;i<nb_centers;i++)
	  {
	    input1[j][i]=floor(randTobound(max_input1));
	    std::cout<<"input"<<input[j][i]<<"\n";
	  }
      }
  int primes1[3]={707801,707813,707827};
  int w1=139;
  long double basew1=1.03063272124421509141114964188;
  int power_rand1=94;
  int power_frac1=76;
  int split1=3549;
  timespec average_time_hosp1,average_time_comp_server1,average_time_decr1;
  average_time_hosp1.tv_sec=0.0;
  average_time_hosp1.tv_nsec=0.0;
  average_time_comp_server1.tv_sec=0.0;
  average_time_comp_server1.tv_nsec=0.0;
  average_time_decr1.tv_sec=0.0;
  average_time_decr1.tv_nsec=0.0;

   for (size_t i=0; i<10; i++)
    {
      significance_function(input1,nb_centers,primes1,3,w1, basew1,power_rand1,power_frac1,split1,average_time_hosp1,average_time_comp_server1,average_time_decr1);
     }
  std::cout<<"average_time_hosp: "<<average_time_hosp1.tv_nsec/10<<"\n";
  std::cout<<"average_time_comp_server: "<<average_time_comp_server1.tv_nsec/10<<"\n";
  std::cout<<"average_time_decr: "<<average_time_decr1.tv_nsec/10<<"\n";

  int max_input2=2500;
  long double input2[4][100];
  for (size_t j=0; j<4;j++)
      {
      for (size_t i=0;i<nb_centers;i++)
	  {
	    input2[j][i]=floor(randTobound(max_input2));
	  }
      } 
  int primes2[3]={707801,707813,707827};
  int w2=108;
  long double basew2=1.03764767257676312608489462749;
  int power_rand2=114;
  int power_frac2=62;
  int split2=3651;
  timespec average_time_hosp2,average_time_comp_server2,average_time_decr2;
  average_time_hosp2.tv_sec=0.0;
  average_time_hosp2.tv_nsec=0.0;
  average_time_comp_server2.tv_sec=0.0;
  average_time_comp_server2.tv_nsec=0.0;
  average_time_decr2.tv_sec=0.0;
  average_time_decr2.tv_nsec=0.0;
   for (size_t i=0; i<10; i++)
    {
      significance_function(input2,nb_centers,primes2,3,w2, basew2,power_rand2,power_frac2,split2,average_time_hosp2,average_time_comp_server2,average_time_decr2);
     }
   std::cout<<"average_time_hosp: "<<average_time_hosp2.tv_sec<<" sec "<<average_time_hosp2.tv_nsec/10<<"\n";
   std::cout<<"average_time_comp_server: "<<average_time_hosp2.tv_sec<<" sec "<<average_time_comp_server2.tv_nsec/10<<"\n";
   std::cout<<"average_time_decr: "<<average_time_hosp2.tv_sec<<" sec "<<average_time_decr2.tv_nsec/10<<"\n";


  int max_input3=25000;
  long double input3[4][100];
  for (size_t j=0; j<4;j++)
      {
      for (size_t i=0;i<nb_centers;i++)
	  {
	    input3[j][i]=floor(randTobound(max_input3));
	  }
      }
  int primes3[3]={707801,707813,707827};
  int w3=87;
  long double basew3=1.04487642192444582300698356032;
  int power_rand3=134;
  int power_frac3=52;
  int split3=3710;
  timespec average_time_hosp3,average_time_comp_server3,average_time_decr3;
  average_time_hosp3.tv_sec=0.0;
  average_time_hosp3.tv_nsec=0.0;
  average_time_comp_server3.tv_sec=0.0;
  average_time_comp_server3.tv_nsec=0.0;
  average_time_decr3.tv_sec=0.0;
  average_time_decr3.tv_nsec=0.0;
   for (size_t i=0; i<10; i++)
    {
      significance_function(input3,nb_centers,primes3,3,w3, basew3,power_rand3,power_frac3,split3,average_time_hosp3,average_time_comp_server3,average_time_decr3);
     }
  std::cout<<"average_time_hosp: "<<average_time_hosp3.tv_nsec/10<<"\n";
  std::cout<<"average_time_comp_server: "<<average_time_comp_server3.tv_nsec/10<<"\n";
  std::cout<<"average_time_decr: "<<average_time_decr3.tv_nsec/10<<"\n";

  return 0;
}

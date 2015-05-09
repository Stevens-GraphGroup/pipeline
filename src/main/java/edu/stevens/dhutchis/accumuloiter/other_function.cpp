#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header_def.cuh"

/*****************************************************************
 * 2. Standard iid null model ("null1")
 *****************************************************************/

/* Function:  NullOne()
 *
 * Purpose:   Calculate the null1 lod score, for sequence <dsq>
 *            of length <L> "aligned" to the base null model <bg>. 
 * 
 */
int NullOne(int L, float *ret_sc)
{
  float p1 = 0.0f;
  
  p1 = (float) L / (float) (L+1);

  *ret_sc = (float) L * logf(p1) + logf(1.-p1);

  return 1;
}

/* Function:  esl_gumbel_surv()
 * Synopsis:  Returns right tail mass above $x$, $P(S > x)$.
 *
 * Purpose:   Calculates the survivor function, $P(X>x)$ for a Gumbel 
 *            (that is, 1-cdf), the right tail's probability mass.
 * 
 *            Let $y=\lambda(x-\mu)$; for 64-bit doubles, 
 *            useful dynamic range for $y$ is $-3.6 <= y <= 708$.
 *            Returns 1.0 for $y$ below lower limit, and 0.0
 *            for $y$ above upper limit.
 */
double esl_gumbel_surv(double x, double mu, double lambda)
{
  double y  = lambda*(x-mu);
  double ey = -exp(-y);

  /* Use 1-e^x ~ -x approximation here when e^-y is small. */
  if (fabs(ey) < eslSMALLX1) return -ey; 
  else                       return 1 - exp(ey);
}


/* For protein models, default iid background frequencies */
/* Function:  p7_AminoFrequencies()
 *
 * Purpose:   Fills a vector <f> with amino acid background frequencies,
 *            in [A..Y] alphabetic order, same order that Easel digital
 *            alphabet uses. Caller must provide <f> allocated for at
 *            least 20 floats.
 *            
 *            These were updated 4 Sept 2007, from SwissProt 50.8,
 *            (Oct 2006), counting over 85956127 (86.0M) residues.
 *
 * Returns:   <eslOK> on success.
 */
int p7_AminoFrequencies(float *f)
{
  f[0] = 0.0787945;   /* A */
  f[1] = 0.0151600;   /* C */
  f[2] = 0.0535222;   /* D */
  f[3] = 0.0668298;   /* E */
  f[4] = 0.0397062;   /* F */
  f[5] = 0.0695071;   /* G */
  f[6] = 0.0229198;   /* H */
  f[7] = 0.0590092;   /* I */
  f[8] = 0.0594422;   /* K */
  f[9] = 0.0963728;   /* L */
  f[10]= 0.0237718;   /* M */
  f[11]= 0.0414386;   /* N */
  f[12]= 0.0482904;   /* P */
  f[13]= 0.0395639;   /* Q */
  f[14]= 0.0540978;   /* R */
  f[15]= 0.0683364;   /* S */
  f[16]= 0.0540687;   /* T */
  f[17]= 0.0673417;   /* V */
  f[18]= 0.0114135;   /* W */
  f[19]= 0.0304133;   /* Y */
  return 1;
}

/* FREE mem of HMM model */
void freeHMM(HMMER_PROFILE* hmm)
{
  free(hmm->f);         /* 20 standard residues */
  free(hmm->MU);        /* 3 floats */
  free(hmm->LAMBDA);    /* 3 floats */

  //dynamic allocation free (temprary solution)
  free(hmm->tran_16bits);
  free(hmm->mat_16bits);
  free(hmm->mat_8bits);

  int i = 0;
  for(i = 0; i < hmm->M; i++)
    free(hmm->log_tran_32bits[i]);
  for(i = 0; i < (hmm->M + 1); i++)
  {
    free(hmm->tran_32bits[i]);
    free(hmm->ins_32bits[i]);
    free(hmm->mat_32bits[i]);
  }
  free(hmm->log_tran_32bits);
  free(hmm->tran_32bits);
  free(hmm->ins_32bits);
  free(hmm->mat_32bits);

  free(hmm);            /* free 'om', not include above extra space */
                        /* only include their pointer address: 3 float pointes: MU,LAMBDA,f */
}

/* Implement roundf() function
 * Example: roundf(2.1) = 2; roundf(2.6) = 3
 *          roundf(-2.1) = -2; roundf(-2.6) = -3
 *
 * Using existing functions: floor(); ceil();
 *
 * Example: floor(2.6) = 2; ceil(2.1) = 3
 */
int round_DIY(float input)
{ 
  int r;
  r = (input - floor(input) < 0.5 ) ? floor(input) : ceil(input);
  return r;
}

/* Get P-values */
void seq_Pass(int begin, int end, int* len, HMMER_PROFILE* hmm, float* sc,bool* passed_msv)
{
  float* nullsc = NULL;
  float  seq_score = 0.0f;         /* Final score of each seq for obtain P-values   */
  double P_value = 0.0f;           /* Used to filter seqs                           */
  long int n_past_msv = 0;         /* # of seqs pass MSV filter     */
  int span = end - begin;          /* null score for each seq with different length */

  nullsc = (float*)malloc(span * sizeof(float));     /* get Null score for each seq */
  for(int i = 0; i < span; i++)
    NullOne(len[i], &nullsc[i]);              /* len has been measured in void GPUs::get_Real_Length(int* L), it is real_L of private member */

  /* get # of seq pass through filter */
  for(int i = 0; i < span; i++)
  {
    seq_score = (sc[i] - nullsc[i]) / eslCONST_LOG2;
    P_value = esl_gumbel_surv(seq_score, hmm->MU[0], hmm->LAMBDA[0]);
    if(P_value <= F1)
       {
        passed_msv[begin + i] = 1;
       }
       else
       {
        passed_msv[begin + i] = 0;
       }
  }
  
  free(sc);
  free (nullsc);
  nullsc = NULL;
  sc = NULL;
}

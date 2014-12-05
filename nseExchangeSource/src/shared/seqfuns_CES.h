#ifndef CES_SEQFUNS_H_
#define CES_SEQFUNS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>		//gettimeofday()
#include <gsl/gsl_rng.h>	//uniform rng
#include <gsl/gsl_randist.h>	//rng from distns
#include <math.h>
#include <R.h>

#include "CES_struct.h"

/*======================*
 * External definitions *
 *======================*/

extern struct amino_acid AA[22];
extern struct codon_struct Codon[64];
extern struct seq_struct Sequence;
extern struct global_variable_struct G;

/*=====================*
 * Function Prototypes *
 *=====================*/

void Convert_Codon_Index_to_AA_Index(int *codon_index_vec,
				          int *aa_index_vec, int aa_count);
void Generate_Random_Codon_Seq_and_Index (char codon_seq[][4],
				              int *codon_index,
				              int *aa_index, int aa_count);
void Generate_Random_Codon_Index(int *codon_index,
				    int *aa_index, int aa_count);

void Calc_NSE_Odds_Vec(int *codon_index_vec, double *NSE_odds_vec,
		  int aa_count);
void Codon_Counts(int * codon_index_vec ,int * codon_cts_vec ,int aa_count);
void Convert_Codon_Index_to_Sigma_Vec(int *index_vec,
	                                 double *sigma_vec, int aa_count);
double Calc_Xi_NSE (double *sigma_vec, int aa_count);
double Calc_Eta_NSE (struct seq_struct * seq);

double Calc_Seq_Mu(int * codon_counts_vec, int aa_count);

void Convert_Sigma_Vec_to_Sigma_Ratio_Vec (double *sigma_vec,
				      double *sigma_ratio_vec, int aa_count);
void Recalc_Sigma_Vec (double *sigma_vec, int mut_codon_pos, double factor,
		  int aa_count);
void eta_min_max(int *aa_index,int aa_count, double * answer);
 void flip_codons(struct seq_struct *currSeq, struct seq_struct *propSeq, double *randNums, int num_flips);

#endif //CES_SEQFUNS_H_

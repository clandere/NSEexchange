#ifndef CES_SYNFUNS_H_
#define CES_SYNFUNS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>		//gettimeofday()
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

int Calc_Dimensionality(int *codon_index_vec, int aa_count);
int Generate_D_Array(int *codon_index_vec,
		         int *D_array[][2], int aa_count);
void Calc_Delta_Eta_NSE_First_and_Second_Term_Vecs(double *sigma_ratio_vec,
		         double *NSE_odds_vec, double *delta_eta_first_term_vec, 
		         double *delta_eta_second_term_vec, int aa_count);
void Calc_Delta_Eta_NSE_Vec(int *codon_index_vec,
		         double *NSE_odds_vec, double *first_term_vec,
		         double *second_term_vec, double *delta_eta_vec,
		         int aa_count, int D_dim);
void Calc_Pi_Vec(double *delta_eta_vec, double *pi_vec, double *ptr_pi_total,
	             double phi, double qPhi, double two_qNePhi, int aa_count, 
	             int D_dim);
void Calc_Replacement_Pr_Vec(double *mutation_vec, double *pi_vec, double *pr_vec,
			     double *ptr_pr_total, int D_dim);



#endif //CES_SYNFUNS_H_

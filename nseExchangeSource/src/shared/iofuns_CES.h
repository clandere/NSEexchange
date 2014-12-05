#ifndef CES_IOFUNS_H_
#define CES_IOFUNS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
 
 int Read_tRNA_File (char *filename);
 int Process_AA_Information ();
 void Process_R_input(int *CODON_INDEX, double *ELONG_PR, 
             double *MUTATION_RATES, int *AA_COUNT, double *PHI, 
             double *POP, double *A_1, double *A_2,
             double *AT_BIAS, double *BEE,
             int *IGNORE, double *GAMMA,
             double *SCALE_FACTOR, char ** AA_VEC,
             char ** CODON_VEC, int *N_CODON);
 void Rprint_codon_sequence( int *codon_index );
 void Generate_Codon_Structures (); 
 void Calc_AA_Moments();
 int wrong();
 
#endif //CES_IOFUNS_H_

#include "CES.h"

void calc_eta_NSE_R_wrapper (int * CODON_INDEX, int * AA_COUNT,  double * ELONG_PR, 
							 int * N_CODONS,    double * A_1,    double * A_2,
							 double * ETA){
	int i;
	
	//1) Set up sequence structure
	//1.a) codon_index
	for(i=0;i<*AA_COUNT;i++){
		Sequence.codon_index[i] = CODON_INDEX[i];
	}
	//1.b) aa_count
	Sequence.aa_count = *AA_COUNT;
	
	//2) Set up codon elong_pr's
	for(i=0;i<*N_CODONS;i++){
		Codon[i].elong_pr = ELONG_PR[i];
	}
	
	//3) Set up initiation and elongation cost 
	G.A1 = *A_1;
	G.A2 = *A_2;
	
	*ETA = Calc_Eta_NSE(&Sequence);
	
}

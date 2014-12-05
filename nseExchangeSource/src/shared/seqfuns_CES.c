#include "seqfuns_CES.h"

/*======================================================*
 * FUNCTION: Convert_Codon_Index_to_AA_Index			*
 * 		Takes the codon_index_vec and assigns each AA	*
 * 		a unique index (0-20) and stores it in a new 	*
 * 		vector aa_index_vec								*
 * =====================================================*/

void
Convert_Codon_Index_to_AA_Index(int *codon_index_vec,
				int *aa_index_vec, int aa_count)
{
	int i;
	for (i = 0; i < aa_count; i++) {
		aa_index_vec[i] = Codon[codon_index_vec[i]].aa_index;
	}

}

/*======================================================*
 * FUNCTION:Generate_Random_Codon_Seq_and_Index			*
 * 		Takes aa_index and chooses a random codon from	*
 * 		the group of synonyms based on relative 		*
 * 		mutation rates.									*
 * =====================================================*/

void
Generate_Random_Codon_Seq_and_Index (char codon_seq[][4],
				     int *codon_index,
				     int *aa_index,
				     int aa_count) {
      int i, j, k, num_codons, test;
      double rnd;

      //things related to rng
      struct timeval curr_time;

      const gsl_rng_type *rng_type;
      gsl_rng *rn;


      //Set up GSL RNG
      gsl_rng_env_setup ();



      //Create a generator chosen by the 
      //environment variable GSL_RNG_TYPE
      //use default gsl for generating uniform rn 
      //from which other rn functions are derived 

      rng_type = gsl_rng_default;
      rn = gsl_rng_alloc (rng_type);

      //seed rng using clock
      gettimeofday (&curr_time, NULL);

      gsl_rng_set (rn, (curr_time.tv_sec + curr_time.tv_usec));	// (const gsl_rng * r, unsigned long int s)   

      //end RNG set up


      for (i = 1; i < aa_count; i++) {
	      k = aa_index[i];
	      num_codons = AA[k].num_codons;
	      if (num_codons > 1) {
		      rnd = gsl_rng_uniform (rn);
		      if (G.benchmark) rnd = 0.5;
		      test = 0;
		      for (j = 0; j < num_codons && test == 0; j++) {
			      if (rnd < AA[k].neutral_obs_pr_cum[j]){ //Picks a codon based on its neutral probability and takes into account AT bias
			         //copy codon index and seq
			         codon_index[i] = AA[k].codon_index[j];
			         strcpy (codon_seq[i], AA[k].codon[j]);
			         test = 1;	//break;
			      }
		      }
	      }
      }
}

/*======================================================*
 * FUNCTION:Generate_Random_Codon_Index	         		*
 * 		Takes aa_index and chooses a random codon from	*
 * 		the group of synonyms based on relative 		*
 * 		mutation rates.	Used to be 						*
 *      Generate_Random_Codon_Seq_and_Index				*
 * =====================================================*/

void
Generate_Random_Codon_Index(int *codon_index,
				    int *aa_index, int aa_count)
{
	int i, j, k, num_codons, test;
	double rnd;

	//things related to rng
	struct timeval curr_time;

	const gsl_rng_type *rng_type;
	gsl_rng *rn;

	//Set up GSL RNG
	gsl_rng_env_setup();

	//Create a generator chosen by the 
	//environment variable GSL_RNG_TYPE
	//use default gsl for generating uniform rn 
	//from which other rn functions are derived 

	rng_type = gsl_rng_default;
	rn = gsl_rng_alloc(rng_type);

	//seed rng using clock
	gettimeofday(&curr_time, NULL);
	gsl_rng_set(rn, (curr_time.tv_sec + curr_time.tv_usec));	// (const gsl_rng * r, unsigned long int s)   

	//end RNG set up

	for (i = 1; i < aa_count; i++) {
		k = aa_index[i];
		num_codons = AA[k].num_codons;
		if (num_codons > 1) {
			rnd = gsl_rng_uniform(rn);
			if (G.benchmark)
				rnd = 0.5;

			test = 0;
			for (j = 0; j < num_codons && test == 0; j++) {
				if (rnd < AA[k].neutral_obs_pr_cum[j])	//Picks a codon based on its neutral probability and takes into account AT bias
				{
					//copy codon index and seq
					codon_index[i] = AA[k].codon_index[j];
					test = 1;	//break;
				}
			}
		}
	}
}

/*======================================================*
 * FUNCTION: Convert_Codon_Index_to_Sigma_Vec			*
 * 		Calculates the probability of successful 		*
 * 		translation up to and including each codon. All	*
 * 		values are stored in a new vector sigma_vec		*
 * =====================================================*/

void
Convert_Codon_Index_to_Sigma_Vec(int *index_vec,
				 double *sigma_vec, int aa_count)
{
	int i;
	double sigma_i;

	sigma_vec[0] = 1.0;
	sigma_i = 1.0;

	for (i = 1; i < aa_count; i++) {	//skip first amino acid
		sigma_i *= Codon[index_vec[i]].elong_pr;
		sigma_vec[i] = sigma_i;
	}

}

/*======================================================*
 * FUNCTION: Calc_B_over_C_Vec							*
 * 		Calculates b/c (NSE rate/elongation rate) for 	*
 * 		all codons in codon_index_vec and stores values	*
 * 		in new vector b_over_c_vec 						*
 * =====================================================*/

void
Calc_NSE_Odds_Vec(int *codon_index_vec, double *NSE_odds_vec,
		  int aa_count)
{
	int i;

	for (i = 0; i < aa_count; i++) {
		NSE_odds_vec[i] = (1-Codon[codon_index_vec[i]].elong_pr) / Codon[codon_index_vec[i]].elong_pr;
	}

}

/*======================================================*
 * FUNCTION: Convert_Sigma_Vec_to_Sigma_Ratio_Vec		*
 * 		Calculates the ratio of sigma_vec[i]/sigma_n	*
 * 		for all values in sigma_vec and stores values	*
 * 		in new vector sigma_ratio_vec					*
 * =====================================================*/

void
Convert_Sigma_Vec_to_Sigma_Ratio_Vec (double *sigma_vec,
				      double *sigma_ratio_vec, int aa_count) {
      int i;
      double sigma_n;
      
      sigma_n = sigma_vec[aa_count - 1];
      for (i = 0; i < aa_count; i++) {	//don't think the first is every used
	      sigma_ratio_vec[i] = sigma_vec[i] / sigma_n;
      }

}

/*======================================================*
 * FUNCTION: Recalc_Sigma_Vec							*
 * 		Recalcs sigma vec given that elong_pr of codon	*
 * 		at position i changes by a factor				*
 * =====================================================*/

void Recalc_Sigma_Vec (double *sigma_vec, int mut_codon_pos, 
                  double factor, int aa_count) 
{
  int i = mut_codon_pos;
  for (; i < aa_count; i++) 
  {
    sigma_vec[i] *= factor;
  }
}


/*======================================================*
 * FUNCTION: Calc_Xi_NSE								*
 * 		Calculates Xi based on a vector of sigma values.*
 * =====================================================*/
 
double Calc_Xi_NSE (double *sigma_vec, int aa_count) {
  int i = 0;
  double xi = 0.0;
  double sigma_n;
  
  
  sigma_n = *(sigma_vec + (aa_count - 1));
  sigma_vec++;
  
  //xi=(G.A1+G.A2)*(1-*sigma_vec);//should be zero
  for (i = 1; i < aa_count; i++) 
  {
    xi += (G.A1 + G.A2 * (i)) * (*(sigma_vec - 1) - *(sigma_vec++));
    //\sigma_{i-1} b/(b+c_i) = \sigma_{i-1}-\sigma_i 
    //modified from G.A2*(i+1)
    //rationale: should be i-1, but C indexing offsets things by 1, so should just be i
  }
  
  xi *= (1 / (1 - sigma_n));
  
  return xi;
}


/*======================================================*
 * FUNCTION: Calc_Eta_NSE								*
 * 		Calculate eta based on xi and sigma_n			*
 * =====================================================*/

double Calc_Eta_NSE (struct seq_struct * seq) 
{
  double eta;
  int aa_count = seq->aa_count;
  
  Convert_Codon_Index_to_Sigma_Vec(seq->codon_index,
                 seq->sigma_vec,
                 aa_count);
  seq->sigma_obs = seq->sigma_vec[aa_count - 1];
  seq->xi_obs = Calc_Xi_NSE(seq->sigma_vec, aa_count);
  
  eta = (1 - seq->sigma_obs) * seq->xi_obs / seq->sigma_obs + (G.A1 + G.A2 * aa_count);
  return eta;
}

void eta_min_max(int *aa_index,int aa_count, double * answer){
    // Purpose: calculate min and max eta, corresponding to optimal and pessimal sequences
    // Inputs: aa_index: memory address of first amino acid index 
    //         aa_count: number of amino acids in sequence
    // outputs: answer: memory address pointing to value of eta_min and eta_max
    //                  --note: *(answer) is eta min
    //                          *(answer + 1) is eta max
    
    int i, j;

    double t_min_elpr,t_max_elpr; //temporary variables for min and max elong pr
    int t_aa_index,t_codon_index; //temp var for aa_index
    double eta_min,eta_max;
    
    struct seq_struct optSeq;
    struct seq_struct pessSeq;
    
    optSeq.aa_count = aa_count;
    pessSeq.aa_count = aa_count;
    
    //*============================================*
    //* 1) Find optimal and pessimal codon indices *
    //*============================================*
    
    for(i=0;i<aa_count;i++){
        t_aa_index = *(aa_index+i);
        t_min_elpr = AA[t_aa_index].min_elong_pr;
        t_max_elpr = AA[t_aa_index].max_elong_pr;
        
        //* Search through codons for t_min_elr and t_max_elr
        for(j=0;j<AA[t_aa_index].num_codons;j++){           
            t_codon_index = AA[t_aa_index].codon_index[j];
            
            //* If this codon has AA's min elong pr, store in pessimal
            //* sequence
            if(Codon[t_codon_index].elong_pr == t_min_elpr){
                pessSeq.codon_index[i] = t_codon_index;
            }
            
            //* If this codon has max elong pr, store in optimal
            if(Codon[t_codon_index].elong_pr == t_max_elpr){
                optSeq.codon_index[i] = t_codon_index;
            }
        }
        
        //* Check to see if this aa position was assigned an optimal
        //* and pessimal codon index
        
        if(pessSeq.codon_index[i]==25443||optSeq.codon_index[i]==25443){ //25443 is the integer equivalent of the character 99
            Rprintf("Error in eta_min_max(): AA position %d was not given an optimal or pessimal codon index",i);
        }
    }
    
    
    //*==========================*
    //* 2) Calculate minimum eta *
    //*==========================*
    
    eta_min = Calc_Eta_NSE(&optSeq);
        
    //*==========================*
    //* 3) Calculate maximum eta *
    //*==========================*
        
    eta_max = Calc_Eta_NSE(&pessSeq);
        
    //*===================================*
    //* Store eta min and max in answer[] *
    //*===================================*
    
    answer[0] = eta_min;
    answer[1] = eta_max;
        
}

 void flip_codons(struct seq_struct *currSeq, struct seq_struct *propSeq, double *randNums, int num_flips){
	 int i,j;
	 int * currCodonIndex = currSeq->codon_index;
	 int currCodon,currAA;
	 int seqLen = currSeq->aa_count; //length of current sequence
	 int flip_count = 0;
	 int newCodon;
	 
	 while(flip_count < num_flips){
		 
		 i = j = 0;
		 //1) Find the codon index that needs flipped
		 i = (int) (*(randNums + 2*flip_count)*seqLen);
		 
		 currCodon = *(currCodonIndex+i);
		 currAA = Codon[currCodon].aa_index;
		 
		 //2) Find the codon it will be flipped to
		 j = (int) (*(randNums + 2*flip_count+1)*AA[currAA].num_codons);
		 
		 newCodon = AA[currAA].codon_index[j];
		 
		 //3) Change codon in propSeq
		 propSeq->codon_index[i] = newCodon;
		 
		 
	
		 flip_count++;
	} 
	 
 }

void Codon_Counts(int * codon_index_vec ,int * codon_cts_vec ,int aa_count){
    int j;
    int codon_indx;
    
    //* Initialize
    for(j=0;j<61;j++){
        *(codon_cts_vec+j)=0;
    }
    
    //* Add up codon counts
    for(j=0;j<aa_count;j++){
        codon_indx = *(codon_index_vec + j);
        codon_cts_vec[codon_indx]++;
    }
}



/*======================================================*
 * FUNCTION: Calc_Seq_Mu 								*
 * 		Calculates mu of sequence 						*
 * 		\mu(\vec{c})=\prod_{i=1}^{61}{\mu_i^{x_i}}		*
 * =====================================================*/
 
 double Calc_Seq_Mu(int * codon_counts_vec, int aa_count) {
	 
	 int i;
	 int codon_counts;
	 double mu, seq_mu;
	 
	 seq_mu = 1;
	 
	 for(i = 0; i < 61; i++) {
		 codon_counts = *(codon_counts_vec++);
		 mu = Codon[i].incoming_mu;
		 
		 seq_mu *= pow(mu,codon_counts);		 
	 }
	 
	 return seq_mu;
 }

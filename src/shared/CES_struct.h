/*==========================*
 * Preprocessor directives  *
 *==========================*/
#ifndef CES_STRUCT_H_
#define CES_STRUCT_H_
 
#define MAX_DDIM 30000
#define MAX_AA 6500
#define MAX_NTS 19500
#define MAX_LOCI 6000		//set max loci to read in and simulate
#define MAX_TIME -5		//Max number of evolutionary steps to simulate,
		    // if set to a negative value it will calculate the appropriate max time to have X substitutions/nt.
		    //no longer used.  Defined via coommand line
#define NE 1.36e7		//effective population size.
		  //#define Q 1.519e-5 //fitness scaling coefficient
#define MU 1E-9			//per generation mutation rate


/*=======================================================*
 *  STRUCTURE: amino_acid								 *
 * 														 *
 *	   Defining structure for a single amino acid        *
 *	   keep track of E(1/c), E(1/c^2), E((c+b)/c)... 	 *
 *	   which are expectations wrt the translation rates  *
 *	   for the aa's set of codons 						 *
 * ======================================================*/
 
struct amino_acid {
      char codon[6][4];
      char aa;
      int num_codons;		//was cc
      int codon_index[6];	//used when generating random sequences

      //expected values for elongation related rates
      double e_omega;				// \omega
      double e_invp;				// 1/p
      double e_omega2;				// \omega^2
      double e_invp2;				// 1/p^2
      double e_omegainvp;			// \omega/p --used in calc cov b/w codons
      double elong_rate[6];			//elongation rate of codons for AA, was tr_rate
      double elong_pr[6];			//Pr of successful elongation p
      double NSE_pr[6];				//Pr of NSE 
      double omega[6];				//NSE odds = NSE_pr/elong_pr
      double max_elong_pr;			//max elongation rate for the amino acid
      double min_elong_pr;			//min Pr of successful elongation for the amino acid
      double neutral_obs_pr[6];		//Pr of observing this codon based solely on its nt sequence and any AT bias.
      double neutral_obs_pr_cum[6];	//Cumulative probability of observing this codon. Useful for when choosing a random codon when there is AT bias
      double mu[6][6];				//hold mutation values
      double incoming_mu[6]; 					  		//sum of all mutation rates to codon
};



/*=======================================================*
 * STRUCTURE: codon_struct								 *
 * 														 *
 * 		Defining structure for a single codon. 			 *
 * 		Keeps track of individual information as well as *
 * 		information relative to one-step neighbors.		 *
 * ======================================================*/

struct codon_struct {
      char codon[4];
      char aa;			//aa letter

      int synonym_index[5];	//Alternative synonyms for the same aa listed by their codon index #
      int one_step_synonym_index[5];	//as above but subset that is one step mutant from codon.

      int aa_index;

      int num_synonym;
      int num_one_step_synonym;
      


      double elong_pr_ratio[6];						
      double delta_omega[6];						
      double relative_mu; 							// relative mutation rate to this codon from all one step synonyms
      double one_step_relative_mutation_rate[5];	//mutation rate to one step neighbors of current codon
													//as above but subset that is one step mutant from codon.
	  double syn_relative_mutation_rate[5];    //mutation rate from this codon to each neighbor (all neighbors,
	                                                //not just one-step neighbors)
	  double incoming_mu;
	 

      double omega;		//(NSE_pr)/(elong_pr) = (1-elong_pr)/(elong_pr)
      //double elong_rate;	//elongation rate of codons for AA, was tr_rate FIXME: this should be removed
      double elong_pr;		
      double NSE_pr;		
};



/*==========================================================*
 * STRUCTURE: seq_struct									*
 * 															*
 * 		Defining structure for one genetic sequence.		*
 * 		Holds information relevant to the whole sequence	*
 * =========================================================*/

struct seq_struct{
      char id[26];		//Gene ID or ORF name
      char codon_seq[MAX_AA][4];	//codon sequence
									
      int codon_index[MAX_AA];	//codon index code
      int aa_index[MAX_AA];		//aa index code

      int aa_count;						//aa count
      double sigma_vec[MAX_AA];			//vector of \sigma(i) values.  
										//Note that sigma_0 repesents the start codon so its value is 1. 
      double sigma_ratio_vec[MAX_AA];	//vector of \sigma(i)/\sigma(n) values
      double max_time;
      int ARR;
      double pA;		//alpha parameter for f(eta)~beta distribution
						//used by gsl
      
      

      double sigma_obs;			//obs sigma value
      double xi_obs;			//obs xi value
      double eta_obs;			//obs eta value
      double mu_obs;            //mutation coefficient for observed sequence = prod_i \mu_i^(x_i(\cvec))
                                //  where mu_i is the sum of relative mutation rates to codon i and x_i
                                //  is the number of codon i in the sequence cvec  
      

	  int codon_cts[64]; //running total of number of codons
	  
      // following added for evolution simulations
      double eta_initial;		//eta value at start of simulation
      double phi_obs;			//assumed production rate
      double z_obs;				//assumed production rate, scaled by Ne and q
      double delta_eta_mean;	//mean of distn of delta_eta values for wt allele
      double delta_eta_var;		//var of distn of delta_eta values for wt allele

      double evol_delta_eta_mean;	//mean diff b/w wt_eta and mut_eta;
      double evol_delta_eta_var;	//var b/w wt_eta and mut_eta;
      double evol_time;				//var b/w wt_eta and mut_eta;

      int evol_steps;		//number of times wt allele is replaced by a mutant
      int num_synon_mut;	//num of mut that were synonymous

};

struct global_variable_struct {
	double A1;          //Cost of initiation of translation in units of ATP. Default val set in command line to  A1 = 2;
	double A2;          //Cost each elongation step in ATPs, default set in command line to A2=4
	double Ne;      //NE;
	double Q;       //4.19e-7;     // Fitness scaling coefficient

	int max_aa;

	int ignore_aa;              //decrease true aa_count by this amount.

	int benchmark;              //Uses predetermined values to generate random sequence and calculate mut_pr and time_step
									//instead of random generated 
};

#endif //CES_STRUCT_H_

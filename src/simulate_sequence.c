#include "CES.h"

void Evolve_Sequence(struct seq_struct *seq, int total_evol_steps, double global_max_time);
void MCMC_Sequence(struct seq_struct *seq, int total_evol_steps);
 
void simulate_sequence(int *CODON_INDEX,                  double *ELONG_PR,         	    double *MUTATION_RATES, 
					   int *AA_COUNT,                     double *PHI,                      double *POP,
					   double *A_1,                       double *A_2,                      double *AT_BIAS,                  
					   int *BENCH,                        double *BEE,                      double *GMT,                      
					   int *IGNORE,                       double *GAMMA,                    double *SCALE_FACTOR, 
					   int *MES,             
					   int *BIS,                          char ** AA_VEC,                   char ** CODON_VEC,
					   int * N_CODON,   				  char ** SIMULATION_METHOD)	
{
	
  int j;
    
  int simulation_steps;
	int burn_in_steps;
	int total_evol_steps;
	
	double global_max_time;
    
  //* Initializing the random number generator
  srand(time(NULL));  //should this be NULL? mikeg
  

  //*============================*
  //*  1) Process input from R   *
  //*============================*

  //* Import Sequence data and genome parameters
  Process_R_input(CODON_INDEX,ELONG_PR,MUTATION_RATES,AA_COUNT,
                  PHI, POP, A_1, A_2, AT_BIAS,
                  BEE, IGNORE, GAMMA, SCALE_FACTOR,
                  AA_VEC, CODON_VEC, N_CODON);
  
  //* Simulation parameters
  G.benchmark = *BENCH;
  global_max_time = *GMT;
  simulation_steps = *MES;
  burn_in_steps = (*BIS) * Sequence.aa_count;
    
    
	Process_AA_Information();
	
  //* adjust aa_counts if we are ignoring aa's on tail end
  //* the idea is that a protein can still function if it is missing
  //* the last few amino acids (i.e. translation has a nonsense error
  //* on the last few codons). We do not use this in the default setting
  if (G.ignore_aa > 0) {
      Sequence.aa_count -= G.ignore_aa;
  }
 
  //*=========================================================*
  //*  2) Calculate length of simulation                      *
  //*     --Default behavior is to use fixed number of        *
  //*       steps, calculated from *MES and *BIS, but we can  *
  //*       also used fixed amount of time, calculated from   *
  //*       *GMT                                              *
  //*=========================================================*
  total_evol_steps = burn_in_steps + simulation_steps; //BIS and MES parameters passed from R
      
  //*==================================*
  //*  3) Calculate initial eta and mu *
  //*==================================*
  Sequence.eta_initial = Sequence.eta_obs = Calc_Eta_NSE(&Sequence);
  Codon_Counts(Sequence.codon_index,Sequence.codon_cts,Sequence.aa_count);
	Sequence.mu_obs = Calc_Seq_Mu(Sequence.codon_cts,Sequence.aa_count);

  //*========================*
  //* 4) Simulate Evolution  *
  //*========================*
  
  switch(**SIMULATION_METHOD){
  case 'E':
	  Evolve_Sequence(&Sequence,total_evol_steps,global_max_time);
	  break;
	
	case 'M':
		MCMC_Sequence(&Sequence,total_evol_steps);
		break;
	}
	
  //*=========================*
  //* 5) Update *CODON_INDEX  *
  //*=========================*
  
  for(j = 0; j < *AA_COUNT; j++) {
      *(CODON_INDEX + j) = Sequence.codon_index[j];
  }
            
}

 

/*======================================================*
 * FUNCTION: Evolve_Sequence                            *
 * =====================================================*/

void Evolve_Sequence(struct seq_struct *seq, int total_evol_steps, double global_max_time)
{
	
	int aa_count = seq->aa_count;
	int aa_position,wt_index,mut_index;

  char wt_codon[4];
  char mut_codon[4];

  int D_array[MAX_DDIM][2];    //array with codon position and codon index of one step neighbors
  int mut_codon_index;
  int wt_codon_index;
  int mut_codon_pos;   //position numbers for the codon and nt that's mutated

  int D_dim;      //number of dimensions
  int i, j;
  int directly_update_eta = 1; //flag to use workaround fix for errors in indirect delta_eta calculation  


  //indices used for Generate_D_Array
  int Dmax;
  int D_index;
  int codon_index;

  double phi = seq -> phi_obs;        //protein production rate
  double qPhi = G.Q * phi;
  double two_qNePhi = 2*G.Ne*qPhi;
  double factor;      //ratio of mut and wt elong_pr, used to rescale sigma_vecs
  double inv_factor;  //inverse of above
  
  /* Parameters related to simulation length*/
  double evol_time_1=0,evol_time_2=0;
  double time_step;
  double wait_parameter;
  double max_time=1,avg_time_step;
  int evol_steps;

  //  double sigma_ratio_vec[MAX_AA];
  double NSE_odds_vec[MAX_AA];
  double delta_eta_first_term_vec[MAX_AA];
  double delta_eta_second_term_vec[MAX_AA];

  double delta_eta_vec[MAX_DDIM]; //array for storing $\Delta \eta_{i,j}$ values for all
  // one step mutants of the resident allele
  double pi_vec[MAX_DDIM];    //array for storing values of $\pi(i->j)$ for all one step mutants
  double pi_total;    //sum over pi_vec values.
  double mutation_vec[MAX_DDIM];  //array for storing values of mutation rates $\mu_{i->j}$ for all one step mutants
  double pr_vec[MAX_DDIM];    //vector of pi * mu values
  double pr_total;    //sum over pr_vec values.
  double mut_pr;      //RV representing replacement allele
  double curr_pr;     //used to find replacement allele

  //double delta_eta_list[MAX_EVOL_STEPS];
  double delta_eta_mean = 0;  //mean effect of synon sub for the current wt seq-- def differs from before
  double delta_eta_var = 0;   //var effect of synon sub for the current wt seq --definition differs from before
  double delta_eta;   //eta_wt-eta_mut

  //variables used for updating eta and calculating delta_eta directly
  double previous_eta;
  double new_eta;
  double direct_calc_delta_eta; //delta_eta based on direct calculations of eta for new and previous sequences.
  double abs_error_delta_eta;
  double rel_error_delta_eta;
 
  struct timeval curr_time;

  const gsl_rng_type *rng_type;
  gsl_rng *rn;
    

  /*========================================================*
   *                1) SET UP GSL RNG                       *
   *      Create a generator chosen by the                  *
   *      environment variable GSL_RNG_TYPE                 *
   *      use default gsl for generating uniform rn         *
   *      from which other rn functions are derived         *
   *========================================================*/

  gsl_rng_env_setup();

  rng_type = gsl_rng_default;
  rn = gsl_rng_alloc(rng_type);

  //seed rng using clock
  gettimeofday(&curr_time, NULL);
  gsl_rng_set(rn, (curr_time.tv_sec + curr_time.tv_usec));    // (const gsl_rng * r, long int s)   

  /*========================================================*
   *    2) CALCULATE SEQUENCE SIGMA_VEC AND B_OVER_C_VEC    *
   *========================================================*/

  Convert_Sigma_Vec_to_Sigma_Ratio_Vec(&(seq->sigma_vec[0]),
                       &(seq->sigma_ratio_vec[0]),
                       aa_count);

  Calc_NSE_Odds_Vec(&(seq->codon_index[0]), NSE_odds_vec, aa_count);
 
  /*========================================================*
   *              3) START SEQ ITERATIONS                   *
   * The first round of simulation steps should burn in the *
   * sequence with a fixed number of steps. The second      *
   * of simulation steps should run for a given amount of   *
   * evolution time.                                        *
   *========================================================*/
    evol_steps = 0;
    time_step = 0;
    
    //* this while loop goes through a prescribed number of burn in steps
    //* total_evol_steps then runs for a prescribed amount of time
    while ((evol_steps < total_evol_steps) || (evol_time_2 < max_time)) 
    {    
        //* figure out time based on average step time duing burn-in
      if(evol_steps==total_evol_steps)
      {
  			if (global_max_time < 0) 
        {
  				avg_time_step = evol_time_1/total_evol_steps;//average time spent at each step
  				max_time = -avg_time_step*aa_count*global_max_time;
  			} else {
  				avg_time_step = evol_time_1/total_evol_steps;//average time spent at each step
  				max_time = global_max_time * avg_time_step;
  			}
		  }
                
      //*=====================================================================*
      //*  4.b) Determine Dimensionality/# one step neighbors of the system.  *
      //*  D_dim DOES change over the simulation!                             *
      //*  The changes in its values are due to the effects of 6 codon aa     *
      //*  For example, for Arginine AGA can mutate to AGG or CGA.            *
      //*  In contrast, CGA can mutate to AGA, CGC, CGG, or CGU.              *
      //*=====================================================================*

      D_dim = Calc_Dimensionality(&(seq->codon_index[0]), aa_count);
              
      //*=============================================================================*
      //*  4.c) Create D_array                                                        *
      //*                                                                             *
      //*  D_array is a 2D matrix of size 2 x Dimensionality of sequence              *
      //*      Dimensionality = # one step syonymous neighbors                        *
      //*      D_array[i][0] = j, where j is the position of an amino acid in the ORF *
      //*      D_array[i][1] = k, where k is a codon_index value for an alternative   *
      //*                            synonymous codon for the wt sequence             *
      //*=============================================================================*
      
      D_index = 0;
        
      for (i = 0; i < aa_count; i++) 
      {
          codon_index = seq->codon_index[i];
          Dmax = Codon[codon_index].num_synonym;

          //go through each synonym
          for (j = 0; j < Dmax; j++) 
          {
              mutation_vec[D_index] = Codon[codon_index].syn_relative_mutation_rate[j];

              (D_array[D_index][0]) = (int)(i);    //amino acid position
              (D_array[D_index++][1]) = (Codon[codon_index].synonym_index[j]);   //synonymous codon for position D_index
          }
      }     
        
		//*==========================================================*
		//*  4.d) Calculate Delta_Eta and fixation probability (pi)  *
		//*  for all one step neighbors.                             *
		//*==========================================================*
		
		Calc_Delta_Eta_NSE_First_and_Second_Term_Vecs(&(seq->sigma_ratio_vec[0]),
													  NSE_odds_vec,
													  delta_eta_first_term_vec,
													  delta_eta_second_term_vec,
													  aa_count);

		Calc_Delta_Eta_NSE_Vec(&(seq->codon_index[0]), NSE_odds_vec,
							   delta_eta_first_term_vec,
							   delta_eta_second_term_vec,
							   delta_eta_vec, aa_count, D_dim);

		
		Calc_Pi_Vec(delta_eta_vec, pi_vec, &pi_total, phi, qPhi,
					two_qNePhi, aa_count, D_dim);

		//*=====================================================================*
		//*  4.e) Calculate replacement probability for all one step neighbors  *
		//*       --> Test to see if there are mutation effects                 *
		//*           -Uses AT bias and transition/transversion bias OR         *
		//*            mutation rates loaded from a mutation file.              *
		//*=====================================================================*

		//* Calculate true pr of replacement, weighing by mutation rates
		Calc_Replacement_Pr_Vec(mutation_vec, pi_vec, pr_vec,
						&pr_total, D_dim);
    
       
        
    //*=============================================================*
    //*  4.f) Calculate expected time until replacement             *
    //*  We are assuming that replacement                           *
    //*  happens quickly. Ideally we would expect the wait          *
    //*  parameter to be >>1, but it likely doesn't matter          *
    //*  The time to replacement should follow a geometric dist     *
    //*  which we approximate with an exponential distribution.     *
    //*  Note the parameter for GSL's RNG argument is the expected  *
    //*  wait time, so it is 1/replacement rate                     *
    //*=============================================================*
    
    wait_parameter = 1 / pr_total;  // Time spent at this allele is inversely proportional
                                    // to the 'failure' or leaving rate of the resident allele.
    
    //* Calculate based on wait_parameter only or random exponential
    
    if (G.benchmark||1)
    {
        time_step = wait_parameter;
    }else{
        time_step = gsl_ran_exponential(rn, wait_parameter);
    }
       
    //*==============================================*
    //*              DEBUGGING TOOLS                 *
    //*==============================================*
    //*  Prints a detailed summary of all one step   *
    //*  neighbors.                                  *
    //*==============================================*
    
    if(1==0){
        
        Rprintf("############# Step Number: %d #############\n",evol_steps);
        Rprintf("\t\t\t\tdelta_eta\tpi_i_j\tpr_i_j\n");
        for(i=0;i<D_dim;i++){
            aa_position = D_array[i][0];
            wt_index = seq->codon_index[aa_position];
            mut_index = D_array[i][1];
            
            //Rprintf("Codon %d: %f  ->  Codon %d: %f\t%g\t%g\t%g\n",wt_index,Codon[wt_index].elong_pr,mut_index,Codon[mut_index].elong_pr,delta_eta_vec[i],pi_vec[i],pr_vec[i]);
            error("Delta eta did not work\n");
        }
    }   
        
    //*==============================================*
    //*  4.g) Choose a RV to determine which allele  *
    //*  replaces the resident.                      * 
    //*==============================================*
    
    mut_pr = gsl_rng_uniform(rn)*pr_total;
    if (G.benchmark) mut_pr = 0.5 * pr_total;

        //* Note: pr_total is not scaled by the constants Ne or mu)
    
    //*===========================================================*
    //*  4.h) Search forward to find substitution allele.         *
    //*  This is likely more efficient than starting at the end   *
    //*  b/c there is weaker selection at the start of the        *
    //*  sequence.                                                *
    //*  Some kind of newton search could be even more efficient  *
    //*===========================================================*
    
    curr_pr = 0.0;
    i = 0;
    while (curr_pr < mut_pr) {
        curr_pr += pr_vec[i++];
    }

    i--;
    D_index = i;

		//***Store This for exchange with bridging***//
    mut_codon_pos = D_array[D_index][0];
    mut_codon_index = D_array[D_index][1];

		//*======================*
		//*  Get wt codon index  *
		//*======================*
        
    wt_codon_index = seq->codon_index[mut_codon_pos];
       
		//*==============*
		//*  Get codons  *
		//*==============*
        
    strcpy(wt_codon, Codon[wt_codon_index].codon);
    strcpy(mut_codon, Codon[mut_codon_index].codon);
    
    if(G.benchmark) {
        Rprintf("Substituting Codon %d %s for Codon %d %s\n",wt_codon_index,wt_codon,mut_codon_index,mut_codon);
    }
           
    //*========================================*
    //*      4.i) Update everything            *
    //*========================================*
        
		//* update codon index
    seq->codon_index[mut_codon_pos] = mut_codon_index;

		//* update codon counts
    seq->codon_cts[wt_codon_index]--;
    seq->codon_cts[mut_codon_index]++;

		//* NOTE: Mike had seen an error in this part of the code at one point. He saw a descrepancy 
        //* in the value calulated by the delta eta function and the value calculated by explicitly 
        //* finding the difference in eta values. However, I never experienced this issue, so 
        //* I stopped exploring it. Below are the original comments from Mike and the code he wrote 
        //* to diagnose the problem. Notice that the if statement is permanently turned off. <whewgley 5-23-14>
		
		//* FIXME: This delta_eta appears to be slightly off!
		//* Using temporary work around <mikeg>
    delta_eta = delta_eta_vec[i];   //should be \eta_i - \eta_j

    if(directly_update_eta&&0)
    {
      //* WORKAROUND: Calculate eta for new sequence from scratch
      //* this line overwrites the updating of sigma_vec done above which might be the cause of our errors in delta_eta
      //* it should be possible to tell since the delta_eta would work for the first evolutionary step but not the later ones.
      //* TODO: Comment out this line to see if there's a problem with how we update sigma_vec above
        
      new_eta = Calc_Eta_NSE (seq);
      //* this overwrites update of seq->eta_obs using delta_eta
      seq->eta_obs = new_eta;

      if (1 == 1) 
      {
        //* check to make sure delta_eta is working right
        direct_calc_delta_eta = (previous_eta - new_eta);
        abs_error_delta_eta = abs(direct_calc_delta_eta - delta_eta);
        rel_error_delta_eta = abs_error_delta_eta/direct_calc_delta_eta;
          
        //Rprintf("Step %d: Position: %d Codon %d %s --> Codon %d %s\td_eta: %f \tdir_d_eta: %f\n",evol_steps,mut_codon_pos,wt_codon_index,wt_codon,mut_codon_index,mut_codon,delta_eta,direct_calc_delta_eta);
        if (rel_error_delta_eta > 1.0E-6) 
        {           
          Rprintf ("ERROR: direct_d_eta: %f != delta_eta: %f \n\tAA: %c Codon %d %s --> Codon %d %s\n", direct_calc_delta_eta, delta_eta,Codon[wt_codon_index].aa,wt_codon_index,Codon[wt_codon_index].codon,mut_codon_index,Codon[mut_codon_index].codon);
          Rprintf("############# Step Number: %d D_dim: %d #############\n",evol_steps,D_dim);
	        Rprintf("\t\t\t\tdelta_eta\tpi_i_j\tpr_i_j\n");
		      for(i=0;i<D_dim;i++)
          {
      			aa_position = D_array[i][0];
      			wt_index = seq->codon_index[aa_position];
      			mut_index = D_array[i][1];
                  
      			Rprintf("Codon %d: %f  ->  Codon %d: %f\t%g\t%g\t%g\n",wt_index,Codon[wt_index].elong_pr,mut_index,Codon[mut_index].elong_pr,delta_eta_vec[i],pi_vec[i],pr_vec[i]);
		      }
		      error("Delta eta did not work\n");
        }
      }
      //Rprintf("delta_eta: %f\tdirect_delta_eta: %f\n",delta_eta,(previous_eta - new_eta));
      delta_eta = (previous_eta - new_eta);
    } else {

      //* Update Sigma_Ratio_Vec by dividing the values up to  
      //* mut_codon_pos by elong_pr_i/elong_pr_j               
      //* values. After this point don't change since they     
      //* involve the factor in both the                       
      //* denominator and numerator.                           
          
      factor = Codon[mut_codon_index].elong_pr / Codon[wt_codon_index].elong_pr;
      inv_factor = Codon[wt_codon_index].elong_pr / Codon[mut_codon_index].elong_pr;

      for (i = 1; i <= mut_codon_pos; i++) 
      {
        seq->sigma_ratio_vec[i] *= inv_factor;
      }
          
      //* update sigma_vec for values at and above mut_codon_pos;
      for (i = mut_codon_pos; i < aa_count; i++) 
      {
        seq->sigma_vec[i] *= factor;
      }
          
      //* update eta_obs by subtracting difference: \eta_wt - (\eta_wt - \eta_mut) = \eta_mut
      seq->eta_obs -= delta_eta;

    } //end else


		//* update NSE_odds_vec
    NSE_odds_vec[mut_codon_pos] = Codon[mut_codon_index].omega;

    //* update evol_delta_mean/var vals
    //evol_delta_mean = delta_eta - delta_eta_mean;
    //evol_delta_eta_mean += evol_delta_mean / (evol_time + time_step);
    //evol_delta_eta_var +=  evol_delta_mean * (delta_eta - evol_delta_eta_mean);

		//* update values for printing
    //seq->evol_delta_eta_mean = evol_delta_eta_mean;
    //seq->evol_delta_eta_var = evol_delta_eta_var;
        
        
        
    if(evol_steps < total_evol_steps){
			evol_time_1 += time_step; //increment evol_time (changed for loop to maximum#steps)
		}else{
			evol_time_2 += time_step;
		}
		
		
		evol_steps++;
  }

  /*========================================================*
   *                 END EVOLUTION ITERATIONS               *
   *========================================================*/
    //update Seq[]
    seq->evol_steps = evol_steps;
    seq->evol_time = evol_time_2;
    seq->delta_eta_mean = delta_eta_mean;
    seq->delta_eta_var = delta_eta_var;
    
    gsl_rng_free(rn);
}

void MCMC_Sequence(struct seq_struct *seq, int total_evol_steps)
{
        
	int i;
	int recent_steps[250]; 	//this is a running buffer of acceptances and rejections. 
							//it is indexed each step as recent_steps[step_num%250] so
							//that each successive step overwrites the step 250 steps ago.
	
	int accept_num=0,step_num=0;
	double accept_ratio;
	double desired_acc_ratio = 0.3;
	double accept_pr;
	
	struct seq_struct tSeq;
	
			
	int num_sub=5;		
	int max_sub = 50;
	double randNums[200];
	double step_eta, step_mu;
	
	struct timeval curr_time;
	
	const gsl_rng_type *rng_type;
  gsl_rng *rn;
    		
	/*========================================================*
	 *                1) SET UP GSL RNG                       *
	 *      Create a generator chosen by the                  *
	 *	    environment variable GSL_RNG_TYPE                 *
	 *      use default gsl for generating uniform rn         *
	 *      from which other rn functions are derived         *
	 *========================================================*/

  gsl_rng_env_setup();

  rng_type = gsl_rng_default;
  rn = gsl_rng_alloc(rng_type);

  //seed rng using clock
  gettimeofday(&curr_time, NULL);
  gsl_rng_set(rn, (curr_time.tv_sec + curr_time.tv_usec));    // (const gsl_rng * r, long int s)   		
        
    
  ///////////////////////////////////
  //1: Initialize tSeq.codon_index //
  ///////////////////////////////////
  
  for(i = 0; i < seq->aa_count; i++)
  {
  	tSeq.codon_index[i] = seq->codon_index[i];
	}
  tSeq.aa_count = seq->aa_count;
    
    
    
  /////////////////////////////////////////////////////
  //2: Start for loop-- start num substitutions at 5 //
  /////////////////////////////////////////////////////
    
	while(step_num < total_evol_steps) 
  {
		
	  //////////////////////////////////	
    //3: Take a step in codon space //
    //////////////////////////////////
    
		for(i = 0; i < 2*num_sub; i++)
    {
			*(randNums + i) = gsl_rng_uniform(rn);
		}
		for(i = 0; i < seq->aa_count; i++)
    {
			tSeq.codon_index[i] = seq->codon_index[i];
		}
		
		flip_codons(seq,&(tSeq),randNums,num_sub);
    
    ///////////////////////////
    //4: Calc new eta and mu //
    ///////////////////////////
		
		tSeq.eta_initial = tSeq.eta_obs = Calc_Eta_NSE(&tSeq);
		
		Codon_Counts(tSeq.codon_index, tSeq.codon_cts, tSeq.aa_count);
		
		tSeq.mu_obs = Calc_Seq_Mu(tSeq.codon_cts, tSeq.aa_count);
	
	  ////////////////////////////////////////	
    //5: accept pr = min(w(new)/w(old),1) //
    ////////////////////////////////////////
		//Assume diploid Wright-Fisher process -- (2*Ne-1)
		accept_pr = fmin(tSeq.mu_obs/seq->mu_obs*exp(-G.Q*(2*G.Ne-1)*seq->phi_obs*(tSeq.eta_obs-seq->eta_obs)), 1);
		//Rprintf("num_sub: %d curr_eta: %f new_eta: %f accept_pr: %f \n",num_sub,seq->eta_obs,tSeq.eta_obs,accept_pr);
    
    /////////////////////////////////
    //6: Draw rn and accept/reject //
    /////////////////////////////////
    
		*randNums = gsl_rng_uniform(rn);
		
		if(*randNums < accept_pr)
    {
			step_eta = tSeq.eta_obs;
			step_mu = tSeq.mu_obs;

      //Replace seq with tSeq - only need to change codon_cts, codon_index, eta_obs, mu_obs
			seq->eta_obs = step_eta;
			seq-> mu_obs = step_mu;
			for(i = 0; i < 64; i++)
      {
				*(seq->codon_cts+i) = tSeq.codon_cts[i];
			}
			for(i = 0; i < seq->aa_count; i++)
      {
				*(seq->codon_index+i) = tSeq.codon_index[i];
			}
			recent_steps[step_num%250] = 1; //keeps a running buffer of accepts and rejects
			step_num++;
			
		}else{ //If this is not a new step, only increment time
			recent_steps[step_num%250] = 0;
			step_num++;
		}
    
    ////////////////////////////////////////////////////////////////////////////
    //7: adjust num_substitutions (analagous to an adaptive variance)         //
    //      -if(accept ratio > desired_acc_ratio): increase num substitutions //
    //      -if(accept ratio < desired_acc_ratio): decrease num substitutions //
    ////////////////////////////////////////////////////////////////////////////
    
		if(step_num > 250 && step_num%50 == 0){
			accept_num = 0;
			//tally up all the acceptances from last 250 steps
			for(i = 0; i < 250; i++)
      {
				accept_num+=recent_steps[i];
			}
			
			accept_ratio = (double)accept_num/250;
			
			if( (accept_ratio > (desired_acc_ratio+0.05)) && (num_sub < max_sub) ){
				num_sub++;
			}else if( (accept_ratio < (desired_acc_ratio-0.05)) && (num_sub > 1) ){
				num_sub--;
			}
		}
		//Rprintf("step_eta: %f ",seq->eta_obs);
		//Rprintf("accept_num: %d step_num: %d ",accept_num,step_num);
		//Rprintf("accept_ratio: %f ",accept_ratio);
		//Rprintf("num_sub: %d\n",num_sub);
	}
  //Rprintf("num_sub: %d\n",num_sub);
	gsl_rng_free(rn);
}


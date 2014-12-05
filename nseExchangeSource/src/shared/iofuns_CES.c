#include "iofuns_CES.h"

/*=======================================================================*
 * DEPRECATED: We now pass parameters from R instead of reading them     *
 *             from a file.                                              *
 * FUNCTION:Read_tRNA_File                                               *
 *      Purpose: This function reads a tRNA file and stores information  *
 *               in amino_acid and codon_struct structures.              *
 *               See readme for details concerning file formats          *
 *      Inputs:  char *filename - string containing tRNA filename        *
 *      Outputs: int aa_processed - number of amino acids in tRNA file   *
 *      Usage:   sprinf(filename,"path/to/tRNA_file");                   *
 *               num_aa = Read_tRNA_File(filename);                      *
 * ======================================================================*/
 
int Read_tRNA_File (char *filename) 
{
      int i, j;
      FILE *file_handle;
      char curr_char;
      char aa_list[] = { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'Z', 'S', 'T', 'V', 'W', 'Y', 'X' };    //'X' represents the 'stop' amino acid, 'Z' represents the smaller set of serine codons: AGT and AGC.  The other ser codons are TCT, TCC, TCG, and TCA.
      //Note program expects 'X' to be at end of list so 'Z' is put in internally.
      int num_codons;
      int aa_processed;

      j = 0;
      //initialize AA structures
      for (i = 0; i < G.max_aa; i++) {
        AA[i].aa = aa_list[i];
        AA[i].num_codons = 0;
      }

      file_handle = fopen (filename, "r");
      if (!file_handle) {
        error("\ntRNA File: %s Doesn't Exist\ntRNA", filename);
      }
      //get aa index 
      curr_char = fgetc(file_handle);

      while (curr_char != EOF) {
        i = 0;      //set AA index to 0
        //flip through AA until you get a match
        while (curr_char != AA[i].aa && i < G.max_aa) {
          i++;
        }

        if (i == G.max_aa) {
          if (curr_char == '"') {
            while (curr_char != '\n'){
                  curr_char = fgetc (file_handle);
            }
            //Get AA index for next line
            curr_char = fgetc (file_handle);
          }else {
            error ("\nAA index '%c' in file did not match any known AA. Exiting...\n",curr_char);
		  }
        } else {
          //get current codon count
          num_codons = AA[i].num_codons;

          if (num_codons >= 6)  //check num_codons value 
          {
            printf ("\nCodon count for AA is greater than maximum possible value of 6.  Exiting...\n");
            exit (1);
          }


          curr_char = fgetc (file_handle);  //get next character. Should be a \t

          //read in codon
          for (j = 0; j < 3; j++)   //load codon sequence and put in codon_index codon index
          {
            AA[i].codon[num_codons][j] = fgetc (file_handle);
          }
          AA[i].codon[num_codons][j] = '\0';


          curr_char = fgetc (file_handle);  //get next char. Should be a \t

          //read in elongation rate
          int exitcode = fscanf(file_handle, "%lf", &AA[i].elong_rate[num_codons]);
			 if(exitcode == 0)
			 {
				printf("No elongation rates read. Exiting...\n");
				exit(1);
			 }
          //elong_rate = AA[i].elong_rate[num_codons];
          //increment codon count for the aa
          AA[i].num_codons++;
          num_codons++;


          //read until end of line or EOF
          while ((curr_char != '\n') && (curr_char != EOF))
            curr_char = fgetc (file_handle);

          //get aa index for next codon if not at EOF
          if (curr_char != EOF)
            curr_char = fgetc (file_handle);
        }
      }
      fclose (file_handle);


      //check to make sure the correct # of AA and codons are defined
      aa_processed = 0;
      for (i = 0; i < G.max_aa; i++) {
        if (AA[i].num_codons > 0) {
          aa_processed++;
        }
      }
      return aa_processed;  //return # of 

}

/*======================================================================================================*
 * FUNCTION: Process_R_input                                                                            *
 *      Purpose: This function is ONLY used in the sequence evolution                                   *
 *               code! It takes inputs from .C() R function (see                                        *
 *               http://users.stat.umn.edu/~geyer/rc/) and imports them                                 *
 *               to their corresponding variable names used in                                          *
 *               simulation_R-ext.c                                                                     *
 *      Inputs: CODON_INDEX - Codon Index Vector - vector containing integer values                     *
 *                respresenting each codon in the given sequence - Use codon index values               *
 *                found in preston/data/codon_summary.tsv                                               *
 *              ELONG_PR - Codon Elongation Vector - vector containing codon elongation         *
 *                rates for each codon in order of codon indices                                        *
 *              MUTATION_RATES - Mutation Rate Vector - vector containing mutation rates                *
 *                for each codon in order of codon indices                                              *
 *              PHI - Phi value for the given codon sequence                                            *
 *              RESIDENCE_TIMES_SIM_GENOTYPE - Initialized vector to fill in the time each simulation   *
 *                step takes- must have length equal to number of evol steps                            *
 *              DELTA_ETA - Initialized vector to fill in (eta_sim - eta_obs) for each step             *
 *              UNSCALED_PR_SIM_GENOTYPES - Initialized vector to fill in \mu*exp{-y(eta_sim - eta_obs} *
 *                at each simulation step                                                               *
 *              UNSCALED_PR_OBS_GENOTYPE - Variable to hold the numerator of the likelihood calculation *
 *                for the original sequence                                                             *
 *              ETA_REF - Variable to hold the reference value for scaling eta to avoid                 *
 *                problems with numerical precision - This should be 0 if first                         *
 *                simulation for sequence, and for subsequent runs, this should                         *
 *                be a value that is passed back from the first simultion.                              *
 *                Now, the first simulation will pass back the first value of eta                       *
 *                as the reference value, so ETA_REF = ETA_OBS                                          *
 *              POP - Effective Population Size- Parameter measuring genetic drift                      *
 *                DEFAULT = 1 (used to be 1.36*10^7 when dealing with phi)                              *
 *              A_1 - Ribosome initation cost in ATPs                                                   *
 *                DEFAULT a1=4                                                                          *
 *              A_2 - Ribosome elongation cost in ATPs, a2 in SEMPPR                                    *
 *                DEFAULT a2=4                                                                          *
 *              AT_BIAS - AT-Bias- Parameter related to mutation-- unused when including                *
 *                mutation file.                                                                        *
 *                DEFAULT = 0.5                                                                         *
 *              BENCH - Benchmark Flag- Flag that replaces rng numbers with predetermined               *
 *                numbers to set repeatable benchmarks                                                  *
 *                DEFAULT = 0                                                                           *
 *              BEE - Background Nonsense Error Rate B                                                  *
 *                DEFAULT = 0.0025                                                                      *
 *              GMT - Global Max Time- When the model is based on fixed time of evolution,              *
 *                this parameter limits the amount of time sequences are allowed to evolve.             *
 *                If the model is based on fixed number of evolution steps, this parameter              *
 *                does nothing.                                                                         *
 *                DEFAULT = -5                                                                          *
 *              IGNORE - When IGNORE is set equal to a positive integer n, nonsense errors can occur    *
 *                on the last n amino acids in a sequence and still produce a functional protein        *
 *                DEFAULT = 0                                                                           *
 *              RANDOM - Random Start Flag- Specifies whether the model should start with a random      *
 *                synonymous sequence or with the given sequence.                                       *
 *                DEFAULT = 0 (no random start)                                                         *
 *              GAMMA -Transition/Transversion Ratio- Parameter related to mutation-- unused when       *
 *                including mutation file.                                                              *
 *                DEFAULT = 1                                                                           *
 *              SCALE_FACTOR - Scaling Factor q- used to scale Effective Population Size Ne and         *
 *                Expression Rate phi.                                                                  *
 *                DEFAULT = 1 (used to be 4.18*10^-7 when dealing with phi)                             *
 *              MES - Minimum Evolution Steps- When the model is based on fixed number of evolution     *
 *                steps (substitutions), this paramter sets the number of evolution steps each          *
 *                sequence will be allowed.                                                             *
 *                DEFAULT = 10000                                                                       *
 *              BIS - Number of burn in steps per amino acid                                            *
 *                DEFAULT = 20                                                                          *
 *      Output: NONE
 *      Usage: refer to source/simulation_R-ext.c
 * =====================================================================================================*/


void Process_R_input(int *CODON_INDEX, double *ELONG_PR, 
             double *MUTATION_RATES, int *AA_COUNT, double *PHI, 
             double *POP, double *A_1, double *A_2,
             double *AT_BIAS, double *BEE,
             int *IGNORE, double *GAMMA,
             double *SCALE_FACTOR, char ** AA_VEC,
             char ** CODON_VEC, int *N_CODON){
    int i, j, counter;
    
    //*=============================================*
    //*               STRUCT SETUP                  *
    //* Initialize AA structure and match codons to *
    //* amino acids                                 *
    //*=============================================*
    
    i=0,j=0;
    char currAA;
    int nCodons=0;
    char aa_list[22] = { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'Z', 'S', 'T', 'V', 'W', 'Y', 'X' };    //'X' represents the 'stop' amino acid, 'Z' represents the smaller set of serine codons: AGT and AGC.  The other ser codons are TCT, TCC, TCG, and TCA.
    
    //* Process info from R vecs *//
    
    
    //* Need to re-initialize num aa and num codons after every call *//
    G.max_aa = 0;
    for(i=0;i<22;i++){
		AA[i].num_codons = 0;
	}
    
    i = 0;
    while(i < *N_CODON){
        j=0;  
        currAA = AA_VEC[i][0];     
        AA[G.max_aa].aa = currAA;
        
        while(i<*N_CODON&&*(AA_VEC+i)[0]==currAA){
            strcpy(AA[G.max_aa].codon[j],CODON_VEC[nCodons++]);     //transfer codon string from R vector
            AA[G.max_aa].codon[j][3] = '\0';                        // must end with this character
            
            AA[G.max_aa].num_codons++;
            
            j++,i++;
        }
        
        G.max_aa++;
    }
    
    //* Check for errors *//

    if(G.max_aa >22) {
        error("Number of amino acids is greater than maximum value of 22\n");
    }
    
    for(i=0;i<G.max_aa;i++){
        j=0;
        
        if(AA[i].num_codons > 6){
            Rprintf("Number of codons for amino acid %d: %c exceeds maximum number of 6\n",i,AA[i].aa);
            exit(EXIT_FAILURE);
        }
        
        while(AA[i].aa!=aa_list[j]&&j<22){
            j++;
        }
        if(j == 22){
            Rprintf("Amino Acid %d: %c did not match any know amino acids\n",i,AA[i].aa);
        }
    }
    
    //*=============================================*
    //*             GENOME PARAMTERS                *
    //* These parameters apply to the entire genome *
    //*=============================================*
    
    counter = 0;
	for(i=0;i<G.max_aa;i++) {
		for(j=0;j<AA[i].num_codons;j++){
			AA[i].elong_pr[j] = *(ELONG_PR + counter);
			AA[i].incoming_mu[j] = *(MUTATION_RATES + counter);
			counter++;
		}
	}

    
    
    
    G.Ne = *POP;
    G.A1 = *A_1;
    G.A2 = *A_2;
    G.Q = *SCALE_FACTOR;
 
    G.ignore_aa = *IGNORE;
    
    //*==============================================*
    //*             SEQUENCE PARAMETERS              *
    //* These parameters only apply to this sequence *
    //*==============================================*
    
    Sequence.aa_count = *AA_COUNT;
    Sequence.phi_obs = *PHI;
    for(i = 0;i < *AA_COUNT;i++) {
        Sequence.codon_index[i] = *(CODON_INDEX + i);
    }

    /* Trace code for debugging
       Rprintf("POP=\t%f\n",*POP);
       Rprintf("A_1=\t%f\n",*A_1);
       Rprintf("A_2=\t%f\n",*A_2);
       Rprintf("AT_BIAS=\t%f\n",*AT_BIAS);
       //Rprintf("BENCH = %d\n",*BENCH);
       Rprintf("BEE=\t%f\n",*BEE);
       Rprintf("IGNORE=\t%d\n",*IGNORE);
       Rprintf("GAMMA=\t%f\n",*GAMMA);
       Rprintf("SCALE_FACTOR=\t%fe-7\n",*SCALE_FACTOR*10000000);
     */

}

void Rprint_codon_sequence( int *codon_index ){
    int i;
    Rprintf("Codon Sequence Start:\n");
    for(i=0;i<Sequence.aa_count;i++){
        Rprintf("%d\t",*(codon_index + i));
    }
    Rprintf("\n\n");
}


/*==============================================================*
 * FUNCTION: Process_AA_Information                             *
 *      Purpose: Find max/min elongation rates, calc elong_pr,  *
 *               fail_pr                                        *
 *      Inputs: None                                            *
 *      Outputs: return(1) if the fcn successfully exits        *
 *      Usage: Process_AA_Information();                        *
 * =============================================================*/
 
int Process_AA_Information () {

      int i, j, k;
      int num_codons;
      double mu_in = 0;
      double sum_neutral_obs_pr;
      double min_elong_pr, max_elong_pr;

      //Added 3/2/09
      //check to make sure the correct # of AA and codons are defined
      j = 0;
      k = 0;
      
      for (i = 0; i < G.max_aa; i++) {
        k += AA[i].num_codons;
      }

      if (k != 64) {
          error("Only %d codons defined. Expecting 64. Exiting", k);
      }

      for (i = 0; i < G.max_aa - 1; i++) {
        num_codons = AA[i].num_codons;

        sum_neutral_obs_pr = 0; //added line 11/21/08

        for (j = 0; j < num_codons; j++) {
          
          AA[i].NSE_pr[j] = 1-AA[i].elong_pr[j];
          AA[i].omega[j]=AA[i].NSE_pr[j]/AA[i].elong_pr[j];
          //first part of calculating the pr of obs codon under neutral model but with mutation
          //set to 1 and then scale based on mutation bias
          AA[i].neutral_obs_pr[j] = 1;
	    }
	  }
      
      /*=============================================*
       * Function Call: Generate_Codon_Structures()  *
       * We need to generate these structures before *
       * calculating neutral observation pr b/c      *
       * they hold relevant mutation information     *
       *=============================================*/

      Generate_Codon_Structures();
      
      for (i = 0; i < G.max_aa - 1; i++) {
        num_codons = AA[i].num_codons;

        sum_neutral_obs_pr = 0; //added line 11/21/08

        for (j = 0; j < num_codons; j++) {
			
			if(AA[i].num_codons == 1){
				AA[i].neutral_obs_pr[j] = 1;
				sum_neutral_obs_pr=1;
			}else{
				mu_in = 0;
				for (k = 0; k < num_codons; k++) {
					  mu_in += AA[i].mu[k][j];
				}
				AA[i].neutral_obs_pr[j] = mu_in;
				sum_neutral_obs_pr+=mu_in;
			}  
        }

        //scale each codon by sum_neutral_obs_pr to ensure they sum to 1
        for (j = 0; j < num_codons; j++) {
          AA[i].neutral_obs_pr[j] /= sum_neutral_obs_pr;
          AA[i].neutral_obs_pr_cum[j] /= sum_neutral_obs_pr;
        }
      }
      
      //Assign min/max_elong_rate/pr
      for (i = 0; i < G.max_aa - 1; i++) {
        num_codons = AA[i].num_codons;

        min_elong_pr = AA[i].elong_pr[0];
        max_elong_pr = min_elong_pr;


        //start at 1 b/c already processed first codon
        for (j = 1; j < num_codons; j++) {

          if (min_elong_pr > AA[i].elong_pr[j]) {
            min_elong_pr = AA[i].elong_pr[j];
          }

          if (max_elong_pr < AA[i].elong_pr[j]) {
            max_elong_pr = AA[i].elong_pr[j];
          }

        }

        AA[i].min_elong_pr = min_elong_pr;
        AA[i].max_elong_pr = max_elong_pr;

      }


      return (1);       //exited successfully
}

void Calc_AA_Moments(){
	// Calculating the first and second moments for each amino acid
	int i,j;
      for (i = 0; i < G.max_aa - 1; i++) {

        AA[i].e_omega = 0;
        AA[i].e_invp = 0;
        AA[i].e_omega2 = 0;
        AA[i].e_invp2 = 0;
        AA[i].e_omegainvp = 0;
        for (j = 0; j < AA[i].num_codons; j++) {	
          AA[i].e_omega +=
            AA[i].omega[j] * AA[i].neutral_obs_pr[j];
          AA[i].e_invp +=
            1/AA[i].elong_pr[j] * AA[i].neutral_obs_pr[j];
          AA[i].e_omega2 +=
            (pow (AA[i].omega[j], 2)) * AA[i].neutral_obs_pr[j];
          AA[i].e_invp2 +=
            1/pow(AA[i].elong_pr[j],2) *
            AA[i].neutral_obs_pr[j];
          AA[i].e_omegainvp +=
            AA[i].omega[j]/AA[i].elong_pr[j] *
            AA[i].neutral_obs_pr[j];
            
        }
        
      }		
}

/*==============================================================*
 * FUNCTION: Generate_Codon_Structures                          *
 *      Purpose: Store aa_index,elong_rate/pr, and info about   *
 *               relative mutation in codon structure--         *
 *      Inputs: None                                            *
 *      Outputs: return(1) if fcn exits correctly               *
 * =============================================================*/

void Generate_Codon_Structures()
{
    int num_syn, num_one_step_syn;
    int num_two_step_synonym;
    double tsn_penalty = 100; //two step neighbor mutation rate will be smaller than osn by this factor
    int i, j, k, l, m;
    int num_matches;
    int codon_index, synon_index;
    int c_index;
    int codon_match[3]; //use to determine which position and nts involved in differences b/w codons.
    double mu_rate = 1;
    int to_codon, to_codon_matcher, from_codon;
    int osn_flag;
    int first_codon;
    int from_matrix_index, to_matrix_index;

    char codon[4];
    char synonym[4];

    codon_index = 0;
    to_codon = to_codon_matcher = from_codon = 0;
    for (i = 0; i < G.max_aa; i++) {
        for (j = 0; j < AA[i].num_codons; j++) {
            AA[i].codon_index[j] = codon_index;
            strcpy(Codon[codon_index].codon, AA[i].codon[j]);
            Codon[codon_index].aa = AA[i].aa;
            Codon[codon_index].aa_index = i;
            Codon[codon_index].incoming_mu = AA[i].incoming_mu[j];
            //Codon[codon_index].elong_rate = AA[i].elong_rate[j]; no longer using elong rates in C code
            Codon[codon_index].elong_pr = AA[i].elong_pr[j];
            Codon[codon_index].NSE_pr = AA[i].NSE_pr[j];
            Codon[codon_index++].omega = AA[i].omega[j];
            
        }
    }

    //Define synonyms
    for (i = 0; i < 64; i++) {
        num_syn = 0;
        for (j = 0; j < 64; j++) {
            if (j != i) {
                if (Codon[i].aa_index == Codon[j].aa_index) {
                    Codon[i].synonym_index[num_syn++] = j;
                }
            }
        }
        Codon[i].num_synonym = num_syn;
    }

    //Define one step synonyms
    for (i = 0; i < 64; i++) {
        num_syn = Codon[i].num_synonym;
        strcpy(codon, Codon[i].codon);
        num_one_step_syn = 0;

        //go through each syn
        for (j = 0; j < num_syn; j++) {
            num_matches = 0;

            synon_index = Codon[i].synonym_index[j];
            strcpy(synonym, Codon[synon_index].codon);
            //go through each nt in codon
            for (k = 0; k < 3; k++) {
                if (codon[k] == synonym[k]) {
                    num_matches++;
                    codon_match[k] = 1;
                } else {
                    codon_match[k] = 0;
                }
            }
            if (num_matches == 2) {
                Codon[i].one_step_synonym_index
                    [num_one_step_syn] = synon_index;

                //find mismatch
                k = 0;
                while (codon_match[k] == 1) {
                    k++;
                }

				//whewgley: We don't need this part of the code anymore 
				//since we are dealing with codon-specific mutation rates. 
				/*
                //mikeg: following code does two things: 
                // 1) Calculate mtuation rates under simple 
                //     AT bias and Transition vs. Transversion bias
                // 2) Assign relative mutation values for one step neighbors
                //Function 1) should be put into a separate function that's
                // used when no mutation file is read in and calculations should
                // should be stored in AA.mu[][] structures as is done when reading
                // mutation information in from a file
                //Function 2) should be done here and it should use the info 
                // from AA.mu

                //set initial relative mutation rate

                relative_mut_rate = 1;

                switch (mut_nt) {
                case 'A':
                case 'a':
                case 'T':
                case 't':
                    relative_mut_rate *= mu_bias;
                    break;

                default:
                    //Do nothing
                    break;
                    
                }   //end case

                //determine if rate is scaled by transition/transversion ratio: gamma_ratio
                if ((orig_nt == 'G' && mut_nt == 'A')
                    || (orig_nt == 'A' && mut_nt == 'G')
                    || (orig_nt == 'C' && mut_nt == 'T')
                    || (orig_nt == 'T' && mut_nt == 'C')) {
                    relative_mut_rate *= gamma_ratio;
                }

                Codon[i].one_step_relative_mutation_rate
                    [num_one_step_syn] = relative_mut_rate;
				*/
                num_one_step_syn++; //increase count

            }   //end match for one step synonym

        }
        Codon[i].num_one_step_synonym = num_one_step_syn;
    }

	
	/* ================================================================*
	 * Since we cannot find individual codon mutation rates, we assume *
	 * all incoming mutation rates from synonymous codons are equal.   *
	 * That is, for example, given a three codon case c_i,c_j,c_k each *
	 * with mu_i=(Sum of synonymous codon mutation rates  to codon i)= *
	 * mu(j to i)+mu(k to i), we say mu(j to i)=mu(k to i)=mu_i/2.     *
	 * Since stationary state depends only on mu_i, choosing arbitrary *
	 * values for synonymous codon mutation rates should conserve      *
	 * equilibrium. <whewgley> 7-26-12                                 *
	 * ================================================================*/
	
	//Initialize mu matrix in AA struct
	for(i=0;i<G.max_aa;i++){
		for(j=0;j<AA[i].num_codons;j++){
			for(l=0;l<AA[i].num_codons;l++){
				AA[i].mu[j][l]=0;
			}
		}
	}
	for (i = 0; i < G.max_aa; i++) {
		
		first_codon = AA[i].codon_index[0];
		
		for (j = 0; j < AA[i].num_codons; j++) {
			
			c_index = AA[i].codon_index[j];
			num_syn = Codon[c_index].num_synonym;
			num_one_step_syn = Codon[c_index].num_one_step_synonym;
			
			if(AA[i].num_codons > 1){
				if(num_syn==num_one_step_syn){
					//This is the mu rate FROM one-step synonyms TO Codon[c_index]
					mu_rate = Codon[c_index].incoming_mu/num_one_step_syn;
					to_codon = c_index;
				}else{
					num_two_step_synonym = num_syn - num_one_step_syn;
					
					//This is the mu rate FROM one-step synonyms TO Codon[c_index]
					//mu rate FROM two-step synonyms TO Codon[c_index] = mu_rate/tsn_penalty
					mu_rate = Codon[c_index].incoming_mu/(num_one_step_syn + num_two_step_synonym/tsn_penalty);
					to_codon = c_index;
				}
			}
			
			for(k = 0; k < Codon[c_index].num_synonym;k++){
				l=0;
				
				from_codon = Codon[c_index].synonym_index[k];
				
				//* Find out if codon is a one step neighbor *//
				osn_flag = 0;
				for(m=0;m<Codon[c_index].num_one_step_synonym;m++){
					if(Codon[c_index].one_step_synonym_index[m] == from_codon) osn_flag = 1;
				}
								
				to_codon_matcher = 0; //* This is what we will use to search for the "to" codon in the "from" codon structure
			   
				do{ to_codon_matcher = Codon[from_codon].synonym_index[l++];}while(to_codon_matcher != to_codon);
				l--; //counteract l++ in while loop
									
				//* from_codon = Codon[from_codon]
				//* to_codon = Codon[to_codon]
				//* from_matrix_index = from_codon minus first_to_codon
				//* to_matrix_index = codon_index minus first_codon_index
				
				from_matrix_index = from_codon - first_codon;
				to_matrix_index = to_codon - first_codon;
				
				if(osn_flag){                    
					Codon[from_codon].syn_relative_mutation_rate[l] = mu_rate; 
					AA[i].mu[from_matrix_index][to_matrix_index] = mu_rate;
				}else{
					Codon[from_codon].syn_relative_mutation_rate[l] = mu_rate/tsn_penalty;
					AA[i].mu[from_matrix_index][to_matrix_index] = mu_rate/tsn_penalty;
				}
			}
		}
	}


    
    //Define factors of elong_pr_ratio for each one step synonym
    for (i = 0; i < 64; i++) {
        num_syn = Codon[i].num_synonym;
        //go through each one step syn
        for (j = 0; j < num_syn; j++) {
            synon_index = Codon[i].synonym_index[j];
            Codon[i].elong_pr_ratio[j] =
                Codon[i].elong_pr / Codon[synon_index].elong_pr;
        }
    }

    //Define differences of (b/c_i - b/c_j) for each one step synonym
    for (i = 0; i < 64; i++) {
        num_syn = Codon[i].num_synonym;
        //go through each one step syn
        for (j = 0; j < num_syn; j++) {
            synon_index = Codon[i].synonym_index[j];
            
            Codon[i].delta_omega[j] =
                Codon[i].omega -
                 Codon[synon_index].omega;
        }
    }
}

/*========================================================*
 * FUNCTION: wrong()                                      *
 *      Purpose: This function prints an error message    *
 *               and exits the program. Used in           *
 *               Read_Commandline_Args and Read_tRNA_File *
 *      Inputs: NONE                                      *
 *      Outputs: NONE                                     *
 *      Usage: printf("Error...");                        *
 *             wrong();                                   *
 *========================================================*/
 
int
wrong () {
   error("\nError in command line:\nCES Exiting...\n");
}

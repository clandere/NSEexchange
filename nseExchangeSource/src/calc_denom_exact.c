#include "CES.h"
 
 void bin_eta(double eta, double * bins, double * bin_lims, double value, int num_bins);
 
 void intervals(double from, double to, int length, double * int_lims);


/*========================================*
 *      CALCULATE ENTIRE DENOMINATOR      *
 * This function will allow us to compare *
 * our estimates of the lik_denom to real *
 * value.                                 *
 *========================================*/

 void calc_denom_exact(int *CODON_INDEX,                  int *AA_COUNT,                    double *PHI,
                 double *ELONG_PR,          double *MUTATION_RATES,           double *POP,                      
                 double *A_1,                       double *A_2,                      double *AT_BIAS,
                 double *BEE,                       int *IGNORE,                      double *GAMMA,                    
                 double *SCALE_FACTOR,              double *LIK,                      double *ETA_MEAN,
                 double *ETA_VAR,                   double *ETA_OBS,                  double *ETA_MIN,
                 double *ETA_MAX,                   double *BINS,                     int *NUM_BINS,
                 char **AA_VEC,                      char **CODON_VEC,                 int *N_CODON){

    
    int i;
    int increment_next=0, counter=0;
    int aa_index,aa_position,aa_count;
    int num_codons;
    int size_codon_space = 1;
    int t_codon_cts[61] = {0} ,t_codon_cts_rel[61] = {0};
    int t_codon_index[MAX_AA]; //holds temporary codon_index
    
    double t_eta,t_mu,t_fitness,t_mu_fitness; 
    double fitness_total = 0, mu_fitness_total=0, eta_total=0, eta_sq_total=0;

    
    double pcnt=0;
    
    double answer[2];
    
    double * bin_lims;
    bin_lims = malloc(sizeof(double)*(*(NUM_BINS)+1));

    //*=================================*
    //* 1) Process input passed from R  *
    //*=================================*
    
    Process_R_input(CODON_INDEX,ELONG_PR,MUTATION_RATES,AA_COUNT,
                        PHI, POP, A_1, A_2, AT_BIAS,
                        BEE, IGNORE, GAMMA, SCALE_FACTOR, AA_VEC,
                        CODON_VEC, N_CODON);

    Process_AA_Information();

    //* initialize codon structures
    Generate_Codon_Structures();
    
    aa_count=Sequence.aa_count;
    
    //*=================================*
    //* 2) Find size of codon space     *
    //*=================================*

    Convert_Codon_Index_to_AA_Index(Sequence.codon_index,
                            Sequence.aa_index, aa_count);

    for(i=0;i<aa_count;i++){
        aa_index = Sequence.aa_index[i];
        num_codons = AA[aa_index].num_codons;
        size_codon_space *= num_codons;
    }
    
    //*=====================================*
    //* 3) Find eta min/max and define bin  *
    //*    limits for eta hist              *
    //*=====================================*
        
    eta_min_max(Sequence.aa_index,Sequence.aa_count,answer);
    *ETA_MIN = *(answer);
    *ETA_MAX = *(answer + 1);

    intervals(*ETA_MIN,*ETA_MAX,*(NUM_BINS)+1,bin_lims);
    
    //*=====================================*
    //* 4) Initialize observed sequence and *
    //*    create temporary codon index     *
    //*=====================================*
        
        //* 4.1) Calculate eta_obs
    Sequence.eta_obs = Calc_Eta_NSE(&Sequence);
    
    //Rprintf("Eta_obs: %f\tSigma_obs: %f\txi_obs: %f\taa_count: %d\n",Sequence.eta_obs,Sequence.sigma_obs,Sequence.xi_obs,aa_count);
        
        //* 4.2) Find initial codon counts
    Codon_Counts(Sequence.codon_index,Sequence.codon_cts,aa_count);
    
    //Rprintf("Sequence.sigma_obs: %f\nSequence.xi_obs %f\nSequence.eta_obs %f\n",Sequence.sigma_obs,Sequence.xi_obs,Sequence.eta_obs);
    
        //* 4.3) Initialize temporary codon index (this will be incremented in the for loop to 
        //*      represent each permutation of the codon sequence). 
    for(i=0;i<aa_count;i++){
        aa_index = Sequence.aa_index[i];
        t_codon_index[i] = AA[aa_index].codon_index[0];
    }
    
        //* 4.4) Find codon counts of t_codon_index
    Codon_Counts(t_codon_index,t_codon_cts,aa_count);
    
        //* 4.5) Find (t_codon_cts - Sequence.codon_cts): This will give the codon counts relative
        //*      to the original seqeunce. When we calculate Mu from these counts, it will give us
        //*      mu_sim/mu_obs. 
    for(i=0;i<61;i++){ //assume no stop codons
        t_codon_cts_rel[i] = t_codon_cts[i] - Sequence.codon_cts[i];
    }
    
    //*=====================================*
    //* 5) Big ass for loop to calculate    *
    //*    fitness of every codon sequence  *
    //*    and keep running total           *
    //*=====================================*
    
    //Rprintf("Calculating denominator with %d synonymous sequences...\n",size_codon_space);
    

    for(i=0;i<size_codon_space;i++) {
        
        //*5.1) Calc Eta

          t_eta = Calc_Eta_NSE (&Sequence);

          
        //* 5.3) Calculate Mu*Fitness
        //*      Fitness = exp(-Q*Ne*phi*eta) or exp(-y*eta)
        //*      Mu=\prod_{i=1}^n(\mu_i)
        
          t_mu = Calc_Seq_Mu(t_codon_cts_rel,aa_count); //This gives mu_sim/mu_obs
        
          t_fitness = exp(-G.Q*G.Ne*Sequence.phi_obs*(t_eta-Sequence.eta_obs)); //calculate relative to eta_obs
                                                                            //to avoid very small number                                                                 
          t_mu_fitness = t_mu*t_fitness;
          
          //* Bin eta
          bin_eta(t_eta,BINS,bin_lims,t_mu_fitness,*NUM_BINS);
          
        //* 5.4) Add fitness to total fitness
        //*      Denominator of liklihood function is the sum of all 
        //*      fitnesses in codon space
        //*      Also add eta to total eta to find avg eta in codon space
          
          fitness_total += t_fitness;
          mu_fitness_total += t_mu_fitness;
          eta_total += t_eta;
          eta_sq_total += t_eta*t_eta;
        
        
        //* 5.5) Increment codon sequence
        //*      This algorithm increments the codon sequence by incrementing the codon index
        //*      of the first position. If the first position is at the highest codon index
        //*      for that amino acid, it resets the first position to the lowest codon index
        //*      for that amino acid and increments the second position. If the second position 
        //*      is at the highest codon index for that amino acid... and so on. This works 
        //*      similar to an old-fashioned rolling wheel counter.
        //*      We also need to update codon counts here...
        
          aa_position=0;
          
          do {
              increment_next=0; //assume that we do not increment the next position
              aa_index = Sequence.aa_index[aa_position];
              num_codons = AA[aa_index].num_codons;
              if(t_codon_index[aa_position] < AA[aa_index].codon_index[num_codons-1]){
              t_codon_cts_rel[t_codon_index[aa_position]++]--; //This line subtracts one from the codon count of the current codon
                                                               //and adds one to the current codon index
              t_codon_cts_rel[t_codon_index[aa_position]]++;   //This line adds one to the codon count of the next codon
              }else{
                    t_codon_cts_rel[t_codon_index[aa_position]]--; //Subtract 1 from codon count of current codon
                    t_codon_index[aa_position] = AA[aa_index].codon_index[0];
                    t_codon_cts_rel[t_codon_index[aa_position]]++; //Add 1 to codon count of next codon
                    increment_next = 1; //turn on increment flag 
              }
              aa_position++; //
          }while(increment_next && aa_position < aa_count);

          counter++;
          
          
          if(counter>size_codon_space/(pcnt*(double)size_codon_space/100)){
              //Rprintf('%d Percent Complete.\n',pcnt);
              pcnt++;
          }
    }
    
    
    //*==========================*
    //* 6) Calculate likelihood  *
    //*==========================*
         
    double likelihood,eta_mean,eta_var;
     
    likelihood = 1/(mu_fitness_total/size_codon_space);
    eta_mean = eta_total/size_codon_space;
    eta_var =  eta_sq_total/size_codon_space - pow(eta_mean,2); //E[x^2] - E[x]^2
          
    *LIK = likelihood;
    *ETA_MEAN = eta_mean;
    *ETA_VAR = eta_var;
    *ETA_OBS = Sequence.eta_obs;
    
    free(bin_lims);
}

/*=======================================================================*
 *                      User Defined Functions                           *
 *=======================================================================*/
 
void bin_eta(double eta, double * bins, double * bin_lims, double value, int num_bins){
    //Purpose: Place eta into bin for histogram. This function will look for the bin interval
    //         (determined by bin_lims) that eta fits into and add the desired number 
    //         (value) to that bin.
    //Inputs: eta: eta value of interest
    //        bin_lims: interval limits for bins. must be in increasing order
    //        value: value added to bin 
    //        num_bins: number of bins, equal to length of bin_lims - 1
    //Outputs: None
    
    int i=0;
    
    
    //*=============================*
    //* Search for eta in bin_lims  *
    //*=============================*
    
    if(eta < *(bin_lims)){
        error("Observed Eta: %f less than min bin limit: %f\n",eta,*(bin_lims));
        
    }
    if(eta > *(bin_lims + num_bins)){
        error("Observed Eta: %f greater than max bin limit %f\n",eta, *(bin_lims+num_bins));
    }
    
      //* Eta should be between two values of bin_lims
    while((!(eta >= *(bin_lims + i-1) && eta < *(bin_lims + i)))&&i<=num_bins) i++;
    
    i--; //counteract i++
    
    //*===================*
    //* Add value to bin  *
    //*===================*
    //Rprintf("eta_bin_low: %f eta: %f eta_bin_high: %f\n",*(bin_lims + i),eta,*(bin_lims + i+1));
    *(bins + i) += value; 
        
}

void intervals(double from, double to, int length, double * int_lims){
    //Purpose: create array sequence of length 'length' from 'from' to 'to' 
    
    int i;
    double interval;
    
    interval = (to - from)/(length-1);
    
    *(int_lims) = from;
    
    for(i=1;i<length;i++){
        *(int_lims + i) = from + interval*i;
    }
    
    //check if last entry equal to 'to'
    if(*(int_lims + length-1)!=to) error("intervals function did not work\nLast entry: %f\tto: %f",*(int_lims+length-1),to);
    
    
}




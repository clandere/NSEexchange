#include "CES.h"

void calc_eta_moments(int *CODON_INDEX,                  int *AA_COUNT,                    double *PHI,
                      double *ELONG_PR,                  double *MUTATION_RATES,           double *POP,                      
                      double *A_1,                       double *A_2,                      double *AT_BIAS,
                      double *BEE,                       int *IGNORE,                      double *GAMMA,                    
                      double *SCALE_FACTOR,              double *ETA_MEAN,                 double *ETA_VAR, 
                      double *ETA_MIN,                   double *ETA_MAX,                  double *ETA_OBS,
                      char ** AA_VEC,                    char ** CODON_VEC,                int * N_CODON){
    //*
    //* Purpose: Calculate the moments (mean and variance) of the eta distribution analytically
    
    int i;
    int aa_count, aa_index_i, aa_index_ip1;
    double prod_e[MAX_AA]={0}, sq_prod_e[MAX_AA]={0}, prod_e_sq[MAX_AA]={0}; // prod_e[i] := \prod_{k=i+1}^n\frac{c_k+b}{c_k}
    double e_T[MAX_AA]={0}, varT[MAX_AA]={0}, H[MAX_AA]={0};
    double inner_loop_sum;
    double var_eta,mean_eta=0;
    double a1pa2i;
    double answer[2];


    //*=================================*
    //* 1) Process input passed from R  *
    //*=================================*
    
    Process_R_input(CODON_INDEX,ELONG_PR,MUTATION_RATES,AA_COUNT,
                        PHI, POP, A_1, A_2, AT_BIAS,
                        BEE, IGNORE, GAMMA, SCALE_FACTOR,AA_VEC,
                        CODON_VEC,N_CODON);


    Process_AA_Information(); //NOTE: this function call now also includes Generate_Codon_Structures()
    
    Calc_AA_Moments();
    
    aa_count = Sequence.aa_count;
        
    Sequence.aa_index[0]=1;
    //* Only need AA sequence now since overall eta distribution is independent of 
    //* observed sequence
    Convert_Codon_Index_to_AA_Index(Sequence.codon_index,
                            Sequence.aa_index, aa_count);
    
    /*======================*
     * 2) Calculate Eta obs *
     *======================*/
     
     Sequence.eta_obs = Calc_Eta_NSE(&Sequence);
                            
    //*========================================*
    //* 2) Iterate backwards to calculate ish  *
    //*========================================* 
    
    i = aa_count-1; //Not sure why we start at aa_count - 2, but this gives us the same eta variance
                  //calculated in msaum's semppr_2 code (originally aa_count-1 in msaum's code)
                  //8-14-13  whewgley - Changed back to aa_count - 1... the values calculated here deviated from
                  //msaums' semppr_2 code for some reason, but I changed this back to aa_count - 1 and now they are
                  //equivalent again. This seems to rely on whether the fasta read by semppr.data.load contains 
                  //a stop codon or not. 
                  
    aa_index_i = Sequence.aa_index[i];
    
    prod_e[i] = sq_prod_e[i] = prod_e_sq[i] = 1; //by definition
    
    e_T[i] = prod_e[i]*AA[aa_index_i].e_omega;
    varT[i] = AA[aa_index_i].e_omega2*prod_e_sq[i] - pow(AA[aa_index_i].e_omega*prod_e[i],2);
    H[i] = AA[aa_index_i].e_omegainvp*prod_e_sq[i]/(prod_e[i]*AA[aa_index_i].e_invp);
    
    i--;
    
    while(i > 0){
        aa_index_ip1 = aa_index_i;
        aa_index_i = Sequence.aa_index[i];
        
        prod_e[i] = prod_e[i+1]*AA[aa_index_ip1].e_invp;
        prod_e_sq[i] = prod_e_sq[i+1]*AA[aa_index_ip1].e_invp2;
        sq_prod_e[i] = sq_prod_e[i+1]*AA[aa_index_ip1].e_invp*AA[aa_index_ip1].e_invp;
        
        e_T[i] = prod_e[i]*AA[aa_index_i].e_omega;
        
        varT[i] = AA[aa_index_i].e_omega2*prod_e_sq[i] - pow(AA[aa_index_i].e_omega*prod_e[i],2);
        H[i] = AA[aa_index_i].e_omegainvp*prod_e_sq[i]/(prod_e[i]*AA[aa_index_i].e_invp);
        
        i--;
    }
    
    //Assume \sum_{j=n+1}^n(a_1+a_2j)Cov(Y_i,Y_j) = 0
    
    //*===================================================*
    //* 3) Iterate backwards again to calculate more ish  *
    //*===================================================*
    
    inner_loop_sum = 0;
    var_eta = 0;
    
    i = Sequence.aa_count - 1;
    
    a1pa2i = G.A1 + G.A2*(i);
      
    var_eta += a1pa2i*a1pa2i*varT[i]+2*a1pa2i*e_T[i]*inner_loop_sum;
    mean_eta += a1pa2i*e_T[i];
        
    inner_loop_sum += a1pa2i*(H[i]-e_T[i]);
    
    i--;
    
    double first_term=0, second_term=0, sum_first_term=0, sum_second_term=0;
    
    while(i > 0)
    {
        a1pa2i -= G.A2;
        first_term = a1pa2i*a1pa2i*varT[i];
        sum_first_term += first_term;
        
        second_term = 2*a1pa2i*e_T[i]*inner_loop_sum;
        sum_second_term += second_term;
        
        var_eta += a1pa2i*a1pa2i*varT[i]+2*a1pa2i*e_T[i]*inner_loop_sum;
        mean_eta += a1pa2i*e_T[i];
        
        inner_loop_sum += a1pa2i*(H[i]-e_T[i]);
        
        i--;
    }
    
    mean_eta = mean_eta + G.A1+G.A2*(aa_count);
    
    //*==========================*
    //* 4) Find min and max eta  *
    //*==========================*
    
    eta_min_max(Sequence.aa_index,Sequence.aa_count,answer);
    
    //*========================================*
    //* 5) Fill in variables to pass back to R *
    //*========================================*
    
    *ETA_VAR = var_eta;
    *ETA_MEAN = mean_eta;
    *ETA_MIN = *(answer);
    *ETA_MAX = *(answer + 1);
    *ETA_OBS = Sequence.eta_obs;
                         
}

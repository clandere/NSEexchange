#ifndef CES_Header_
#define CES_Header_

/*===========================*
 *  Declaring Header Files   *
 *===========================*/

#include <stdio.h>
#include <stdlib.h>
//#include <sys/stat.h>       //for mkdir
//#include <sys/types.h>      //for types w/in mkdir
#include <math.h>
#include <time.h>
#include <sys/time.h>       //gettimeofday()
#include <string.h>
//#include <omp.h>
//#include <gsl/gsl_errno.h>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_integration.h>
//#include <gsl/gsl_min.h>
//#include <gsl/gsl_multimin.h>
//#include <gsl/gsl_roots.h>
#include <gsl/gsl_rng.h>    //uniform rng
#include <gsl/gsl_randist.h>    //rng from distns
#include <R.h>

#include "shared/CES_struct.h"
#include "shared/iofuns_CES.h"
#include "shared/synfuns_CES.h"
#include "shared/seqfuns_CES.h"

/*================*
 * Define Structs *
 *================*/

struct amino_acid AA[22];
struct codon_struct Codon[64];
struct seq_struct Sequence;
struct global_variable_struct G;

#endif

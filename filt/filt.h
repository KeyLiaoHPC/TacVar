/**
 * @file filtv3.h 
 * @author Key Liao
 * 
 */

#include <stdint.h>

#define i64 int64_t
#define u64 uint64_t
#define i32 int32_t
#define u32 uint32_t

/*=== BEGIN: Global Variables ===*/


typedef struct FilT_Param_T {
    i64 ref_time;           // Reference running time (ns) of timing fluctuations.
    i64 rnd;                // Rounding bin edge to an integer multiple of rnd.
    u64 nsamp, nsim;        // Number of samples in a simulation, and the number of simulation for a bin
    i64 step;
    int met_col, tf_col;
    double p_low;           // Lower threshold of probability ofa single bin
    double p_xcut, p_ycut;  // Cut the highest fraction.
    double hz_ns;           // tick per nanosecond
    char *in_tm_file;      // Input measurement results.
    char *in_tf_file;       // Input timing fluctuation samples.
    char *out_tr_file;     // Output real time estimation.
} filt_param_t;

typedef struct Prob_Bin_T {
    i64 t;
    double p;
} prob_bin_t;

typedef struct Prob_Hist_T {
    u64 nbin;
    prob_bin_t *pbin;
} prob_hist_t;

//extern int errno;

/*=== END: Global Variables ===*/

/*=== BEGIN: Interfaces ===*/

/**
 * Parsing input command-line arguments.
 * @param argc
 * @para argv
 */
// int arg_parse(int argc, char **argv);

/**
 * Slicing a rounded array to bins.
 * The min width is larger than args->rnd, the min probability is larger than args->p_low
 * @param arr Input array being sliced.
 * @param len The length of the input array.
 * @param rnd.
 * @param p_low.
 * @param p_cut.
 * @param phist A pointer to the histogram of sliced running times.
 */
int slice(i64 *arr, u64 len, i64 rnd, double p_low, double p_cut, prob_hist_t *phist);

/**
 * Parsing csv file, extracting timing number from csv and slicing the results.
 * @param fpath File path.
 * @param col The target column in csv file started from 0.
 * @param adj Adjustment of array
 * @param rnd Round the edge of bins.
 * @param hz_ns tick per ns.
 * @param p_low The lower threshold of probability per bin.
 * @param p_cut The highest probability being cut.
 * @param arr The sorted array.
 * @param len The length of returned array.
 * @param hist The parsed and sliced histogram.
 */
int parse_csv(char *fpath, int col, i64 adj, i64 rnd, double hz_ns, double p_low, double p_hcut, 
    u64 *len, i64 **arr,  prob_hist_t *hist);


/**
 * Performing filtering with FilT, enclosing subroutines
 * @param args FilT parameters.
 * @param tr_hist Pointer to the histogram of filtered running times.
 */
int run_filt(filt_param_t *args, prob_hist_t *tr_hist);

/**
 * Comparing function for qsort.
 */
int cmp(const void *a, const void *b);

/**
 * @brief 
 * @param tmh   Binned meausrement times.
 * @param tfh   Binned timing fluctuation samples.
 * @param trh   Binned real time estimations.
 * @param tf_arr    Sorted timing fluctuation raw array after cut.
 * @param tf_len    The length of tf_arr.
 * @param args      
 * @return int 
 */
int calc_tr(prob_hist_t *tmh, prob_hist_t *tfh, prob_hist_t *trh, i64 *tf_arr, u64 tf_len, filt_param_t *args);

/**
 * @brief 
 * @param tmh 
 * @param trh 
 * @param tf_arr 
 * @param tf_len 
 * @param psim 
 * @param args 
 * @return int 
 */
int sim_met(prob_hist_t *tmh, prob_hist_t *trh, i64 *tf_arr, u64 tf_len, double *psim, filt_param_t *args);

/*=== END: Interfaces ===*/

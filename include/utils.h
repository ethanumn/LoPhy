#pragma once 

// preprocessor directives
#define NO_CHILDREN -1
#define ROOT 0
#define NO_PARENT -1

#include <iostream>
#include <string>
#include <vector>
#include <span>
#include <cstdint>
#include <stdexcept>
#include <numeric>
#include <algorithm>
#include <random>
#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include <stack>
#include <tuple>
#include <queue>
#include <span>
#include "custom_unordered_map.h"


using std::queue;
using std::stack;
using std::deque;
using std::mt19937;
using std::stable_sort;
using std::iota;
using std::cout; 
using std::endl; 
using std::string;
using std::vector;
using std::invalid_argument;
using std::unordered_set;
using std::unordered_map;
using std::max_element;
using std::accumulate;
using std::numeric_limits;
using std::sort;
using std::span;

using std::tgamma; 
using std::lgamma;
using std::tuple;
using std::get;
using std::pair; 
using std::make_tuple;
using std::span;
using std::hash;

/*------------------*/
// DATA STRUCTURES //
/*-----------------*/
enum CNA_TYPE {
  NA,
  LOSS_ALT,
  LOSS_REF,
  GAIN_ALT,
  GAIN_REF,
  CNLOH_ALT,
  CNLOH_REF
};


typedef pair<int, CNA_TYPE> CNA; // <region, CNA_TYPE> 
typedef vector<vector<CNA>> CNAs; // <clone, <region, CNA_TYPE>>
typedef custom_unordered_map<double, vector<double>> CACHE;
const double NINF = -numeric_limits<double>::infinity();

// hashing function for pair
namespace std {
    template<>
    struct hash<std::pair<int, CNA_TYPE>> {
        std::size_t operator()(const std::pair<int, CNA_TYPE>& p) const noexcept {
            std::size_t h1 = std::hash<int>{}(p.first);
            std::size_t h2 = std::hash<CNA_TYPE>{}(p.second);
            return h1 ^ (h2 << 1); // combine the two
        }
    };
}

/*-----------------*/
// PRINT FUNCTIONS //
/*-----------------*/

/* Print functions source */
template<typename T>
void print(T & v, string comment)
{
    cout << comment << endl;
    cout << "[";
    for(int i = 0; i < static_cast<int>(v.size()); ++i)
    {
        cout << v[i];
        if (i < v.size() - 1)
        {
            cout << ",";
        }
    }
    cout << "]" << endl;
}

template<typename T>
void print_vector(const vector<T> & v, string comment)
{
    print(v, comment);
}

template<typename T>
void print_set(const unordered_set<T> & v)
{
    for (auto it = v.begin(); it != v.end(); ++it)
    {
        cout << *it;
        if (it != v.end())
            cout << ", ";
    }
    cout << endl;
}

template<typename T>
void print_matrix(const vector<vector<T>> & matrix, string comment)
{
    cout << comment << endl;
    for(const vector<T> & v : matrix)
        print_vector(v, "");
}

/*-----------*/
// UTILITIES //
/*-----------*/

/* Sorting functions source */
template<typename T>
vector<int> argsort(const vector<T> & arr) {

    vector<int> sorted_indices(arr.size());
    iota(sorted_indices.begin(), sorted_indices.end(), 0);

    stable_sort(
        sorted_indices.begin(), sorted_indices.end(),
        [&arr](int i, int j) {return arr[i] > arr[j];});

    return sorted_indices;
}

/* Log-sum-exp */
inline double logsumexp(const vector<double> & x) 
{
    double max_val = *max_element(x.begin(), x.end());
    // avoid returning nan
    if (max_val == NINF) 
        return NINF; 
    else 
    {
        double sum = 0.0;
        for (double val : x)
            sum += exp(val - max_val);
        return max_val + log(sum);
    }
}

/* Calculates posterior probabilities */
inline vector<double> posterior(const vector<double> & vec) 
{
    const int size = vec.size();
    vector<double> result(size, 0.0);

    // find max
    double logprob = logsumexp(vec);

    if(logprob == NINF)
    {
        for(double & val : result)
            val = 1.0 / size;
    }

    for (int i = 0; i < size; ++i) {
        result[i] = exp(vec[i] - logprob);  
    }

    return result;
}

/* softmax function */
inline vector<double> softmax(const vector<double> & vec) 
{
    const int size = vec.size();
    vector<double> result(size);

    // find max
    double max_val = *max_element(vec.begin(), vec.end());

    if(max_val == NINF)
    {
        for(double & val : result)
            val = 1.0 / size;
    }

    // compute denominator
    double sum_exp = 0.0;
    for (int i = 0; i < size; ++i) {
        result[i] = exp(vec[i] - max_val);  
        sum_exp += result[i];
    }

    // compute numerator
    for (double & val : result) {
        val /= sum_exp;
    }

    return result;
}

/* Returns true if vec1 is a subset of vec2 */
inline bool vector_subset(const vector<int> & vec1, const vector<int> & vec2) 
{
    for (const auto val : vec1) {
        if (find(vec2.begin(), vec2.end(), val) == vec2.end()) {
            return false; 
        }
    }
    return true; 
}

// beta function 
inline double log_beta(double x, double y) 
{
    return lgamma(x) + lgamma(y) - lgamma(x + y);
}

// for computing log binomial coefficient
inline double log_binom_coeff(double n, double k) noexcept
{
    return lgamma(n+1)-lgamma(n-k+1)-lgamma(k+1);
}

/* Discretizes vector<double> to vector<int> */
inline vector<int> discretize_vector(const vector<double> & vec) 
{
    vector<int> discretized;
    discretized.reserve(vec.size()); // Preallocate memory for efficiency

    for (double val : vec) {
        discretized.push_back(static_cast<int>(round(val * 1000)));
    }

    return discretized;
}


/*----------------*/
// HASH FUNCTIONS //
/*---------------*/
inline int hash_function1(const int Cr, 
                               const int Ca, 
                               const double dropout_rate, 
                               const int locus)
{
    return Cr + 20*Ca + 400*dropout_rate + 8000*locus;
}

inline int hash_function2(const int region, 
                               const double region_proportion)
{
    return region + 500*round(region_proportion*1000);
}

inline int hash_function3(const vector<int> & parents, 
                               const int Cr, 
                               const int Ca, 
                               const double dropout_rate, 
                               const int locus)
{
    uint64_t hash = 0;

    // Hash the original components
    hash += Cr;
    hash += 20ull * Ca;
    hash += 400ull * static_cast<int>(round(dropout_rate * 1000)); // scale for precision
    hash += 8000ull * locus;

    // hash parents vector
    for (int i = 0; i < static_cast<int>(parents.size()); ++i) 
    {
        int val = parents[i];
        uint64_t mapped_val = static_cast<uint64_t>(static_cast<int64_t>(val) + 100000); 
        hash = hash * 31 + mapped_val;
    }

    return static_cast<int>(hash % UINT32_MAX);
}


/* Hash function 4 */
inline int hash_function4(const vector<int> & parents, 
                               const int Ca, 
                               const int Cr, 
                               const double alt_dropout_rate,
                               const double ref_dropout_rate, 
                               const int locus)
{
    uint64_t hash = 0;

    hash += Cr;
    hash += 20ull * Ca;
    hash += 400ull * static_cast<int>(round(ref_dropout_rate * 1000));  
    hash += 8000ull * static_cast<int>(round(alt_dropout_rate * 1000));
    hash += 16000000ull * locus;  

    for (int i = 0; i < static_cast<int>(parents.size()); ++i) 
    {
        int val = parents[i];
        uint64_t mapped_val = static_cast<uint64_t>(static_cast<int64_t>(val) + 100000); 
        hash = hash * 31 + mapped_val;
    }

    return static_cast<int>(hash % UINT32_MAX);
}


/*----------------------*/
// LIKELIHOOD FUNCTIONS //
/*----------------------*/

/* Calculates unnormalized likelihood that a region has a specific proportion of reads fall into it */
inline const vector<double> & cn_llh(CACHE & cn_llh_cache,
                                     const int & region, 
                                     const double & region_proportion,
                                     span<const vector<int>> region_reads, 
                                     span<const int> cell_sums,
                                     const double & theta)
{
    int hash = hash_function2(region, region_proportion);

    // return cached log likelihoods if they exist
    if(cn_llh_cache.has_key(hash))
        return cn_llh_cache[hash];
    else
    {
        // compute copy number log likelihoods from scratch
        int n_cells = cell_sums.size();

        vector<double> llhs(n_cells, NINF);
        for (int cell = 0; cell < n_cells; ++cell)
        {
            double mu = cell_sums[cell] * region_proportion;

            llhs[cell] = lgamma(region_reads[cell][region] + theta)  - lgamma(theta) - lgamma(region_reads[cell][region] + 1) 
                                + theta * log(theta / (theta + mu)) 
                                + region_reads[cell][region] * log(mu / (theta + mu));
        }
        cn_llh_cache[hash] = llhs;
        return cn_llh_cache[hash];
    }
}


/* Evaluate likelihood that cells have locus Ca alternate alleles and Cr reference alleles  */
inline const vector<double> & calculate_locus_llh(CACHE & locus_llh_cache,
                                                  CACHE & dropped_alt_alleles_cache,
                                                  CACHE & dropped_ref_alleles_cache,
                                                  const vector<int> & parents,
                                                  const int & locus,
                                                  const int & _Ca,
                                                  const int & _Cr,
                                                  const double & alt_dropout_rate,
                                                  const double & ref_dropout_rate,
                                                  span<const vector<int>> variant_reads, 
                                                  span<const vector<int>> total_reads, 
                                                  const double & fp, 
                                                  const double & fn,
                                                  const double & hom_precision,
                                                  const double & het_precision) 
{
    int Ca = _Ca, Cr = _Cr;
    if (Cr == 0) Ca = 1;
    else if (Ca == 0) Cr = 1;

    int hash = hash_function4(parents, Ca, Cr, alt_dropout_rate, ref_dropout_rate, locus);

    // return cached log likelihoods if they exist
    if(locus_llh_cache.has_key(hash) && dropped_alt_alleles_cache.has_key(hash)  && dropped_ref_alleles_cache.has_key(hash))
        return locus_llh_cache[hash];
    else 
    {
        // otherwise compute log likelihoods from scratch
        const int n_cells = variant_reads.size();
        double alpha, beta, dropout_prob, read_count_prob, precision, pCa, pCr, f;
        int Ca_amp, Cr_amp, Ct_amp, Ca_lost, Cr_lost, term;

        // compute only once
        const double & log_all_alleles_dropped_prob = log(1 - pow(alt_dropout_rate, Ca)*pow(ref_dropout_rate, Cr));

        // reset and resize term_llhs 
        const int num_terms = (Ca + 1) * (Cr + 1) - 1;
        vector<vector<double>> term_llhs(n_cells, vector<double>(num_terms, NINF));
        term = 0;
        
        for(Ca_amp = 0; Ca_amp <= Ca; ++Ca_amp)
        {
            for(Cr_amp = 0; Cr_amp <= Cr; ++Cr_amp)
            {
                if(Ca_amp == 0 && Cr_amp == 0) continue;
                else 
                {
                    // set precision accordingly
                    if(Ca_amp == 0 || Cr_amp == 0)
                        precision = hom_precision;
                    else    
                        precision = het_precision;

                    // dropout probability
                    Ca_lost = Ca - Ca_amp;
                    Cr_lost = Cr - Cr_amp;
                    dropout_prob = log_binom_coeff(Ca, Ca_amp) + log_binom_coeff(Cr, Cr_amp) + Ca_lost*log(alt_dropout_rate)
                                    + Cr_lost*log(ref_dropout_rate) + Ca_amp*log(1-alt_dropout_rate) + Cr_amp*log(1-ref_dropout_rate) - log_all_alleles_dropped_prob;   
                    
                    // parameters for read count probability
                    Ct_amp = Ca_amp + Cr_amp;
                    pCa = (double)Ca_amp / (double)Ct_amp;
                    pCr = (double)Cr_amp / (double)Ct_amp;
                    f = pCa*(1-fn) + pCr*fp;
                    alpha = f * precision;
                    beta = (1-f) * precision;

                    for (int cell = 0; cell < n_cells; ++cell) 
                    {
                        // read counts for cell i, mutation j
                        const int & k = variant_reads[cell][locus]; 
                        const int & n = total_reads[cell][locus];  

                        // log likelihood of dropout + read counts 
                        read_count_prob = lgamma(k + alpha) + lgamma(n - k + beta) - lgamma(n + precision) - lgamma(alpha) - lgamma(beta) + lgamma(precision);
                        term_llhs[cell][term] = dropout_prob + read_count_prob; 
                    }
                    term++;
                }
            }
        }

        // log likelihood for each cell that considers possible dropouts
        vector<double> llhs(n_cells, NINF);
        for(int cell = 0; cell < n_cells; ++cell)
            llhs[cell] = logsumexp(term_llhs[cell]);

        vector<double> E_dropped_alt_alleles(n_cells, 0.0);
        vector<double> E_dropped_ref_alleles(n_cells, 0.0);
        term = 0;
        for(Ca_amp = 0; Ca_amp <= Ca; ++Ca_amp)
        {
            for(Cr_amp = 0; Cr_amp <= Cr; ++Cr_amp)
            {
                if(Ca_amp == 0 && Cr_amp == 0) continue;
                else 
                {
                    for(int cell = 0; cell < n_cells; ++cell)
                    {
                        // calculate likelihood for cell attachment to clone
                        E_dropped_alt_alleles[cell] += exp(term_llhs[cell][term] - llhs[cell]) * (double)(Ca - Ca_amp);
                        E_dropped_ref_alleles[cell] += exp(term_llhs[cell][term] - llhs[cell]) * (double)(Cr - Cr_amp);
                    }
                    term++;
                }
            }
        }

        locus_llh_cache[hash] = llhs;
        dropped_alt_alleles_cache[hash] = E_dropped_alt_alleles;
        dropped_ref_alleles_cache[hash] = E_dropped_ref_alleles;
        return locus_llh_cache[hash];
    }
}

/*-------------------------*/
// PREPROCESSING FUNCTIONS //
/*-------------------------*/


/* Sum function source */
template<typename T>
vector<T> sum(vector<vector<T>> matrix, int16_t axis, T exclude = -1) 
{
    // throw error if there's no data to sum over
    if(matrix.size() == 0)
        throw invalid_argument("Matrix is empty!");

    vector<T> axis_sum;
    int n_cols = matrix[0].size();

    // sum all entries
    if(axis == 0)
    {
        axis_sum.resize(matrix[0].size(), 0);
    }
    else if (axis == 1 || axis == -1)
    {
        axis_sum.resize(matrix.size(), 0);
    }
    else 
    {
        throw invalid_argument("axis parameter is invalid, must be one of {-1, 0, 1}");
    }

    // sum over values
    for(int i = 0; i < static_cast<int>(matrix.size()); i++)
    {
        if(axis == 1)
        {
            for(int j = 0; j < n_cols; j++)
            {
                if(matrix[i][j] != exclude)
                    axis_sum[i] += matrix[i][j];
            }
        }
        else 
        {
            for(int j = 0; j < n_cols; j++)
            {
                if(matrix[i][j] != exclude)
                    axis_sum[j] += matrix[i][j];
            }

        }
    }

    if(axis == -1)
    {
        T result = reduce(axis_sum.begin(), axis_sum.end());
        return {result};
    }
    else 
    {
        return axis_sum;
    }

}

inline unordered_map<int, vector<int>> order_variants(const vector<vector<int>> & character_matrix, 
                                                                const vector<int> & variant_samples)
{
    const int m = static_cast<int>(variant_samples.size());
    vector<int> order(m);
    iota(order.begin(), order.end(), 0);
    vector<int> row_sum = sum(character_matrix, 1);

    // Sort by sample index ascending, then by row_sum in descending order
    sort(order.begin(), order.end(), [&](size_t i, size_t j) {
        if (variant_samples[i] != variant_samples[j])
            return variant_samples[i] < variant_samples[j]; // loci that appear in earlier samples first
        return row_sum[i] > row_sum[j]; // loci that are mutated in more cells appear first
    });

    // create vectors to contains order of variants for each sample
    unordered_map<int, vector<int>> variant_orders;
    for (int i = 0; i < m; ++i) 
        variant_orders[variant_samples[order[i]]].push_back(order[i]);

    return variant_orders;
}

/*----------------*/
// TREE FUCNTIONS //
/*---------------*/

/* Computes combinations of elements */
void comb(vector<vector<int>> & choices,
          const vector<int> & arr, 
          const int & k);

/* Gets all descendants of a node in a parents vector */
void get_descendants(unordered_set<int> & descendant_list, 
                     const vector<int> & parents, 
                     const int & node);

/* Gets all ancestors of a node in a parents vector */
void get_ancestors(unordered_set<int> &     ancestors_list, 
                   const vector<int> & parents, 
                   const int & node);

/* breadth-first topological sort from parents vector */
inline vector<int> topological_sort(const vector<int> & parents, const int & m) 
{
    int num_clones = parents.size() + 1;  
    vector<int> topo_order = {ROOT}; 
    unordered_map<int, vector<int>> children_map; 

    // children map from parents vector
    for (int i = 0; i < num_clones - 1; ++i) 
    {
        const int p = parents[i];
        if (p != ROOT && p <= 2*m && p != NO_PARENT)
            children_map[p].push_back(i+1); 
    }

    // find nodes that are children of the root
    deque<int> node_queue;
    for (int i = 0; i < num_clones - 1; ++i) {
        if (parents[i] == ROOT) 
            node_queue.push_back(i+1);  
    }

    // perform BFS 
    while (!node_queue.empty()) 
    {
        int node = node_queue.front();
        node_queue.pop_front();
        topo_order.push_back(node);  

        // add all children to queue
        for (int child : children_map[node]) 
            node_queue.push_back(child);

    }
    return topo_order;
}

/*---------*/
// STRUCTS //
/*--------*/

/* Struct for caching computations */
struct SCORE_CACHE
{
    const vector<double> & get_dropped_alleles(const bool & use_alt_cache,
                                               const vector<int> & parents,
                                               const int & locus,
                                               const int & _Ca,
                                               const int & _Cr,
                                               const double & _alt_dropout_rate,
                                               const double & _ref_dropout_rate,
                                               span<const vector<int>> variant_reads, 
                                               span<const vector<int>> total_reads, 
                                               const double & fp, 
                                               const double & fn,
                                               const double & hom_precision,
                                               const double & het_precision)
    {
        // if homozygous the likelihood doesn't change regardless of the number of other alleles
        int Ca = _Ca, Cr = _Cr;
        double alt_dropout_rate = _alt_dropout_rate, ref_dropout_rate = _ref_dropout_rate;
        if(Cr == 0) Ca = 1;
        else if(Ca == 0) Cr = 1;

        if(Ca == 0 || Cr == 0)
        {
            alt_dropout_rate = 0.1;
            ref_dropout_rate = 0.1;
        }

        int hash = hash_function4(parents, Ca, Cr, alt_dropout_rate, ref_dropout_rate, locus);
        if(use_alt_cache)
        {
            if(!this->E_dropped_alt_alleles.has_key(hash))
            {
                calculate_locus_llh(this->locus_llh,
                                    this->E_dropped_alt_alleles,
                                    this->E_dropped_ref_alleles,
                                    parents,
                                    locus,
                                    Ca,
                                    Cr,
                                    alt_dropout_rate,
                                    ref_dropout_rate,
                                    variant_reads, 
                                    total_reads, 
                                    fp, 
                                    fn,
                                    hom_precision,
                                    het_precision);
            }
            return this->E_dropped_alt_alleles[hash];

        }
        else 
        {
            if(!this->E_dropped_ref_alleles.has_key(hash))
            {
                calculate_locus_llh(this->locus_llh,
                                    this->E_dropped_alt_alleles,
                                    this->E_dropped_ref_alleles,
                                    parents,
                                    locus,
                                    Ca,
                                    Cr,
                                    alt_dropout_rate,
                                    ref_dropout_rate,
                                    variant_reads, 
                                    total_reads, 
                                    fp, 
                                    fn,
                                    hom_precision,
                                    het_precision);
            }
            return this->E_dropped_ref_alleles[hash];

        }

    }


    const vector<double> & get_dropped_alt_alleles(const vector<int> & parents,
                                                   const int & locus,
                                                   const int & Ca,
                                                   const int & Cr,
                                                   const double & alt_dropout_rate,
                                                   const double & ref_dropout_rate,
                                                   span<const vector<int>> variant_reads, 
                                                   span<const vector<int>> total_reads, 
                                                   const double & fp, 
                                                   const double & fn,
                                                   const double & hom_precision,
                                                   const double & het_precision)
    {
        return get_dropped_alleles(true,
                                   parents,
                                   locus,
                                   Ca,
                                   Cr,
                                   alt_dropout_rate,
                                   ref_dropout_rate,
                                   variant_reads, 
                                   total_reads, 
                                   fp,
                                   fn,
                                   hom_precision,
                                   het_precision);

                                   
    }


    const vector<double> & get_dropped_ref_alleles(const vector<int> & parents,
                                                   const int & locus,
                                                   const int & Ca,
                                                   const int & Cr,
                                                   const double & alt_dropout_rate,
                                                   const double & ref_dropout_rate,
                                                   span<const vector<int>> variant_reads, 
                                                   span<const vector<int>> total_reads, 
                                                   const double & fp, 
                                                   const double & fn,
                                                   const double & hom_precision,
                                                   const double & het_precision)
    {
        return get_dropped_alleles(false,
                                   parents,
                                   locus,
                                   Ca,
                                   Cr,
                                   alt_dropout_rate,
                                   ref_dropout_rate,
                                   variant_reads, 
                                   total_reads, 
                                   fp,
                                   fn,
                                   hom_precision,
                                   het_precision);
    }
            
    void clear(void) 
    {
        this->locus_llh.clear();
        this->cn_llh.clear();
        this->E_dropped_alt_alleles.clear();
        this->E_dropped_ref_alleles.clear();
    }
    CACHE locus_llh;
    CACHE cn_llh;
    CACHE E_dropped_alt_alleles;
    CACHE E_dropped_ref_alleles;

};

/* Struct for packaging parameters */
struct Params
{
    Params() : character_matrix_file(), 
               meta_file(),
               out_path(), 
               out_prefix(),
               variant_reads_file(),
               total_reads_file(),
               region_reads_file(),
               cell_samples_file(),
               fp(0.02), 
               fn(0.02), 
               fn_indel(0.06),
               hom_precision(50.0),
               het_precision(8.0),
               hom_precision_indel(15.0),
               het_precision_indel(4.0),
               dropout_concentration(100.0),
               dropout_rate_prior(0.05),
               theta(6.0),
               iters(200),
               delta(20),
               CNA_penalty(85.0),
               tau1(0.5),
               tau2(0.1),
               seed(time(NULL)) {}

    string character_matrix_file; // character matrix file where rows are cells and columns are mutations
    string meta_file; // meta data file
    string out_path;
    string out_prefix;
    string variant_reads_file;
    string total_reads_file;
    string region_reads_file;
    string cell_samples_file; // line separated file indicating which sample each cell comes from
    double fp; // false positive rate
    double fn; // false negative rate (SNVs)
    double fn_indel; // false negative rate (indels)
    double hom_precision; // the precision for the beta binomial model given a homozygous locus (SNVs)
    double het_precision; // the precision for the beta binomial model given a heterozygous locus (SNVs)
    double hom_precision_indel; // the precision for the beta binomial model given a homozygous locus (indels)
    double het_precision_indel; // the precision for the beta binomial model given a homozygous locus (indels)
    double dropout_concentration; // concentration parameter used when computing the allelic dropout rate
    double dropout_rate_prior; // prior value for dropout rate 
    double theta; // theta parameter for negative binomial
    int iters; // number of iterations to run hill climbing
    int delta; // maximum number of hill climbing iterations to perform without likelihood improvement
    double CNA_penalty; // penalty for including disjoint CNAs 
    double tau1; // probability of performing an SNV relocation during hill climbing
    double tau2; // probability adding a CNA clone during hill climbing
    int seed;
};

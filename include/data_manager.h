#pragma once 

#include <vector>
#include <unordered_map>
#include <cstdint>
#include <span>
#include <tuple>

using std::vector;
using std::unordered_map;
using std::span;
using std::tuple;

/* Manages data and parameters */
class data_manager {      
    public:         
        data_manager(const double & fp, 
                     const double & fn, 
                     const double & fn_indel,
                     const double & hom_precision, 
                     const double & het_precision,
                     const double & hom_precision_indel, 
                     const double & het_precision_indel,
                     const double & dropout_concentration,
                     const double & dropout_rate_prior,
                     const double & theta,
                     const int & iters,
                     const int & delta,
                     const double & CNA_penalty,
                     const double & tau1,
                     const double & tau2); 
        double score(const int & cell, const int & mutation, const int observed_value, const int value) const;
        void fill(const vector<vector<int>> & character_matrix,
                  const vector<vector<int>> & variant_reads, 
                  const vector<vector<int>> & total_reads,
                  const vector<vector<int>> & region_reads,
                  const vector<int> & cell_samples,
                  const vector<int> & locus_regions,
                  const vector<bool> & is_germline,
                  const vector<bool> & is_sbs,
                  const vector<int> & locus_samples,
                  const vector<int> & chromosomes);
        template<typename T>
        const span<const vector<T>> get_sample_span(const int & sample, const vector<vector<T>> & matrix) const
        {
            const tuple<int, int> indices = this->sample_indices.at(sample);
            const int sample_start_index = get<0>(indices);
            const int sample_end_index = get<1>(indices);
            const int num_elements = sample_end_index - sample_start_index;
            return span<const vector<T>>(matrix.data() + sample_start_index, num_elements);
        }
        const span<const int> get_sample_span(const int & sample, const vector<int> & data) const;
        void compute_region_weights(const int sample, const vector<int> & assignments, const vector<int> & parents);
        void print(void) const;
        span<const vector<int>> get_character_matrix(const int & sample) const {return this->get_sample_span(sample, this->character_matrix);}
        span<const vector<int>> get_variant_reads(const int & sample) const {return this->get_sample_span(sample, this->variant_reads);}
        span<const vector<int>> get_total_reads(const int & sample) const {return this->get_sample_span(sample, this->total_reads);}
        span<const vector<int>> get_region_reads(const int & sample) const {return this->get_sample_span(sample, this->region_reads);}
        const vector<int> & get_locus_regions(void) const {return this->locus_regions;}
        const vector<bool> & get_is_germline(void) const {return this->is_germline;}
        const vector<bool> & get_is_sbs(void) const {return this->is_sbs;}
        const vector<int> & get_locus_samples(void) const {return this->locus_samples;}
        const vector<int> & get_chromosomes(void) const {return this->chromosomes;}
        const unordered_map<int, int> & get_region_to_chr(void) const {return this->region_to_chr;}
        const vector<double> & get_region_weights(const int & sample) const {return this->region_weights[sample];}
        const vector<bool> & get_reliable_regions(const int & sample) const {return this->reliable_regions[sample];}
        span<const int>  get_cell_read_sums(const int & sample) const {return this->get_sample_span(sample, this->cell_read_sums);}
        const tuple<int, int> get_sample_indices(const int & sample) const {return this->sample_indices.at(sample);}
        const vector<int> & get_variant_order(const int & sample) const {return this->variant_orders.at(sample);}

        const vector<int> & get_samples(void) const {return this->samples;}
        int get_num_samples(void) const {return this->samples.size();}
        double get_fp(void) const {return this->fp;}
        double get_fn(void) const {return this->fn;}
        double get_fn_indel(void) const {return this->fn_indel;}
        double get_hom_precision(void) const {return this->hom_precision;}
        double get_het_precision(void) const {return this->het_precision;}
        double get_hom_precision_indel(void) const {return this->hom_precision_indel;}
        double get_het_precision_indel(void) const {return this->het_precision_indel;}
        double get_dropout_concentration(void) const {return this->dropout_concentration;}
        double get_dropout_rate_prior(void) const {return this->dropout_rate_prior;}
        double get_theta(void) const {return this->theta;}
        int get_num_cells(const int sample) const;
        int get_num_loci(void) const {return this->num_loci;}
        int get_num_regions(void) const {return this->num_regions;}
        int get_iters(void) const {return this->iters;}
        int get_delta(void) const {return this->delta;}
        double get_CNA_penalty(void) const {return this->CNA_penalty;}
        double get_tau1(void) const {return this->tau1;}
        double get_tau2(void) const {return this->tau2;}

        // other functions
        void define_variant_orders(void);

    private:

        void find_sample_indices(void);

        double fp;
        double fn;
        double fn_indel;
        int num_cells;
        int num_loci;
        int num_regions;
        double hom_precision;
        double het_precision;
        double hom_precision_indel;
        double het_precision_indel;
        double dropout_concentration;
        double  dropout_rate_prior;
        double theta;
        vector<vector<int>> character_matrix;
        vector<vector<int>> variant_reads;
        vector<vector<int>> total_reads;
        vector<vector<int>> region_reads;
        vector<vector<double>> region_weights;
        vector<int> locus_regions;
        vector<int> locus_samples;
        vector<bool> is_germline;
        vector<bool> is_sbs;
        unordered_map<int, int> region_to_chr;
        vector<vector<bool>> reliable_regions;
        vector<int> chromosomes;
        vector<int> cell_read_sums;
        vector<int> cell_samples;
        vector<int> samples;
        unordered_map<int, int> sample_to_index;
        unordered_map<int, tuple<int, int>> sample_indices;
        unordered_map<int, vector<int>> variant_orders;
        int iters;
        int delta;
        double CNA_penalty;
        double tau1;
        double tau2;
};
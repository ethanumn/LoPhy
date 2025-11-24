
#pragma once 

#include "utils.h"
#include "data_manager.h"

/* Tree class */
class Tree
{
    public:
        Tree() : phi_hat(-1e-30), llh(NINF) {}
        Tree(int num_loci,
             int total_num_cells,
             double dropout_rate_prior) : phi_hat(-1e-30), 
                                          llh(NINF),
                                          num_loci(num_loci),
                                          parents(num_loci, NO_PARENT),
                                          clones({ROOT}),
                                          cnas(num_loci+1, vector<CNA>()),
                                          variants(num_loci+1, vector<int>()),
                                          loci(),
                                          priors(num_loci+1, 0.0),
                                          dropout_rates(num_loci, dropout_rate_prior),
                                          ref_dropout_rates(num_loci, dropout_rate_prior),
                                          alt_dropout_rates(num_loci, dropout_rate_prior),
                                          cell_assignments(total_num_cells, ROOT),
                                          dropout_rates_collection(),
                                          ref_dropout_rates_collection(),
                                          alt_dropout_rates_collection(),
                                          priors_collection(),
                                          dummy_clone_count(0) {}

        // update parents
        void update(const data_manager & manager, const int & sample);

        // add/check for CNAs
        void _add_cna(const int & clone, const CNA & cna) {this->cnas[clone].push_back(cna); sort(this->cnas[clone].begin(), this->cnas[clone].end());}
        bool has_CNA_in_region(const int & clone, const int & region) const;
        void remove_CNA_from_region(const int & clone, const int & region);
        
        // remove variant
        void remove_variant(const int & variant, const vector<int> & variants_in_sample, const data_manager & manager);

        // getters
        int get_parent(int clone_index) const {if(clone_index < static_cast<int>(this->parents.size()) && clone_index >= 0) return this->parents[clone_index]; else return NO_PARENT;}
        const vector<int> & get_parents(void) const {return this->parents;}
        vector<int> get_parents(void) {return this->parents;}
        const vector<int> & get_clones(void) const {return this->clones;}
        const vector<CNA> & get_clone_cnas(int clone) const {return this->cnas[clone];}
        const vector<vector<int>> & get_variants(void) const {return this->variants;}
        const vector<int> & get_loci(void) const {return this->loci;}
        const vector<double> & get_dropout_rates(void) const {return this->dropout_rates;}
        const vector<double> & get_ref_dropout_rates(void) const {return this->ref_dropout_rates;}
        const vector<double> & get_alt_dropout_rates(void) const {return this->alt_dropout_rates;}
        const double & get_dropout_rate(int locus) const {return this->dropout_rates[locus];}
        const double & get_ref_dropout_rate(int locus) const {return this->ref_dropout_rates[locus];}
        const double & get_alt_dropout_rate(int locus) const {return this->alt_dropout_rates[locus];}
        const vector<double> & get_priors(void) {return this->priors;}
        double get_llh(void) const {return this->llh;}
        double get_phi_hat(void) {return this->phi_hat;}
        const vector<int> & get_cell_assignments(void) const {return this->cell_assignments;}
        const vector<vector<double>> & get_priors_collection(void) const {return this->priors_collection;}
        const vector<vector<double>> & get_dropout_rates_collection(void) const {return this->dropout_rates_collection;}
        const vector<vector<double>> & get_ref_dropout_rates_collection(void) const {return this->ref_dropout_rates_collection;}
        const vector<vector<double>> & get_alt_dropout_rates_collection(void) const {return this->alt_dropout_rates_collection;}
        const unordered_set<int> & get_regions(void) const {return this->regions;}
        const unordered_set<int> & get_candidate_regions(const int & clone) const {return this->clone_candidate_regions.at(clone);}
        const unordered_set<int> & get_candidate_regions(void) const {return this->all_candidate_regions;}
        const CNAs & get_cnas(void) const {return this->cnas;}

        // setters
        void set_parent(int clone_index, int clone_parent);
        void set_priors(vector<double> priors) {this->priors = priors;}
        void set_cnas(CNAs cnas) {this->cnas = cnas;}
        void set_phi_hat(double phi_hat) {this->phi_hat = phi_hat;}
        void set_llh(double llh) {this->llh = llh;}
        void set_cell_assignment(const int index, const int clone) {if(index < static_cast<int>(this->cell_assignments.size())) this->cell_assignments[index] = clone;}

        // functions for estimating parameters
        void reset_priors(void);
        void reset_dropout_rates(const double dropout_rate_prior);
        void save_dropout_rates(void) {this->dropout_rates_collection.push_back(this->dropout_rates); 
                                        this->ref_dropout_rates_collection.push_back(this->ref_dropout_rates);
                                        this->alt_dropout_rates_collection.push_back(this->alt_dropout_rates);}
        void find_candidate_regions(const data_manager & manager, 
                                    const vector<int> & clones_in_sample, 
                                    const int & sample);
        void save_priors(void) {this->priors_collection.push_back(this->priors);}
        void EM(const vector<vector<int>> & A,
                const vector<vector<int>> & B,
                const vector<vector<double>> & log_likelihoods,
                const data_manager & manager,
                const int & sample,
                SCORE_CACHE & cache);
        void add_cna(const int & clone,
                     const int & region,
                     const CNA_TYPE & cna_type,
                     const vector<int> & clones_in_sample);
        void score(const data_manager & manager,     
                    const int & sample,
                    SCORE_CACHE & cache);
        void CNAs_to_matrices(const data_manager & manager,
                              vector<vector<int>> & A,
                              vector<vector<int>> & B) const;
        bool has_CNA_overfitting(void);
        bool has_min_cells_attached(const data_manager & manager, const int sample, int min_cells = 1);    

        // dummy clone functions
        int add_dummy_clone(const int & sample);
        bool remove_empty_dummy_clones(const int & sample);
        const vector<int> get_dummy_clones_in_sample(const int & sample);
        int get_num_empty_clones(void);
        int get_num_dummy_clones_wo_cnas(void);
        int get_num_cnas_in_dummy_clones(void);

        // operator overloads for comparing Tree objects
        bool operator<(const Tree & rhs) const
        {
            if(abs(llh - rhs.llh) < 1e-300)
            {
                return false;
            }
            else 
            {
                return llh < rhs.llh;
            }
        }

        bool operator> (const Tree & rhs) const
        {
            if(abs(llh - rhs.llh) < 1e-300)
            {
                return false;
            }
            else 
            {
                return llh > rhs.llh;
            }
        }

        bool operator==(const Tree & rhs) const
        {
            if(abs(llh - rhs.llh) < 1e-300)
            {
                return true;
            }
            else 
            {
                return false;
            }
        }
    private:
        void _score(const vector<vector<int>> & A,
                    const vector<vector<int>> & B,
                    vector<vector<double>> & log_likelihoods,
                    const data_manager & manager,     
                    const int & sample,
                    SCORE_CACHE & cache);
        double phi_hat; // sum of log likelihood for each tree of size j = 0, 1, ...
        double llh; // log likelihood of current tree
        int num_loci; // number of loci
        vector<int> parents; // parents vector, parents[i] = j means that i is a child of j
        unordered_map<int, vector<int>> children;
        vector<int> clones; // vector of clones in tree in breadth-first order
        unordered_set<int> regions; // regions that can be impacted by CNAs
        CNAs cnas; // copy numbers impacting each clone
        vector<vector<int>> variants; // loci that are mutated in each clone
        vector<int> loci; // loci that are mutated in this tree
        vector<double> priors; // prior probabilities for each clone
        vector<double> dropout_rates; // dropout rates for current setting of cnas and variants assuming a fixed tree
        vector<double> ref_dropout_rates; // ref dropout rates for current setting of cnas and variants assuming a fixed tree
        vector<double> alt_dropout_rates; // alt dropout rates for current setting of cnas and variants assuming a fixed tree
        vector<int> cell_assignments; // clone assignments for each cell
        unordered_map<int,unordered_set<int>> clone_candidate_regions; // possible regions to introduce CNAs per clone
        unordered_set<int> all_candidate_regions; // possible regions across all clones where CNAs may be possible

        // functions that allow us to collect priors and dropout rates for each longitudinal sample
        vector<vector<double>> dropout_rates_collection; // collection of dropout rates for each sample
        vector<vector<double>> ref_dropout_rates_collection; // collection of ref dropout rates for each sample
        vector<vector<double>> alt_dropout_rates_collection; // collection of alt dropout rates for each sample
        vector<vector<double>> priors_collection; // collection of priors for each sample

        // for adding/removing dummy clones
        int dummy_clone_count;
        unordered_map<int, vector<int>> sample_dummy_clones;

};
#include "tree.h"

#include <numeric>
#include <tuple>
#include <map>
#include <cmath> 
#include <algorithm>

using std::iota;
using std::tuple;
using std::min;
using std::pair;
using std::make_pair;
using std::get;
using std::map;
using std::find;


void Tree::set_parent(int clone_index, int clone_parent) 
{
    if(clone_index < static_cast<int>(this->parents.size())) this->parents[clone_index] = clone_parent; 
    else 
        throw std::runtime_error("Invalid clone index passed");
}

bool Tree::has_CNA_in_region(const int & clone, const int & region) const
{
    const auto & clone_cnas = this->cnas.at(clone);
    return std::any_of(clone_cnas.begin(), clone_cnas.end(), [region](const CNA & cna) {
        return cna.first == region;
    });
}
void Tree::remove_CNA_from_region(const int & clone, const int & region) 
{
    this->cnas[clone].erase(remove_if(this->cnas[clone].begin(), this->cnas[clone].end(), [region](CNA cna) {return cna.first == region;}), this->cnas[clone].end());
}

void Tree::reset_priors(void)
{
    fill(this->priors.begin(), this->priors.end(), 0.0);
    const auto c = this->clones.size();
    for(const auto & clone : this->clones)
        this->priors[clone] = 1.0/(double)c;
}

void Tree::reset_dropout_rates(const double dropout_rate_prior)
{
    fill(this->dropout_rates.begin(), this->dropout_rates.end(), dropout_rate_prior);
    fill(this->ref_dropout_rates.begin(), this->ref_dropout_rates.end(), dropout_rate_prior);
    fill(this->alt_dropout_rates.begin(), this->alt_dropout_rates.end(), dropout_rate_prior);
}

void Tree::find_candidate_regions(const data_manager & manager, 
                                  const vector<int> & clones_in_sample, 
                                  const int & sample)
{
    // obtain necessary data
    tuple<int, int> sample_indices = manager.get_sample_indices(sample);
    const auto sample_start_index = get<0>(sample_indices);
    const auto num_clones = this->parents.size() + 1;
    int num_cells = manager.get_num_cells(sample);
    const auto & num_regions = manager.get_num_regions();
    const auto & region_reads = manager.get_region_reads(sample);
    const auto & region_weights = manager.get_region_weights(sample);
    const auto & cell_read_sums = manager.get_cell_read_sums(sample);

    // clear candidate regions
    this->clone_candidate_regions.clear();
    this->all_candidate_regions.clear();

    // initialize new candidate region keys
    for(const auto & clone : this->clones)
        this->clone_candidate_regions[clone] = unordered_set<int>();

    // record average fraction that fall in each region for each clone based on hard cell assignments
    vector<unordered_map<int,vector<double>>> node_region_weights(num_clones, unordered_map<int,vector<double>>());
    for(const auto & r : this->regions)
    {
        for(auto cell = 0; cell < num_cells; ++cell)
        {
            const auto & clone = this->cell_assignments[sample_start_index + cell];
            node_region_weights[clone][r].push_back(
                region_reads[cell][r] / cell_read_sums[cell]
            );
        }
    }
    // Select regions whose probability is different between the root and another node
    for(const auto & r : this->regions)
    {
        for(const auto & clone : clones_in_sample)
        {
            if(clone != ROOT)
            {
                if (node_region_weights[clone][r].size() >= std::max(40.0, 0.03*num_cells))
                {
                    double weight = 0;
                    for (double w : node_region_weights[clone][r]) 
                        weight += w / node_region_weights[clone][r].size();
                    if ((weight > region_weights[r]*1.275 || region_weights[r] > weight*1.35) && (region_weights[r] > (0.05 / num_regions)))
                    { 
                        this->clone_candidate_regions[clone].insert(r);
                        this->all_candidate_regions.insert(r);
                        break;
                    }
                }
            }
        }
    }

}


bool Tree::has_min_cells_attached(const data_manager & manager, const int sample, int min_cells)
{
    unordered_map<int, int> counts;
    const int num_loci = manager.get_num_loci();
    tuple<int, int> sample_indices = manager.get_sample_indices(sample);
    const auto sample_start_index = get<0>(sample_indices);
    const vector<int> & variants_in_sample = manager.get_variant_order(sample);
    const auto & num_cells = manager.get_num_cells(sample);

    counts[ROOT] = 0;
    // make sure we initialize counts properly
    for(const auto & v : variants_in_sample)
    {
        if(this->get_parent(v) != NO_PARENT)
            if((int)this->get_parent(v) <= 2*num_loci)
                counts[v+1] = 0;
    }

    for(const auto & v : this->sample_dummy_clones[sample])
        counts[v+1] = 0;

    // count number of cells attached to each clone for this sample
    for (auto cell = 0; cell < num_cells; ++cell) 
    {
        const auto & val = this->cell_assignments[sample_start_index + cell];
        counts[val]++;
    }
    
    bool clones_have_min_cells = true;
    for (const auto & pair : counts) 
    {
        if((pair.second < min_cells) && (pair.first != ROOT))
        {
            clones_have_min_cells = false;
            break;
        }
    }

    return clones_have_min_cells;
}

void Tree::update(const data_manager & manager, const int & sample)
{

    const vector<int> & variants_in_sample = manager.get_variant_order(sample);
    const auto & locus_regions = manager.get_locus_regions();
    const auto & is_germline = manager.get_is_germline();
    const auto & region_is_reliable = manager.get_reliable_regions(sample);
    const auto & num_regions = manager.get_num_regions();
    const int num_clones = this->parents.size() + 1;
    const int num_loci = static_cast<int>(manager.get_num_loci());
    this->variants.clear();
    this->variants.resize(num_clones, vector<int>());
    this->loci.clear();
    this->regions.clear();
    this->children.clear();

    // Collect variants in each clone
    for (int j = 0; j < static_cast<int>(this->parents.size()); ++j)
    {
        const int p = this->parents[j];

        const bool is_locus = (j < num_loci);

        if (p == NO_PARENT)
        {
            if (is_locus && is_germline[j])
            {
                // Variant originates at root
                this->loci.push_back(j);
                this->variants[ROOT].push_back(j);
            }
            continue;
        }

        if (is_locus)
        {
            this->loci.push_back(j);
            if (region_is_reliable[locus_regions[j]])
                this->regions.insert(locus_regions[j]);
        }

        if (p <= 2 * num_loci)
        {
            const int clone = j + 1;
            this->children[p].push_back(clone);
            if (is_locus)
                this->variants[clone].push_back(j);
        }
        else
        {
            const int clone = p - 2 * num_loci;
            this->variants[clone].push_back(j);
        }
    }

    // if all variants from this sample are incorporated into the tree utilize all regions
    if(vector_subset(variants_in_sample, this->loci))
    {
        for(auto r = 0; r < num_regions; ++r)
            if(region_is_reliable[r])
                this->regions.insert(r);
    }

    // collect clones in tree in breadth-first order
    this->clones = topological_sort(this->get_parents(), num_loci);

    // reset priors
    this->reset_priors();
}


void Tree::EM(const vector<vector<int>> & A,
              const vector<vector<int>> & B,
              const vector<vector<double>> & log_likelihoods,
              const data_manager & manager,
              const int & sample,
              SCORE_CACHE & cache)
{
    // initialize variables
    const auto & num_loci = manager.get_num_loci();
    int num_cells = manager.get_num_cells(sample);
    int num_clones = this->parents.size() + 1;

    // grab constant references to variables from manager
    const auto & dropout_rate_prior = manager.get_dropout_rate_prior();
    const auto & dropout_concentration = manager.get_dropout_concentration();
    const auto & variant_reads = manager.get_variant_reads(sample);
    const auto & total_reads = manager.get_total_reads(sample);
    const auto & locus_regions = manager.get_locus_regions();
    const auto & is_sbs = manager.get_is_sbs();

    // create matrix of posterior probabilities
    vector<vector<double>> cell_attachment_probs(num_cells, (vector<double>(num_clones, 0.0)));
    vector<double> clone_attachment_probs(num_clones, 0.0);
    vector<double> cell_llhs(num_cells, 0.0);
    for(auto cell = 0; cell < num_cells; ++cell)
        cell_llhs[cell] = logsumexp(log_likelihoods[cell]);

    for(auto cell = 0; cell < num_cells; ++cell)
    {
        for(const auto & clone : this->clones)
        {
            cell_attachment_probs[cell][clone] = exp(log_likelihoods[cell][clone] - cell_llhs[cell]);
            clone_attachment_probs[clone] += cell_attachment_probs[cell][clone];
        }
    }

    ////////////////
    //// E-step ////
    ////////////////
    int Ca, Cr, Ct;
    vector<double> E_dropped_ref_alleles, E_dropped_alt_alleles;
    vector<double> dropped_ref_alleles(num_loci, 0.0);
    vector<double> dropped_alt_alleles(num_loci, 0.0);
    vector<double> total_ref_alleles(num_loci, 0.0);
    vector<double> total_alt_alleles(num_loci, 0.0);

    // collect probability of each genotype for each cell
    for(const auto & locus : this->loci)
    {
        const auto & alt_dropout_rate = this->get_alt_dropout_rate(locus);
        const auto & ref_dropout_rate = this->get_ref_dropout_rate(locus);

        // choose precision and fn value based on whether or not locus has single base substitution or indel
        auto hom_precision = manager.get_hom_precision();
        auto het_precision = manager.get_het_precision();
        auto fn = manager.get_fn();
        if(!is_sbs[locus])
        {
            hom_precision = manager.get_hom_precision_indel();
            het_precision = manager.get_het_precision_indel();
            fn = manager.get_fn_indel();
        }

        map<pair<int, int>, vector<double>> genotype_probs;
        for (const auto & clone : this->clones)
        {
            // copy numbers for clone
            Ca = B[clone][locus];
            Ct = A[clone][locus_regions[locus]];  
            Cr = Ct - Ca;
            pair<int, int> genotype = make_pair(Ca, Cr); 

            // initialize genotype in map
            if (!genotype_probs.count(genotype))
                genotype_probs[genotype] = vector<double>(num_cells, 0.0);

            // collect the posteriors for each cell for this particular genotype
            for (auto cell = 0; cell < num_cells; ++cell)
                genotype_probs[genotype][cell] += cell_attachment_probs[cell][clone];
        }

        // compute expected number of lost alleles for each cell
        for(const auto & gp : genotype_probs)
        {
            const auto & genotype = gp.first;
            const auto & probs = gp.second;
            Ca = genotype.first;
            Cr = genotype.second;

            if(Cr > 0 && Ca > 0)
            {
                E_dropped_ref_alleles = cache.get_dropped_ref_alleles(this->get_parents(),
                                                                      locus, 
                                                                      Ca, 
                                                                      Cr, 
                                                                      alt_dropout_rate,
                                                                      ref_dropout_rate,
                                                                      variant_reads, 
                                                                      total_reads, 
                                                                      manager.get_fp(), 
                                                                      fn,
                                                                      hom_precision,
                                                                      het_precision);
                E_dropped_alt_alleles = cache.get_dropped_alt_alleles(this->get_parents(),
                                                                      locus, 
                                                                      Ca, 
                                                                      Cr, 
                                                                      alt_dropout_rate,
                                                                      ref_dropout_rate,
                                                                      variant_reads, 
                                                                      total_reads, 
                                                                      manager.get_fp(), 
                                                                      fn,
                                                                      hom_precision,
                                                                      het_precision);
                for(auto cell = 0; cell < num_cells; ++cell)
                {
                    if(total_reads[cell][locus] > 4)
                    {
                        dropped_alt_alleles[locus] += probs[cell] * E_dropped_alt_alleles[cell];
                        dropped_ref_alleles[locus] += probs[cell] * E_dropped_ref_alleles[cell];
                        total_alt_alleles[locus] += probs[cell] * (double)Ca;
                        total_ref_alleles[locus] += probs[cell] * (double)Cr;
                    }
                }

            }
        }
    }

    ////////////////
    //// M-step ////
    ////////////////

    // maximum likelihood estimate for prior probabilities
    vector<double> priors(num_clones, 0.0);
    for (auto clone = 0; clone < num_clones; ++clone) 
        priors[clone] += clone_attachment_probs[clone] / (double)num_cells;
    this->set_priors(priors);

    // maximum likelihood estimate for dropout rates
    for(const auto & locus : this->loci)
    {
        const auto dropout_rate = ((dropout_rate_prior * dropout_concentration - 1)*2 + dropped_ref_alleles[locus] + dropped_alt_alleles[locus]) 
                                    / ((dropout_concentration-2)*2 + total_ref_alleles[locus] + total_alt_alleles[locus]);
        const auto alt_dropout_rate = (dropout_rate_prior * dropout_concentration - 1 + dropped_alt_alleles[locus]) 
                                    / (dropout_concentration - 2 + total_alt_alleles[locus]);
        const auto ref_dropout_rate = (dropout_rate_prior * dropout_concentration - 1 + dropped_ref_alleles[locus]) 
                                    / (dropout_concentration - 2 + total_ref_alleles[locus]);
        const auto diff_dropout_rates_llh = dropped_ref_alleles[locus] * log(ref_dropout_rate) + (total_ref_alleles[locus] - dropped_ref_alleles[locus])*log(1-ref_dropout_rate)
                                            + dropped_alt_alleles[locus] * log(alt_dropout_rate) + (total_alt_alleles[locus] - dropped_alt_alleles[locus])*log(1-alt_dropout_rate)
                                            + (dropout_rate_prior * dropout_concentration - 1) * (log(ref_dropout_rate) + log(alt_dropout_rate))
                                            + ((1-dropout_rate_prior) * dropout_concentration - 1) * (log(1-ref_dropout_rate) + log(1-alt_dropout_rate));
        const auto same_dropout_rates_llh = (dropped_ref_alleles[locus] + dropped_alt_alleles[locus]) * log(dropout_rate)
                                            + (total_alt_alleles[locus] + total_ref_alleles[locus] - dropped_alt_alleles[locus] - dropped_ref_alleles[locus]) * log(1-dropout_rate)
                                            + (dropout_rate_prior * dropout_concentration - 1) * 2 * log(dropout_rate)
                                            + ((1 - dropout_rate_prior) * dropout_concentration - 1) * 2 * log(1 - dropout_rate);
        this->dropout_rates[locus] = fmin(0.5, fmax(dropout_rate, 0.01));
        if(same_dropout_rates_llh > diff_dropout_rates_llh - 60)
        {
            this->ref_dropout_rates[locus] = fmin(0.5, fmax(dropout_rate, 0.01));
            this->alt_dropout_rates[locus] = fmin(0.5, fmax(dropout_rate, 0.01));
        }
        else 
        {
            this->ref_dropout_rates[locus] = fmin(0.5, fmax(ref_dropout_rate, 0.01));
            this->alt_dropout_rates[locus] = fmin(0.5, fmax(alt_dropout_rate, 0.01));
        }

    }
}

void Tree::add_cna(const int & clone,
                   const int & region,
                   const CNA_TYPE & cna_type,
                   const vector<int> & clones_in_sample)
{
    // remove CNAs from region for all clones in this sample that are ancestral or descendant
    unordered_set<int> ancestors;
    unordered_set<int> descendants;
    get_ancestors(ancestors, this->get_parents(), clone);
    get_descendants(descendants, this->get_parents(), clone);

    for(const auto & c : clones_in_sample)
        if(find(ancestors.begin(), ancestors.end(), c) != ancestors.end() || find(descendants.begin(), descendants.end(), c) != descendants.end() || c == clone)
            this->remove_CNA_from_region(c, region);

    // add cna to clone
    if(cna_type != NA)
        this->_add_cna(clone, make_pair(region, cna_type));    
}

void Tree::score(const data_manager & manager,
                 const int & sample,
                 SCORE_CACHE & cache)
{
    // variables
    double avg_diff_dropout_rates = 10.0, avg_diff_node_probs = 10.0;

    // basic information about data
    const auto & num_regions = manager.get_num_regions();
    const auto & num_loci = manager.get_num_loci();
    const auto & num_cells = manager.get_num_cells(sample);
    const int num_clones = static_cast<int>(this->parents.size()) + 1;

    // matrices for computing likelihoods
    vector<vector<int>> A(num_clones, vector<int>(num_regions, 0));
    vector<vector<int>> B(num_clones, vector<int>(num_loci, 0));   
    vector<vector<double>> log_likelihoods(num_cells, vector<double>(num_clones, NINF));

    // fill copy number matrices based on tree
    this->CNAs_to_matrices(manager, A, B);

    // run EM with each CNA to determine which tree is optimal
    int epoch = 0;
    while((avg_diff_dropout_rates > 0.0001 || avg_diff_node_probs > 0.0005) && (epoch < 100))
    {

        const auto prev_ref_dropout_rates = this->ref_dropout_rates;
        const auto prev_alt_dropout_rates = this->alt_dropout_rates;
        const auto prev_priors = this->priors;
        this->_score(A, B, log_likelihoods, manager, sample, cache); // score tree
        EM(A, B, log_likelihoods, manager, sample, cache);

        avg_diff_dropout_rates = 0.0;
        for(int i = 0; i < num_loci; ++i)
            avg_diff_dropout_rates += (abs(this->ref_dropout_rates[i] - prev_ref_dropout_rates[i]) + abs(this->alt_dropout_rates[i] - prev_alt_dropout_rates[i]))
                                      / (double)this->loci.size();

        avg_diff_node_probs = 0.0;
        for(int i = 0; i < num_clones; ++i)
            avg_diff_node_probs += abs(this->priors[i] - prev_priors[i]) / (double)this->clones.size();
        epoch += 1;
    }
}

/* Assign cells to clones using the optimal total copy number (A) and mutant copy number (B) matrices */
void  Tree::_score(const vector<vector<int>> & A,
                   const vector<vector<int>> & B,
                   vector<vector<double>> & log_likelihoods,
                   const data_manager & manager,     
                   const int & sample,
                   SCORE_CACHE & cache) 
{
    // initialize variables
    double normalizing_constant;
    vector<double> mu;
    vector<double> term_probs;
 
    // get relevant data from data_manager
    const auto & num_cells = manager.get_num_cells(sample);
    const auto & num_regions = manager.get_num_regions();
    const auto & cell_read_sums = manager.get_cell_read_sums(sample);
    const auto & variant_reads = manager.get_variant_reads(sample);
    const auto & total_reads = manager.get_total_reads(sample);
    const auto & region_reads = manager.get_region_reads(sample);
    const auto & region_weights = manager.get_region_weights(sample);
    const auto & region_is_reliable = manager.get_reliable_regions(sample);
    const auto & locus_regions = manager.get_locus_regions();
    const auto & chromosomes = manager.get_chromosomes();
    const auto & dropout_concentration = manager.get_dropout_concentration();
    const auto & dropout_rate_prior = manager.get_dropout_rate_prior();
    const auto & region_chromosomes = manager.get_chromosomes();
    const auto & theta = manager.get_theta();
    const auto & is_sbs = manager.get_is_sbs();
    tuple<int, int> sample_indices = manager.get_sample_indices(sample);
    const auto sample_start_index = get<0>(sample_indices);
    vector<double> region_proportions(num_regions, 0.0);
    const double coeff = 0.2 + (double)num_cells / 8000.0; // for penalties
    const double CNA_penalty = manager.get_CNA_penalty();


    // reset log likelihoods
    for(auto & llhs : log_likelihoods)
        fill(llhs.begin(), llhs.end(), NINF);

    // Calculate likelihood of each clone assignment for each cell
    for (const auto & clone : this->clones) 
    {          
        // log prior 
        for(auto cell = 0; cell < num_cells; ++cell)
            log_likelihoods[cell][clone] = log(this->priors[clone]);

        // calculate region proportion for clone
        normalizing_constant = 0.0;
        fill(region_proportions.begin(), region_proportions.end(), 0.0);
        for(auto r = 0; r < num_regions; ++r)
        {
            if(region_is_reliable[r])
            {
                region_proportions[r] = (double)A[clone][r] * region_weights[r];
                normalizing_constant += region_proportions[r];
            }
        }

        // region read count likelihood
        for(const auto & r : this->regions)
        {
            const auto & llhs = cn_llh(cache.cn_llh,
                                       r, 
                                       region_proportions[r] / normalizing_constant,
                                       region_reads, 
                                       cell_read_sums,
                                       theta);

            for(auto cell = 0; cell < num_cells; ++cell)
                log_likelihoods[cell][clone] += llhs[cell];
        }

        // allelic read count likelihood
        for(const auto & locus : this->loci)
        {
            // choose precision and fn value based on whether or not locus has single base substitution or indel
            auto hom_precision = manager.get_hom_precision();
            auto het_precision = manager.get_het_precision();
            auto fn = manager.get_fn();
            if(!is_sbs[locus])
            {
                hom_precision = manager.get_hom_precision_indel();
                het_precision = manager.get_het_precision_indel();
                fn = manager.get_fn_indel();
            }

            const auto & ref_dropout_rate = this->get_ref_dropout_rate(locus);
            const auto & alt_dropout_rate = this->get_alt_dropout_rate(locus);
            const vector<double> & llhs = calculate_locus_llh(cache.locus_llh,
                                                              cache.E_dropped_alt_alleles, 
                                                              cache.E_dropped_ref_alleles, 
                                                              this->get_parents(),
                                                              locus, 
                                                              B[clone][locus], 
                                                              A[clone][locus_regions[locus]] - B[clone][locus], 
                                                              alt_dropout_rate,
                                                              ref_dropout_rate,
                                                              variant_reads, 
                                                              total_reads, 
                                                              manager.get_fp(), 
                                                              fn,
                                                              hom_precision,
                                                              het_precision);

            for(auto cell = 0; cell < num_cells; ++cell)
                log_likelihoods[cell][clone] += llhs[cell];
        }
    }

    // Find the maximum likelihood clone assignment, and calculate the complete likelihood
    double llh = 0.0;
    for (auto cell = 0; cell < num_cells; ++cell) 
    {
        const auto & llhs = log_likelihoods[cell];
        auto clone = distance(llhs.begin(), max_element(llhs.begin(), llhs.end()));
        this->cell_assignments[sample_start_index + cell] = clone;
        llh += logsumexp(llhs);
    }

    // PRIORS
    // collect all mutations in each region
    unordered_map<int, vector<int>> loci_in_region;
    for(const auto & locus : this->loci)
        loci_in_region[locus_regions[locus]].push_back(locus);

    // lambda function to determine whether two regions are adjacent if you remove all unreliable 
    // regions between them in the list of regions
    auto reliable_adjacent = [&](int a, int b) -> bool {
        if (!region_is_reliable[a] || !region_is_reliable[b]) return false;
        if (a > b) std::swap(a, b);
        for (int i = a + 1; i < b; ++i) {
            if (region_is_reliable[i]) return false;
        }
        return true;
    };

    // go through each clone's cnas and record which chromosome they are on
    for(const auto & clone : this->clones)
    {
        // get the number of descendant clones to penalize all CNAs
        unordered_set<int> descendants;
        get_descendants(descendants, this->get_parents(), clone);

        // initalize CNA_events, and collect all CNAs on each chromosome
        unordered_map<int, vector<CNA>> CNA_events;
        for (const auto & chrom : chromosomes) 
            CNA_events[chrom] = vector<CNA>(); // Initialize with an empty vector

        for(const auto & cna : this->get_clone_cnas(clone))
        {
            const auto & region = cna.first;
            CNA_events[region_chromosomes.at(region)].push_back(cna);
        }

        // penalize disjoint CNAs
        int prev_region, region;
        CNA_TYPE type, prev_type;

        for(const auto & chrom : chromosomes)
        {
            // sort all events based on region
            sort(CNA_events[chrom].begin(), CNA_events[chrom].end(),
                [](const CNA & a, const CNA & b) {
                    return a.first < b.first;
            });

            prev_region = -10;
            region = -10;
            prev_type = NA;
            // penalize disjoint CNAs
            for(const CNA & e : CNA_events[chrom])
            {
                region = e.first;
                type = e.second;
                bool disjoint = true;
                if(prev_region >= 0)
                {
                    if(type == CNLOH_ALT || type == CNLOH_REF)
                        disjoint = prev_region+1 != region;
                    else 
                        disjoint = reliable_adjacent(prev_region, region);
                }

                // penalize CNA if it's disjoint
                if(disjoint || type != prev_type)
                {
                    if(type == CNLOH_ALT || type == CNLOH_REF)
                        llh -= 2*coeff*CNA_penalty;
                    else
                        llh -= coeff*CNA_penalty;
                }
                prev_region = region;
                prev_type = type;

            }
        }
    }

    // penalize dropout rate
    for(const auto & locus : this->loci)
    {
        const auto & ref_dropout_rate = this->get_ref_dropout_rate(locus);
        const auto & alt_dropout_rate = this->get_alt_dropout_rate(locus);
        if(abs(ref_dropout_rate - alt_dropout_rate) > 0.001)
            llh -= 70;
        llh += (dropout_rate_prior*dropout_concentration-1)*(log(ref_dropout_rate) + log(alt_dropout_rate)) 
                    + ((1.0-dropout_rate_prior) * dropout_concentration-1)*(log(1.0-ref_dropout_rate) + log(1.0-alt_dropout_rate));
    }

    // penalty for number of non-CNA defined clones
    llh -= (this->clones.size() - this->dummy_clone_count) * (4 + this->loci.size());

    // // penalty for CNA clones
    // llh -= 50 * coeff * this->dummy_clone_count * (8 + this->loci.size());

    // heavily penalize trees with nodes that do not have any cells attached
    // if(!this->has_min_cells_attached(manager, sample, 1))
    //     llh -= 1e8;
    
    // heavily penalize if the same CNAs occur in multiple children of a clone
    if(this->has_CNA_overfitting())
        llh -= 1e8;

    this->set_llh(llh);
}

void Tree::CNAs_to_matrices(const data_manager & manager,
                            vector<vector<int>> & A,
                            vector<vector<int>> & B) const
{
    // clear A and B
    for(auto & row : A)
        fill(row.begin(), row.end(), 0);
    for(auto & row : B)
        fill(row.begin(), row.end(), 0);

    // get relevant information about sequencing data
    const auto & locus_regions = manager.get_locus_regions();
    const auto & is_germline = manager.get_is_germline();

    // fill in copy number matrices from the root down
    for(const int & clone : this->clones)
    {
        if(clone == ROOT)
        {
            fill(A[ROOT].begin(), A[ROOT].end(), 2);

            for(const auto & locus : this->loci)
            {
                if(is_germline[locus])
                    B[ROOT][locus] = 1;
                else 
                    B[ROOT][locus] = 0;
            }
        }
        else
        {
            
            // initialize copy number for each clone with its parent's copy numbers
            const auto & clone_parent = this->parents[clone-1];
            A[clone] = A[clone_parent];
            B[clone] = B[clone_parent];

            // add new variants unique to this clone
            for(const auto & locus : this->variants[clone])
                B[clone][locus] = 1;

            // collect locus that are mutated in this clone
            vector<int> variants_in_clone;
            for(const auto & locus : this->loci)
            {
                if(B[clone][locus] == 1)
                    variants_in_clone.push_back(locus);
            }
            
            // add cnas unique to this clone
            for(const CNA & cna : this->get_clone_cnas(clone))
            {
                const auto & region = cna.first;
                const auto & type = cna.second;
                switch(type)
                {
                    case LOSS_ALT:
                        A[clone][region] = 1;
                        for(const auto & locus : variants_in_clone)
                            if(locus_regions[locus] == region)
                                B[clone][locus] = 0;
                        break;
                    case LOSS_REF:
                        A[clone][region] = 1;
                        break;
                    case GAIN_ALT:
                        A[clone][region] = 3;
                        for(const auto & locus :variants_in_clone)
                            if(locus_regions[locus] == region)
                                B[clone][locus] = 2;
                        break;
                    case GAIN_REF:
                        A[clone][region] = 3;
                        break;
                    case CNLOH_ALT:
                        A[clone][region] = 2;
                        for(const auto & locus : variants_in_clone)
                            if(locus_regions[locus] == region)
                                B[clone][locus] = 0;
                        break;
                    case CNLOH_REF:
                        A[clone][region] = 2;
                        for(const auto & locus : variants_in_clone)
                            if(locus_regions[locus] == region)
                                B[clone][locus] = 2;
                        break;
                    default:
                        break;
                }
            }
        }
    }
}  
bool Tree::has_CNA_overfitting(void)
{
    for (const auto& it : this->children)
    {
        const auto & children = it.second;

        if (children.size() < 2)
            continue;

        std::unordered_map<CNA, int> cna_count;
        for (int child : children)
        {
            for (const CNA& cna : this->cnas[child])
            {
                ++cna_count[cna];
            }
        }

        for (const auto& c : cna_count)
        {
            if (c.second >= 2)
                return true;
        }
    }

    return false;
}


int Tree::add_dummy_clone(const int & sample)
{
    this->parents.push_back(NO_PARENT);
    this->cnas.push_back(vector<CNA>());
    this->variants.push_back(vector<int>());
    this->priors.push_back(0.0);
    const auto dummy_clone_index = this->num_loci + this->dummy_clone_count;
    this->dummy_clone_count++;
    if(this->sample_dummy_clones.find(sample) != this->sample_dummy_clones.end())
        this->sample_dummy_clones[sample].push_back(dummy_clone_index);
    else    
        this->sample_dummy_clones[sample] = {dummy_clone_index};
    return dummy_clone_index;
}

void Tree::remove_variant(const int & variant, const vector<int> & variants_in_sample, const data_manager & manager)
{
    const int parent_node = this->get_parent(variant);
    const int num_loci = manager.get_num_loci();

    if(parent_node != NO_PARENT)

    {
        this->set_parent(variant, NO_PARENT); // remove variant
        if(parent_node <= 2*num_loci)
        {
            // collect all children of variant
            vector<int> children;
            for(int clone = 0; clone < static_cast<int>(this->parents.size()); ++clone)
                if(this->get_parent(clone) == variant+1)
                    children.push_back(clone);

            // collect all SNVs in clone defined by variant
            vector<int> SNVs_in_clone;
            for(const auto & v : variants_in_sample)
                if(this->get_parent(v) == (2*num_loci + variant))
                    SNVs_in_clone.push_back(v);

            if(SNVs_in_clone.size() > 0)
            {
                // make new node with one of the variants and set its parent to the original node's parent
                int new_node = SNVs_in_clone[0];
                this->set_parent(new_node,parent_node);

                // place all SNVs in this new node
                for(int i = 1; i < static_cast<int>(SNVs_in_clone.size()); ++i)
                    this->set_parent(SNVs_in_clone[i],2*num_loci + new_node);

                // set all children to be parented by this new node
                for(const auto & v : children)
                    this->set_parent(v,new_node+1);
            }
            else // set children to be parented by node's parent 
            { 
                for(const auto v : children)
                    this->set_parent(v,parent_node);
            }
            // remove all CNAs
            this->cnas[variant+1].clear();
        }
    }
}

bool Tree::remove_empty_dummy_clones(const int & sample)
{
    vector<int> indices_to_remove;

    // Identify dummy clones to remove
    for (int idx = this->parents.size() - 1; idx >= this->num_loci; --idx)
    {
        if ((this->parents[idx] == NO_PARENT) || this->cnas[idx + 1].empty())
            indices_to_remove.push_back(idx);
    }

    // Remove in reverse to avoid index shifting
    for (const auto & idx : indices_to_remove)
    {
        this->parents.erase(this->parents.begin() + idx);
        this->cnas.erase(this->cnas.begin() + idx);
        this->variants.erase(this->variants.begin() + idx);
        this->priors.erase(this->priors.begin() + idx);

        // Also remove from sample_dummy_clones
        auto & clones = this->sample_dummy_clones[sample];
        clones.erase(std::remove(clones.begin(), clones.end(), idx), clones.end());
    }

    // Adjust dummy count
    this->dummy_clone_count = this->parents.size() - this->num_loci;

    return !indices_to_remove.empty();
}

const vector<int> Tree::get_dummy_clones_in_sample(const int & sample) 
{
    if(this->sample_dummy_clones.find(sample) != this->sample_dummy_clones.end())
        return this->sample_dummy_clones[sample];
    else
        return vector<int>();
}

int Tree::get_num_empty_clones(void)
{
    int count = 0;
    for(const auto & clone : this->clones)
    {   
        if(clone == ROOT)
            continue;

        if(this->variants[clone].empty() && this->cnas[clone].empty())
            count++;
    }
    return count;
}

int Tree::get_num_cnas_in_dummy_clones(void)
{
    int count = 0;
    for(const auto & clone : this->clones)
    {
        if((this->variants[clone].size() == 0) && (this->cnas[clone].size() > 0))
            count += this->cnas[clone].size();
    }
    return count;
}

int Tree::get_num_dummy_clones_wo_cnas(void)
{
    int count = 0;
    for(const auto & clone : this->clones)
    {
        bool found_cnloh = false;
        if((this->variants[clone].size() == 0) && (this->cnas[clone].size() > 0))
        {
            for(const CNA & cna : this->get_clone_cnas(clone))
            {
                const auto & type = cna.second;
                if(type == CNLOH_REF || type == CNLOH_ALT)
                {
                    found_cnloh = true;
                    break;
                }
            }
            if(!found_cnloh)
                count++;
        }
    }
    return count;
}


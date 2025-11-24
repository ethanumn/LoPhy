#include "data_manager.h"
#include "utils.h"

#include <iostream>

using std::accumulate;

data_manager::data_manager(const double & fp, 
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
                           const double & tau2) : fp(fp),
                                                  fn(fn), 
                                                  fn_indel(fn_indel),
                                                  num_cells(0), 
                                                  num_loci(0), 
                                                  num_regions(0),
                                                  hom_precision(hom_precision),
                                                  het_precision(het_precision),
                                                  hom_precision_indel(hom_precision_indel),
                                                  het_precision_indel(het_precision_indel),
                                                  dropout_concentration(dropout_concentration),
                                                  dropout_rate_prior(dropout_rate_prior),
                                                  theta(theta),
                                                  variant_reads(),
                                                  total_reads(),
                                                  region_reads(),
                                                  region_weights(),
                                                  locus_regions(),
                                                  locus_samples(),
                                                  is_germline(),
                                                  is_sbs(),
                                                  region_to_chr(),
                                                  reliable_regions(),
                                                  chromosomes(),
                                                  cell_read_sums(),
                                                  cell_samples(),
                                                  samples(),
                                                  variant_orders(),
                                                  iters(iters),
                                                  delta(delta),
                                                  CNA_penalty(CNA_penalty),
                                                  tau1(tau1),
                                                  tau2(tau2) {}

const span<const int> data_manager::get_sample_span(const int & sample, const vector<int> & vec) const
{
    const tuple<int, int> indices = this->sample_indices.at(sample);
    const int sample_start_index = get<0>(indices);
    const int sample_end_index = get<1>(indices);
    const int num_elements = sample_end_index - sample_start_index;
    return span<const int>(vec.data() + sample_start_index, num_elements);
}

int data_manager::get_num_cells(const int sample) const
{
    const tuple<int, int> indices = this->sample_indices.at(sample);
    const int sample_start_index = get<0>(indices);
    const int sample_end_index = get<1>(indices);
    return sample_end_index - sample_start_index;
}

void data_manager::fill(const vector<vector<int>> & character_matrix, 
                        const vector<vector<int>> & variant_reads, 
                        const vector<vector<int>> & total_reads,
                        const vector<vector<int>> & region_reads,
                        const vector<int> & cell_samples,
                        const vector<int> & locus_regions,
                        const vector<bool> & is_germline,
                        const vector<bool> & is_sbs,
                        const vector<int> & locus_samples,
                        const vector<int> & chromosomes)
{
    // assign member variables
    this->character_matrix = character_matrix;
    this->variant_reads = variant_reads;
    this->total_reads = total_reads;
    this->region_reads = region_reads;
    this->cell_samples = cell_samples;
    this->locus_regions = locus_regions;
    this->is_germline = is_germline;
    this->is_sbs = is_sbs;
    this->locus_samples = locus_samples;
    this->chromosomes = chromosomes;

    // fill in basic data
    this->num_cells = character_matrix.size();
    this->num_loci = character_matrix[0].size();
    this->num_regions = this->region_reads[0].size();

    // find list of unique sample ids in ascending order
    this->samples = cell_samples;
    sort(this->samples.begin(), this->samples.end());
    this->samples.erase(unique(this->samples.begin(), this->samples.end()), this->samples.end());\
    this->find_sample_indices();

    // create map from region index to chromosome
    for(int i = 0; i < static_cast<int>(this->chromosomes.size()); ++i)
        this->region_to_chr[i] = this->chromosomes[i];

    // initialize sum of reads for each cell, and the reliabilty of each region
    this->cell_read_sums = sum(this->region_reads, 1);
    this->reliable_regions.resize(this->samples.size(), vector<bool>(this->num_regions, false));
    this->region_weights.resize(this->samples.size(), vector<double>(this->num_regions, 0.0));

    // define variant orders for each sample
    for(int i = 0; i < static_cast<int>(this->samples.size()); ++i)
        this->variant_orders[this->samples[i]] = vector<int>();
    this->define_variant_orders();
}

void data_manager::compute_region_weights(const int sample, 
                                          const vector<int> & cell_assignments, 
                                          const vector<int> & parents)
{
    auto threshold = 1.0 / (double)this->num_regions / 15.0;
    vector<int> region_sums(this->num_regions, 0);
    tuple<int, int> indices = this->sample_indices.at(sample);
    const int sample_start_index = get<0>(indices);
    const int sample_end_index = get<1>(indices);
    const int num_cells_sample = sample_end_index - sample_start_index;

    // determine which regions are reliable (same as COMPASS)
    for(auto r_i = 0; r_i < this->num_regions; ++r_i)
    {
        double mean_frac = 0.0;
        int under_threshold = 0;
        for(int i = sample_start_index; i < sample_start_index + num_cells_sample; ++i)
        {
            double frac = (double)this->region_reads[i][r_i] / (double)this->cell_read_sums[i];
            mean_frac += frac / (double)num_cells_sample;
            if(frac <= threshold)
                under_threshold += 1;
        }

        if(((double)under_threshold / (double)num_cells_sample) <= 0.04 && mean_frac >= 0.2 / (double)this->num_regions)
            this->reliable_regions[sample][r_i] = true;
    }

    // if root has no cells assigned to it, then the sample must consist of 100% cancerous cells
    // and if there are two distinct lineages coming off of the root, then we cannot appl
    unordered_map<int,int> assignment_counts;
    assignment_counts[ROOT] = 0;
    for(auto i = sample_start_index; i < sample_start_index + num_cells_sample; ++i)
    {
        if(assignment_counts.find(cell_assignments[i]) != assignment_counts.end())
            assignment_counts[cell_assignments[i]] += 1;
        else
            assignment_counts[cell_assignments[i]] = 1;
    }

    int root_node = ROOT;
    int min_cells = 5;
    if(assignment_counts[ROOT] < min_cells)
    {
        int child_of_root = -1;
        for(int i = 0; i < static_cast<int>(parents.size()); ++i)
        {
            if(parents[i] == ROOT)
            {
                if(child_of_root == -1 && assignment_counts[i+1] >= min_cells)
                    child_of_root = i+1;
                else 
                    throw std::runtime_error("Data likely has 100% purity and some cells likely have no common ancestor");

            }
        }
        root_node = child_of_root;
    }

    // compute probability a read falls in a region
    for(auto i = sample_start_index; i < sample_start_index + num_cells_sample; ++i)
    {
        if(cell_assignments[i] == root_node)
            for(auto r_i = 0; r_i < this->num_regions; ++r_i)
                region_sums[r_i] += this->region_reads[i][r_i];
    }
    
    int total = accumulate(region_sums.begin(), region_sums.end(), 0);
    for(auto r_i = 0; r_i < this->num_regions; ++r_i)
        this->region_weights[sample][r_i] = (double)region_sums[r_i] / (double)total;
}

void data_manager::find_sample_indices(void)
{
    for(const auto & sample : this->samples)
    {
        int sample_start_index, sample_end_index = this->cell_samples.size();
        bool found_start = false;
        for(int i = 0; i < static_cast<int>(this->cell_samples.size()); ++i)
        {
            if(this->cell_samples[i] == sample)
            {
                if(found_start)
                {
                    continue;
                }
                else
                {
                    found_start = true;
                    sample_start_index = i;
                    sample_end_index = cell_samples.size();
                    continue;
                }
            }
            else if(found_start)
            {
                sample_end_index = i;
                break;
            }

        }
        this->sample_indices[sample] = tuple<int, int>(sample_start_index, sample_end_index);
    }
}

void data_manager::define_variant_orders(void)
{
    for(const auto & sample : this->samples)
    {
        const tuple<int, int> indices = this->sample_indices.at(sample);
        const int sample_start_index = get<0>(indices);
        const int sample_end_index = get<1>(indices);

        // select loci that were first mutated in this sample
        vector<int> sample_loci;
        for (int i = 0; i < static_cast<int>(this->locus_samples.size()); ++i) 
            if (this->locus_samples[i] == sample)
                sample_loci.push_back(i);

        // sum rows for cells in this sample for the columns corresponding to these loci
        vector<int> locus_sums(sample_loci.size(), 0);
        for (int idx = 0; idx < static_cast<int>(sample_loci.size()); ++idx) 
        {
            const int & locus = sample_loci[idx];
            for (auto cell = sample_start_index; cell < sample_end_index; ++cell) 
            {
                if(this->character_matrix[cell][locus] == 1)
                    locus_sums[idx]++;
            }
        }

        // Sort sample_loci by descending order of locus_sums
        vector<int> order(sample_loci.size());
        iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(), [&](int i, int j) {
            return locus_sums[i] > locus_sums[j];
        });

        // Append ordered loci to variant_orders if they are SNPs
        for (const auto & idx : order) {
            const auto & locus = sample_loci[idx];
            if (!this->is_germline[locus]) {
                this->variant_orders[sample].push_back(locus);
            }
        }
    }
}
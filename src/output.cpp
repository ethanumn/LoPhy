#include "output.h"
#include "csv.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <limits>
#include <filesystem>
#include <unordered_map>

using std::stringstream;
using std::ofstream;
using std::ifstream;
using std::ostringstream;
using std::to_string;
using std::numeric_limits;
using std::to_string;
using std::endl;
using std::distance;
using std::sort;
using std::unique;
using std::unordered_map;
using std::find;

/* Adjust a vector containing integer values into one that contains unique contiguous integers, e.g., [4, 2, 5, 9] -> [2,1,3,4] */
vector<int> adjust_to_contiguous(const vector<int> & vec, const vector<int> & unique_values) 
{
    
    int start_val = 0;
    vector<int> new_vec = vec;
    auto it = find(unique_values.begin(), unique_values.end(), ROOT);
    if(it == unique_values.end())
        start_val=1;

    // Create a mapping from original value to contiguous value
    unordered_map<int, int> value_map;
    for (int i = 0; i < static_cast<int>(unique_values.size()); ++i) 
        value_map[unique_values[i]] = start_val + i;
    
    // Adjust the original vector using the new mapping
    for (int & num : new_vec)
        num = value_map[num];
    return new_vec;
}

template <typename T>
vector<vector<T>> extract_rows(const vector<vector<T>> & matrix, const vector<int> & indices) 
{
    // Check if the matrix is empty
    if (matrix.empty()) 
        throw invalid_argument("Matrix is empty.");

    // matrix to store data
    vector<vector<T>> new_matrix;

    // Iterate through the list of indices
    for (int index : indices) 
    {
        // Check if the index is within the bounds of the matrix
        if (index < 0 || index >= static_cast<int>(matrix.size())) 
            throw out_of_range("Index out of range.");

        // Add the corresponding row to the matrix
        new_matrix.push_back(matrix[index]);
    }

    return new_matrix;
}


/* Joins elements in a vector with a separator */
template <typename Iter>
string join(Iter begin, Iter end, string const & sep)
{
  ostringstream result;
  if (begin != end)
    result << "'" << *begin++ << "'";
  while (begin != end)
    result << sep << "'" << *begin++ << "'";
  return result.str();
}

/* Turns a matrix into a string */
template <typename T>
string matrix_to_graph_string(string key, const vector<vector<T>> & matrix) 
{
    ostringstream output_ss;
    output_ss << key << "=\"[";
    const int matrix_size = static_cast<int>(matrix.size());

    for (int i = 0; i < matrix_size; ++i) {
        output_ss << "[";
        // Convert each row (vector<double>) to a comma-separated string
        output_ss << join(matrix[i].begin(), matrix[i].end(), ",");
        output_ss << "]";
        if (i < matrix_size - 1) {
            output_ss << ","; // Separate rows with commas
        }
    }

    output_ss << "]\",\n";
    return output_ss.str();
}


/* Determines how to name output files */
string out_path_prefix(const string & in_file, 
                       const string & out_path, 
                       const string & prefix)
{
    string new_out_path = out_path; // create new outpath depending on input

    // make sure out_path ends in a backslash
    if(new_out_path.back() != '/')
        new_out_path += '/';

	if(prefix.empty())
    {
        std::filesystem::path p(in_file);
		return new_out_path + p.stem().string();
	} 
    else 
    {
    	return new_out_path + prefix;
    }
}


/* Reads weights from a file indicating the probability a read falls into a genomic region */
vector<double> read_weights(const string & weights_file, 
                            const int & num_regions)
{
    vector<double> region_weights;
	ifstream file(weights_file.c_str());

    if(file) 
    {
        // read weights
        float w;
        for (auto i = 0; i < num_regions; i++) {
            file >> w;
            region_weights.push_back(w);
        }
    }
    return region_weights;
}

/* Creates a string indicating each cell's attachment point in the clone tree */
string compute_attachment_strings(const vector<int> cell_assignments,
                                  const vector<string> & cell_names)
{
	stringstream attachment_ss;

    for(int i = 0; i < static_cast<int>(cell_assignments.size()); ++i)
    {
        auto index = cell_assignments[i];
        string node_id = "Clone_" + to_string(index);
        attachment_ss << "\"" << node_id << "\"" << " -> " << "\"" << cell_names[i] << "\""  << ";\n";
    }

	return attachment_ss.str();
}

/* Write params to string */
string params_to_string(const Params & params)
{
    stringstream ss;
    ss << "      fp=" << params.fp << ",\n";
    ss << "      fn=" << params.fn << ",\n";
    ss << "      fn_indel=" << params.fn_indel << ",\n";
    ss << "      hom_precision=" << params.hom_precision << ",\n";
    ss << "      het_precision=" << params.het_precision << ",\n";
    ss << "      hom_precision_indel=" << params.hom_precision_indel << ",\n";
    ss << "      het_precision_indel=" << params.het_precision_indel << ",\n";
    ss << "      dropout_concentration=" << params.dropout_concentration << ",\n";
    ss << "      dropout_rate_prior=" << params.dropout_rate_prior << ",\n";
    ss << "      theta_param=" << params.theta << ",\n";
    ss << "      iters=" << params.iters << ",\n";
    ss << "      delta=" << params.delta << ",\n";
    ss << "      CNA_penalty=" << params.CNA_penalty << ",\n";
    ss << "      tau1=" << params.tau1 << ",\n";
    ss << "      tau2=" << params.tau2 << ",\n";
    ss << "      seed=" << params.seed << ",\n";
    ss << "      character_matrix_file=\"" << params.character_matrix_file << "\",\n";
    ss << "      meta_file=\"" << params.meta_file << "\",\n";
    ss << "      out_path=\"" << params.out_path << "\",\n";
    ss << "      out_prefix=\"" << params.out_prefix << "\",\n";
    ss << "      variant_reads_file=\"" << params.variant_reads_file << "\",\n";
    ss << "      total_reads_file=\"" << params.total_reads_file << "\",\n";
    ss << "      region_reads_file=\"" << params.region_reads_file << "\",\n";
    ss << "      cell_samples_file=\"" << params.cell_samples_file << "\",\n";
    return ss.str();
}

/* writes the out file for a particular parents vector */
void write_out_file(const string out_file, 
                    const string graph_string, 
                    const string attachment_string,
                    const Tree & tree,
                    const data_manager & manager,
                    const Params & params,
                    const vector<vector<int>> & alt_cns,
                    const vector<vector<int>> & total_cns,
                    const vector<vector<int>> & genotypes,
                    const vector<int> & cell_assignments,
                    const vector<string> & cell_names,
                    const vector<string> & variants)
{
    const vector<int> & samples = manager.get_samples();
    const int n_samples = static_cast<int>(samples.size());
    const double & theta = manager.get_theta();
    const vector<vector<double>> & dropout_rates_collection = tree.get_dropout_rates_collection();
    const vector<vector<double>> & alt_dropout_rates_collection = tree.get_alt_dropout_rates_collection();
    const vector<vector<double>> & ref_dropout_rates_collection = tree.get_ref_dropout_rates_collection();
    const vector<vector<double>> & priors_collection = tree.get_priors_collection(); 
    vector<vector<double>> region_weights(n_samples, vector<double>()); 
    
    for(int i = 0; i < n_samples; ++i)
        region_weights[i] = manager.get_region_weights(samples[i]);

    // build string content of file
	stringstream output_ss;
	output_ss << "digraph G {\n";
    output_ss << "  graph [cells=\"[" << join(cell_names.begin(), cell_names.end(), ",") << "]\",\n";
    output_ss << "      variants=\"[" << join(variants.begin(), variants.end(), ",") << "]\",\n";
    output_ss << matrix_to_graph_string("mutant_copy_numbers", alt_cns);
    output_ss << matrix_to_graph_string("total_copy_numbers", total_cns);
    output_ss << matrix_to_graph_string("genotypes", genotypes);
    output_ss << matrix_to_graph_string("dropout_rates", dropout_rates_collection);
    output_ss << matrix_to_graph_string("alt_dropout_rates", alt_dropout_rates_collection);
    output_ss << matrix_to_graph_string("ref_dropout_rates", ref_dropout_rates_collection);
    output_ss << matrix_to_graph_string("region_weights", region_weights);
    output_ss << matrix_to_graph_string("priors", priors_collection);
    output_ss << "      cell_assignments=\"[" << join(cell_assignments.begin(), cell_assignments.end(), ",") << "]\",\n";
    output_ss << "      theta=" << theta << ",\n";
    output_ss << "      root_id=0,\n";
    output_ss << "      root_name=\"" << "Clone_0" << "\",\n";
    output_ss << "      type=cell_tree,\n";
    output_ss << params_to_string(params);
    output_ss << "\n];";
	output_ss << "node [color=deeppink4, style=filled, fontcolor=white];\n";
    output_ss << graph_string;
	output_ss << "node [color=lightgrey, style=filled, fontcolor=black];\n";
    output_ss << attachment_string;
	output_ss <<  "}\n";

    // write content to file
	ofstream output;
	output.open(out_file);
	output << output_ss.str();
	output.close();
}

/* Function for writing all clone trees and their corresponding data to file. */
void write_all_results(string path_with_prefix, 
                       const vector<string> & cell_names,
                       const vector<string> & gene_names,
                       const vector<string> & region_names,
                       const Tree & tree, 
                       const data_manager & manager,
                       const Params & params)
{
    // initialize variables 
    string attachment_string;
    const int m = manager.get_num_loci();
    const int & num_regions = manager.get_num_regions();
    const vector<int> & parents = tree.get_parents();
    const int num_clones = parents.size() + 1;
    vector<vector<int>> A(num_clones, vector<int>(num_regions, 0));
    vector<vector<int>> B(num_clones, vector<int>(m, 0));  
    const auto & chromosomes = manager.get_chromosomes();
    const auto & locus_regions = manager.get_locus_regions();

    // get relevant data from data_manager
    const vector<bool> & is_germline = manager.get_is_germline();
    const CNAs cnas = tree.get_cnas();

    stringstream out_file_ss;
    out_file_ss << path_with_prefix << "_ml" << 0 << ".gv";
    stringstream graph_ss;
    vector<string> variants;

    // Each vector contains the SNPs/SNVs that make up a clone
    vector<vector<string>> clone_variants(num_clones, vector<string>());

    // If some of the variants are marked as SNPs, then these will be used to name the root clone
    vector<string> root_strs;
    for(auto j = 0; j < m; ++j)
        if(is_germline[j])
            root_strs.push_back(gene_names[j+1]);
    if(root_strs.size() == 0)
        clone_variants[ROOT].push_back("root");

    // add variants initially to their own clones
    for(auto clone = 1; clone <= m; ++clone)
    {
        ostringstream oss;
        oss << gene_names[clone] << " (chr" << chromosomes[locus_regions[clone-1]] << ")";
        clone_variants[clone].push_back(oss.str());
    }
    
    // merge variants basd on their clone assignment
    // we use a little trick here -- if parents[j] > 2*m, then mutation j is actually part of clone parents[j] - 2*m
    for(auto j = 0; j < m; ++j)
    {   
        if(parents[j] == NO_PARENT && is_germline[j])
        {
            ostringstream oss;
            oss << gene_names[j+1] << " (chr" << chromosomes[locus_regions[j]] << ")";
            clone_variants[ROOT].push_back(oss.str());
            clone_variants[j+1].clear();
        }
        else if(parents[j] > 2*m)
        {
            int p = parents[j] - 2*m;
            ostringstream oss;
            oss << gene_names[j+1] << " (chr" << chromosomes[locus_regions[j]] << ")";
            clone_variants[p].push_back(oss.str());
            clone_variants[j+1].clear();
        }
    }

    // now that we have all mutation names assigned to each clone, join then into one string
    vector<string> graph_strings(num_clones);
    for(auto clone = 0; clone < num_clones; ++clone)
    {
        auto & v = clone_variants[clone];
        if(v.empty() && clone <= m)
            continue;
        else 
        {
            bool is_cna_clone = true;
            stringstream ss;

            if(clone <= m)
            {
                is_cna_clone = false; // this clone has at least one SNV

                // add all SNVs/SNPs
                sort(v.begin(), v.end());
                for (int vi = 0; vi < static_cast<int>(v.size()); ++vi) 
                {
                    if (vi != 0) {
                        ss << "\\n";
                    }
                    ss << v[vi];
                }
            }

            bool added_first_cna = false;

            for(const auto & cna : tree.get_clone_cnas(clone))
            {
                const auto & region = get<0>(cna);
                const auto & type = get<1>(cna);
                
                // add a new line if (1) the clone contains variants or 
                // (2) we've already added at least one CNA
                if(!is_cna_clone || added_first_cna)
                {
                    ss << "\\n";
                }
                else 
                {
                    added_first_cna = true;
                }

                switch(type)
                {
                    case LOSS_ALT:
                        ss << "Loss " << region_names[region] << " ALT";
                        break;
                    case LOSS_REF:
                        ss << "Loss " << region_names[region] << " REF";
                        break;
                    case GAIN_ALT:
                        ss << "Gain " << region_names[region] << " ALT";
                        break;
                    case GAIN_REF:
                        ss << "Gain " << region_names[region] << " REF";
                        break;
                    case CNLOH_ALT:
                        ss << "CNLOH " << region_names[region] << " ALT";
                        break;
                    case CNLOH_REF:
                        ss << "CNLOH " << region_names[region] << " REF";
                    default:
                        break;
                }
            }
            graph_strings[clone] = ss.str();
        }
    }

    // create the graphviz string for each cell's attachment in the clone tree
    attachment_string = compute_attachment_strings(tree.get_cell_assignments(), cell_names);

    // First, define all nodes with unique IDs and their labels
    for (int clone = 0; clone < num_clones; ++clone) 
    {
        if (parents[clone - 1] == NO_PARENT || parents[clone - 1] > 2 * m)
            continue;
        string node_id = "Clone_" + to_string(clone);
        graph_ss << node_id << " [label=\"" << graph_strings[clone] << "\"];\n";
    }

    // create the graphviz string for the clone tree
    for (int clone = 1; clone < num_clones; ++clone) {
        if (parents[clone - 1] == NO_PARENT || parents[clone - 1] > 2 * m)
            continue;

        variants.push_back(gene_names[clone]);
        int index = parents[clone-1];
        string parent_id = "Clone_" + to_string(index);
        string child_id = "Clone_" + to_string(clone);
        graph_ss << parent_id << " -> " << child_id << ";\n";
    }

    // remove data for empty clones, i.e., cleaning up data structures before writing it to file
    vector<int> indices_to_keep;
    for(int idx = 0; idx < static_cast<int>(graph_strings.size()); ++idx)
    {
        if(!graph_strings[idx].empty())
            indices_to_keep.push_back(idx);
    }

    graph_strings.erase(graph_strings.begin()); // remove root node from graph_strings

    graph_strings.erase(std::remove_if(graph_strings.begin(), graph_strings.end(), [](const string & s) {
        return s.empty();
    }), graph_strings.end());

    // fill copy number matrices based on tree
    for(auto & row : A)
        fill(row.begin(), row.end(), 0);
    for(auto & row : B)
        fill(row.begin(), row.end(), 0);
    tree.CNAs_to_matrices(manager, A, B);
    vector<vector<int>> A_prime = extract_rows(A, indices_to_keep);
    vector<vector<int>> B_prime = extract_rows(B, indices_to_keep);
    vector<vector<int>> genotypes = B_prime;
    for(auto & g : genotypes)
    {
        for(auto & g_i : g)
        {
            if(g_i > 0)
                g_i = 1; 
        }
    }

    // write all data to a file 
    write_out_file(out_file_ss.str(), 
                    graph_ss.str(), 
                    attachment_string, 
                    tree,
                    manager,
                    params,
                    B_prime,
                    A_prime,
                    genotypes,
                    adjust_to_contiguous(tree.get_cell_assignments(), indices_to_keep),
                    cell_names, 
                    graph_strings);

}

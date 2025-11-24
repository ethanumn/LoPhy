#pragma once 

#include "utils.h"
#include "data_manager.h"
#include "tree.h"

#include <string.h>
#include <vector>
#include <map>

using std::vector;
using std::string;
using std::map;

string out_path_prefix(const string & in_file, 
                       const string & out_path, 
                       const string & prefix);

vector<double> read_weights(const string & weights_file, 
                            const int & num_regions);

vector<double> get_region_read_probs(const string & region_probs_file, 
                                     const int & num_regions);

void write_all_results(string prefix, 
                       const vector<string> & cell_names,
                       const vector<string> & gene_names,
                       const vector<string> & region_names,
                       const Tree & tree, 
                       const data_manager & manager,
                       const Params & params);



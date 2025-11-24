#pragma once

#include "utils.h"
#include "tree.h"
#include "data_manager.h"
#include <vector>
#include <random>
#include <map>

using std::mt19937;
using std::vector;
using std::map; 

/* Function for computing 1-node extensions */
Tree find_extension(Tree tree, 
                    const data_manager & manager,
                    SCORE_CACHE & cache,
                    const int & sample,
                    const int & v,
                    bool infer_cnas,
                    int seed);

/* Main search function */
Tree search(Tree tree,
            data_manager & manager,
            SCORE_CACHE & cache,
            const int & sample,
            bool infer_cnas,
            bool include_prev_variants,
            int seed);
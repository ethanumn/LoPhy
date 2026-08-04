#include "search.h"
#include "tree.h"

#include <set>
#include <algorithm>
#include <limits>
#include <iterator>
#include <utility>

using std::vector;
using std::min;
using std::numeric_limits;
using std::sort;
using std::set;
using std::swap; 
using std::unordered_map;
using std::make_pair;

/* Add or removes CNAs from the tree */
Tree add_remove_CNAs(Tree tree,
                     const data_manager & manager,                       
                     SCORE_CACHE & cache,
                     const int & sample,
                     mt19937 & gen,
                     bool evaluate_all = false)
{
    // get relevant data from manager
    const auto & num_loci = manager.get_num_loci();
    const auto & num_cells = manager.get_num_cells(sample);
    const auto & locus_regions = manager.get_locus_regions();
    const auto & parents = tree.get_parents();
    const auto & variants = tree.get_variants();

    // find clones first identified in this sample
    vector<int> clones_in_sample;
    const vector<int> & variants_in_sample = manager.get_variant_order(sample);
    const vector<int> dummy_clones_in_sample = tree.get_dummy_clones_in_sample(sample);
    for(const auto locus : variants_in_sample)
        if((tree.get_parent(locus) != NO_PARENT) && (tree.get_parent(locus) <= 2*num_loci))
            clones_in_sample.push_back(locus+1);
    for(const auto index : dummy_clones_in_sample)
        if(tree.get_parent(index) != NO_PARENT)
            clones_in_sample.push_back(index+1);

    tree.find_candidate_regions(manager, clones_in_sample, sample); // update candidate regions

    unordered_set<int> ancestors;
    vector<pair<int, int>> cna_possibilities;
    unordered_set<int> candidate_regions;
    // for each clone in this sample, see what regions can be impacted by CNAs
    for(const auto & clone : clones_in_sample)
    {
        // get ancestors of this clone
        get_ancestors(ancestors, parents, clone);
        vector<int> ancestors_vec(ancestors.begin(), ancestors.end());
        sort(ancestors_vec.begin(), ancestors_vec.end());

        if(clone <= num_loci)
            candidate_regions = tree.get_candidate_regions(clone);
        else 
            candidate_regions = tree.get_candidate_regions();

        for(const auto & r : candidate_regions)
        {
            bool already_impacted = false;
            // if this region is already impacted by a CNA in an ancestor that's not in this sample, do not modify it
            for(const auto & ancestor : ancestors_vec)
            {
                if(tree.has_CNA_in_region(ancestor, r) && find(clones_in_sample.begin(), clones_in_sample.end(), ancestor) == clones_in_sample.end())
                {
                    already_impacted = true;
                    break;
                }
            }
                    
            // if this region hasn't already been impacted we can try adding CNAs to it
            if(!already_impacted)
                cna_possibilities.push_back(make_pair(clone, r));
        }
    }

    // sample some of the (clone, region) pairs without replacement
    int n_possibilities = cna_possibilities.size();
    vector<int> indices(n_possibilities);
    iota(indices.begin(), indices.end(), 0); 
    shuffle(indices.begin(), indices.end(), gen);

    const double original_tree_likelihood = tree.get_llh();
    Tree best_tree = tree;

    int n_to_eval = std::min(1,n_possibilities);
    if(evaluate_all)
        n_to_eval = n_possibilities;

    // evaluate possible CNA's to see which ones are most likely
    for (auto idx = 0; idx < n_to_eval; ++idx) 
    {
        // get the clone/region to evaluate
        const auto clone =  cna_possibilities[indices[idx]].first;
        const auto region = cna_possibilities[indices[idx]].second;

        // determine which types of CNAs are valid for this clone
        vector<CNA_TYPE> valid_cnas;
        const int max_cn = manager.get_max_cn();
        switch(max_cn) {
            case 5:
                valid_cnas = {NA, LOSS_REF, GAIN_REF, GAIN_REF_2, GAIN_REF_3};
                break;
            case 4:
                valid_cnas = {NA, LOSS_REF, GAIN_REF, GAIN_REF_2};
                break;
            case 3:
                valid_cnas = {NA, LOSS_REF, GAIN_REF};
                break;
            case 2:
                valid_cnas = {NA, LOSS_REF};
                break;
            default:
                throw std::runtime_error("max_cn must be at least 2 for there to be any CNAs");
        }

        vector<int> loci_in_region;
        bool at_least_one_variant = false;

        // check if there is a locus mutated in at least one clone (in this clone or ancestral) in this lineage
        get_ancestors(ancestors, parents, clone);
        vector<int> clones_to_check(ancestors.begin(), ancestors.end());
        clones_to_check.push_back(clone); // SNVs in this clone can also be impacted by CNAs
        clones_to_check.push_back(ROOT); // need to add this to check if SNPs are impacted by CNAs in descendant clones
        for(const auto & c : clones_to_check)
        {
            for(const auto & locus : variants[c])
            {
                if(locus_regions[locus] == region)
                {
                    loci_in_region.push_back(locus);

                    // if there's at least one locus that has a variant allele 
                    // then we'll allow CNLOH and GAIN/LOSS of the variant allele
                    if(!at_least_one_variant)
                    {   
                        at_least_one_variant = true;
                        valid_cnas.push_back(LOSS_ALT);
                        valid_cnas.push_back(CNLOH_ALT);
                        valid_cnas.push_back(CNLOH_REF);
                        switch(max_cn) {
                            case 5:
                                valid_cnas.push_back(GAIN_ALT);
                                valid_cnas.push_back(GAIN_ALT_2);
                                valid_cnas.push_back(GAIN_ALT_3);
                                break;
                            case 4:
                                valid_cnas.push_back(GAIN_ALT);
                                valid_cnas.push_back(GAIN_ALT_2);
                                break;
                            case 3:
                                valid_cnas.push_back(GAIN_ALT);
                                break;
                            case 2:
                                break;
                            default:
                                throw std::runtime_error("max_cn must be at least 2 for there to be any CNAs");
                        }
                    }
                }
            }
        }

        for(const CNA_TYPE & cna_type : valid_cnas)
        {
            Tree tree_prime = best_tree;
            tree_prime.add_cna(clone,
                               region,
                               cna_type,
                               clones_in_sample);

            tree_prime.update(manager, sample); // update tree with new data structures
            tree_prime.score(manager, sample, cache);
            if(cna_type == GAIN_REF_2 || cna_type == GAIN_REF_3 || cna_type == GAIN_ALT_2 || cna_type == GAIN_ALT_3)
                cout << "Evaluating CNA: " << cna_type << " for clone " << clone << " in region " << region << " with likelihood: " << tree_prime.get_llh() << endl;

            const double delta_best = tree_prime.get_llh() - best_tree.get_llh();

            const double delta_original = tree_prime.get_llh() - original_tree_likelihood;

            // only add the CNA if it improves the likelihood of the tree and the likelihood of the tree is better than the original tree by a significant margin (to avoid overfitting)
            if ((delta_best > 0) && (delta_original > num_cells * log(manager.get_psi())))
            {
                best_tree = tree_prime;
            }
        }
    }
    return best_tree;
}

/* Find all possible modifications of the tree */
Tree find_extension(Tree tree, 
                    const data_manager & manager,
                    SCORE_CACHE & cache,
                    const int & sample,
                    const int & v,
                    bool infer_cnas,
                    int seed)
{
    // get data
    const auto phi_hat_prev = tree.get_phi_hat(); // MAKE SURE THIS IS A COPY
    const auto num_loci = manager.get_num_loci();
    const auto & samples = manager.get_samples();

    // find clones in first identified in this sample
    vector<int> clones_in_sample;
    if(sample == samples[0])
        clones_in_sample.push_back(ROOT);
    const vector<int> & variants_in_sample = manager.get_variant_order(sample);
    const vector<int> dummy_clones_in_sample = tree.get_dummy_clones_in_sample(sample);
    for(const auto locus : variants_in_sample)
        if((tree.get_parent(locus) != NO_PARENT) && (tree.get_parent(locus) <= 2*num_loci))
            clones_in_sample.push_back(locus+1);
    for(const auto index : dummy_clones_in_sample)
        if(tree.get_parent(index) != NO_PARENT)
            clones_in_sample.push_back(index+1);

    // initialize variables for search
    vector<int> children;
    vector<vector<int>> choices;
    mt19937 gen(seed);

    // keeping track of best tree
    bool found_tree = false;
    Tree best_tree = tree;
    double best_llh = -std::numeric_limits<double>::infinity();

    // find possible parents for clone v 
    vector<int> possible_parents = {ROOT};
    for(int i = 0; i < static_cast<int>(tree.get_parents().size()); ++i)
        if(tree.get_parent(i) != NO_PARENT)
            if((int)tree.get_parent(i) <= 2*num_loci)
                possible_parents.push_back(i + 1); 

    // try all ways to incorporate v into the tree
    for(const auto & u : possible_parents) 
    {
        // fix u -> v
        Tree tree_prime = tree;
        tree_prime.set_parent(v, u);

        // (1) score u -> v
        tree_prime.update(manager, sample); // update tree with new data structures
        tree_prime.score(manager, sample, cache);
        if(infer_cnas)
            tree_prime = add_remove_CNAs(tree_prime, manager, cache, sample, gen, true); // sample CNAs
        if(!found_tree || tree_prime.get_llh() > best_llh)
        {
            best_llh = tree_prime.get_llh();
            best_tree = std::move(tree_prime);
            found_tree = true;
        }

        // (2) score v is added to clone u
        auto it = find(clones_in_sample.begin(), clones_in_sample.end(), u);
        if((u != ROOT) && (u <= num_loci) && (v < num_loci) && it != clones_in_sample.end()) // only add SNVs to clones with at least one SNVs
        {
            Tree tree_prime = tree;
            tree_prime.set_parent(v, 2*num_loci+u);
            tree_prime.update(manager, sample); // update tree with new data structures
            tree_prime.score(manager, sample, cache);
            if(infer_cnas)
                tree_prime = add_remove_CNAs(tree_prime, manager, cache, sample, gen, true); // sample CNAs
            if(!found_tree || tree_prime.get_llh() > best_llh)
            {
                best_llh = tree_prime.get_llh();
                best_tree = std::move(tree_prime);
                found_tree = true;
            }
        }

        // (3) score u -> v -> ch, where ch is some subset of child clones first identified in this sample
        children.clear();
        for(const auto & c : clones_in_sample)
        {
            if(c == ROOT)
                continue;
            else if(tree.get_parent(c-1) == u)
                children.push_back(c);
        }

        choices.clear(); // clear all ways to incorporate v into the tree when it's a parent of one or more clones in the tree

        // iterate through all possible ways to incorporate v into the tree when it's a child of u 
        // and a parent of one or more clones
        for(int ck = 1; ck <= static_cast<int>(children.size()); ++ck)
            comb(choices, children, ck); // fill with all choices
        
        // now evaluate all possible u -> v -> ch
        for(const vector<int> & ch : choices)
        {
            Tree tree_prime = tree;
            tree_prime.set_parent(v, u);
            // make v the parent of the nodes in ch
            for(const int & c : ch)
                tree_prime.set_parent(c-1, v+1);

            // score u -> v -> ch
            tree_prime.update(manager, sample); // update tree with new data structures
            tree_prime.score(manager, sample, cache);
            if(infer_cnas)
                tree_prime = add_remove_CNAs(tree_prime, manager, cache, sample, gen, true); // sample CNAs
            if(!found_tree || tree_prime.get_llh() > best_llh)
            {
                best_llh = tree_prime.get_llh();
                best_tree = std::move(tree_prime);
                found_tree = true;
            }
        }
    }

    if(!found_tree)
        throw std::runtime_error("No valid extensions for tree found!");

    best_tree.set_phi_hat(best_tree.get_llh() + phi_hat_prev);

    return best_tree;
}

/* Main search function */
Tree search(Tree tree,
            data_manager & manager,
            SCORE_CACHE & cache,
            const int & sample,
            bool infer_cnas,
            bool include_prev_variants,
            int seed)
{

    // collect variants to include in tre
    vector<int> variants_in_sample;
    if(include_prev_variants)
    {
        // samples should already be in ascending order from earliest to latest
        for(const auto prev_sample : manager.get_samples())
        {
            if(prev_sample <= sample)
            {
                const vector<int> variants_in_prev_sample = manager.get_variant_order(prev_sample);
                variants_in_sample.insert(variants_in_sample.end(), variants_in_prev_sample.begin(), variants_in_prev_sample.end());
            }
        }
    }
    else 
        variants_in_sample = manager.get_variant_order(sample);

    // add all variants from this sample to the tree to initialize tree structure 
    for(const auto & v : variants_in_sample)
    {
        seed = (seed + 1) % numeric_limits<int>::max();
        tree = find_extension(tree, 
                              manager,
                              cache,
                              sample,
                              v,
                              false,
                              seed);
    } 

    // even if sample does not have any new variants, we still want to assign cells to clones
    if(variants_in_sample.size() == 0)
    {
        tree.update(manager, sample);
        tree.score(manager, sample, cache);
    }

    // hill climbing
    // the base case is: if we don't have any new variants in this sample, we'll just assign cells to clones and score the tree 
    // otherwise, perform a hill climbing search to find a better tree:
    // (1) draw a random number alpha ~ U(0,1)
    // (2) if alpha < tau1, randomly select a variant v and remove it from the tree and find the best way to re-add it to the tree
    // (3) if alpha - tau1 < tau2, add a dummy clone to the tree and find the best way to add it and modify the CNA inheritance
    // (4) otherwise, find the best way to add/remove CNAs from the tree
    if(variants_in_sample.size() > 0)
    {
        mt19937 gen(seed);

        std::uniform_real_distribution<> U(0.0, 1.0); // uniform distribution over reals (0,1)
        std::uniform_int_distribution<int> v_dist(0, variants_in_sample.size() - 1); // uniform distribution over variants

        int i = 0;
        int n_since_improvement = 0;
        
        seed = (seed + 1) % numeric_limits<int>::max();
        const double tau1 = manager.get_tau1();
        const double tau2 = manager.get_tau2();
        const int iters = manager.get_iters();
        const int delta = manager.get_delta();
        while(i < iters && n_since_improvement < delta)
        {

            // make sure cache doesn't get too large
            if(cache.size() > 100000)
            {
                cache.clear();
            }
            
            double alpha = U(gen); // (1)
            bool improving = true;
            Tree tree_prime = tree;

            
            if(alpha < tau1 || !infer_cnas) // (2)
            {
                const int v = variants_in_sample[v_dist(gen)];
                tree_prime.remove_variant(v, variants_in_sample, manager);
                tree_prime = find_extension(tree_prime, 
                                            manager,
                                            cache,
                                            sample,
                                            v,
                                            false,
                                            seed);
            }
            else
            {

                if(alpha-tau1 < tau2) // (3)
                {
                    const auto v = tree_prime.add_dummy_clone(sample);

                    tree_prime = find_extension(tree_prime, 
                                                manager,
                                                cache,
                                                sample,
                                                v,
                                                true,
                                                seed);

                    improving = !tree_prime.remove_empty_dummy_clones(sample);
                }
                else // (4)
                {
                    tree_prime = add_remove_CNAs(tree_prime, manager, cache, sample, gen); 
                }
            }

            if((tree_prime.get_llh() > tree.get_llh()) && improving)
            {
                n_since_improvement = 0;
                tree = tree_prime;
            }
            n_since_improvement++;
            i++;
        }
    }

    tree.save_priors();
    tree.save_dropout_rates();
    return tree;
}
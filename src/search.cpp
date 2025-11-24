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
        vector<CNA_TYPE> valid_cnas = {NA, LOSS_REF, GAIN_REF};
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
                        valid_cnas.push_back(GAIN_ALT);
                        valid_cnas.push_back(CNLOH_REF);
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
            if(tree_prime.get_llh() - best_tree.get_llh() > num_cells * log(1.05)) 
                best_tree = tree_prime;
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
    vector<Tree> S;
    vector<int> children;
    vector<vector<int>> choices;
    mt19937 gen(seed);

    // find possible parents for clone v 
    vector<int> possible_parents = {ROOT};
    for(int i = 0; i < static_cast<int>(tree.get_parents().size()); ++i)
        if(tree.get_parent(i) != NO_PARENT)
            if((int)tree.get_parent(i) <= 2*num_loci)
                possible_parents.push_back(i + 1); 

    // reserve an estimated capacity to avoid repeated reallocations and copies
    S.reserve(possible_parents.size() * 4);

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
        S.emplace_back(std::move(tree_prime));

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
            S.emplace_back(std::move(tree_prime));
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
            S.emplace_back(std::move(tree_prime));


        }
        
    }

    // verify we found at least one valid extension
    if(S.size() == 0)
        throw std::runtime_error("No valid extensions for tree found!");
    
    // sort all extensions by log likelihood in descending order
    sort(S.begin(), S.end(), std::greater<Tree>());

    S[0].set_phi_hat(S[0].get_llh() + phi_hat_prev); // sum of log likelihoods for approximate 

    return S[0];
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

    // add all variants from this sample to the tree
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
    if(infer_cnas && variants_in_sample.size() > 0)
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
            double alpha = U(gen); // draw move randomly 
            bool improving = true;
            Tree tree_prime = tree;

            if(alpha < tau1)
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

                if(alpha-tau1 < tau2)
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
                else 
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
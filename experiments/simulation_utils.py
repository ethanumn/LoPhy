# imports
import random
from collections import defaultdict, deque

import os
import numpy as np
import pandas as pd
import networkx as nx
import pydot
import anndata as ad
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import nbinom, betabinom

import oncophylo as op

# CONSTANTS
MCN_MAE = "MCN-MAE"
TCN_MAE = "TCN_MAE"
FER = "FER"
LOPHY = "LoPhy"
COMPASS = "COMPASS"
SCITE = "SCITE"
LACE = "LACE"
SEED = "seed"

# metrics
def mae(A1, A2):
    "Mean absolute error"
    d1,d2 = A1.shape
    return (1.0/(d1*d2))*np.sum(np.abs(A1-A2))

def fer(A_pred=None, A_true=None, B_pred=None, B_true=None, cell_samples=None):
    "False emergence rate"
    if cell_samples is None:
        raise ValueError("cell_samples must be provided.")

    use_total = A_pred is not None and A_true is not None
    use_alternate = B_pred is not None and B_true is not None

    if not (use_total or use_alternate):
        raise ValueError("At least one of (A_pred, A_true) or (B_pred, B_true) must be provided.")

    n = len(cell_samples)
    unique_samples = np.unique(cell_samples)
    num_samples = len(unique_samples)

    false_emergence_cells = 0

    # Precompute observed copy numbers per sample
    if use_total:
        _, r = A_true.shape
        A_obs_per_sample = []
        for s in unique_samples:
            sample_mask = cell_samples == s
            A_obs = [set(A_true[sample_mask, k]) for k in range(r)]
            A_obs_per_sample.append(A_obs)

    if use_alternate:
        _, m = B_true.shape
        B_obs_per_sample = []
        for s in unique_samples:
            sample_mask = cell_samples == s
            B_obs = [set(B_true[sample_mask, i]) for i in range(m)]
            B_obs_per_sample.append(B_obs)

    # Check each cell for false emergence
    for i, s in enumerate(unique_samples):
        sample_mask = cell_samples == s
        cells_in_sample = np.where(sample_mask)[0]

        # Observed in future samples
        if use_total:
            A_future_obs = [set() for _ in range(r)]
            for j in range(i + 1, num_samples):
                for k in range(r):
                    A_future_obs[k].update(A_obs_per_sample[j][k])

        if use_alternate:
            B_future_obs = [set() for _ in range(m)]
            for j in range(i + 1, num_samples):
                for l in range(m):
                    B_future_obs[l].update(B_obs_per_sample[j][l])

        # Evaluate each cell
        for j in cells_in_sample:
            fe = False

            if use_total:
                for k in range(r):
                    val = A_pred[j, k]
                    if val not in A_obs_per_sample[i][k] and val in A_future_obs[k]:
                        fe = True
                        break

            if not fe and use_alternate:
                for l in range(m):
                    val = B_pred[j, l]
                    if val not in B_obs_per_sample[i][l] and val in B_future_obs[l]:
                        fe = True
                        break

            if fe:
                false_emergence_cells += 1

    return false_emergence_cells / n

# EVALUATION HELPERS
def load_tree(adata, name):
    return nx.node_link_graph(json.loads(adata.uns[name]))

def infer_copy_number_tree(parents, CNAs):
    
    # create empty parents/CNAs for copy number tree
    cn_tree_parents = [-1]
    cn_tree_CNAs = [[]]
    
    # perform BFS on parents, removing any nodes without CNAs 
    # and updating the CNAs/parents for the copy number tree
    node_map = {0:0}
    node_id = 1
    q = [0]
    while len(q) > 0:
        n = q.pop()
        p = parents[n]
        for c in np.flatnonzero(parents == n):
            q.append(c)
            if len(CNAs[c]) > 0:
                cn_tree_parents.append(node_map[n])
                cn_tree_CNAs.append(CNAs[c])
                node_map[c] = node_id
                node_id += 1
            else:
                node_map[c] = node_map[n]
                        
    dot_string = 'digraph G{\n'
    dot_string += 'node [color=dimgray fontsize=24 fontcolor=black fontname=Helvetica penwidth=5];\n'
    for n in range(1, len(cn_tree_parents)):
        dot_string += f'"{cn_tree_parents[n]}" -> "{n}" [color="dimgray" penwidth=4 weight=2];\n'
    for n in range(len(cn_tree_parents)):
        label = f'"{n}" [label=<'
        for CNA in cn_tree_CNAs[n]:
            if CNA[1] > 0:
                label += "Gain"
            elif CNA[1] < 0:
                label += "Loss"
            else:
                label += "CNLOH"
            label += f" {CNA[0]}-{','.join(map(str, CNA[2]))}<br/>"

        if len(CNAs[n]) == 0:
            label += " "
        label += ">];\n"
        dot_string += label
    dot_string += "}"
    return cn_tree_parents, cn_tree_CNAs, dot_string

# SIMULATION HELPER FUNCTIONS
def build_tree(parents):
    children = defaultdict(list)
    for i, p in enumerate(parents):
        if p != -1:
            children[p].append(i)
    return children

def choose_total_clones_per_sample(rng, num_clones, k, d):
    """Determines the number of new clones appearing in each sample. First sample has at least 2 clones (root and cancerous),
    then the remaining samples are initialized with d new clones. Then the remaining clones are assigned with uniform probability
    to each sample."""
    num_clones_to_be_assigned = num_clones - 2 - d*(k-1)
    print(num_clones_to_be_assigned)
    assert num_clones_to_be_assigned >= 0, f"Cannot have at least {d} new clones per sample when there {num_clones} clones and {k} samples"  
    new_clones_per_sample = [2] + [d]*(k-1)
    for _ in range(num_clones_to_be_assigned):
        sample = rng.choice(k)
        new_clones_per_sample[sample] += 1
    total_clones_per_sample = [sum(new_clones_per_sample[:i+1]) for i in range(k)]
    return total_clones_per_sample
    
def random_traversal(rng, parents, k, d=1, p=0.7):
    children = build_tree(parents)
    root = parents.index(-1)
    num_clones = len(parents)

    visited_sets = []
    visited = set()
    queue = deque([root])

    total_clones_per_sample = choose_total_clones_per_sample(rng, num_clones, k, d)
    sample = 0
    
    while queue and len(visited_sets) < (k - 1):  # Collect k-1 snapshots
        node = queue.popleft()
        visited.add(node)
        for child in children.get(node, []):
            if rng.random() < p:
                queue.append(child)
            else:
                queue.appendleft(child)
        
        current_snapshot = set(visited)

        if len(current_snapshot) == total_clones_per_sample[sample]:
            visited_sets.append(sorted(list(current_snapshot.copy())))
            sample += 1
    
    # Final snapshot with full tree
    full_tree = list(range(len(parents)))
    visited_sets.append(full_tree)
    
    return visited_sets


# MORAN PROCESS
def moran_process(subtrees, 
                  num_clones, 
                  N=2000,
                  steps_per_epoch=5000,
                  init_clonal_frac=0.01,
                  fitness_std=0.3,
                  drift_mean=0.0,
                  drift_std=0.01,
                  record_every=5,
                  seed=None,
                  fig_path=""):
    
    # initialize 
    population = []
    types_seen = []
    fitness = {}
    history = []
    time_points = []
    sample_freqs = []
    sample_markers = []

    rng = np.random.default_rng(seed)

    # Moran step
    def moran_step(population, fitness, rng):
        unique_types, counts = np.unique(population, return_counts=True)
        freqs = dict(zip(unique_types, counts))
        weights = np.array([fitness[t] * freqs[t] for t in unique_types])
        probs = weights / weights.sum()
    
        reproducer_type = rng.choice(unique_types, p=probs)
        death_index = rng.integers(0, len(population))
        population[death_index] = reproducer_type
    
        return population
    
    # birth, death, selection
    for epoch, new_clones in enumerate(subtrees):
        if len(set(new_clones)) != len(new_clones):
            raise ValueError("subtree contains duplicate nodes")
        elif epoch == 0:
            population = rng.choice(new_clones, size=N)
        else:
            for clone in new_clones:
                if clone not in types_seen:
                    num_to_assign = int(N * init_clonal_frac)
                    replace_idxs = rng.choice(N, size=num_to_assign, replace=False)
                    for idx in replace_idxs:
                        population[idx] = clone
    
        # Add new clones with random initial fitness
        if len(fitness) > 0:
            max_fitness = np.max(list(fitness.values()))
        else:
            max_fitness = 1.0
        for clone in new_clones:
            if clone not in types_seen:
                types_seen.append(clone)
                fitness[clone] = rng.lognormal(mean=max_fitness, sigma=fitness_std)
    
        # Run Moran process for this epoch
        for step in range(steps_per_epoch):
            population = moran_step(population, fitness, rng)
    
            # Slowly drift clone fitness values
            for clone in fitness:
                fitness[clone] *= rng.lognormal(mean=drift_mean, sigma=drift_std)
                fitness[clone] = max(0.1, fitness[clone])
                fitness[clone] = min(10.0, fitness[clone])
    
            # Record population composition
            freqs = [np.sum(np.array(population) == t) / N for t in types_seen]
            freqs += [0.0] * (num_clones - len(freqs))
            if step % record_every == 0 or step == steps_per_epoch - 1:
                history.append(freqs)
                time_point = epoch * steps_per_epoch + step
                time_points.append(time_point)
    
        # Record final timepoint of this epoch as sample time
        sample_freqs.append(freqs)
        sample_markers.append((epoch + 1) * steps_per_epoch)

    # plot
    cmap=plt.get_cmap('tab10')
    colors = [cmap(i) for i in range(1,num_clones+1)]
    history = np.array(history)
    plt.figure(figsize=(10, 5))
    plt.stackplot(time_points, history.T, colors=colors, labels=types_seen, alpha=0.85)
    
    for i, marker in enumerate(sample_markers):
        plt.axvline(marker, linestyle="--", color="gray", linewidth=1)
        plt.text(marker + 10, 1.02, f"S{i+1}", rotation=0, verticalalignment='bottom', color='gray')
    
    plt.title(f"Moran Process Clonal Frequencies (N={N}, epochs={steps_per_epoch})")
    plt.xlabel("Time")
    plt.ylabel("Frequency")
    plt.ylim(0, 1.1)
    plt.legend(loc='upper left', bbox_to_anchor=(1.01, 1.0))
    plt.tight_layout()
    plt.savefig(fig_path, bbox_inches='tight')

    return sample_freqs

# SIMULATE FUNCTION
import json 

def load_tree(adata, name):
    return nx.node_link_graph(json.loads(adata.uns[name]))
    
# simulator
def simulate(num_mutations,
             num_cells, 
             num_CNAs,
             num_nodes,
             num_samples=1,
             num_regions=20,
             fp_rate=0.02,
             fn_alpha=1.0,
             fn_beta=19.0,
             allow_empty_regions=False,
             allow_CNLOH=True, 
             allow_CNA_clones=True,
             max_cn=3,
             concentration=1,
             mean_coverage=20,
             region_concentration=1.5,
             doublet_rate=0.0,
             frac_cells_per_sample=[],
             n_clones_different=1,
             correlation=-1,
             theta=10000,
             min_perc_normal=0.05,
             pop_size=2000,
             steps_per_epoch=5000,
             record_every=5,
             init_clone_frac=0.01,
             fitness_std=0.03,
             drift_mean=0.0,
             drift_std=0.01,
             rho_alpha=9.0,
             rho_beta=1.0,
             save_path="",
             seed=0,
             plot_tree=False):
    """Simulate Tapestri data for longitudinal cancer phylogenies.

    This simulator is a modification of the code found here: 
        https://github.com/cbg-ethz/COMPASS/blob/master/Experiments/simulations/generate_synthetic_data.py

    Parameters
    -----------
    num_mutations: int
        The number of SNVs/SNPs in the data.
    num_cells: int
        The number of cells in total across all samples.
    num_CNAs: int
        The number of CNAs that occurred during evolution.
    num_nodes: int
        The number of nodes (clones) that appear during evolution.
    num_samples: int
        The number of samples taken during evolution. Default is 1.
    num_regions: int
        The number of target regions sequenced. Default is 20.
    fp_rate: float
        The false positive rate during sequencing (sequencing error). Default is 0.02.
    fn_alpha: float
        The alpha parameter for the Beta distribution used to sample dropout rates. Default is 1.
    fn_beta: float
        The beta parameter for the Beta distribution used to sample dropout rates. Default is 20.
    allow_empty_regions: bool
        Flag to enable a target region to not contain an SNV/SNP. Default is False.
    allow_CNA_clones: bool
        Flag to allow clones to be solely distinguished by CNAs. Default is True.
    concentration: int
        Concentration parameter for Dirichlet distribution used to draw node probabilities. Default is 1.
    mean_coverage: int
        The average number of reads covering each locus. Default is 20.
    rho_mean: float
        The mean for the Normal distribution which region probabilities are sample from. Default is 1.0.
    rho_std: float
        Standard deviation for the Normal distribution from which region probabilities are sample from. Default is 5.0
    doublet_rate: float
        The fraction of cells that are doublet during sequencing. Default is 0.0.
    frac_cells_per_sample: list
        A list fractions specifying what percentage of num_cells should be from each sample. The list should sum to 1 and have as many 
        entries as there are samples. Default is [], which uniformly distributes cells among the samples.
    n_clones_different: int
        The minimum number of new clones to appear in each subsequent longitudinal sample. Default is 1.
    theta: int
        Shape parameter for Gamma distribution used to introduce correlation between regions. Default is 10000.
    min_perc_normal: float
        The minimum percentage of each sample that must be "normal" or non-cancerous. Default is 0.05.
    max_clonal_sweeps: int
        The maximum number of clonal sweeps to occur throughout evolution. Default is 0.
    clonal_sweep_probability:
        The probably that for any sample (after the first one) a clonal sweep occurs. Default is 0.9.
    clonal_sweep_concentration: float
        The concentration parameter used to sample node probabilities from the Dirichlet if a clonal sweep is occuring.
        Default is 0.5. Smaller values will result in sparser population frequencies of clones in a sample where a sweep occurs.
    seed: int
        Seed value used to reproduce results

    Returns
    --------
    AnnData
        An annotated DataFrame object with all of the simulated data

    """

    assert (len(frac_cells_per_sample) == 0) or len(frac_cells_per_sample) == num_samples, \
        "The number of longitudinal samples must match the fraction of cells per sample"
    assert (len(frac_cells_per_sample) == 0) or (sum(frac_cells_per_sample) == 1.0), "Fraction of cells per sample must sum to 1.0"
    assert max_cn >= 2, "max_cn must be at least 2 for there to be any CNAs"

    if not os.path.exists(save_path) and save_path != "":
        os.mkdir(save_path)
        
    all_muts_present = False
    rng = np.random.default_rng(seed)

    # all_muts_present is a flag that makes sure that each mutation occurs in at least one clone
    while not all_muts_present:

        ###########################################################
        # 1) Generate random tree structure using prufer sequence #
        ###########################################################
        prufer_code = rng.choice(num_nodes,size=num_nodes-2,replace=True)
        parents = [-1] * num_nodes
        nodes_to_add = set(range(num_nodes-1))
        for i in range(len(prufer_code)):
            n = 0
            while (not n in nodes_to_add) or (n in prufer_code[i:]): n+=1
            parents[n] = prufer_code[i]
            nodes_to_add.remove(n)

        # root is num_nodes-1 from prufer, change it so root is 0 and those 
        # parented by 0 are now parents by num_nodes-1
        parents[list(nodes_to_add)[0]] = num_nodes-1
        parents[0],parents[num_nodes-1] = parents[num_nodes-1],parents[0]
        for i in range(len(parents)):
            if parents[i] == 0: 
                parents[i] = num_nodes-1
            elif parents[i] == num_nodes-1: 
                parents[i] = 0

        # collect children
        children = [[] for x in range(num_nodes)]
        for i in range(1,num_nodes): children[parents[i]].append(i)
        DFT_order=[]
        stack = [0]
        while stack!=[]:
            top = stack.pop()
            DFT_order.append(top)
            for child in children[top]: stack.append(child)

        ######################################
        # 2) Randomly assign SNVs to regions #
        ######################################
        if allow_empty_regions:
            SNV_to_region = sorted(list(rng.integers(0,num_regions,num_mutations)))
        else: # Force each region to contain at least one SNV.
            SNV_to_region = sorted([x for x in range(num_regions)] + list(rng.integers(0,num_regions,num_mutations-num_regions)))
            SNV_to_region = sorted(SNV_to_region)
        region_to_SNVs = [ [] for i in range(num_regions)]
        for i in range(num_mutations):
            region_to_SNVs[SNV_to_region[i]].append(i)

        #########################################
        # 3) Randomly assign SNVs/CNAs to nodes #
        #########################################
        SNVs = [[] for x in range(num_nodes)]
        CNAs = [[] for x in range(num_nodes)] # each CNA is a triplet (region,type,alleles)

        # Assign SNVs/CNAs to nodes (allow_CNA_clones == True permits clones to be solely defined by CNAs)
        if allow_CNA_clones:
            nodes = list(range(1,num_nodes)) + list(rng.choice(range(1,num_nodes),num_mutations+num_CNAs-num_nodes+1,replace=True))
            rng.shuffle(nodes)
        else:
            snv_nodes = list(range(1,num_nodes)) + list(rng.choice(range(1,num_nodes),num_mutations-num_nodes+1,replace=True))
            rng.shuffle(snv_nodes)
            cna_nodes = list(rng.choice(range(1,num_nodes), num_CNAs, replace=True))
            nodes = snv_nodes+cna_nodes


        # collect SNVs
        for i in range(num_mutations):
            SNVs[nodes[i]].append(i)
        for n in range(num_nodes):
            SNVs[n] = sorted(SNVs[n])

        # check if mutation is already in an ancestor
        def mut_in_ancestors(SNV,node):
            if node < 0: return False
            elif SNV in SNVs[node]: return True
            elif node<=0: return False
            else: return mut_in_ancestors(SNV,parents[node])

        # select the type of CNAs to add to each node
        regions_with_CNA = sorted(rng.choice(num_regions,num_CNAs,replace=False))
        t_max_values = list(range(1,max_cn-2+1)) 
        for i in range(num_CNAs):
            n = nodes[num_mutations+i]
            region = regions_with_CNA[i]
            node_has_hetSNVs = False # Can only have a CNLOH if there are heterozygous SNVs in the region
            for mut in region_to_SNVs[region]:
                if mut_in_ancestors(mut,n): node_has_hetSNVs= True
            if node_has_hetSNVs and allow_CNLOH: 
                t = rng.choice([-1,0]+t_max_values)
            else:
                t = rng.choice([-1]+t_max_values)
            alleles = []
            # CNAs can only affect alleles present in this node. In case of CNLOH, we cannot lose the ALT allele in the same node as the SNV. 
            # Also, a SNV followed by a CNLOH in the next node can also be explained by two parallel branches without CNLOH, so we exclude this case
            for snv in region_to_SNVs[region]:
                if (mut_in_ancestors(snv,n) and t!=0) or (mut_in_ancestors(snv,parents[n]) and (not snv in SNVs[parents[n]] or len(children[parents[n]])==1)): 
                    alleles.append(rng.integers(0,2))
                else: 
                    alleles.append(0)

            CNAs[n].append((region, t, alleles))


        # Compute the ground truth copy number profiles of each node
        cnum_regions = 2*np.ones((num_regions,num_nodes),dtype=int)
        n_ref_allele = 2*np.ones((num_mutations,num_nodes),dtype=int)
        n_alt_allele =  np.zeros((num_mutations,num_nodes),dtype=int)

        for n in DFT_order:
            if n!=0:
                cnum_regions[:,n] = cnum_regions[:,parents[n]]
                n_ref_allele[:,n] = n_ref_allele[:,parents[n]]
                n_alt_allele[:,n] = n_alt_allele[:,parents[n]]
            for SNV in SNVs[n]:
                if n_ref_allele[SNV,n]>0:
                    n_ref_allele[SNV,n]-=1
                    n_alt_allele[SNV,n]+=1
            for CNA in CNAs[n]:
                region,t,alleles = CNA
                cnum_regions[region,n] += t
                for i in range(len(alleles)):
                    if t!=0: # gain or loss
                        if alleles[i]==0:
                            n_ref_allele[region_to_SNVs[region][i],n]+=t
                        else:
                            if n_alt_allele[region_to_SNVs[region][i],n]==0:
                                raise Exception("Invalid CNA: affects the alt allele, but its copy number before CNA is 0")
                            n_alt_allele[region_to_SNVs[region][i],n]+=t
                    else: #CNLOH
                        if alleles[i]==0:
                            #if lose the ref allele but do not have the alt allele, do nothing
                            if n_alt_allele[region_to_SNVs[region][i],n]==0:
                                pass
                                #raise Exception("Invalid CNLOH: affects the ref allele, but is homozygous ref")
                            else:
                                n_ref_allele[region_to_SNVs[region][i],n]-=1
                                n_alt_allele[region_to_SNVs[region][i],n]+=1
                        else:
                            if n_alt_allele[region_to_SNVs[region][i],n]==0:
                                pass
                                #raise Exception("Invalid CNA: affects the alt allele, but its copy number before CNA is 0")
                            else:
                                n_ref_allele[region_to_SNVs[region][i],n]+=1
                                n_alt_allele[region_to_SNVs[region][i],n]-=1


                if len(alleles)>0 and (n_ref_allele[region_to_SNVs[region][0],n]<0 or n_alt_allele[region_to_SNVs[region][0],n]<0):
                    raise Exception("Invalid CNA: copy number of one allele < 0.")
                if len(alleles)>0 and (n_ref_allele[region_to_SNVs[region][0],n]+n_alt_allele[region_to_SNVs[region][0],n]!=cnum_regions[region,n]):
                    raise Exception("Copy number of region does not match copy number of locus.")
        
        if np.all(n_alt_allele.sum(axis=1) >= 1): 
            all_muts_present = True
        else:
            seed += 1

    # cell data from sequencing experiment
    depths = np.zeros((num_regions,num_cells),dtype=int)
    ref_reads = np.zeros((num_mutations,num_cells),dtype=int)
    alt_reads = np.zeros((num_mutations,num_cells),dtype=int)
    genotypes = np.zeros((num_mutations,num_cells),dtype=int)

    # ground truth cell data
    true_ref_reads = np.zeros((num_mutations,num_cells),dtype=int)
    true_alt_reads = np.zeros((num_mutations,num_cells),dtype=int)

    # calculate number of cells in each sample
    cell_samples = []
    if len(frac_cells_per_sample) == 0:
        num_cells_in_samples = [int((1.0/num_samples)*num_cells)]*num_samples
    else:
        num_cells_in_samples = [int(num_cells * frac) for frac in frac_cells_per_sample]

    leftover = num_cells - sum(num_cells_in_samples)
    num_cells_in_samples[0] += leftover

    ################################################################
    # 4) Randomized traversal to generate subtrees for each sample #
    ################################################################
    subtrees = random_traversal(rng, parents, num_samples, n_clones_different)
    labeled_subtrees = [[f"Clone {_id}" for _id in subtree if _id != 0] for subtree in subtrees]
    variant_sample = np.array([-1]*num_mutations)
    for s,clones in enumerate(subtrees):
        for c in clones:
            for SNV in SNVs[c]:
                if variant_sample[SNV] == -1:
                    variant_sample[SNV] = s

    fig_path = ""
    if os.path.exists(save_path):
        fig_path = os.path.join(save_path, "fig.pdf")

    cancer_clone_freqs = moran_process(labeled_subtrees, 
                                       num_nodes-1, 
                                       N=pop_size,
                                       steps_per_epoch=steps_per_epoch,
                                       init_clonal_frac=init_clone_frac,
                                       fitness_std=fitness_std,
                                       drift_mean=drift_mean,
                                       drift_std=drift_std,
                                       record_every=record_every,
                                       seed=seed,
                                       fig_path=fig_path)
        
    ###################################################
    # 5) Main loop to generate data for each sample s #
    ###################################################
    cell_attachments = []
    sample_region_probabilities = []
    sample_dropout_rates = []
    sweep_count = 0
    for s in range(num_samples):

        # data for this sample
        clones = subtrees[s]
        num_clones = len(clones)
        SNVs_in_sample = []
        for c in clones:
            SNVs_in_sample += SNVs[c]

        ################################
        # 6) Sample node probabilities #
        ################################
        freqs = cancer_clone_freqs[s] 
        rho = np.clip(rng.beta(rho_alpha, rho_beta), 0.2, 0.95, dtype=float)  # Tumor purity
        node_probabilities = [1.0-rho] + [rho*f for f in freqs]
        node_probabilities = node_probabilities[:num_clones]

        ###########################
        # 7) Sample dropout rates #
        ###########################
        dropout_rates = np.abs(rng.beta(fn_alpha,fn_beta,num_mutations))
        dropout_rates = np.minimum(0.2, np.maximum(1e-3, dropout_rates))
        for i in range(len(dropout_rates)):
            if i in SNVs_in_sample:
                continue
            else:
                dropout_rates[i] = 0.0
        sample_dropout_rates.append(dropout_rates)

        ##################################
        # 8) Sample region probabilities #
        ##################################
        region_probabilities = rng.dirichlet([region_concentration]*num_regions)
        region_probabilities = [max(1.0/(num_regions*5),x) for x in region_probabilities]
        sample_region_probabilities.append(region_probabilities)
        

        ############################################################
        # 9/10) Obtain number of cells in sample and generate data #
        ############################################################
        cell_count = num_cells_in_samples[s] 
        start_index = sum(num_cells_in_samples[:s])
        for j in range(start_index, start_index+cell_count):
            cell_samples.append(s)

            ###############################
            # 10a) Sample node attachment #
            ###############################
            if rng.random() < doublet_rate:
                # doublet code matches COMPASS simulator
                node1 = clones[int(rng.choice(num_clones, p=node_probabilities))]
                node2 = clones[int(rng.choice(num_clones, p=node_probabilities))]
                cell_attachments.append(node1)
                cnum_regions_cell = (np.copy(cnum_regions[:,node1]) + np.copy(cnum_regions[:,node2]) ) / 2
                n_ref_allele_cell = np.copy(n_ref_allele[:,node1]) + np.copy(n_ref_allele[:,node2])
                n_alt_allele_cell = np.copy(n_alt_allele[:,node1]) + np.copy(n_alt_allele[:,node2])
    
            else:
                # singlet code matches COMPASS simulator
                node = clones[int(rng.choice(num_clones, p=node_probabilities))]
                cell_attachments.append(node)
                cnum_regions_cell = np.copy(cnum_regions[:,node])
                n_ref_allele_cell = np.copy(n_ref_allele[:,node]) 
                n_alt_allele_cell = np.copy(n_alt_allele[:,node])
            if correlation<0:
                for k in range(num_regions):
                    expected_depth = mean_coverage * (num_regions * region_probabilities[k]) * cnum_regions_cell[k] / 2.0
                    rate = rng.gamma(shape=theta,scale=expected_depth/theta)
                    depths[k,j] = rng.poisson(rate)
            else:
                means = np.zeros(num_regions)
                for k in range(num_regions):
                    means[k] = mean_coverage * (num_regions * region_probabilities[k]) * cnum_regions_cell[k] / 2.0 
                depths[:,j] = rng.multivariate_normal(means,covariance_regions)

            ###########################################################
            # 10b) Generate read count data for each region/variant i #
            ###########################################################
            for i in range(num_mutations):
                depth = depths[SNV_to_region[i],j]
                c_r=0 # number of copies of the ref allele that did not get dropped out
                c_a=0 # number of copies of the alt allele that did not get dropped out
                for x in range(n_ref_allele_cell[i]):
                    if rng.random()>=dropout_rates[i]: c_r+=1
                for x in range(n_alt_allele_cell[i]):
                    if (rng.random())>= dropout_rates[i]: c_a+=1
                if (c_r==0 and c_a==0):
                    if n_ref_allele_cell[i]>0 and n_alt_allele_cell[i]:
                        if rng.choice(2)==0: c_r+=1
                        else: c_a+=1
                    elif n_ref_allele_cell[i]>0: c_r+=1
                    else: c_a+=1
                f = c_a / (c_r+c_a) * (1-fp_rate) + c_r/(c_r+c_a) * fp_rate
                omega=1000 # high omega --> beta binomial ~= binomial
                a = f*omega
                b = (1-f)*omega
                alt_reads[i,j] = betabinom.rvs(depth,a,b)
                ref_reads[i,j] = depth-alt_reads[i,j]
    
                if ref_reads[i,j] + alt_reads[i,j]<8:
                    genotypes[i,j]=3
                elif min(ref_reads[i,j],alt_reads[i,j])>0.17 * (ref_reads[i,j] + alt_reads[i,j]): 
                    genotypes[i,j]=1
                elif alt_reads[i,j]>0.88 * (alt_reads[i,j] + ref_reads[i,j]):
                    genotypes[i,j]=2
                elif ref_reads[i,j]>0.88 * (alt_reads[i,j] + ref_reads[i,j]):
                    genotypes[i,j]=0
                else:
                    genotypes[i,j]=3
    
    regions = np.array(depths)
    regions_df = pd.DataFrame(regions,index = [str(i)+"_Region"+str(i) for i in range(num_regions)]).astype(int)

    colors = ["lightcoral","skyblue3","sandybrown","paleturquoise3","thistle","darkolivegreen3","lightpink","mediumpurple","darkseagreen3","navajowhite","gold"]

    # Create the DOT string with all labels quoted
    dot_string = 'digraph G{\n'
    dot_string += 'node [color=dimgray fontsize=24 fontcolor=black fontname=Helvetica penwidth=5];\n'
    for n in range(1, num_nodes):
        dot_string += f'"{parents[n]}" -> "{n}" [color="dimgray" penwidth=4 weight=2];\n'
    for n in range(num_nodes):
        label = f'"{n}" [label=<'
        for SNV in SNVs[n]:
            label += f"{SNV}-SNV{SNV}<br/>"
        for CNA in CNAs[n]:
            if CNA[1] == 1:
                label += "Gain"
            elif CNA[1] == 2:
                label += "2-Gains"
            elif CNA[1] == 3:
                label += "3-Gains"
            elif CNA[1] == -1:
                label += "Loss"
            else:
                label += "CNLOH"
            label += f" {CNA[0]}-{','.join(map(str, CNA[2]))}<br/>"

        if len(SNVs[n]) == 0 and len(CNAs[n]) == 0:
            label += " "
        label += ">];\n"
        dot_string += label

    # Convert DOT string to NetworkX DiGraph
    dot_graph = pydot.graph_from_dot_data(dot_string + "}")[0]
    T_mut = nx.nx_pydot.from_pydot(dot_graph)

    # plot the simulated mutation tree
    if plot_tree:
        op.pl.show_tree(T_mut)
    
    for n in range(num_nodes):
        dot_string += f'"{n}" -> "{n + num_nodes}" [dir="none" style="dashed" weight=1 penwidth=5 color="{colors[n % len(colors)]}"];\n'
    for n in range(num_nodes):
        size = np.sqrt(100.0 * node_probabilities[n]) / 3.0
        dot_string += f'"{n + num_nodes}" [label="{100 * node_probabilities[n]:.2f}%" style="filled" width="{size}" height="{size}" color="{colors[n % len(colors)]}"];\n'
    dot_string += '}'

    # Convert DOT string to NetworkX DiGraph
    dot_graph = pydot.graph_from_dot_data(dot_string)[0]
    T_clonal = nx.nx_pydot.from_pydot(dot_graph)
    
    # save copy number tree
    cn_tree_parents, cn_tree_CNAs, cn_tree_dot_string = infer_copy_number_tree(np.array(parents), CNAs)
    cn_tree_dot_graph = pydot.graph_from_dot_data(cn_tree_dot_string)[0]
    T_copy_number = nx.nx_pydot.from_pydot(cn_tree_dot_graph)
    T_copy_number.graph["parents"] = cn_tree_parents
    T_copy_number.graph[f"CNAs"] = json.dumps([[(int(cna[0]), int(cna[1]), tuple(map(int,cna[2]))) for cna in node_CNAs] for node_CNAs in cn_tree_CNAs])

    # create DataFrames that contain all necessary information about data set
    cells = [i for i in range(1,num_cells+1)]
    SNVs = ["SNV%d" % i for i in range(1,num_mutations+1)]
    character_matrix = pd.DataFrame(n_alt_allele.T[cell_attachments], index=cells, columns=SNVs).replace(2,1)
    observed_character_matrix = pd.DataFrame(genotypes.T, index=cells, columns=SNVs).replace(2,1)
    variant_reads_df = pd.DataFrame(true_alt_reads.T, index=cells, columns=SNVs)
    total_reads_df = pd.DataFrame((true_ref_reads + true_alt_reads).T, index=cells, columns=SNVs)
    noisy_variant_reads_df = pd.DataFrame(alt_reads.T, index=cells, columns=SNVs)
    noisy_total_reads_df = pd.DataFrame((ref_reads + alt_reads).T, index=cells, columns=SNVs)
    
    adata = ad.AnnData(pd.DataFrame(genotypes.T, index=cells, columns=SNVs).replace(2,1))
    adata.var_names = ["SNV"+str(snv) for snv in range(0,num_mutations)]
    adata.var["CHR"] = [str(SNV_to_region[snv]) for snv in range(num_mutations)]
    adata.var["REGION"] = ["Region"+str(SNV_to_region[snv]) for snv in range(num_mutations)]
    adata.var["NAME"] = adata.var_names
    adata.var["VARIANT_TYPE"] = ["SNV"] * num_mutations
    adata.var["SAMPLE"] = variant_sample

    # add variables to tree graph
    T_mut.graph["root_name"] = "root"
    T_mut.graph["variants"] = SNVs
    T_mut.graph["cell_assignments"] = cell_attachments

    T_clonal.graph["root_name"] = "root"
    T_clonal.graph["variants"] = SNVs
    T_clonal.graph["cell_assignments"] = cell_attachments

    # collect all trees
    adata.uns[op.ul.DATA.CELL_TREE] = None
    adata.uns[op.ul.DATA.MUTATION_TREE] = str(json.dumps(nx.node_link_data(T_mut)))
    adata.uns[op.ul.DATA.CN_TREE] = str(json.dumps(nx.node_link_data(T_copy_number)))
    adata.uns[op.ul.DATA.CLONAL_TREE] = str(json.dumps(nx.node_link_data(T_clonal)))
    adata.uns[op.ul.DATA.DROPOUT_RATES] = sample_dropout_rates
    
    # collect all cell/mutation data 
    adata.layers[op.ul.DATA.TRUE_DATA] = character_matrix
    adata.layers[op.ul.DATA.OBS_DATA] = observed_character_matrix
    adata.layers[op.ul.DATA.VARIANT_READS] = variant_reads_df
    adata.layers[op.ul.DATA.TOTAL_READS] = total_reads_df
    adata.layers[op.ul.DATA.VARIANT_READS_CORRUPT] = noisy_variant_reads_df
    adata.layers[op.ul.DATA.TOTAL_READS_CORRUPT] = noisy_total_reads_df
    
    # collect cell specific data
    adata.obs[op.ul.DATA.CLUSTER_ID] = cell_attachments
    adata.obs[op.ul.DATA.CELL_SAMPLE] = cell_samples
    
    # compute FPR and FNR and missing rate
    obs_values = observed_character_matrix.values
    true_values = character_matrix.values

    # make sure columns and index are strings for regions_df
    regions_df.index = regions_df.index.map(str)
    regions_df.columns = regions_df.columns.map(str)

    adata.uns[op.ul.DATA.REGION_PROBABILITIES] = sample_region_probabilities
    adata.uns[op.ul.DATA.MUTANT_COPY_NUMBERS] = n_alt_allele.T
    adata.uns[op.ul.DATA.TOTAL_COPY_NUMBERS] = cnum_regions
    adata.uns[op.ul.SIM_KEYS.FPR] = np.maximum(1e-6, np.sum((obs_values == 1) & (true_values == 0)) / np.sum(true_values == 0)) # percentage of 0's flipped to 1's
    adata.uns[op.ul.SIM_KEYS.FNR] = np.maximum(1e-6, np.sum((obs_values == 0) & (true_values == 1)) / np.sum(true_values == 1)) # percentage of 1's flipped to 0's
    adata.uns[op.ul.SIM_KEYS.MISSING_RATE] = np.sum(obs_values == 3) / obs_values.size # percentage of entries flipped to 3
    adata.uns[op.ul.DATA.REGION_READS] = regions_df
    adata.uns["seed"] = seed
    
    if os.path.exists(save_path):
        adata.write_h5ad(os.path.join(save_path, "adata.h5ad"))
            
    return adata

# RUN MULTIPLE SIMULATIONS
def run_simulations(n_trials,
                    num_mutations,
                    num_cells, 
                    num_CNAs,
                    num_nodes,
                    num_samples,
                    num_regions=20,
                    fp_rate=0.02,
                    fn_alpha=1.0,
                    fn_beta=20.0,
                    allow_empty_regions=False,
                    allow_CNLOH=True, 
                    allow_CNA_clones=True,
                    concentration=1,
                    mean_coverage=20,
                    region_concentration=1.5,
                    doublet_rate=0.0,
                    frac_cells_per_sample=[],
                    n_clones_different=1,
                    theta=10000,
                    min_perc_normal=0.05,
                    pop_size=2000,
                    steps_per_epoch=5000,
                    record_every=5,
                    init_clone_frac=0.01,
                    fitness_std=0.03,
                    drift_mean=0.0,
                    drift_std=0.01,
                    rho_alpha=9.0,
                    rho_beta=1.0,
                    max_cn=3,
                    save_path="",
                    include_LoPhy=True,
                    include_COMPASS=True,
                    include_SCITE=True,
                    include_LACE=True,
                    LoPhy_kwargs={"hom_precision":15.0, "het_precision":4.0, "fp":0.02, "fn":0.05},
                    start_from=0,
                    seed=0,
                    plot_tree=False):

    # make sure path exists for saving results
    if save_path == "":
        raise ValueError("save_path must be a valid directory")
    if not os.path.exists(save_path):
        os.mkdir(save_path)
    
    # set up random number generator and seeds
    rng = np.random.default_rng(seed)
    seeds = rng.integers(low=0, high=2**32, size=n_trials)

    # collect results in dictionary
    results = {
        op.ul.EVAL_KEYS.MODEL:[],
        op.ul.EVAL_KEYS.RUNTIME:[],
        MCN_MAE:[],
        TCN_MAE:[],
        FER: [],
        SEED: []
    }

    # run simulation trials
    for i in range(n_trials):
        if i < start_from:
            continue
        print("Sim %d (seed %d)" % (i, seeds[i]))
        trial_folder = f"n{num_cells}_m{num_mutations}_CNA{num_CNAs}_clones{num_nodes}_r{num_regions}_s{num_samples}_trial{i}"
        adata = simulate(num_mutations=num_mutations,
                         num_cells=num_cells,
                         num_CNAs=num_CNAs,
                         num_nodes=num_nodes,
                         num_regions=num_regions,
                         num_samples=num_samples,
                         fp_rate=fp_rate,
                         fn_alpha=fn_alpha,
                         fn_beta=fn_beta,
                         allow_empty_regions=allow_empty_regions,
                         allow_CNLOH=allow_CNLOH, 
                         allow_CNA_clones=allow_CNA_clones,
                         concentration=concentration,
                         mean_coverage=mean_coverage,
                         region_concentration=region_concentration,
                         doublet_rate=doublet_rate,
                         frac_cells_per_sample=frac_cells_per_sample,
                         n_clones_different=n_clones_different,
                         theta=theta,
                         min_perc_normal=min_perc_normal,
                         pop_size=pop_size,
                         steps_per_epoch=steps_per_epoch,
                         record_every=record_every,
                         init_clone_frac=init_clone_frac,
                         fitness_std=fitness_std,
                         drift_mean=drift_mean,
                         drift_std=drift_std,
                         rho_alpha=rho_alpha,
                         rho_beta=rho_beta,
                         max_cn=max_cn,
                         save_path=os.path.join(save_path, trial_folder),
                         seed=seeds[i],
                         plot_tree=plot_tree)

        cells = adata.obs.index
        mutations = adata.var.index
        character_matrix = pd.DataFrame(adata.X, index=cells, columns=mutations)
        variant_reads_df = pd.DataFrame(adata.layers[op.ul.DATA.VARIANT_READS_CORRUPT], index=cells, columns=mutations)
        total_reads_df = pd.DataFrame(adata.layers[op.ul.DATA.TOTAL_READS_CORRUPT], index=cells, columns=mutations)
        regions_df = adata.uns[op.ul.DATA.REGION_READS]
        regions_df.columns = regions_df.columns.astype(int)
        meta_df = adata.var.set_index("CHR").drop(columns="VARIANT_TYPE")
        cell_samples = adata.obs[op.ul.DATA.CELL_SAMPLE]
    
        true_genotypes = pd.DataFrame(adata.layers[op.ul.DATA.TRUE_DATA], index=cells, columns=mutations)
        B_true = pd.DataFrame(adata.uns[op.ul.DATA.MUTANT_COPY_NUMBERS][adata.obs[op.ul.DATA.CLUSTER_ID]], index=cells)
        A_true = pd.DataFrame(adata.uns[op.ul.DATA.TOTAL_COPY_NUMBERS][:,adata.obs[op.ul.DATA.CLUSTER_ID]].T, index=cells)
        
        if include_LoPhy:
            print("Running LoPhy")
            
            # run LoPhy
            LoPhy_folder = os.path.join(save_path, trial_folder, "LoPhy")
        
            sol_LoPhy = op.tl.solver.LoPhy(character_matrix.replace(2,1).replace(3,-1), 
                                           variant_reads_df, 
                                           total_reads_df,                        
                                           regions_df, 
                                           meta_df,
                                           cell_samples=cell_samples,
                                           destination_dir=LoPhy_folder,
                                           **LoPhy_kwargs)

            # process LoPhy's results
            T_LoPhy = sol_LoPhy[op.ul.DATA.CELL_TREE]
            cell_assignments_LoPhy = np.array(T_LoPhy.graph["cell_assignments"], dtype=int)
            alt_cns = np.array(T_LoPhy.graph[op.ul.DATA.MUTANT_COPY_NUMBERS], dtype=int)
            total_cns = np.array(T_LoPhy.graph[op.ul.DATA.TOTAL_COPY_NUMBERS], dtype=int)
            B_pred_LoPhy = pd.DataFrame(alt_cns[cell_assignments_LoPhy], index=cells)
            A_pred_LoPhy = pd.DataFrame(total_cns[cell_assignments_LoPhy], index=cells)
            predicted_genotypes_LoPhy = sol_LoPhy[op.ul.DATA.PRED_DATA].astype(int)
        
            results[op.ul.EVAL_KEYS.MODEL].append(LOPHY)
            results[op.ul.EVAL_KEYS.RUNTIME].append(sol_LoPhy[op.ul.EVAL_KEYS.RUNTIME])
            results[MCN_MAE].append(mae(B_true.values, B_pred_LoPhy.values))
            results[TCN_MAE].append(mae(A_true.values, A_pred_LoPhy.values))
            results[FER].append(fer(A_true=A_true.values, B_true=B_true.values, A_pred=A_pred_LoPhy.values, B_pred=B_pred_LoPhy.values, cell_samples=cell_samples))
            results[SEED].append(seeds[i])
        
        if include_COMPASS:
            print("Running COMPASS")
            # run COMPASS
            COMPASS_folder = os.path.join(save_path, trial_folder, "COMPASS")
            
            sol_COMPASS = op.tl.solver.COMPASS(character_matrix, 
                                               variant_reads_df, 
                                               total_reads_df, 
                                               regions_df, 
                                               meta_df, 
                                               remove_temp_dir=True,
                                               destination_dir=COMPASS_folder)
            
            # process COMPASS's results
            T_compass = sol_COMPASS[op.ul.DATA.CELL_TREE]
            cell_assignments_COMPASS = np.array(T_compass.graph["cell_assignments"], dtype=int)
            A_pred_COMPASS = T_compass.graph[op.ul.DATA.TOTAL_COPY_NUMBERS][cell_assignments_COMPASS]
            B_pred_COMPASS = sol_COMPASS[op.ul.DATA.PRED_DATA].astype(int)
        
            results[op.ul.EVAL_KEYS.MODEL].append(COMPASS)
            results[op.ul.EVAL_KEYS.RUNTIME].append(sol_COMPASS[op.ul.EVAL_KEYS.RUNTIME])
            results[MCN_MAE].append(mae(B_true.values, B_pred_COMPASS.values))
            results[TCN_MAE].append(mae(A_true.values, A_pred_COMPASS))
            results[FER].append(fer(A_pred=A_pred_COMPASS, B_pred=B_pred_COMPASS.values, A_true=A_true.values, B_true=B_true.values, cell_samples=cell_samples))
            results[SEED].append(seeds[i])
    
        if include_SCITE:
            print("Running SCITE")
            SCITE_folder = os.path.join(save_path, trial_folder, "SCITE")
            sol_SCITE = op.tl.solver.SCITE(character_matrix.replace(-1,3),
                                           chain_length=90000,
                                           repetitions=3,
                                           fp=fp_rate,
                                           fn1=0.05,
                                           fn2=0.05,
                                           cc=0.01,
                                           destination_dir=SCITE_folder)
            
            B_pred_SCITE = sol_SCITE[op.ul.DATA.PRED_DATA].astype(int)

            true_genotypes_binary = true_genotypes.replace(2,1).values
            results[op.ul.EVAL_KEYS.MODEL].append(SCITE)
            results[op.ul.EVAL_KEYS.RUNTIME].append(sol_SCITE[op.ul.EVAL_KEYS.RUNTIME])
            results[MCN_MAE].append(mae(true_genotypes_binary, B_pred_SCITE.values))
            results[TCN_MAE].append(np.nan)
            results[FER].append(fer(B_pred=B_pred_SCITE.values, B_true=true_genotypes_binary, cell_samples=cell_samples))
            results[SEED].append(seeds[i])
    
        if include_LACE:
            print("Running LACE")
            LACE_folder = os.path.join(save_path, trial_folder, "LACE")
            sol_LACE = op.tl.solver.LACE(character_matrix,
                                         cell_samples,
                                         fp=fp_rate,
                                         fn=0.05,
                                         destination_dir=LACE_folder)
            
            B_pred_LACE = sol_LACE[op.ul.DATA.PRED_DATA].astype(int)

            true_genotypes_binary = true_genotypes.replace(2,1).values
            results[op.ul.EVAL_KEYS.MODEL].append(LACE)
            results[op.ul.EVAL_KEYS.RUNTIME].append(sol_LACE[op.ul.EVAL_KEYS.RUNTIME])
            results[MCN_MAE].append(mae(true_genotypes_binary, B_pred_LACE.values))
            results[TCN_MAE].append(np.nan)
            results[FER].append(fer(B_pred=B_pred_LACE.values, B_true=true_genotypes_binary, cell_samples=cell_samples))
            results[SEED].append(seeds[i])


    return pd.DataFrame(results)

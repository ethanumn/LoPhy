# imports
import random
from collections import defaultdict, deque

import os, json
import numpy as np
import pandas as pd
from scipy.stats import nbinom, betabinom
import pydot
import anndata as ad
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx
from networkx.exception import NetworkXError

import oncophylo as op


def load_tree(adata, name):
    return nx.node_link_graph(json.loads(adata.uns[name]))

def add_population_nodes(tree, regions_df, meta_df, cell_samples, use_labels=False):
    colors = ["lightcoral", "#6CA6CD", "sandybrown", "paleturquoise3", "thistle",
              "#A2CD5A", "lightpink", "mediumpurple", "darkseagreen3",
              "navajowhite", "gold", "#c3b091", "lightsteelblue", "rosybrown",
              "mediumaquamarine", "palevioletred", "burlywood", "orchid3",
              "cadetblue3", "lightsalmon3"]

    if use_labels:
        root_name = tree.graph["root_name"]
        variants = tree.graph["variants"]
        nodes = [root_name] + variants
    else:
        nodes = sorted(list(tree.nodes))
    cell_assignments = np.array(tree.graph["cell_assignments"], dtype=int)

    samples = np.unique(cell_samples)
    clones = np.unique(cell_assignments)
    node_to_int = {n: i for i, n in enumerate(clones)}

    # Track population nodes per clone and sample
    pop_nodes_by_clone = {n: [] for n in clones}
    pop_nodes_by_sample = {s: [] for s in samples}

    for s in samples:
        sample_mask = (cell_samples == s)
        sample_cells = np.flatnonzero(sample_mask)
        sample_assignments = cell_assignments[sample_cells]
        total_cells = len(sample_cells)

        for n in clones:
            count = np.sum(sample_assignments == n)
            freq = round(100 * count / total_cells) if total_cells > 0 else 0

            internal_name = f"pop_clone{n}_sample{s}"  # unique internal name
            label = f"{count} cells\n{freq}\\%"
            node_color = colors[node_to_int[n] % len(colors)]
            size = np.sqrt(freq) / 3.0 if freq > 0 else 0.15

            tree.add_node(internal_name,
                          label=label,
                          color=node_color,
                          style="filled",
                          height=size,
                          width=size,
                          fillcolor="white" if count == 0 else node_color)

            pop_nodes_by_clone[n].append((s, internal_name))
            pop_nodes_by_sample[s].append(internal_name)

    # Connect population nodes
    for n, sample_nodes in pop_nodes_by_clone.items():
        sample_nodes.sort()
        clone_node = nodes[nodes.index(n)] if isinstance(n, str) else nodes[n]
        if sample_nodes:
            tree.add_edge(clone_node, sample_nodes[0][1],
                          style="dashed", penwidth=5, color="gray", dir="none")
            for (_, src), (_, tgt) in zip(sample_nodes[:-1], sample_nodes[1:]):
                tree.add_edge(src, tgt,
                              style="dashed", penwidth=3, color="gray", dir="none")

    # Add rank alignment info
    tree.graph["rank_same"] = list(pop_nodes_by_sample.values())
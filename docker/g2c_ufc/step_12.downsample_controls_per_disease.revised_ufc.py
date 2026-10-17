#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Adapted from Ryan L. Collins, Riaz Gillani, Jett Crowdis and the Van Allen Laboratory
# Copyright (c) 2025 Noah Fields and the Van Allen Laboratory
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Noah D. Fields <Noah_Fields@dfci.harvard.edu>

"""
Downsample control samples to ancestry-match individual pediatric cancer subtypes
"""


import argparse
import numpy as np
import pandas as pd
import random
import string
import math
import re
import networkx as nx
from scipy.stats import chisquare, fisher_exact

#### Noah's Version ####

import pandas as pd

def select_unrelated_samples(kinship_file, df):
    """
    Returns a set of person_id values to KEEP.
    """

    # ---- Normalize IDs in meta ----
    df = df.copy()
    df['person_id'] = df['person_id'].apply(normalize_id)

    df = df.set_index('person_id')

    # ---- Load kinship (only first two columns matter) ----
    kin = pd.read_csv(
        kinship_file,
        sep="\t",
        usecols=[0, 1],
        names=["ID1", "ID2"],
        header=0
    )

    kin['ID1'] = kin['ID1'].apply(normalize_id)
    kin['ID2'] = kin['ID2'].apply(normalize_id)

    related_ids = set(kin['ID1']).union(set(kin['ID2']))

    # ---- Start with everyone who is NOT related ----
    keep = set(df.index) - related_ids

    # ---- Function to check exact cancer match ----
    def has_preferred_cancer(dx):
        if pd.isna(dx):
            return False
        return preferred_cancer_type in dx.split(';')

    # ---- Resolve each related pair ----
    for _, row in kin.iterrows():
        id1, id2 = row['ID1'], row['ID2']

        if id1 not in df.index and id2 not in df.index:
            continue
        elif id1 in df.index and id2 not in df.index:
            keep.add(id1)
            continue
        elif id1 not in df.index and id2 in df.index:
            keep.add(id2)
            continue

        s1 = df.loc[id1]
        s2 = df.loc[id2]

        # Age-based decision
        if s1['age'] <= s2['age']:
            keep.add(id1)
        else:
            keep.add(id2)

    return keep

def normalize_id(x):
    try:
        return str(int(float(x)))
    except (ValueError, TypeError):
        return str(x)


def extract_non_familial_set(samples, kinship_file):
    """
    Given a set of sample IDs and a kinship file (with columns #ID1 and ID2),
    return the subset of sample IDs that do NOT appear in either column.
    """
    kinship = pd.read_csv(kinship_file, delim_whitespace=True)

    # Collect all IDs that appear in the kinship pairs
    related_ids = set(kinship['#ID1']).union(set(kinship['ID2']))

    # Return the samples that are not related to anyone
    return samples - related_ids
    
def extract_family_units(kinship_file):
    """
    Given a kinship file with columns #ID1 and ID2,
    returns a list of sets, where each set represents a family unit
    (i.e., a connected component in the kinship graph).
    """
    kinship = pd.read_csv(kinship_file, delim_whitespace=True)
    
    G = nx.Graph()
    
    # Add edges from kinship file
    for _, row in kinship.iterrows():
        G.add_edge(str(int(row['#ID1'])), str(int(row['ID2'])))

    # Get connected components (family units)
    family_units = [set(component) for component in nx.connected_components(G)]
    
    return family_units


def dfs_independent_set(G, node_list, selected, best, cancer_status):
    # Base Case if the node list is empty
    if not node_list:
        selected_cases = {n for n in selected if cancer_status.get(n, 'control') != 'control'}
        best_cases = {n for n in best[0] if cancer_status.get(n, 'control') != 'control'}

        # Prefer more cancer cases; break ties with total size
        if (len(selected_cases) > len(best_cases)) or \
           (len(selected_cases) == len(best_cases) and len(selected) > len(best[0])):
            best[0] = selected.copy()
        return

    # Recursive Case where node list is not empty
    node = node_list[0]

    # Option 1: include node
    if all(neigh not in selected for neigh in G.neighbors(node)):
        selected.add(node)
        remaining = [n for n in node_list[1:] if n not in G.neighbors(node)]
        dfs_independent_set(G, remaining, selected, best, cancer_status)
        selected.remove(node)

    # Option 2: exclude node
    dfs_independent_set(G, node_list[1:], selected, best, cancer_status)



def maximal_non_related_subset_dfs(family, kinship_file, meta):
    """
    DFS-based search to find a maximal independent set within a family group.
    Only includes individuals present in the metadata.
    """

    # Normalize person_id to clean strings (handles floats like 1234.0)
    meta = meta.copy()
    meta['person_id'] = (
        pd.to_numeric(meta['person_id'], errors='coerce')
          .astype('Int64')
          .astype(str)
    )

    # Normalize family IDs to clean strings
    family = {
        str(int(x)) if isinstance(x, float) and x.is_integer() else str(x)
        for x in family
    }

    # Filter family members to only those present in metadata
    available_ids = set(meta['person_id'])
    family = family & available_ids  # intersection

    # Load kinship
    kinship = pd.read_csv(kinship_file, delim_whitespace=True)
    edges = kinship[
        kinship['#ID1'].isin(family) & kinship['ID2'].isin(family)
    ][['#ID1', 'ID2']].values


    if not family:
        return set()  # nothing to do if no valid individuals

    # Build graph
    G = nx.Graph()
    G.add_nodes_from(family)
    G.add_edges_from(edges)

    # Cancer info
    cancer_status = meta.set_index('person_id').loc[list(family)]['cancer'].to_dict()
    #cancer_status = meta.set_index('original_id').loc[family]['cancer'].to_dict()
    nx.set_node_attributes(G, cancer_status, name='cancer')

    # Sort nodes to prioritize cases first
    node_list = sorted(
        G.nodes, key=lambda n: (cancer_status.get(n, 'control') == 'control', G.degree[n])
    )
    # Sort nodes: preferred cancer type first, other cancers next, then controls
    # node_list = sorted(
    #     G.nodes,
    #     key=lambda n: (
    #         0 if cancer_status.get(n) == preferred_cancer_type else
    #         1 if cancer_status.get(n, 'control') != 'control' else
    #         2,
    #         G.degree[n]  # within each group, sort by degree
    #     )
    # )

    # DFS setup
    selected = set()
    best = [set()]  # use list as mutable container to hold best set
    dfs_independent_set(G, node_list, selected, best, cancer_status)

    return best[0]

def dfs_independent_set_family_style(G, node_list, selected, best, cancer_status, preferred_cancer_type):

    def has_preferred_cancer(node): c = cancer_status.get(node, ''); return isinstance(c, str) and preferred_cancer_type in [x.strip() for x in c.split(';')]

    def score(nodes): return (sum(has_preferred_cancer(n) for n in nodes), sum(cancer_status.get(n, 'control') != 'control' for n in nodes), len(nodes))

    if not node_list:
        if score(selected) > score(best[0]): best[0] = selected.copy()
        return

    node = node_list[0]

    if all(neigh not in selected for neigh in G.neighbors(node)):
        selected.add(node)
        remaining = [n for n in node_list[1:] if n not in G.neighbors(node)]
        dfs_independent_set_family_style(G, remaining, selected, best, cancer_status, preferred_cancer_type)
        selected.remove(node)

    dfs_independent_set_family_style(G, node_list[1:], selected, best, cancer_status, preferred_cancer_type)


def maximal_non_related_subset_dfs_family_style(family, kinship_file, meta, preferred_cancer_type):
    """
    DFS-based search to find a maximal independent set within a family group.
    Only includes individuals present in the metadata.
    """

    # Filter family members to only those present in metadata
    available_ids = set(meta['person_id'])
    family = family & available_ids  # intersection

    # Load kinship
    kinship = pd.read_csv(kinship_file, delim_whitespace=True)
    edges = kinship[
        kinship['#ID1'].isin(family) & kinship['ID2'].isin(family)
    ][['#ID1', 'ID2']].values


    if not family:
        return set()  # nothing to do if no valid individuals

    # Build graph
    G = nx.Graph()
    G.add_nodes_from(family)
    G.add_edges_from(edges)

    # Cancer info
    cancer_status = meta.set_index('person_id').loc[list(family)]['cancer'].to_dict()
    #cancer_status = meta.set_index('original_id').loc[family]['cancer'].to_dict()
    nx.set_node_attributes(G, cancer_status, name='cancer')

    # Sort nodes to prioritize cases first
    # node_list = sorted(
    #     G.nodes, key=lambda n: (cancer_status.get(n, 'control') == 'control', G.degree[n])
    # )
    # Sort nodes: preferred cancer type first, other cancers next, then controls

    def has_preferred_cancer(cancer_str, preferred):
        if not isinstance(cancer_str, str):
            return False
        return preferred in [x.strip() for x in cancer_str.split(';')]

    node_list = sorted(
        G.nodes,
        key=lambda n: (
            0 if has_preferred_cancer(cancer_status.get(n, ''), preferred_cancer_type) else
            1 if cancer_status.get(n, 'control') != 'control' else
            2,
            G.degree[n]  # within each group, sort by degree
        )
    )

    # DFS setup
    selected = set()
    best = [set()]  # use list as mutable container to hold best set
    dfs_independent_set(G, node_list, selected, best, cancer_status,preferred_cancer_type)

    return best[0]

def match_controls_by_ancestry(meta: pd.DataFrame, max_controls_per_case: int = 3, seed: int = 42) -> pd.DataFrame:
    """
    Match controls to cases based on ancestry (ancestry_pred), sex, and cohort if available.
    Try to balance age within each matched group.
    
    Parameters:
        meta: DataFrame with at least ['cancer', 'ancestry_pred', 'sex', 'age', 'cohort']
        max_controls_per_case: Max number of controls to retain per case (default: 3)
        seed: Seed for reproducible sampling

    Returns:
        DataFrame of all cases and matched controls.
    """
    np.random.seed(seed)

    already_matched_control_ids = set()

    # Define control and case status
    is_control = meta['cancer'].str.lower().str.contains('control')
    is_case = ~meta['cancer'].isin(['control', 'unknown'])

    cases = meta[is_case].copy()
    controls = meta[is_control].copy()
    matched_controls = []

    # Iterate over ancestry + sex combinations in cases
    for (pop, sex), case_group in cases.groupby(['ancestry_pred', 'sex_karyotype']):
        case_count = len(case_group)

        # To account for ancestry groups with disproportionately few cases, 
        # we cap the control-to-case ratio at 5:1 when such imbalances are present. 
        if (max_controls_per_case < 5):
            max_controls = case_count * 5
        else:
            max_controls = case_count * max_controls_per_case

        median_age = case_group['age'].median()

        # Find matching controls for ancestry and sex
        potential_controls = controls[
            (controls['ancestry_pred'] == pop) &
            (controls['sex_karyotype'] == sex) &
            (~controls['person_id'].isin(already_matched_control_ids))
        ].copy()

        if potential_controls.empty:
            continue

        # Prioritize controls by age difference
        potential_controls['age_diff'] = (potential_controls['age'] - median_age).abs()
        potential_controls = potential_controls.drop_duplicates(subset='person_id')
        selected = potential_controls.sort_values('age_diff').head(max_controls)

        already_matched_control_ids.update(selected['person_id'])
        matched_controls.append(selected)

    # Combine everything
    matched_controls_df = pd.concat(matched_controls, ignore_index=True)
    final_meta = pd.concat([cases, matched_controls_df], ignore_index=True)

    return final_meta


def get_min_case_control_ratio(df: pd.DataFrame):
    """
    Returns the case:control ratio for each ancestry group and identifies the smallest ratio.
    
    Parameters:
        df: DataFrame with 'cancer' and 'ancestry_pred' columns
    
    Returns:
        A Series with case/control ratios per ancestry, and a printed note about the smallest one.
    """
    # Define cases and controls
    df = df.copy()
    df['is_case'] = df['pancancer_cancer'].eq(1)

    # Group by ancestry
    grouped = df.groupby('ancestry_pred')['is_case']

    # Count cases and controls
    counts = grouped.agg(['sum', 'count'])  # sum = cases, count - sum = controls
    counts['controls'] = counts['count'] - counts['sum']
    counts['control_case_ratio'] =  counts['controls'] / counts['sum']

    # Handle div-by-zero (if any population has 0 controls)
    counts['control_case_ratio'] = counts['control_case_ratio'].replace([float('inf')], float('nan'))

    # Drop intermediate columns for clarity
    ratio_series = counts['control_case_ratio']
    
    # Print info on the ancestry with the lowest ratio
    min_pop = ratio_series.idxmin()
    min_val = ratio_series.min()
    print(f"Smallest control:case ratio is {min_val:.3f} in population: {min_pop}")
    print(f"Ratio Series: {ratio_series}")
    return min_val


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--phenotype-data', required=True, help='data describing phenotypes')
    parser.add_argument('--metadata',required=True)
    parser.add_argument('--sample-list', required=True, help='list of AoU CDRv7 UFC samples (n = 3218)')
    parser.add_argument('--cancer-type', required=True, help='cancer subtype being analyzed')
    parser.add_argument('--pca', required = True, help='.eigenvec file')
    parser.add_argument('--kinship', required = True, help='plink king file represtenting sparse matrix of kinship')
    parser.add_argument('--outfile', required=True, help='output.tsv')
    parser.add_argument('--sex_karyotypes',required=False)
    parser.add_argument('--exclude-samples',required=False, help='Samples to Exclude')
    parser.add_argument('--log-file', required=False, help='Path to log file')
    parser.add_argument('--familial', required=False, type=int, choices=[0, 1, 2], default=0)
    args = parser.parse_args()

    # Make the log file
    if not args.log_file:
        if args.cancer_type:
            args.log_file = f"{args.cancer_type}.cohort.log"
        else:
            args.log_file = "cancer.cohort.log"

    with open(args.log_file, "w") as f:
        f.write("Size\tNum_Filtered\tPercent_Filtered\tExclusion_Criteria\n")
        

    # Load phenotype data
    phenotype_data = pd.read_csv(args.phenotype_data,sep='\t',index_col=False)[['original_id','inferred_sex','intake_qc_pop','sex_karyotype']]
    phenotype_data['person_id'] = phenotype_data['original_id'].astype(str)

    metadata = pd.read_csv(args.metadata, sep='\t',index_col=False)
    metadata['person_id'] = metadata['person_id'].astype(str)
    metadata = metadata.merge(phenotype_data, how = "left", on="person_id")

    # Remove any straggler controls with >= 1 reported related with cancer
    metadata = metadata[~((metadata['pancancer_cancer'].fillna(0) == 0) & (metadata['num_relatives_w_cancer'] >= 1))]

    # Remove any straggler case with <3 reported relatives with cancer
    metadata = metadata[~((metadata['pancancer_cancer'] == 1) & (metadata['num_relatives_w_cancer'] <= 3))]

    # Read in PCAs
    pca = pd.read_csv(args.pca,sep='\t',index_col=False)
    pca['person_id'] = pca['#IID'].astype(str)

    metadata = metadata.merge(pca, how="left",on="person_id")

    # Load list of samples to exclude
    if args.exclude_samples:
        with open(args.exclude_samples) as f3:
            exclude_samples = set(patient.strip() for patient in f3)

        # If these samples do exist, exclude them
        metadata = metadata[~metadata['person_id'].isin(exclude_samples)]

    # Filter to just cases in our study as well as cases in the specific subtype
    metadata = metadata[
        (metadata[f"{args.cancer_type}_cancer"] == 1) | 
        (metadata["pancancer_cancer"].fillna(0) == 0)
        ]
    
    if args.familial == 1:
        metadata = metadata[
            (metadata[f"{args.cancer_type}_cancer"] == 1) & (metadata[f"fdrs_w_{args.cancer_type}_cancer"] >= 1) | 
            (metadata["pancancer_cancer"].fillna(0) == 0)
        ]
    elif args.familial == 2:
        metadata = metadata[
            (metadata[f"{args.cancer_type}_cancer"] == 1) & (metadata[f"fdrs_w_{args.cancer_type}_cancer"] >= 2) | 
            (metadata["pancancer_cancer"].fillna(0) == 0)
        ]

    metadata = metadata[metadata['sex_karyotype'].isin(args.sex_karyotypes.split(","))]
    # with open(args.log_file, "a") as f:
    #     f.write(f"{len(phenotype_data)}\t0\t0\tInitial Samples\n")


    # Filter to Cancer Subtypes of Interest
    # print("Args Variables")

    # --- Filtering logic ---

    # Remove samples with irrelevant cancer diagnosis for this study
    # sample_size4 = len(phenotype_data)

    # with open(args.log_file, "a") as f:
    #     f.write(f"{sample_size4}\t{(sample_size3 - sample_size4)}\t{round(((sample_size3 - sample_size4)/sample_size3),3) * 100}\tRemove samples that are not {args.sex_karyotypes} and not in {args.cohorts} cohort(s).\n")

    ## Grab maximally unrelated set ##

    # Grab all cases not involved in a family
    non_familial_set = extract_non_familial_set(samples=set(metadata['person_id']),kinship_file=args.kinship)

    # Grab family units
    family_units = extract_family_units(kinship_file=args.kinship)

    keep_ids = select_unrelated_samples(args.kinship,metadata)

    metadata = metadata[metadata['person_id'].isin(keep_ids)]

    # Downsample controls to match on cases ancestry
    min_ancestry_case_control_ratio = math.ceil(get_min_case_control_ratio(metadata))
    metadata = match_controls_by_ancestry(metadata, max_controls_per_case = min_ancestry_case_control_ratio)

    metadata['cancer_status'] = metadata['pancancer_cancer'].eq(1)

    ## Print Summary Statistics
    with open(args.log_file, 'a') as f:
        f.write("===== Summary Report =====\n\n")
        
        case_df = metadata[metadata['pancancer_cancer'].fillna(0) > 0]
        control_df = metadata[metadata['pancancer_cancer'].fillna(0) == 0]

        # Total number of individuals
        total_count = len(metadata)
        f.write(f"Total individuals: {total_count}\n\n")

        # Case and Control Counts
        f.write(f"Case Count: {len(case_df)}\n")
        f.write(f"Control Count: {len(control_df)}\n")
        f.write("\n\n")

        # Counts of each sex
        f.write("Sex Stratification:\n")
        f.write(f"Case Count:\n{case_df['sex_at_birth'].value_counts(dropna=False).to_string()}\n")
        f.write(f"Control Count:\n{control_df['sex_at_birth'].value_counts(dropna=False).to_string()}\n")
        f.write("\n\n")

        # Mean age for cases and controls
        f.write("Mean age by cancer type:\n")
        # mean_case_age = round(metadata[metadata['pancancer_cancer'] == 1][f'{args.cancer_type}_age'].mean(), 2)
        # mean_control_age = round(metadata[metadata['pancancer_cancer'].fillna(0) == 0]['age'].mean(), 2)
        # f.write(f"Mean Case Age\t{mean_case_age}\n")
        # f.write(f"Mean Control Age\t{mean_control_age}\n")
        f.write("\n\n")

        # ancestry_pred count per cancer
        f.write("Counts of ancestry_pred by cancer status:\n")
        case_cohort = metadata[metadata['pancancer_cancer'].fillna(0) > 0]
        control_cohort = metadata[metadata['pancancer_cancer'].fillna(0) == 0]
        case_ancestry = case_cohort['ancestry_pred'].value_counts(dropna=False).to_string()
        control_ancestry = control_cohort['ancestry_pred'].value_counts(dropna=False).to_string()
        f.write("Cases:\n")
        f.write(case_ancestry + "\n\n")
        f.write("Controls:\n")
        f.write(control_ancestry + "\n\n")

        f.write("==========================\n\n")




    # Write to outfile
    metadata.to_csv(args.outfile, sep='\t', index=False, na_rep="NA")

if __name__ == '__main__':
    main()

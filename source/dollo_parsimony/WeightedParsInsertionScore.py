#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 16:05:35 2021

@author: claraiglhaut
"""

import numpy as np


def WParsScoreLeaf(leaf, cost_matrix):
    
    '''
    Initializes weighted parsimony scores at the leaves with zero for the 
    observed residue and infinity everywhere else.

    Parameters
    ----------
    leaf : PhyloNode
        leaf nodes of a PhlyoNode.
    
    cost_matrix : double dictionary
        cost matrix specifying the cost of substitutions, insertions and 
        deletions.

    Returns
    -------
    None.

    '''
    
    characters = [key for key in cost_matrix.keys()]
    
    for i in range(len(leaf.sequence)):
        w_site_scores = {}
        for character in characters:
            if leaf.sequence[i] == character:
                w_site_scores[character] = 0
            else:
                w_site_scores[character] = np.inf
            
            leaf.w_parsimony_scores[i] = w_site_scores


def WParsScoreInternal(tree, cost_matrix, i):
    
    '''
    Calculates the weightes parsimony score at site i for the given trees 
    internal nodes.

    Parameters
    ----------
    tree : PhyloNode.
        Internal nodes of a tree structure.
    cost_matrix : double dictionary
        cost matrix specifying the cost of substitutions, insertions and 
        deletions.
    i : int.
        Index for the associated sequence of the node.

    Returns
    -------
    None.

    '''
    
    characters = [key for key in cost_matrix.keys()]
    left_scores = tree.children[0].w_parsimony_scores[i]
    
    right_scores = tree.children[1].w_parsimony_scores[i]
    
    w_site_scores = {}
    
    for character in characters:
        min_score_left = np.inf
        min_score_right = np.inf
        for key in characters:
            if left_scores[key] + cost_matrix[character][key] < min_score_left:
                min_score_left = left_scores[key] + cost_matrix[character][key]
                
            if right_scores[key] + cost_matrix[character][key] < min_score_right:
                min_score_right = right_scores[key] + cost_matrix[character][key]
            
        w_site_scores[character] = min_score_left + min_score_right
    
    tree.w_parsimony_scores[i] = w_site_scores
    print(tree.w_parsimony_scores[i])
        
        
def WeightedParsWithInsertionScore(tree, cost_matrix):
    '''
    Calculates the weighted parsimony score for the whole tree while accounting 
    for insertions and deletions.

    Parameters
    ----------
    tree : PhyloNode
        Input tree with an associated alignment.
    cost_matrix : double dictionary
        cost matrix specifying the cost of substitutions, insertions and 
        deletions.

    Returns
    -------
    tree_score : float
        Weighted parsimony score of the whole tree.

    '''
    
    tree_score = 0
    length_MSA = len(tree.get_leaves()[0].sequence)
    
    #initilaize scores for every node, insertion flags are set to zero
    for node in tree.traverse('postorder'):
        w_pars_scores = [{}]*length_MSA
        ins_flags = [False]*length_MSA
        
        node.add_features(w_parsimony_scores = w_pars_scores)
        node.add_features(insertion_flags = ins_flags)
        
    #scores for leaves
    for leaf in tree.iter_leaves():
        WParsScoreLeaf(leaf, cost_matrix)
        print(leaf.w_parsimony_scores)
        
        
    #find most parsimonious insertion points
    for i in range(length_MSA):
        leaf_res = []
        for leaf in tree.iter_leaves():
            if leaf.sequence[i] != '-':
                leaf_res.append(leaf)
                
        if len(leaf_res) == 1:
            ancestor = leaf_res[0]
        else:
            ancestor = tree.get_common_ancestor(leaf_res)
            
        #calculate internal weightes parsimony scores for nodes with residues   
        for node in ancestor.traverse('postorder'):
            if not node.is_leaf():
                WParsScoreInternal(node, cost_matrix, i)
             
        #set weighted parsimony scores for nodes without residues and 
        #mark them as insertions
        for node in tree.traverse('postorder'):
            if node not in ancestor.get_descendants() and (node != ancestor):
                node.w_parsimony_scores[i] = ancestor.w_parsimony_scores[i]
                node.insertion_flags[i] = True
                print(i, node)
        
        
        min_key = min(tree.w_parsimony_scores[i].keys(), 
                      key=(lambda k: tree.w_parsimony_scores[i][k]))
        tree_score += tree.w_parsimony_scores[i][min_key]
            
        
    return tree_score
        
        




#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 14:06:20 2021

@author: claraiglhaut
"""


characters = ['A', 'T', 'C', 'G', '-']

def ParsInsertionsLeaf(leaf):
    '''
    Initializes the parsimony sets and scores at the leaf nodes.

    Parameters
    ----------
    leaf : PhlyoNode or PhyloTree
        tree leaves with the aligned sequences.

    Returns
    -------
    None.

    '''
    
    pars_sets = []
    for character in leaf.sequence:
        pars_sets.append(set(character))
    leaf.add_features(parsimony_sets = pars_sets)
    
    pars_scores = [0]*len(leaf.sequence)
    leaf.add_features(parsimony_scores = pars_scores)    

  
def ParsInsertionsInternal(tree, i):
    '''
    Creates the parsimony sets and scores for the internal nodes. 

    Parameters
    ----------
    tree : PhyloNode or PhlyoTree
        Internal nodes a tree structue.
    i : int
        Indey for the sequence.

    Returns
    -------
    None.

    '''
    
    left_set = tree.children[0].parsimony_sets[i]
    left_score = tree.children[0].parsimony_scores[i]
    
    right_set = tree.children[1].parsimony_sets[i]
    right_score  = tree.children[1].parsimony_scores[i]
    
    if left_set == set('-') and right_set == set('-'):
        tree.parsimony_sets[i] = set('-')
        tree.parsimony_scores[i] = left_score + right_score
        
    elif (left_set == set('-') and right_set != set('-')):
        if tree.insertion_flags[i]:
            tree.parsimony_sets[i] = set('-')
            tree.parsimony_scores[i] = left_score + right_score
        else: 
            tree.parsimony_sets[i] = right_set
            tree.parsimony_scores[i] = left_score + right_score + 1
        
    elif (left_set != set('-') and right_set == set('-')):
        if tree.insertion_flags[i]:
            tree.parsimony_sets[i] = set('-')
            tree.parsimony_scores[i] = left_score + right_score
        else: 
            tree.parsimony_sets[i] = left_set
            tree.parsimony_scores[i] = left_score + right_score + 1
        
    elif not left_set.intersection(right_set):
        tree.parsimony_sets[i] = left_set.union(right_set)
        tree.parsimony_scores[i] = left_score + right_score + 1
    else:
        tree.parsimony_sets[i] = left_set.intersection(right_set)
        tree.parsimony_scores[i] = left_score + right_score
       

def ParsInsertionsScore(tree):
    '''
    Calculates the parsimony score for the whole tree while accounting for 
    insertions and deletions.

    Parameters
    ----------
    tree : PhyloNode and PhyloTree
        Input tree with the alignment.

    Returns
    -------
    tree_score : int
        parsimony score with insertions for the whole tree.

    '''
    
    length_MSA = len(tree.get_leaves()[0].sequence)
    
    #empty sets and scores for every node, insertion flags are set to False
    for node in tree.traverse('postorder'):
        pars_sets = [set()] * length_MSA
        pars_scores = [0] * length_MSA
        ins_flags = [False] * length_MSA
        node.add_features(parsimony_sets = pars_sets)
        node.add_features(parsimony_scores = pars_scores)
        node.add_features(insertion_flags = ins_flags)
    
    #sets and scores for leaves
    for leaf in tree.iter_leaves():
        ParsInsertionsLeaf(leaf)
    
    # find insertion points and mark them with an insertion flag set to True
    for i in range(length_MSA):
        leaf_res = []
        for leaf in tree.iter_leaves():
            if leaf.sequence[i] != '-':
                leaf_res.append(leaf)
               
        if len(leaf_res) == 1:
            ancestor = leaf_res[0]
        else:
            ancestor = tree.get_common_ancestor(leaf_res)
    
        if not ancestor.is_root():
            ancestor.up.insertion_flags[i] = True    
    
    #find internal sets and scores
    for i in range(length_MSA):
        for node in tree.traverse('postorder'):
            if not node.is_leaf():
                ParsInsertionsInternal(node, i)
    
    # sum the parsimony scores at the root over the whole sequence
    tree_score = sum(tree.parsimony_scores)
    
    return tree_score
            

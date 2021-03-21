#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 12:47:33 2021

@author: claraiglhaut
"""

from ete3 import PhyloTree, PhyloTree
import random


characters = ['A', 'T', 'C', 'G', '-']

def parsimony_leaves(leaf):
    
    pars_sets = []
    for character in leaf.sequence:
        pars_sets.append(set(character))
    leaf.add_features(parsimony_sets = pars_sets)
    
    pars_scores = [0]*len(leaf.sequence)
    leaf.add_features(parsimony_scores = pars_scores)    
    
def parsimony_internal_per_site(tree, i):
    
    left_set = tree.children[0].parsimony_sets[i]
    left_score = tree.children[0].parsimony_scores[i]
    
    right_set = tree.children[1].parsimony_sets[i]
    right_score  = tree.children[1].parsimony_scores[i]
    
    if not left_set.intersection(right_set):
        tree.parsimony_sets[i] = left_set.union(right_set)
        tree.parsimony_scores[i] = left_score + right_score + 1
    else:
        tree.parsimony_sets[i] = left_set.intersection(right_set)
        tree.parsimony_scores[i] = left_score + right_score
            


def DolloParsimony(tree):
    ''' Input: tree of type PhyloTree
        Output: parsimony score of the tree with Dollo's law
    '''
    
    leaves = [] #get all leaves
    for leaf in tree.iter_leaves():
        leaves.append(leaf)
    
    length_MSA = len(leaves[0].sequence)
    
    #sets and scores for every node
    for node in tree.traverse('postorder'):
        pars_sets = [set()] * length_MSA
        pars_scores = [0] * length_MSA
        node.add_features(parsimony_sets = pars_sets)
        node.add_features(parsimony_scores = pars_scores)
    
    #sets and scores for leaves
    for leaf in tree.iter_leaves():
        parsimony_leaves(leaf)

    
    for i in range(length_MSA):
        #find all leaves with no gap in this position
        no_gap_leaves = [leaf for leaf in leaves if leaf.sequence[i] != '-']
   
        #find the insertion event
        subtree = tree.get_common_ancestor(no_gap_leaves)
        print(subtree)
        # calculate parsimony score
        for node in subtree.traverse('postorder'):
            if not node.is_leaf():
                parsimony_internal_per_site(node, i)
                
        # move the parsimony score to the root node
        tree.parsimony_sets[i] = subtree.parsimony_sets[i]
        tree.parsimony_scores[i] = subtree.parsimony_scores[i]
        
        
    # calculate the total tree score 
    tree_score = sum(tree.parsimony_scores)
    
    return tree_score

def construct_pars_tree_Dollo(tree):
    DolloParsimony(tree)
    sequence_constr = []
    for node in tree.traverse('levelorder'):
        if not node.is_leaf():
            for i in range(len(node.parsimony_sets)):
                if node.parsimony_sets[i] == set():
                    sequence_constr.append('-')
                else:    
                    char = random.sample(node.parsimony_sets[i], 1)[0]
                    sequence_constr.append(char)
    
        node.add_features(sequence = sequence_constr)

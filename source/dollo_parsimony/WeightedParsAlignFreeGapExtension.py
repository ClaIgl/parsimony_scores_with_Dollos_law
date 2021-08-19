#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 13:26:32 2021

@author: claraiglhaut
"""

import numpy as np

cost_matrix = {'T':{'T':0, 'C':1, 'A':1.5, 'G':1.5, '-':10},
               'C':{'T':1, 'C':0, 'A':1.5, 'G':1.5, '-':10},
               'A':{'T':1.5, 'C':1.5, 'A':0, 'G':1, '-':10},
               'G':{'T':1.5, 'C':1.5, 'A':1, 'G':0, '-':10},
               '-':{'T':10, 'C':10, 'A':10, 'G':10, '-':0}}

def InitalizeSetsAndAlignment(leaf):
    '''
    Initializes the nucleotide sets and the alignments at the leaf nodes

    Parameters
    ----------
    leaf : PhlyoNode or PhyloTree
        Tree leaves with ungapped sequences

    Returns
    -------
    None.

    '''

    pars_sets = []
    align = np.empty((1, 0), dtype=str)
    for i in range(len(leaf.sequence)):
        if leaf.sequence[i] != '-':
            character = np.array(leaf.sequence[i]).reshape((1,1))
            align = np.concatenate((align, character), axis=1)
            pars_sets.append(set(leaf.sequence[i]))
                     
    leaf.add_features(parsimony_sets = pars_sets)
    leaf.add_features(alignment = align)
    
   
def GenerateMatricesFreeGapE(tree, cost_matrix):
    '''
    Forward phase of the progressive algorithm. Generates the matrix S with
    the parsimony score and the trace back matrix T, which records the moves. 
    Returns the parsimony score and the trace back matrix. Adds the nucleotide 
    sets and the optimal alignment to the tree.

    Parameters
    ----------
    tree : PhlyoTree or PhyloNode
        Current (sub-)tree

    Returns
    -------
    parsimony_score : numpy.float64
        parsimony score of the alignment for the given tree
    T : numpy.ndarray
        Trace back matrix 

    '''
   
    
    left_sets = tree.children[0].parsimony_sets
    right_sets = tree.children[1].parsimony_sets
    
    S = np.zeros((len(left_sets)+1, len(right_sets)+1), dtype=float)
    T = np.zeros((len(left_sets)+1, len(right_sets)+1), dtype=float)
    
    #fill the first column of the scoring and trace back matrix 
    #only gaps for the right alignment - move vertical
    for i in range(1, len(left_sets)+1):
        costs = [cost_matrix['-'][character] for character in left_sets[i-1]]
        S[i][0] = min(costs)
        T[i][0] = 3
    
    #fill the first row of the scoring and trace back matrix
    #only gaps for the left alignment - move horizontal
    for j in range(1,len(right_sets)+1):
        costs = [cost_matrix['-'][character] for character in right_sets[j-1]]
        S[0][j] = min(costs)
        T[0][j] = 2
    
    #fill the score matrix S
    #matching sets with a non empty intersection +0
    #matching sets with an empty intersection +1
    #gap penalty +1
    for i in range(1, len(left_sets)+1):
        for j in range(1, len(right_sets)+1):
            if left_sets[i-1].intersection(right_sets[j-1]):
                intersection = left_sets[i-1].intersection(right_sets[j-1])
                min_score = np.inf
                for left_character in intersection:
                    for right_character in intersection:
                        score = cost_matrix[left_character][right_character]
                        
                        if score < min_score:
                            min_score = score
            
                score_intersection =  S[i-1][j-1] + min_score
                
            elif not left_sets[i-1].intersection(right_sets[j-1]):
                min_score = np.inf
                for left_character in left_sets[i-1]:
                    for right_character in right_sets[j-1]:
                        score = cost_matrix[left_character][right_character]
                        
                        if score < min_score:
                            min_score = score
                
                score_intersection = S[i-1][j-1] + min_score
            
                
            if T[i][j-1] == 2:
                score_gap_left = S[i][j-1]
                
            else:
                min_score = np.inf
                for left_character in left_sets[i-1]:
                    score = cost_matrix[left_character]['-']
                    
                    if score < min_score:
                            min_score = score
                
                score_gap_left = S[i][j-1] + min_score
           
                
            if T[i-1][j] == 3:
                score_gap_right = S[i-1][j]
            else:
                min_score = np.inf
                for right_character in right_sets[j-1]:
                    score = cost_matrix[right_character]['-']
                    
                    if score < min_score:
                            min_score = score
                            
                score_gap_right = S[i-1][j] + min_score
            
                
            S[i][j] = min(score_intersection, score_gap_left, score_gap_right)
            
            # if path is not unique favor the diagonal move
            if S[i][j] == score_intersection:
                T[i][j] = 1 #move diagonal
            elif S[i][j] == score_gap_right:
                T[i][j] = 3 #move vertical
            elif S[i][j] == score_gap_left:
                T[i][j] = 2 #move horizontal
    
    parsimony_score = S[len(left_sets)][len(right_sets)]
    
    return parsimony_score, T


def TraceBack(T, tree):
    '''
    Finds the alignment for the (sub-)tree and adds it to the (sub-)tree root. 
    Adds the nucleotide sets to the (sub-)tree root.

    Parameters
    ----------
    T : numpy.ndarray
        Trace back matrix for the alognment
    tree : PhyloTree or PhyloNode
        Current (sub-)tree

    Returns
    -------
    None.

    '''
    
    #get the alignemts from the left and right child
    left_alignment = tree.children[0].alignment
    right_alignment = tree.children[1].alignment
    
    #get the nucleotide sets from the left and right child
    left_sets = tree.children[0].parsimony_sets
    right_sets = tree.children[1].parsimony_sets

    number_of_rows = len(left_alignment) + len(right_alignment)
    
    align = np.empty((number_of_rows, 0), dtype=str)
    
    
    i = left_alignment.shape[1]
    j = right_alignment.shape[1]
    
    pars_sets = []
    
    while i > 0 or j > 0:
        
        #move diagonal - match two colums with residues 
        if T[i][j] == 1:
            new_col = np.empty((number_of_rows,1), dtype=str)
            for n in range(len(left_alignment)):
                new_col[n] = left_alignment[n][i-1]
            for m in range(len(right_alignment)):
                new_col[len(left_alignment)+m] = right_alignment[m][j-1]
                
            if not left_sets[i-1].intersection(right_sets[j-1]):
                pars_sets.insert(0, left_sets[i-1].union(right_sets[j-1]))
            else:
                pars_sets.insert(0, left_sets[i-1].intersection(right_sets[j-1]))

            i = i-1
            j = j-1 
        
        #move vertical - put a gap column in the left alignment and match 
        #with the column of the right alignment
        elif T[i][j] == 2:
            new_col = np.empty((number_of_rows,1), dtype=str)
            for n in range(len(left_alignment)):
                new_col[n] = '-'
            for m in range(len(right_alignment)):
                new_col[len(left_alignment)+m] = right_alignment[m][j-1]
           
            pars_sets.insert(0, right_sets[j-1])
            
            j = j-1    
        
        #move horizontal - put a gap column in the right alignment and match 
        #with the column of the left alignment
        elif T[i][j] == 3:
            new_col = np.empty((number_of_rows,1), dtype=str)
            for n in range(len(left_alignment)):
                new_col[n] = left_alignment[n][i-1]
                
            for m in range(len(right_alignment)):
                new_col[len(left_alignment)+m] = '-' 
      
            pars_sets.insert(0, left_sets[i-1])
 
            i = i-1   

        align = np.concatenate((new_col, align), axis=1)
        
    tree.add_features(alignment = align)
    tree.add_features(parsimony_sets = pars_sets)  
 
    
def WeightedParsAlignFreeGapE(tree):
    '''
    Finds the Multiple Sequence Alignment for the given tree.
    
    Parameters
    ----------
    tree : PhyloNode or PhyloTree
        Phylogenetic Tree with ungapped sequences at the leaves

    Returns
    -------
    parsimony_score : numpy.float64
        parsimony score for the alignment on the tree
    alignment : numpy.ndarray
        Multiple sequence alignment for the given tree 

    '''
    
    parsimony_score = 0   
    for node in tree.traverse('postorder'):
        if node.is_leaf():
            InitalizeSetsAndAlignment(node)    
        else:
            pars_score, T = GenerateMatricesFreeGapE(node, cost_matrix)
            TraceBack(T, node)
            parsimony_score = parsimony_score + pars_score
    alignment = tree.alignment 
    
    return parsimony_score, alignment


   






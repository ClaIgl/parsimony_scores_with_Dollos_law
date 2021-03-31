#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 16:19:26 2021

@author: claraiglhaut
"""
from ete3 import PhyloTree
import numpy as np
import random

sequences = '/Users/claraiglhaut/Desktop/ZHAW/TrackModule2/test_data_MSA/test_MSA_sequence1'
newick = '/Users/claraiglhaut/Desktop/ZHAW/TrackModule2/test_data_MSA/test_MSA_tree1'

tree = PhyloTree(newick=newick, alignment=sequences)
print(tree)
  

def ParsimonySetsLeaves(leaf):
    pars_sets = [set(character) for character in leaf.sequence]
    align = np.empty((1, len(leaf.sequence)), dtype=str)
    for i in range(len(leaf.sequence)):
        align[0][i] = leaf.sequence[i]
    
                     
    leaf.add_features(parsimony_sets = pars_sets)
    leaf.add_features(alignment = align)

#%%    
def GenerateMatrices(tree):
    left_sets = tree.children[0].parsimony_sets
    right_sets = tree.children[1].parsimony_sets
    
    S = np.zeros((len(left_sets)+1, len(right_sets)+1))
    T = np.zeros((len(left_sets)+1, len(right_sets)+1))
    
    for i in range(len(left_sets)+1):
        S[i][0] = i
        if i == 0:
            T[i][0] = 0
        else:
            T[i][0] = 3
            
    for j in range(len(right_sets)+1):
        S[0][j] = j
        if j == 0:
            T[0][j] = 0
        else:
            T[0][j] = 2
        
    for i in range(1, len(left_sets)+1):
        for j in range(1, len(right_sets)+1):
            if left_sets[i-1].intersection(right_sets[j-1]):
                score_intersection =  S[i-1][j-1] 
            elif not left_sets[i-1].intersection(right_sets[j-1]):
                score_intersection = S[i-1][j-1] + 1
            score_gap_left = S[i][j-1] + 1
            score_gap_right = S[i-1][j] + 1
            
            S[i][j] = min(score_intersection, score_gap_left, score_gap_right)
            
            moves = list()
            if S[i][j] == score_intersection:
                moves.append(1) #move diagonal
            if S[i][j] == score_gap_left:
                moves.append(2) #move horizontal
            if S[i][j] == score_gap_right:
                moves.append(3) #move vertical
            
            T[i][j] = random.choice(moves)
    
    return S, T


def TraceBack(T, tree):

    left_alignment = tree.children[0].alignment
    right_alignment = tree.children[1].alignment

    number_of_rows = len(left_alignment) + len(right_alignment)
    
    align = np.empty((number_of_rows, 0), dtype=str)
    
    i = left_alignment.shape[1]
    j = right_alignment.shape[1]
    
    
    while i > 0 or j > 0:
        
        if T[i][j] == 1:
            new_col = np.empty((number_of_rows,1), dtype=str)
            for n in range(len(left_alignment)):
                new_col[n] = left_alignment[n][i-1]
            for m in range(len(right_alignment)):
                new_col[len(left_alignment)+m] = right_alignment[m][j-1]  
                
            i = i-1
            j = j-1 
            
        elif T[i][j] == 2:
            new_col = np.empty((number_of_rows,1), dtype=str)
            for n in range(len(left_alignment)):
                new_col[n] = '-'
            for m in range(len(right_alignment)):
                new_col[len(left_alignment)+m] = right_alignment[m][j-1]
            
            j = j-1    
            
        elif T[i][j] == 3:
            new_col = np.empty((number_of_rows,1), dtype=str)
            for n in range(len(left_alignment)):
                new_col[n] = left_alignment[n][i-1]
                
            for m in range(len(right_alignment)):
                new_col[len(left_alignment)+m] = '-' 
  
            i = i-1   

        align = np.concatenate((new_col, align), axis=1)
        
    tree.add_features(alignment = align)
        
    
    
    
    
#%%    
for leaf in tree.iter_leaves():
    ParsimonySetsLeaves(leaf)
    print(leaf.parsimony_sets)
    t = leaf.alignment
print(tree.children[0])
print(tree.children[1])
S, T  = GenerateMatrices(tree)
TraceBack(T, tree)
print(tree.alignment)
print(S)
print(T)
    

                
                
                
            
            
        
    
    
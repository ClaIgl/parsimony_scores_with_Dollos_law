#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 11 17:31:16 2021

@author: claraiglhaut
"""
from ete3 import PhyloNode
import pytest
from dollo_parsimony.dollo_parsimony_MSA import GenerateMatrices

characters = characters = ['A', 'T', 'C', 'G']

@pytest.mark.parametrize('child0_pars_set, child1_pars_set, expected_score, expected_T, message',
             [([set(characters[0])], [set(characters[1])], 1, [[0,2],[3, 1]], 'one character per set'),
              ([set(characters[0])], [set(characters[0]), set(characters[1])], 1, [[0,2,2],[3,1,2]], 'one character left set, two characters right set'),
              ([set(characters[2]), set(characters[3])], [set(characters[3])], 1, [[0,2],[3,1],[3,1]],'two characters left set, one character right set')])


def test_GenerateMatrices(child0_pars_set, child1_pars_set, expected_score, expected_T, message):
    newick = '(A:1,B:1):1;'
    node = PhyloNode(newick=newick)
    node.children[0].parsimony_sets = child0_pars_set
    node.children[1].parsimony_sets = child1_pars_set
    pars_score, T = GenerateMatrices(node)
    
    assert pars_score == expected_score, 'wrong score for' + message
    assert (T == expected_T).all(), 'wrong trace back matrix for' + message
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 20:30:54 2021

@author: claraiglhaut
"""

import numpy as np
from ete3 import PhyloNode
import pytest
from dollo_parsimony.dollo_parsimony_MSA import TraceBack

characters = characters = ['A', 'T', 'C', 'G']

@pytest.mark.parametrize(
    '''child0_pars_set, child1_pars_set, child0_alignment,  child1_alignment, T, 
    expected_alignment, expected_pars_sets, message''',
    [([set(characters[0])], [set(characters[1])], np.array([characters[0]], dtype=str, ndmin=2), np.array([characters[1]], dtype=str, ndmin=2), 
      [[0,2],[3, 1]], [[characters[0]], [characters[1]]], [set(characters[0]).union(characters[1])], 'one character per set'),
     ([set(characters[0])], [set(characters[0]), set(characters[1])], np.array([characters[0]], dtype=str, ndmin=2), np.array([characters[0], characters[1]],dtype=str, ndmin=2), [[0,2,2],[3,1,2]], 
      [[characters[0], '-'], [characters[0], characters[1]]], [set(characters[0]), set(characters[1])], 'one character left set, two characters right set'),
     ([set(characters[2]), set(characters[3])], [set(characters[3])], np.array([characters[2], characters[3]], dtype=str, ndmin=2), np.array([characters[3]], dtype=str, ndmin=2), 
      [[0,2],[3,1],[3,1]], [[characters[2], characters[3]], ['-', characters[3]]], [set(characters[2]), set(characters[3])], 'two characters left set, one character right set')])


def test_TraceBack(child0_pars_set, child1_pars_set, child0_alignment,  child1_alignment, T, 
    expected_alignment, expected_pars_sets, message):
    newick = '(A:1,B:1):1;'
    node = PhyloNode(newick=newick)
    node.children[0].parsimony_sets = child0_pars_set
    node.children[1].parsimony_sets = child1_pars_set
    node.children[0].alignment = child0_alignment
    node.children[1].alignment = child1_alignment
    
    TraceBack(T, node)
    
    assert (node.alignment == expected_alignment).all(), 'wrong alignment for ' + message
    assert node.parsimony_sets == expected_pars_sets, 'wrong sets for ' + message
    

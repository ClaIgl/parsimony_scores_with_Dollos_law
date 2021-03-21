#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 12:07:03 2021

@author: claraiglhaut
"""

import pytest

from ete3 import PhyloTree

from dollo_parsimony.dollo_parsimony_score import DolloParsimony

@pytest.mark.parametrize("newick,alignment,score,message",
    [('../test_data/test_tree','../test_data/test_sequence.txt',2,"wrong score for test_tree and sequences.txt"),
     ('../test_data/test_tree','../test_data/test_sequence1',5,"wrong score for test_tree and sequences1"),
     ('../test_data/test_tree1','../test_data/test_sequence2',3,"wrong score for test_tree and sequences1")])

def test_example_trees_sequences(newick, alignment, score, message):
    tree = PhyloTree(newick=newick, alignment=alignment)
    assert DolloParsimony(tree) == score, message


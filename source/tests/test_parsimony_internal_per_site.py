#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 11:38:55 2021

@author: claraiglhaut
"""

import pytest

from ete3 import PhyloNode

from dollo_parsimony.ParsInsertionsScore import ParsInsertionsInternal
from dollo_parsimony.ParsInsertionsScore import characters

@pytest.mark.parametrize("ch0_score,ch0_set,ch1_score,ch1_set,exp_score,exp_set, insertion, message",
    [(0,set('-'),0,set('-'),0,set('-'),False,"two gap sets with 0 score"),
    (1,set('-'),0,set('-'),1,set('-'),False,"two gap sets with non-zero score on first"),
    (0,set('-'),2,set('-'),2,set('-'),False,"two gap sets with non-zero score on second"),
    (1,set('-'),2,set('-'),3,set('-'),False,"two gap sets with non-zero scores"),
    (1,set(characters[0]),0,set('-'),2,set(characters[0]),False, "one non-empty set, single element"),
    (1,set(characters[0]),0,set('-'),1,set('-'),True, "one non-empty set, single element"),
    (1,set(characters[0:2]),0,set('-'),2,set(characters[0:2]),False, "one non-empty set, two elements"),
    (1,set(characters[0:2]),0,set('-'),1,set('-'),True, "one non-empty set, two elements"),
    (1,set('-'),0,set(characters[1:3]),2,set(characters[1:3]),False,"second non-empty set, two elements"),
    (1,set('-'),0,set(characters[1:3]),1,set('-'),True,"second non-empty set, two elements"),
    (1,set(characters[0]),2,set(characters[0:2]),3,set(characters[0]),False,"intersection not empty, single element"),
    (3,set(characters[1]),2,set(characters[0:2]),5,set(characters[1]),False,"intersection not empty, single element"),
    (3,set(characters[0:]),2,set(characters[0:]),5,set(characters[0:]),False,"intersection not empty, multiple elements"),
    (1,set(characters[0]),1,set(characters[1]),3,set(characters[0:2]),False,"intersection empty, two elements"),
    (1,set(characters[0]),2,set(characters[1:]),4,set(characters[0:]),False,"intersection empty, all elements")])

def test(ch0_score, ch0_set, ch1_score, ch1_set, exp_score, exp_set, insertion, message):
    newick = '(A:1,B:1):1;'
    tree = PhyloNode(newick=newick)
    
    tree.parsimony_scores = [0]
    tree.parsimony_sets = [set()]
    tree.insertion_flags = [insertion]
    tree.children[0].parsimony_scores = [ch0_score]
    tree.children[0].parsimony_sets = [ch0_set]
    tree.children[1].parsimony_scores = [ch1_score]
    tree.children[1].parsimony_sets = [ch1_set]
    
    ParsInsertionsInternal(tree, 0)
    
    assert tree.parsimony_scores[0] == exp_score, "wrong score for " + message
    assert tree.parsimony_sets[0] == exp_set, "wrong char set for " + message
    

    
    
    

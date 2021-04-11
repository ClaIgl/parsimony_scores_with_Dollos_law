#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 14:50:52 2021

@author: claraiglhaut
"""
import pytest

from ete3 import PhyloNode
import numpy as np
from dollo_parsimony.dollo_parsimony_MSA import InitalizeSetsAndAlignment


characters = ['A', 'T', 'C', 'G']

@pytest.mark.parametrize('sequence, expected_set, expected_alignment, message',
                         [(characters[0], [set(characters[0])], [[characters[0]]], 'one character'),
                          (characters[1], [set(characters[1])], [[characters[1]]], 'another character'),
                          ([characters[0], characters[1]], [set(characters[0]), set(characters[1])], [[characters[0], characters[1]]], 'two characters')])


def test_initial(sequence, expected_set, expected_alignment, message):
    node = PhyloNode()
    node.sequence = sequence
    InitalizeSetsAndAlignment(node)
    assert len(node.sequence) == len(node.parsimony_sets), 'wrong number of sets for' + message
    assert len(node.alignment) == 1, 'wrong size of alignment for one sequence with' + message
    assert node.parsimony_sets == expected_set, 'wrong sets for' + message 
    assert (node.alignment == expected_alignment).all()

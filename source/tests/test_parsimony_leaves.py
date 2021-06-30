import pytest

from ete3 import PhyloNode

from dollo_parsimony.ParsInsertionsScore import ParsInsertionsLeaf
from dollo_parsimony.ParsInsertionsScore import characters

@pytest.mark.parametrize("sequence,expected_set,message",
    [(characters[0],set(characters[0]),"one char"),
    (characters[3],set(characters[3]),"another char"),
    (characters[4],set(characters[4]),"gap")])

def test_one_char(sequence, expected_set, message):
    node = PhyloNode()
    node.sequence = sequence
    ParsInsertionsLeaf(node)
    assert len(node.parsimony_scores) == len(node.sequence), "wrong number of scores for " + message
    assert len(node.parsimony_sets) == len(node.sequence), "wrong number of sets for " + message
    assert node.parsimony_scores[0] == 0, "wrong score for " + message
    assert node.parsimony_sets[0] == expected_set, "wrong set for " + message

@pytest.mark.parametrize("sequence,message",
    [(str(characters),"multiple chars"),
    (str(characters[0]*200),"multiple chars try two")])

def test_multiple_chars(sequence, message):
    node = PhyloNode()
    node.sequence = sequence
    ParsInsertionsLeaf(node)
    assert len(node.parsimony_scores) == len(node.sequence), "wrong number of scores for " + message
    assert len(node.parsimony_sets) == len(node.sequence), "wrong number of sets for " + message

import scripts.codon_utils as codon_utils


def test_get_degeneracy():
    assert codon_utils.get_degeneracy('ATT', 2) == 2
    assert codon_utils.get_degeneracy('TTT', 0) == 0
    assert codon_utils.get_degeneracy('TTT', 1) == 0
    assert codon_utils.get_degeneracy('TTT', 2) == 2
    assert codon_utils.get_degeneracy('GTT', 2) == 4
    assert codon_utils.get_degeneracy('ATG', 0) == 0
    assert codon_utils.get_degeneracy('ATG', 1) == 0
    assert codon_utils.get_degeneracy('ATG', 2) == 0
    assert codon_utils.get_degeneracy('TAA', 2) == 2


def test_is_synonymous():
    assert codon_utils.is_synonymous('TTT', 'TTC')
    assert codon_utils.is_synonymous('CCA', 'CCG')
    assert codon_utils.is_synonymous('TAA', 'TAG')
    assert codon_utils.is_synonymous('ATG', 'ATG')
    assert not codon_utils.is_synonymous('AAA', 'CAA')
    assert not codon_utils.is_synonymous('TGA', 'TCA')


def test_is_synonymous():
    assert codon_utils.mutate('AGT', 'C', 0) == 'CGT'
    assert codon_utils.mutate('AGT', 'C', 1) == 'ACT'
    assert codon_utils.mutate('AGT', 'C', 2) == 'AGC'

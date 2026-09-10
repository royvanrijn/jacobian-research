#!/usr/bin/env sage-python
"""Read-only replay of the sealed generic input, without the producer."""
import json
from sage.all import *
from freeze_x1092_class1_arithmetic import check, DEST, sha

check()
data=json.loads((DEST/'generic-arithmetic.json').read_text())
seal=json.loads((DEST/'generic-input-seal.json').read_text())
assert seal['generic_arithmetic_sha256']==sha(DEST/'generic-arithmetic.json')
assert data['parent_sha256']==sha(DEST/'parent.json')
parent=json.loads((DEST/'parent.json').read_text())
sections=json.loads((DEST/'section-certificate.json').read_text())
R=PolynomialRing(QQ,'s');K=R.fraction_field()
def dec(r):return K(R(r['numerator']))/R(r['denominator'])
S=PolynomialRing(K,'z');z=S.gen()
coeff=[dec(c) for c in data['cubic_coefficients_low_to_high']]
f=S(coeff);E=EllipticCurve(K,[dec(a) for a in parent['a_invariants']])
assert f==z**3+E.a4()*z+E.a6() and f.gcd(f.derivative())==1
assert dec(data['cubic_discriminant'])==f.discriminant()
assert len(data['inherited_generators'])==17
assert data['generic_height_gram']==parent['generic_height_gram']
for i,row in enumerate(data['inherited_generators']):
    assert row['basis_index']==i
    x,y=map(dec,parent['basis_weierstrass_coordinates'][i]);E(x,y)
    assert dec(row['x'])==x and dec(row['y'])==y
    gamma=S([dec(c) for c in row['kummer_representative_coefficients']])
    assert gamma==x-z
    M=matrix(K,[[dec(c) for c in r] for r in row['multiplication_matrix']])
    for j in range(3):assert S(list(M.column(j)))==(gamma*z**j)%f
    assert M.det()==x**3+E.a4()*x+E.a6()==y*y==dec(row['norm'])
    assert dec(row['norm_square_root'])==y
    assert row['source_frame_coordinates']==sections['sections'][i]['child_frame_coordinates']
assert data['parameters_selected']==data['prospective_strict_classes']==[]
for key in ('complete_inherited_everywhere_even_half_ideal_image','strict_character_subspace','ordinary_unramified_character_subspace'):
    assert data[key]=='UNKNOWN'
print('PASS sealed generic arithmetic: exact cubic and 17 norm-square representatives; no arithmetic-space completeness or novelty claimed')

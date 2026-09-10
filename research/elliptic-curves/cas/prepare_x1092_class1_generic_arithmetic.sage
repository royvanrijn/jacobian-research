#!/usr/bin/env sage-python
"""Generic arithmetic input only, on the byte-frozen certified equation."""
import json
from pathlib import Path
from sage.all import *
from freeze_x1092_class1_arithmetic import check, DEST, sha, write_once

freeze=check();parent=json.loads((DEST/'parent.json').read_text())
R=PolynomialRing(QQ,'s');K=R.fraction_field();s=R.gen()
def dec(r):return K(R(r['numerator']))/R(r['denominator'])
def rf(v):
    v=K(v)
    return {'numerator':list(map(str,v.numerator().list())) or ['0'],
            'denominator':list(map(str,v.denominator().list())) or ['1']}
E=EllipticCurve(K,[dec(r) for r in parent['a_invariants']])
assert E.a1()==E.a2()==E.a3()==0
S=PolynomialRing(K,'theta');theta=S.gen();f=theta**3+E.a4()*theta+E.a6()
assert f.gcd(f.derivative())==1
G=matrix(ZZ,parent['generic_height_gram'])
assert G.is_positive_definite() and G.det()==1092 and G.nrows()==17
records=[]
for i,record in enumerate(parent['basis_weierstrass_coordinates']):
    x,y=map(dec,record);E(x,y)
    gamma=S(x-theta)
    norm=f.resultant(gamma)
    assert norm==y*y
    # Multiplication in the etale cubic quotient: a second exact norm check.
    columns=[(gamma*theta**j)%f for j in range(3)]
    mult=matrix(K,3,3,lambda a,b:columns[b][a])
    assert mult.det()==norm
    records.append({'basis_index':i,'x':rf(x),'y':rf(y),
      'kummer_representative_coefficients':[rf(x),rf(-1),rf(0)],
      'norm_square_root':rf(y),'norm':rf(norm),
      'multiplication_matrix':[[rf(c) for c in row] for row in mult],
      'source_frame_coordinates':json.loads((DEST/'section-certificate.json').read_text())['sections'][i]['child_frame_coordinates']})
out={'schema':'x1092.class1.generic-arithmetic-inputs.v1',
 'status':'PASS_EXACT_GENERIC_CUBIC_AND_17_INHERITED_NORM_IDENTITIES',
 'parent_sha256':sha(DEST/'parent.json'),'freeze_sha256':sha(DEST/'freeze.json'),
 'producer_sha256':sha(Path(__file__)),
 'base_field':'Q(s)','algebra':'Q(s)[theta]/(theta^3+A(s)*theta+B(s))',
 'cubic_coefficients_low_to_high':[rf(f[i]) for i in range(4)],
 'cubic_discriminant':rf(f.discriminant()),'etale_over_Q_s':True,
 'marked_basis':'1,theta,theta^2','generic_height_gram':parent['generic_height_gram'],
 'inherited_generators':records,
 'complete_inherited_everywhere_even_half_ideal_image':'UNKNOWN',
 'strict_character_subspace':'UNKNOWN','ordinary_unramified_character_subspace':'UNKNOWN',
 'maximal_order':'NOT_COMPUTED','bad_place_support':'NOT_COMPUTED',
 'prospective_strict_classes':[], 'parameters_selected':[],
 'boundary':'These are generic section Kummer representatives and exact principal norm-square identities, not a complete half-ideal image, arithmetic strict classes, soluble-cover certificates, or quotient directions.'}
write_once(DEST/'generic-arithmetic.json',(json.dumps(out,indent=2,sort_keys=True)+'\n').encode())
print(out['status'])

#!/usr/bin/env sage-python
"""Independent algebraic replay for Br(X)/Br(Q)=0 on the fixed parent.

Reuses immutable Fp^2 counts. Checks the reciprocal zeta convention, maximal
arithmetic Picard lattice and good reduction without importing the producer.
The cohomological implication is proved in the canonical note, not by CAS.
"""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,vector,block_diagonal_matrix
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json';OLD=ART/'curve302_parent_geometric_picard19_v1.json'
OUT=ART/'det1092_surface_brauer_triviality_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for path,expected in read(OUT/'protocol.json')['inputs'].items():assert sha(ROOT/path)==expected
d=read(PARENT);old=read(OLD);assert old['status']=='PASS'
for path,expected in old['sources'].items():assert sha(ROOT/path)==expected
G=matrix(ZZ,d['generic_height_gram']);N=block_diagonal_matrix(matrix(ZZ,[[-2,1],[1,0]]),-G)
lat=read(OUT/'arithmetic-lattice.json');assert N.det()==lat['determinant']==1092
assert list(ZZ(1092).factor())==[(ZZ(2),2),(ZZ(3),1),(ZZ(7),1),(ZZ(13),1)]
v=vector(ZZ,lat['unique_order_two_dual_word']);assert any(x%2 for x in v)
assert N.change_ring(GF(2)).right_kernel().dimension()==1
assert all(x%2==0 for x in N*v) and v*N*v==-180
# Thus v/2 is the unique nonzero order-two class in the discriminant group,
# but its square is odd. No even index-two overlattice can exist.
assert abs(N.det())==1092 and str(QQ(v*N*v)/4)==lat['half_word_square']
R=PolynomialRing(QQ,'t');K=R.fraction_field()
def dec(x):return K(R(x['numerator']))/R(x['denominator'])
a1,a2,a3,a4,a6=[dec(x) for x in d['a_invariants']]
b2=a1*a1+4*a2;b4=a1*a3+2*a4;b6=a3*a3+4*a6
A=R(b4/2-b2*b2/48);B=R(b6/4-b2*b4/24+b2**3/864)
assert A.degree()==8 and B.degree()==12
S=PolynomialRing(QQ,'u');u=S.gen();results=[]
reductions=read(OUT/'reductions.json')['reductions']
assert [r['prime'] for r in reductions]==[149,151]
for saved,prior in zip(reductions,old['reductions']):
    p=ZZ(saved['prime']);assert p.is_prime(proof=True) and prior['prime']==p
    Rp=PolynomialRing(GF(p),'t');aa,bb=Rp(A),Rp(B);disc=4*aa**3+27*bb**2
    assert disc.degree()==24 and disc.gcd(disc.derivative()).degree()==0 and disc[24]
    assert disc.gcd(aa).degree()==0
    char=S(prior['full_H2_polynomial']);assert char.degree()==22 and char.is_monic()
    reciprocal=S(list(reversed(char.list())))
    residual,rem=reciprocal.quo_rem((1-p*u)**19);assert rem==0 and residual(QQ(1)/p)!=0
    ap=ZZ(saved['residual_pair_trace']);assert residual==(1+p*u)*(1-ap*u+p*p*u*u)
    n1=1+p*p-reciprocal[1]
    n2=1+p**4+reciprocal[1]**2-2*reciprocal[2]
    assert n1==saved['base_field_surface_count_independent']==prior['counts'][0]['surface_count']
    assert n2==saved['squared_field_count_reused']==prior['counts'][1]['surface_count']
    # Independent Artin--Tate normalization using det(1-u*Frob), not the
    # producer's monic characteristic polynomial evaluation at p.
    value=p*residual(QQ(1)/p)/abs(N.det());assert value==1
    assert value==saved['Brauer_group_order'] and p*p*1092==QQ(saved['Artin_Tate_numerator'])
    results.append({'p':int(p),'Picard_rank_over_Fp':19,'Picard_index':1,
        'Brauer_order':1,'prime_to_p_Kummer_fixed_space':'exactly the specialized global Picard group modulo n, for every n coprime to p'})
assert results[0]['p']!=results[1]['p']
paths=[PARENT,OLD,OUT/'protocol.json',OUT/'arithmetic-lattice.json',OUT/'reductions.json',Path(__file__)]
report={'status':'PASS_ARITHMETIC_INPUTS_FOR_GLOBAL_BRAUER_TRIVIALITY',
    'classification':'new deduction from verified finite-field arithmetic and established cohomology',
    'reductions':results,'global_conclusion':'Br(X)=Br(Q), via the rational zero-section point',
    'proof_gate':'For each torsion prime ell, use one of149,151 different from ell. Kummer plus finite-field Hochschild--Serre identifies Frobenius-fixed H2(mu_ell) with Pic(X_p)/ell because Br(X_p)=0. Smooth proper base change and Picard specialization identify it with global divisor classes. Any global Brauer ell-class is therefore algebraic, and trivial Galois action on the torsion-free Picard group makes it constant. Normalize at a rational point to pass from prime order to every torsion order.',
    'boundary':'Not vanishing of the geometric Brauer group, Br(Xbar)^G, or any specialized elliptic Selmer/Sha group. No spectral-Jacobian torsion computation or seed construction.',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths}}
dest=OUT/'independent-replay.json'
if dest.exists():assert read(dest)==report
else:
    with dest.open('x') as s:json.dump(report,s,indent=2,sort_keys=True);s.write('\n')
print('PASS_ARITHMETIC_INPUTS_FOR_GLOBAL_BRAUER_TRIVIALITY',[(r['p'],r['Brauer_order']) for r in results],flush=True)

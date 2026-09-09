#!/usr/bin/env sage-python
"""Reinterpret the two committed Frobenius certificates by Artin--Tate over Fp.

No new prime, extension-field count, Jacobian/Selmer or class-group calculation.
Direct Fp point counts only; external 25-second cap.
"""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector,block_diagonal_matrix
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
OLD=ART/'curve302_parent_geometric_picard19_v1.json'
OUT=ART/'det1092_surface_brauer_triviality_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    OUT.mkdir(exist_ok=True)
    path=OUT/name
    if path.exists():assert read(path)==data
    else:
        with path.open('x') as s:json.dump(data,s,indent=2,sort_keys=True);s.write('\n')
paths=[PARENT,OLD,Path(__file__)]
save('protocol.json',{'classification':'new arithmetic interpretation of immutable geometric evidence',
    'rule':'Use precisely the completed149 and151 Frobenius records. Recheck good reduction, the integral rank19 arithmetic divisor lattice and Artin--Tate over Fp rather than Fp squared. Recount Fp points directly; reuse the certified Fp-squared counts.',
    'limits':{'seconds':25,'primes':[149,151],'new_extension_field_counts':0,'new_test_parameters':0,'point_searches':0,'Selmer_groups':0,'class_groups':0,'exceptional_points':0},
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths}})
d=read(PARENT);old=read(OLD);assert old['status']=='PASS'
for path,expected in old['sources'].items():assert sha(ROOT/path)==expected
R=PolynomialRing(QQ,'t');t=R.gen();K=R.fraction_field()
def dec(r):return K(R(r['numerator']))/R(r['denominator'])
E=EllipticCurve(K,[dec(r) for r in d['a_invariants']]);A=R(-E.c4()/48);B=R(-E.c6()/864)
delta=-16*(4*A**3+27*B**2);assert [A.degree(),B.degree(),delta.degree()]==[8,12,24]
G=matrix(ZZ,d['generic_height_gram']);L=block_diagonal_matrix(matrix(ZZ,[[-2,1],[1,0]]),-G)
assert L.det()==1092 and G.is_positive_definite()
smith=L.smith_form()[0];assert [abs(smith[i,i]) for i in range(19)]==[1]*18+[1092]
kernel=L.change_ring(GF(2)).right_kernel_matrix();assert kernel.nrows()==1
v=vector(ZZ,[int(x) for x in kernel.row(0)]);half=vector(QQ,v)/2
assert all(x in ZZ for x in L*half) and half*L*half==-45
save('arithmetic-lattice.json',{'status':'NO_PROPER_EVEN_SAME_RANK_OVERLATTICE',
    'rank':19,'determinant':1092,'smith_diagonal':[1]*18+[1092],
    'only_possible_nontrivial_index':2,'unique_order_two_dual_word':list(map(int,v)),
    'half_word_square':'-45','interpretation':'At either fixed prime the Frobenius eigenvalue p has multiplicity19. The inherited rank19 lattice is therefore the whole arithmetic Neron--Severi lattice: the sole possible index2 enlargement has odd square.'})
S=PolynomialRing(QQ,'T');T=S.gen();records=[]
for p,rec in zip([149,151],old['reductions']):
    assert rec['prime']==p
    Rp=PolynomialRing(GF(p),'t');ap,bp=Rp(A),Rp(B);dp=-16*(4*ap**3+27*bp**2)
    assert dp.degree()==24 and dp.gcd(dp.derivative()).degree()==0 and ap.gcd(dp).degree()==0
    assert 4*ap[8]**3+27*bp[12]**2
    counts=[row['surface_count'] for row in rec['counts']]
    total=0
    for addr in list(range(p))+[None]:
        a,b=(int(ap(addr)),int(bp(addr))) if addr is not None else (int(ap[8]),int(bp[12]))
        fibre=1
        for x in range(p):
            value=(x*x*x+a*x+b)%p
            fibre+=1 if value==0 else (2 if pow(value,(p-1)//2,p)==1 else 0)
        total+=fibre
    assert total==counts[0]
    s1=ZZ(counts[0])-1-p*p-19*p;s2=ZZ(counts[1])-1-p**4-19*p*p
    signs=[e for e in [-1,1] if s2==s1*s1-2*e*p*s1];assert signs==[-1]
    trace=s1+p;residual=(T+p)*(T*T-trace*T+p*p)
    full=(T-p)**19*residual;assert full==S(rec['full_H2_polynomial']) and residual(p)!=0
    # K3 Artin--Tate in characteristic-polynomial convention:
    # residual(p) = p^(22-19-1) * #Br(X_p) * |disc NS(X_p)|.
    numerator=QQ(residual(p));denominator=QQ(p*p*1092);order=numerator/denominator
    assert order==1 and (4*p-2*trace)==1092
    records.append({'prime':p,'good_reduction_24I1':True,'arithmetic_Picard_rank':19,
        'geometric_Picard_rank':20,'arithmetic_NS_determinant':1092,
        'base_field_surface_count_independent':total,'squared_field_count_reused':counts[1],
        'residual_pair_trace':int(trace),'residual_cubic':list(map(str,residual.list())),
        'Artin_Tate_numerator':str(numerator),'Artin_Tate_denominator':str(denominator),
        'Brauer_group_order':1})
save('reductions.json',{'status':'TRIVIAL_BRAUER_GROUPS_AT_BOTH_FIXED_REDUCTIONS',
    'classification':'verified application of Tate and Artin--Tate for K3 surfaces',
    'reductions':records,'boundary':'Uses the completed Fp-squared counts, not new independent extension-field counts. Arithmetic NS has rank19; the larger geometric rank20 lattice is not substituted in the formula.'})
print('PASS_REDUCTION_BRAUER_ORDERS',[(r['prime'],r['Brauer_group_order']) for r in records],flush=True)

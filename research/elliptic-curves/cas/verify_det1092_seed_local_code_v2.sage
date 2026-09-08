#!/usr/bin/env sage-python
"""Exact replay using rational and finite-ring arithmetic only.

No number-field maximal-order, prime-ideal, unit, logarithm, or local
squareclass routine is used. A suggested order is certified p-maximal by
the radical multiplier criterion. Square tests use Frobenius (odd p) or
all 64 possible roots modulo4 (p=2). Nonunit nonsquares are certified by
their nontrivial radical principal ideals. Each case has a25-second cap.
"""
import argparse,hashlib,json,signal
from itertools import product
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,vector,identity_matrix,lcm
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
DIR=ART/'det1092_seed_local_code_v5'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True,default=int)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def provenance(d):
    for path,digest in d['inputs'].items():assert sha(ROOT/path)==digest,(path,'HASH_MISMATCH')
def residue(a,m):
    a=QQ(a)
    assert a.denominator().gcd(m)==1
    return int(a.numerator()*a.denominator().inverse_mod(m)%m)
def signs_at_roots(f,x):
    seq=[f,f.derivative()]
    while seq[-1].degree()>0:seq.append(-(seq[-2]%seq[-1]))
    def changes(s):
        s=[a for a in s if a]
        return sum(a!=b for a,b in zip(s,s[1:]))
    lo=changes([int(g.leading_coefficient().sign())*(-1)**g.degree() for g in seq])
    hi=changes([int(g.leading_coefficient().sign()) for g in seq])
    at=changes([int(g(x).sign()) for g in seq])
    assert f(x) and lo-hi in [1,3]
    return [0]*(lo-at)+[1]*(at-hi)

class Algebra:
    def __init__(self,f,basis):
        self.f=f;self.R=f.parent();self.basis=basis
        self.B=matrix(QQ,[self.polycoords(b) for b in basis]).transpose()
        self.Binv=self.B.inverse()
        assert basis[0]==1 and self.B.det()
        assert all(v in ZZ for v in self.Binv.list())
        self.table=[[list(self.coords((a*b)%f)) for b in basis] for a in basis]
        assert all(v in ZZ for row in self.table for col in row for v in col)
        self.mults=[matrix(ZZ,self.table[j]).transpose() for j in range(3)]
        self.one=(1,0,0);self.cache={};self.locals={};self.unit_squares_8=None
    def polycoords(self,a):return list(a)+[QQ(0)]*(3-len(a.list()))
    def coords(self,a):return self.Binv*vector(QQ,self.polycoords(a%self.f))
    def mul(self,a,b,m):
        if m not in self.cache:
            self.cache[m]=[[[residue(v,m) for v in col] for col in row] for row in self.table]
        tab=self.cache[m];out=[0,0,0]
        for i,u in enumerate(a):
            if not u:continue
            for j,v in enumerate(b):
                if not v:continue
                for k,w in enumerate(tab[i][j]):out[k]+=int(u)*int(v)*w
        return tuple(a%m for a in out)
    def power(self,a,n,p):
        result=self.one
        while n:
            if n&1:result=self.mul(result,a,p)
            a=self.mul(a,a,p);n//=2
        return result
    def finite_mult(self,a,p):
        return matrix(GF(p),[self.mul(a,tuple(int(i==j) for i in range(3)),p) for j in range(3)]).transpose()
    def local_order(self,p):
        F=matrix(GF(p),[self.power(tuple(int(i==j) for i in range(3)),p,p) for j in range(3)]).transpose()
        radical=(F*F if p==2 else F).right_kernel().basis_matrix()
        generators=(p*identity_matrix(ZZ,3)).stack(radical.change_ring(ZZ))
        H=generators.hermite_form(include_zero_rows=False).transpose()
        assert H.nrows()==H.ncols()==3
        Hi=H.inverse()
        conjugates=[Hi*M*H for M in self.mults]
        assert all(a in ZZ for M in conjugates for a in M.list())
        constraints=matrix(GF(p),[[conjugates[j].list()[i] for j in range(3)] for i in range(9)])
        assert constraints.rank()==3,'P_MAXIMALITY_NOT_CERTIFIED'
        components=3-(F-identity_matrix(GF(p),3)).rank()
        self.locals[p]=(F,components)
        return {'place':p,'p_maximal_certified':True,'radical_dimension':radical.nrows(),
                'multiplier_constraint_rank':int(constraints.rank()),'prime_count':components}
    def square_test(self,value,certificate,p):
        scale=self.R(list(map(QQ,certificate['scale'])))
        assert scale
        u=(value*scale.inverse_mod(self.f)**2)%self.f
        c=self.coords(u)
        assert all(a.denominator()%p for a in c),'NORMALIZATION_NOT_P_INTEGRAL'
        norm=self.f.resultant(u)
        nv=int(norm.valuation(p));assert nv>=0
        residues=tuple(residue(a,p) for a in c)
        F,components=self.locals[p]
        if nv:
            inv=self.coords((p*u.inverse_mod(self.f))%self.f)
            assert all(a.denominator()%p for a in inv),'PRINCIPAL_IDEAL_DOES_NOT_CONTAIN_P'
            Mu=self.finite_mult(residues,p)
            q=3-Mu.rank()
            assert q==nv and q>0
            assert Mu.augment(F).rank()-Mu.rank()==q,'PRINCIPAL_QUOTIENT_NOT_REDUCED'
            square=False;method='nontrivial_radical_principal_ideal'
        elif p==2:
            if self.unit_squares_8 is None:
                self.unit_squares_8={self.mul(v,v,8) for v in product(range(4),repeat=3)}
            square=tuple(residue(a,8) for a in c) in self.unit_squares_8
            method='all_64_roots_mod4_squared_mod8'
        else:
            twist=self.finite_mult(self.power(residues,(p-1)//2,p),p)*F
            square=3-(twist-identity_matrix(GF(p),3)).rank()==components
            method='frobenius_fixed_space_of_quadratic_unit_algebra'
        assert square==certificate['expected_square'],(p,'SQUARE_TEST_MISMATCH',method)
        return method

def verify(index):
    generic_path=DIR/f'case-{index:02d}-generic.json'
    generic=read(generic_path);provenance(generic)
    protocol=read(DIR/'protocol.json');provenance(protocol)
    inp=read(DIR/f'case-{index:02d}-input.json');provenance(inp)
    assert generic['case']==protocol['cases'][index]==inp['case']
    parent=read(ART/'curve302_recovered_mw17_parent_v1.json')
    arith=read(ART/'rank_jump_curve302_strict_constructor_arithmetic_inputs_v1.json')
    primes=[int(p) for p,e in arith['discriminant_factors']]
    assert protocol['places']==primes+['infinity'] and inp['places']==primes
    assert all(ZZ(p).is_prime(proof=True) for p in primes)
    tau=QQ(inp['case']['parameter']);T=PolynomialRing(QQ,'t')
    def ev(r):return T(r['numerator'])(tau)/T(r['denominator'])(tau)
    ai=list(map(ev,parent['a_invariants']));assert ai[:3]==[1,1,1]
    d=tau.denominator()
    d*=lcm([((16*ai[3]+8)*d**8).denominator(),((64*ai[4]+16)*d**12).denominator()])
    assert str(d)==inp['scale_denominator']
    R=PolynomialRing(QQ,'x');x=R.gen()
    f=x**3+5*d**4*x*x+(16*ai[3]+8)*d**8*x+(64*ai[4]+16)*d**12
    assert f.discriminant() and all(a in ZZ for a in f)
    assert list(map(str,f.list()))==inp['polynomial']==generic['polynomial']
    coordinates=[]
    for P in parent['basis_weierstrass_coordinates']:
        xp,yp=map(ev,P);X,Y=4*d**4*xp,d**6*(8*yp+4*xp+4)
        assert f(X)==Y*Y and Y
        coordinates.append([str(X),str(Y)])
    assert coordinates==inp['generic_cubic_points']
    betas=[QQ(X)-x for X,Y in coordinates]
    algebra=Algebra(f,[R(list(map(QQ,b))) for b in generic['local_order_basis']])
    local_reports=[];matrices=[];counts={}
    def verify_square(value,certificate,p):
        method=algebra.square_test(value,certificate,p)
        counts[method]=counts.get(method,0)+1
    for row,p in zip(generic['local'][:-1],primes):
        assert row['place']==p
        order=algebra.local_order(p)
        rank=row['generic_rank'];pivots=row['basis_indices'];cc=row['generic_coordinates']
        assert len(pivots)==rank and len(set(pivots))==rank and all(0<=j<17 for j in pivots)
        assert len(cc)==17 and all(len(c)==rank for c in cc)
        C=matrix(GF(2),17,rank,[a for c in cc for a in c]).transpose()
        assert C.rank()==rank
        assert C.matrix_from_columns(pivots)==identity_matrix(GF(2),rank)
        for j,beta in enumerate(betas):
            divisor=R(1)
            for k,a in zip(pivots,cc[j]):
                assert a in [0,1]
                if a:divisor=(divisor*betas[k])%f
            value=(beta*divisor.inverse_mod(f))%f
            cert=row['relation_square_tests'][j];assert cert['expected_square'] is True
            verify_square(value,cert,p)
        assert [r['mask'] for r in row['basis_nonsquare_tests']]==list(range(1,1<<rank))
        for cert in row['basis_nonsquare_tests']:
            value=R(1)
            for j,k in enumerate(pivots):
                if cert['mask']>>j&1:value=(value*betas[k])%f
            assert cert['expected_square'] is False
            verify_square(value,cert,p)
        full=order['prime_count']-1+(p==2)
        assert row['prime_count']==order['prime_count'] and full==row['full_local_point_dimension']
        local_reports.append({**order,'generic_rank':rank,'full_local_point_dimension':full})
        matrices.append(C)
    real=generic['local'][-1];assert real['place']=='infinity'
    signs=[signs_at_roots(f,QQ(X)) for X,Y in coordinates]
    assert signs==real['root_signs']
    M=matrix(GF(2),signs).transpose();pivots=real['basis_indices']
    C=matrix(GF(2),17,len(pivots),[a for c in real['generic_coordinates'] for a in c]).transpose()
    assert M.matrix_from_columns(pivots)*C==M and C.rank()==M.rank()==real['generic_rank']
    assert real['full_local_point_dimension']==int(len(signs[0])==3)
    local_reports.append({'place':'infinity','generic_rank':M.rank(),'full_local_point_dimension':int(len(signs[0])==3)})
    matrices.append(C);H=matrices[0]
    for M in matrices[1:]:H=H.stack(M)
    assert [list(map(int,row)) for row in H.rows()]==generic['generic_matrix']
    checks=H.left_kernel().basis_matrix()
    assert [list(map(int,row)) for row in checks.rows()]==generic['compatibility_checks']
    assert H.rank()==generic['generic_joint_rank']
    assert H.nrows()==generic['product_of_generic_local_images_dimension']
    full=all(r['generic_rank']==r['full_local_point_dimension'] for r in local_reports)
    assert full==generic['all_generic_local_images_full']
    result={'status':'PASS_INDEPENDENT_GENERIC_LOCAL_CODE','classification':'verified application',
        'case':generic['case'],'local':local_reports,'generic_joint_rank':H.rank(),
        'product_dimension':H.nrows(),'compatibility_check_dimension':checks.nrows(),
        'all_generic_local_images_full':full,'square_test_counts':counts,
        'limitations':'Finite-place compatibility only; neither a Selmer dimension nor rational solubility of complementary classes.',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [generic_path,DIR/f'case-{index:02d}-input.json',Path(__file__)]}}
    # This independently verified generic code precedes all first-seed evaluation.
    retain(DIR/f'case-{index:02d}-replay.json',result)
    if index==8:
        old=ART/'det1092_seed_local_code_v3'
        evaluation=read(old/'first-seed-evaluation.json');provenance(evaluation)
        old_generic=read(old/'case-08-generic.json')
        assert old_generic['generic_matrix']==generic['generic_matrix']
        assert old_generic['local']==generic['local']
        source=read(ART/'det1092_first_centre_rr_net_replay_v1.json')
        xp,yp=map(QQ,source['reconstructed_point_literal302'])
        X,Y=4*xp,8*yp+4*xp+4
        assert d==1 and f(X)==Y*Y
        assert f.change_ring(GF(31)).is_irreducible(),'NO_RATIONAL_2_TORSION_NOT_CERTIFIED'
        assert evaluation['point_cubic']==[str(X),str(Y)]
        seed=X-x;coords=[]
        for row,cert in zip(generic['local'][:-1],evaluation['local_square_tests']):
            assert cert['place']==row['place'] and cert['expected_square'] is True
            cc=cert['coordinates'];assert len(cc)==row['generic_rank']
            divisor=R(1)
            for k,a in zip(row['basis_indices'],cc):
                assert a in [0,1]
                if a:divisor=(divisor*betas[k])%f
            verify_square((seed*divisor.inverse_mod(f))%f,cert,row['place'])
            coords.extend(cc)
        sr=signs_at_roots(f,X);assert sr==evaluation['real_signs']
        C=matrix(GF(2),signs).transpose().matrix_from_columns(real['basis_indices'])
        coords.extend(map(int,C.solve_right(vector(GF(2),sr))))
        assert coords==evaluation['coordinates']
        col=vector(GF(2),coords);syndrome=list(map(int,checks*col))
        assert syndrome==evaluation['syndrome'] and any(syndrome)
        assert checks.nrows()==5 and H.augment(col).rank()==18
        # Rank17 of17 Kummer columns makes the generic subgroup2-saturated:
        # its rational saturation has odd index (the old torsion-free proof applies).
        # A nonzero syndrome of a rational point therefore proves non-generic rank.
        offsets=[];start=0
        for row in generic['local']:
            offsets.append((row['place'],start,start+row['generic_rank']))
            start+=row['generic_rank']
        separators=[]
        for mask in range(1,32):
            word=sum((checks[j] for j in range(5) if mask>>j&1),checks[0]*0)
            if word*col:
                support=[p for p,lo,hi in offsets if any(word[lo:hi])]
                separators.append({'mask':mask,'word':list(map(int,word)),
                                   'places':support,'place_count':len(support)})
        separators.sort(key=lambda r:(r['place_count'],r['mask']))
        assert separators==evaluation['all_separating_checks']
        assert separators[0]==evaluation['minimal_separator'] and len(separators)==16
        retain(DIR/'first-seed-replay.json',{'status':'PASS_INDEPENDENT_FIRST_SEED_LOCAL_COMPATIBILITY_DEFECT',
            'classification':'verified application and new diagnostic deduction',
            'generic_joint_rank':17,'with_first_seed_rank':18,'syndrome':syndrome,
            'smallest_separating_place_count_in_fixed_set':separators[0]['place_count'],
            'smallest_separator':separators[0],'nonzero_checks_tested':31,
            'input_boundary':'All generic code inputs are equations and17 generic sections only. The historical first point is evaluation only; no later points or V3 inputs.',
            'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [old/'first-seed-evaluation.json',DIR/'case-08-replay.json',Path(__file__)]}})
    print('PASS_INDEPENDENT_LOCAL_CODE',index,H.nrows(),H.rank(),checks.nrows(),full,flush=True)
if __name__=='__main__':
    signal.alarm(25)
    parser=argparse.ArgumentParser();parser.add_argument('--case',type=int,choices=range(9),required=True)
    verify(parser.parse_args().case)

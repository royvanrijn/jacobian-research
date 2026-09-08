#!/usr/bin/env sage-python
"""Independent local image proof, without PARI number-field/unit routines.

Certify2-maximality by the radical multiplier criterion; enumerate squares
of units modulo8; verify p-adic normalization suggestions by exact finite
ring congruences; reject every nonzero generator combination modulo all
eight rational squareclasses. No idealstar/ideallog is used in this replay.
"""
import argparse,hashlib,json,signal
from itertools import product
from pathlib import Path
from sage.all import QQ,ZZ,GF,Qp,PolynomialRing,matrix,vector,identity_matrix
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_rr_dyadic_images_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def residue(a,m):
    a=QQ(a);return int(a.numerator()*a.denominator().inverse_mod(m)%m)
def bits_rank(columns):
    pivots={}
    for values in columns:
        n=sum(int(b)<<j for j,b in enumerate(values))
        while n:
            k=n.bit_length()-1
            if k in pivots:n^=pivots[k]
            else:pivots[k]=n;break
    return len(pivots)
class LocalRing:
    def __init__(self,g,basis):
        self.g=g;self.basis=basis;self.R=g.parent();self.n=6
        def coordinates(p):return list(p)+[QQ(0)]*(6-len(p.list()))
        self.polycoords=coordinates
        self.B=matrix(QQ,[coordinates(b) for b in basis]).transpose();self.Binv=self.B.inverse()
        assert basis[0]==1 and self.B.det()
        assert all(a.valuation(2)>=0 for a in self.Binv.list() if a)
        self.table=[[[QQ(v) for v in self.Binv*vector(QQ,coordinates((a*b)%g))] for b in basis] for a in basis]
        assert all(a.valuation(2)>=0 for row in self.table for col in row for a in col if a)
        self.one=(1,0,0,0,0,0);self.zero=(0,)*6;self.cached={}
        self.units={}
        for v in product(range(2),repeat=6):
            cols=[self.mul(v,tuple(int(i==j) for i in range(6)),2) for j in range(6)]
            self.units[v]=bits_rank(cols)==6
        self.unit_order_mod8=sum(self.units.values())*4**6
    def mul(self,a,b,m):
        if m not in self.cached:
            self.cached[m]=[[[residue(v,m) for v in col] for col in row] for row in self.table]
        t=self.cached[m];out=[0]*6
        for i,u in enumerate(a):
            if not u:continue
            for j,v in enumerate(b):
                if not v:continue
                c=u*v
                for k,w in enumerate(t[i][j]):out[k]+=c*w
        return tuple(x%m for x in out)
    def power(self,a,n,m):
        if n<0:
            assert m==8 and self.units[tuple(x%2 for x in a)]
            a=self.power(a,self.unit_order_mod8-1,8);n=-n
        result=self.one
        while n:
            if n&1:result=self.mul(result,a,m)
            a=self.mul(a,a,m);n//=2
        return result
    def max_order_and_uniformizers(self):
        F=matrix(GF(2),[[residue(self.table[j][j][i],2) for j in range(6)] for i in range(6)])
        radical=(F**3).right_kernel();rb=[tuple(map(int,v)) for v in radical.basis()]
        rset={tuple(map(int,v)) for v in radical};passing=[]
        for y in product(range(2),repeat=6):
            if y not in rset:continue
            good=all(self.mul(y,tuple(int(k==j) for k in range(6)),2) in rset for j in range(6))
            for r in rb:
                z=self.mul(y,r,4)
                if any(a%2 for a in z) or tuple(a//2 for a in z) not in rset:good=False;break
            if good:passing.append(y)
        assert passing==[self.zero], 'CANDIDATE_ORDER_NOT_CERTIFIED_2_MAXIMAL'
        idempotents=[tuple(map(int,v)) for v in (F-identity_matrix(GF(2),6)).right_kernel()]
        primitive=[e for e in idempotents if e!=self.zero and not any(
            f not in [e,self.zero] and self.mul(e,f,2)==f for f in idempotents)]
        primitive.sort();assert tuple(sum(e[j] for e in primitive)%2 for j in range(6))==self.one
        records=[];uniformizers=[]
        for e in primitive:
            comp=[self.mul(e,tuple(int(k==j) for k in range(6)),2) for j in range(6)]
            nil=[self.mul(e,r,2) for r in rb]
            n=bits_rank(comp);f=n-bits_rank(nil);assert f>0 and n%f==0
            lift=e
            for step in range(3):
                sq=self.mul(lift,lift,256);cube=self.mul(sq,lift,256)
                lift=tuple((3*a-2*b)%256 for a,b in zip(sq,cube))
            assert self.mul(lift,lift,256)==lift
            candidate=None;candidate_index=None
            for j,r in enumerate([(2,0,0,0,0,0)]+rb):
                er=self.mul(lift,r,256)
                pi=tuple((a-b+c)%256 for a,b,c in zip(self.one,lift,er))
                columns=[self.mul(pi,tuple(int(k==j) for k in range(6)),256) for j in range(6)]
                norm=matrix(ZZ,columns).det()%256
                if norm and ZZ(norm).valuation(2)==f:
                    candidate=pi;candidate_index=j;break
            assert candidate is not None,'NO_BOUNDED_UNIFORMIZER'
            uniformizers.append(candidate)
            records.append({'idempotent_mod2':e,'idempotent_mod256':lift,
                            'local_degree':n,'residue_degree':f,'ramification_degree':n//f,
                            'uniformizer_basis_coordinates':candidate,'uniformizer_candidate_index':candidate_index,
                            'uniformizer_norm_v2':f})
        self.uniformizers=uniformizers;self.local_records=records
        self.padics=Qp(2,prec=512);self.inverse_matrices=[]
        for pi in uniformizers:
            columns=[]
            for j in range(6):
                columns.append([sum(QQ(pi[i])*self.table[i][j][k] for i in range(6)) for k in range(6)])
            self.inverse_matrices.append(matrix(self.padics,columns).transpose().inverse())
        return {'radical_basis_mod2':rb,'multiplier_candidates_tested':64,
                'passing_multiplier_classes':passing,'primitive_idempotents':records,
                'two_maximal_order_certified':True}
    def normalize(self,P):
        P=P%self.g;norm_v=int(self.g.resultant(P).valuation(2));assert 0<=norm_v<=128
        exact=self.Binv*vector(QQ,self.polycoords(P));assert all(a.valuation(2)>=0 for a in exact if a)
        c=vector(self.padics,exact);valuations=[];total=0
        for info,inverse in zip(self.local_records,self.inverse_matrices):
            count=0
            while True:
                nxt=inverse*c
                assert all(a.precision_absolute()>0 for a in nxt)
                if any(a.valuation()<0 for a in nxt):break
                c=nxt;count+=1;total+=info['residue_degree']
                assert total<=norm_v,'NORMALIZATION_VALUATION_CAP'
            valuations.append(count)
        assert total==norm_v
        L=norm_v+3;m=2**L
        assert all(a.precision_absolute()>=L and a.valuation()>=0 for a in c)
        unit=tuple(int(a.lift())%m for a in c)
        assert self.units[tuple(a%2 for a in unit)]
        principal=self.one
        for pi,e in zip(self.uniformizers,valuations):principal=self.mul(principal,self.power(pi,e,m),m)
        # This exact congruence, together with v2(Norm(principal))=norm_v,
        # certifies the unit part modulo8 independently of Qp arithmetic.
        assert self.mul(principal,unit,m)==tuple(residue(a,m) for a in exact)
        return {'valuations':valuations,'norm_v2':norm_v,'verification_precision':L,
                'unit_basis_residues':unit,'unit_mod8':tuple(a%8 for a in unit),
                'exact_normalization_congruence_verified':True}
def verify(i):
    path=OUT/('case-%02d.json'%i);d=read(path);protocol=read(OUT/'protocol.json')
    for item in [d,protocol]:
        for p,h in item['inputs'].items():assert sha(ROOT/p)==h
    prior=ART/'det1092_rr_real_dyadic_panel_replay_v1.json';dimensions=read(prior)['cases'][i]
    assert d['dimensions']==dimensions
    R=PolynomialRing(QQ,'x');x=R.gen();source=read(ROOT/d['source']);q=R(source['q']).monic()
    s=d['T_equals_two_power_times_x'];f=R(q(2**s*x)/2**(6*s));g=R(d['proxy_polynomial'])
    N=d['coefficient_precision'];D=d['integral_monic_discriminant_v2']
    assert N==512 and f.discriminant().valuation(2)==g.discriminant().valuation(2)==D and N>2*D
    assert f.is_monic() and g.is_monic() and all(a.valuation(2)>=0 for a in f if a)
    assert all((a-b).valuation(2)>=N for a,b in zip(f,g) if a!=b)
    p=d['proxy_irreducibility_prime'];assert ZZ(p).is_prime()
    assert PolynomialRing(GF(p),'x')(g.list()).is_irreducible()
    x0=QQ(source['base_x'])
    def evaluation_check(P,record):
        v=record['input_coefficient_min_v2'];integral=R(P/2**v)
        assert all(a.valuation(2)>=0 for a in integral if a) and any(a.valuation(2)==0 for a in integral if a)
        small=R(record['normalized_residue_polynomial']);assert integral.degree()==small.degree()
        assert all((a-b).valuation(2)>=N for a,b in zip(integral,small) if a!=b)
        vn=int(g.resultant(small).valuation(2));assert vn==record['norm_v2']>=0
        assert N-D>vn+2 and record['root_matching_lower_bound']==N-D
        return small
    anchor=evaluation_check(R([x0,-2**s]),d['anchor']);polynomials=[]
    for j,(saved,row) in enumerate(zip(d['generic_evaluations'],source['generic_divisors'])):
        original=R(row);assert j==saved['generic_divisor_index'] and saved['degree']==original.degree()
        polynomials.append(evaluation_check(R(original(2**s*x)),saved))
    assert len(polynomials)==17
    ring=LocalRing(g,[R(row) for row in d['two_maximal_order_basis']]);order=ring.max_order_and_uniformizers()
    assert sorted((r['ramification_degree'],r['residue_degree']) for r in ring.local_records)==sorted((r['e'],r['f']) for r in d['prime_data'])
    chosen=d['independent_generic_divisor_indices'];target=dimensions['Q2_fake_Kummer_dimension']
    assert len(chosen)==len(set(chosen))==target and all(0<=j<17 for j in chosen)
    normalization={'anchor':ring.normalize(anchor)}
    normalization['scalars']=[ring.normalize(R(a)) for a in [-1,2,5]]
    normalization['generic']=[{'index':j,**ring.normalize(polynomials[j])} for j in chosen]
    # Complete unit squares modulo8: squaring is unchanged by adding4O.
    squares=set();units_tested=0
    for a in product(range(4),repeat=6):
        if ring.units[tuple(v%2 for v in a)]:
            units_tested+=1;squares.add(ring.mul(a,a,8))
    assert ring.one in squares and units_tested*4**6//2**6==ring.unit_order_mod8
    classes=[];base=normalization['anchor']
    for item in normalization['generic']:
        degree=d['generic_evaluations'][item['index']]['degree']
        vals=tuple((a-degree*b)%2 for a,b in zip(item['valuations'],base['valuations']))
        unit=ring.mul(item['unit_mod8'],ring.power(base['unit_mod8'],-degree,8),8)
        if degree%2:unit=tuple(-a%8 for a in unit)
        classes.append((vals,unit))
    scalars=[(tuple(v%2 for v in a['valuations']),a['unit_mod8']) for a in normalization['scalars']]
    def multiply_classes(items,mask):
        vals=[0]*len(ring.uniformizers);unit=ring.one
        for j,(v,u) in enumerate(items):
            if (mask>>j)&1:
                vals=[a^b for a,b in zip(vals,v)];unit=ring.mul(unit,u,8)
        return vals,unit
    scalar_values=[multiply_classes(scalars,mask) for mask in range(8)];tests=[]
    for mask in range(1,2**target):
        v,u=multiply_classes(classes,mask)
        for sm,(w,z) in enumerate(scalar_values):
            parity=[a^b for a,b in zip(v,w)];unit=ring.mul(u,z,8)
            if any(parity):reason='ODD_LOCAL_VALUATION'
            else:
                assert unit not in squares,'DEPENDENT_GENERIC_CLASSES_MODULO_RATIONAL_SCALARS'
                reason='UNIT_NOT_A_SQUARE_MOD_EIGHT'
            tests.append({'generic_mask':mask,'rational_scalar_mask':sm,'valuation_parity':parity,
                          'unit_mod8':unit,'nonsquare_reason':reason})
    kernel=dimensions['Q2_true_Kummer_dimension']-target;assert kernel in [0,1]
    result={'classification':'verified application and independent exact finite-ring replay',
            'case_index':i,'status':'PASS_INDEPENDENT_COMPLETE_DYADIC_KUMMER_IMAGE',
            'generic_fake_image_rank':target,'generic_true_image_rank_including_D0':target+kernel,
            'independent_generic_divisor_indices':chosen,'include_inherited_D0':bool(kernel),
            'order_certificate':order,'normalizations':normalization,
            'unit_square_residues_mod8':sorted(squares),'unit_residues_mod4_tested':units_tested,
            'nonsquare_tests':tests,
            'scope':'Complete true and fake2-adic images are generated by the inherited subgroup. Proxy squareclass transport is certified for these finitely many divisors. No global Selmer group is computed.',
            'limits':{'wall_seconds':25,'proxy_precision':512,'multiplier_candidates':64,
                      'residues_mod4':4096,'nonsquare_tests':248,'normalization_norm_v2_cap':128,
                      'PARI_number_field_calls':0,'class_groups':0,'point_searches':0,'pilot_changes':0},
            'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [path,OUT/'protocol.json',prior]},
            'checker_sha256':sha(Path(__file__))}
    retain(OUT/('case-%02d-finite-replay.json'%i),result)
    print('case',i,result['status'],'true rank',target+kernel,'square residues',len(squares),'tests',len(tests),flush=True)
    return result
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--case',type=int,required=True,choices=range(10));args=ap.parse_args()
    signal.alarm(25);verify(args.case)

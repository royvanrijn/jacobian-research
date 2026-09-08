#!/usr/bin/env sage-python
"""Extend the retained Kummer diagnostic to all17 generic section divisors.

Same curve and same64 primes. The extra fibre-divisor class has trivial fake
Kummer image and is handled separately, not silently omitted from the group.
"""
import hashlib, json, signal
from pathlib import Path
from sage.all import QQ, ZZ, GF, PolynomialRing, EllipticCurve, matrix, vector, power_mod
ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL = ROOT/'artifacts/local/elliptic-curves/det1092-rr-full-inherited-jacobian-v1'
OUT = ART/'det1092_rr_full_inherited_jacobian_v1.json'
PROTOCOL = ART/'det1092_rr_full_inherited_jacobian_protocol_v1.json'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def retain(p,d):
    s=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists(): assert p.read_text()==s
    else: p.write_text(s)
def build():
    paths=[ART/n for n in ['curve302_recovered_mw17_parent_v1.json','det1092_first_centre_rr_net_v1.json',
           'det1092_rr_residual_jacobian_class_v1.json','det1092_rr_residual_jacobian_protocol_v1.json']]
    parent,net,diag,oldprotocol=[json.loads(p.read_text()) for p in paths]
    protocol={'classification':'frozen full-inherited-subgroup retrospective diagnostic',
              'curve_choice':'Same u=u0,v=0 as the retained diagnostic; retrospective, not blind.',
              'generic_generators':'Restrictions of all17 generic sections, plus the fibre divisor; O restricts to the base point.',
              'primes':oldprotocol['primes'],
              'limits':{'wall_seconds':25,'generic_section_divisors':17,'primes':64,
                        'point_searches':0,'global_Selmer_runs':0,'class_group_runs':0,'pilot_changes':0},
              'script_sha256':sha(Path(__file__))}
    retain(PROTOCOL,protocol)
    R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field()
    def rat(r):return K(R(r['numerator']))/R(r['denominator'])
    E=EllipticCurve(K,[rat(r) for r in parent['a_invariants']])
    basis=[E([rat(r) for r in P]) for P in parent['basis_weierstrass_coordinates']]
    w=-vector(QQ,net['trace_word']);G=matrix(QQ,parent['generic_height_gram'])
    C=sum((n*P for n,P in zip(w,basis)),E(0))
    h=R(C[0].denominator().sqrt());cx=C[0]+E.b2()/12
    A,B=[[R(r) for r in net[k]] for k in ['A','B']]
    # Freeze the complete generic rational functions before using u0.
    generic=[]
    for j,P in enumerate(basis):
        r=-(B[0]+B[1]*P[0]+B[2]*P[1])/(A[0]+A[1]*P[0]+A[2]*P[1])
        degree=ZZ(1+G[j,j]-(w*G)[j])
        assert max(r.numerator().degree(),r.denominator().degree())==degree
        generic.append({'basis_index':j,'degree':int(degree),
                        'numerator':list(map(str,r.numerator().list())),
                        'denominator':list(map(str,r.denominator().list()))})
        retain(LOCAL/f'generic-divisor-{j:02d}.json',generic[-1])
    u0=QQ(diag['curve']['u']);q=R(diag['curve']['q']);c=QQ(diag['curve']['scale'])
    x0=QQ(diag['inherited_x_values'][diag['base_pair_index']])
    f1,f2=B[1]+u0*A[1],B[2]+u0*A[2];m=-f1+E.a1()*f2/2
    polynomials=[];divisors=[]
    for j,record in enumerate(generic):
        g=(R(record['numerator'])-u0*R(record['denominator'])).monic()
        assert g.degree()==record['degree'] and g.gcd(q)==1 and g.gcd(g.derivative())==1
        s=((2*(basis[j][0]+E.b2()/12)+cx)*f2**2-m**2)/h**3
        assert s.denominator().gcd(g)==1
        y=(s.numerator()*s.denominator().inverse_mod(g))%g
        assert (y*y-c*q)%g==0
        polynomials.append(g)
        divisors.append({'basis_index':j,'degree':int(g.degree()),'g':list(map(str,g.list())),
                         's_mod_g':list(map(str,y.list()))})
        retain(LOCAL/f'special-divisor-{j:02d}.json',divisors[-1])
    rows=[];trials=[]
    for p in protocol['primes']:
        r={'p':p}
        coeffs=list(q)+[x0]+[a for g in polynomials for a in g]
        if any(a.denominator()%p==0 for a in coeffs):r['status']='SKIP_DENOMINATOR'
        else:
            F=GF(p);S=PolynomialRing(F,'X');X=S.gen();f=S(q.list())
            gs=[S(g.list()) for g in polynomials]
            if f.degree()!=6 or f.gcd(f.derivative())!=1:r['status']='SKIP_BAD_SEXTIC'
            elif f(F(x0))==0 or f(0)==0 or any(g.gcd(f)!=1 for g in gs):r['status']='SKIP_NONUNIT_DIVISOR'
            else:
                assert c.valuation(p)%2==0
                factors=sorted([a.monic() for a,e in f.factor()],key=lambda a:(a.degree(),list(map(int,a))))
                bits=[]
                for a in factors:
                    exp=(ZZ(p)**a.degree()-1)//2
                    anchor=power_mod(S(F(x0))-X,exp,a)
                    vals=[(((-1)**g.degree())*g)%a for g in gs]+[-X]
                    degs=[g.degree() for g in gs]+[1]
                    row=[]
                    for value,degree in zip(vals,degs):
                        sign=(power_mod(value,exp,a)*anchor**degree)%a
                        assert sign in [S(1),S(-1)]
                        row.append(int(sign==S(-1)))
                    bits.append(row)
                raw=matrix(GF(2),bits)
                assert all(sum(raw[:,j].list(),GF(2)(0))==0 for j in range(18))
                quotient=matrix(GF(2),1,len(factors),[a.degree()%2 for a in factors]).right_kernel().basis_matrix()
                block=quotient*raw;rows.extend([list(map(int,z)) for z in block.rows()])
                r.update({'status':'PASS_LOCAL_KUMMER_BLOCK','factors':[list(map(int,a.list())) for a in factors],
                          'raw_rows':[list(map(int,z)) for z in raw.rows()],
                          'quotient_rows':[list(map(int,z)) for z in quotient.rows()],
                          'block_rows':[list(map(int,z)) for z in block.rows()]})
        trials.append(r);retain(LOCAL/f'prime-{p:04d}.json',r)
    M=matrix(GF(2),rows);H=M[:,:-1]
    z=next((z for z in H.left_kernel().basis() if z*M[:,-1]),None)
    result={'classification':'retrospective verified application; independent replay required',
            'status':'PASS_FULL_INHERITED_FAKE_KUMMER_SEPARATION' if z is not None else 'INCONCLUSIVE_FULL_INHERITED_CLASS_NOT_SEPARATED',
            'generic_r_functions':generic,'divisors':divisors,'trials':trials,'matrix_rows':rows,
            'inherited_character_rank':int(H.rank()),'total_character_rank':int(M.rank()),
            'separator':None if z is None else list(map(int,z)),
            'fibre_class':'D0=[K_C-2P0]; its fake Kummer image is zero.',
            'theta_divisibility_gate':'A rational half of D0 is equivalent to a rational theta characteristic. Degree6 irreducibility excludes odd theta characteristics; a1+5 factorization excludes invariant3+3 partitions and hence even theta characteristics.',
            'rank_claim':'NOT_YET_PROMOTED; replay and exact inherited-lattice/theta arguments required.',
            'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+[PROTOCOL,Path(__file__)]}}
    retain(OUT,result)
    print(result['status'],'ranks',H.rank(),M.rank(),flush=True)
if __name__=='__main__':
    signal.alarm(25);LOCAL.mkdir(parents=True,exist_ok=True);build()

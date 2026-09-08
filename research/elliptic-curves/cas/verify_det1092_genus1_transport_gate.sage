#!/usr/bin/env sage-python
"""Independent integer replay of the ten fixed genus1 transport obstructions.

No constructor import, finite-field point counter, factorization, number
field, isogeny or point search. Every old prime is replayed;25-second cap.
"""
import hashlib,json,signal,itertools,math
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,identity_matrix
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
DIR=ART/'det1092_genus1_transport_gate_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def provenance(d):
    for path,digest in d['inputs'].items():assert sha(ROOT/path)==digest,(path,'HASH_MISMATCH')
def retain(p,d):
    s=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==s,('IMMUTABLE_MISMATCH',str(p))
    else:p.write_text(s)
def mod(v,p):
    v=QQ(v);assert v.denominator()%p
    return int(v.numerator())*pow(int(v.denominator()),-1,p)%p
def disc(a,b,c):return a*a*b*b-4*b**3-4*a**3*c-27*c*c+18*a*b*c
def matrix_key(M):return tuple(map(int,M.list()))
def verify():
    protocol=read(DIR/'protocol.json');panel=read(DIR/'panel.json')
    provenance(protocol);provenance(panel)
    parentpath=ART/'curve302_recovered_mw17_parent_v1.json'
    oldpath=ART/'det1092_genus1_picard_controls_v1/protocol.json'
    selectionpath=ART/'det1092_genus1_picard_controls_v1/selection.json'
    parent,old,selection=map(read,[parentpath,oldpath,selectionpath])
    assert protocol['primes']==old['primes'] and len(old['primes'])==64
    assert protocol['limits']['pairs']==10 and protocol['limits']['exceptional_coordinate_inputs']==0
    R=PolynomialRing(QQ,'t');t=R.gen();K=R.fraction_field()
    def dec(d):return K(R(d['numerator']))/R(d['denominator'])
    ai=list(map(dec,parent['a_invariants']));assert ai[:3]==[1,1,1]
    checks=[];paths=[DIR/'protocol.json',DIR/'panel.json',parentpath,oldpath,selectionpath,Path(__file__)]
    for i in range(10):
        eqpath=DIR/f'case-{i:02d}-equations.json';outpath=DIR/f'case-{i:02d}.json'
        equations,result=map(read,[eqpath,outpath]);provenance(result);paths +=[eqpath,outpath]
        sourcepath=ROOT/equations['source'];source=read(sourcepath);paths.append(sourcepath)
        if i==0:
            assert sourcepath==ART/'det1092_genus1_picard_image_v1/case-01-generic.json'
            assert equations['tau']=='0' and equations['label']=='302-first-carrier'
        else:
            sel=selection['cases'][i-1]
            assert sourcepath==ART/f'det1092_genus1_picard_controls_v1/case-{i-1:02d}-generic.json'
            assert equations['label']==sel['label']==source['label']
            assert equations['tau']==sel['tau']==source['tau'] and source['z']==sel['z']
        tau=QQ(equations['tau']);q=R(source['quartic']);t0,s=map(QQ,source['basepoint'])
        assert q.degree()==4 and q.gcd(q.derivative())==1 and q(t0)==s*s and s
        q0,q1,q2,q3,q4=q(t+t0).list()
        J=[q2,q1*q3-4*s*s*q4,s*s*q3*q3+q1*q1*q4-4*s*s*q2*q4]
        E=[QQ(5),16*ai[3](tau)+8,64*ai[4](tau)+16]
        assert J==list(map(QQ,equations['carrier_abc'])) and E==list(map(QQ,equations['original_abc']))
        assert disc(*J) and disc(*E)
        witnesses=dict.fromkeys(['carrier_irreducible','carrier_transposition','original_irreducible',
          'original_transposition','distinct_quadratic_fields','geometric_nonisogeny'])
        assert [r['p'] for r in result['trials']]==protocol['primes']
        good=0
        for row in result['trials']:
            p=row['p'];assert ZZ(p).is_prime()
            if any(v.denominator()%p==0 for abc in [J,E] for v in abc):
                assert row['status']=='SKIP_MODEL_DENOMINATOR';continue
            abc=[[mod(v,p) for v in model] for model in [J,E]]
            ds=[disc(*v)%p for v in abc]
            assert row['abc']==abc and row['discriminants_mod_p']==ds
            if 0 in ds:
                assert row['status']=='SKIP_BAD_REDUCTION';continue
            assert row['status']=='COMPLETE_GOOD_PAIR';good+=1
            squares={y*y%p for y in range(1,p)};types=[];counts=[];traces=[]
            for a,b,c in abc:
                vals=[((x+a)*x*x+b*x+c)%p for x in range(p)]
                roots=vals.count(0);assert roots in [0,1,3]
                typ={0:[3],1:[1,2],3:[1,1,1]}[roots];types.append(typ)
                count=1+sum(1 if v==0 else 2 if v in squares else 0 for v in vals)
                counts.append(count);traces.append(p+1-count)
            assert row['factor_types']==types and row['point_counts']==counts and row['traces']==traces
            for name,typ in zip(['carrier','original'],types):
                key=name+('_irreducible' if typ==[3] else '_transposition' if typ==[1,2] else '_other')
                if key in witnesses and witnesses[key] is None:witnesses[key]=p
            signs=[int(d in squares) for d in ds];assert signs==row['discriminant_is_square']
            if signs[0]!=signs[1] and witnesses['distinct_quadratic_fields'] is None:
                witnesses['distinct_quadratic_fields']=p
            fd=[a*a-4*p for a in traces];assert all(v<0 for v in fd)
            ordinary=all(a%p for a in traces);prod=fd[0]*fd[1]
            different=ordinary and math.isqrt(prod)**2!=prod
            assert row['Frobenius_discriminants']==fd and row['ordinary_pair']==ordinary
            assert row['distinct_Frobenius_fields']==different
            if different and witnesses['geometric_nonisogeny'] is None:witnesses['geometric_nonisogeny']=p
        assert witnesses==result['witnesses'] and all(v is not None for v in witnesses.values())
        assert result['status']=='CERTIFIED_TRANSPORT_OBSTRUCTION'
        wr=next(r for r in result['trials'] if r['p']==witnesses['geometric_nonisogeny'])
        checks.append({'label':equations['label'],'complete_good_pairs':good,'witnesses':witnesses,
          'ordinary_pair':{k:wr[k] for k in ['p','abc','point_counts','traces','Frobenius_discriminants']},
          'geometric_Hom':'ZERO_IN_BOTH_DIRECTIONS',
          'joint_2_torsion_Galois_group':'S3_X_S3',
          'rational_torsion_module_Hom':'ZERO_IN_BOTH_DIRECTIONS'})
    assert len(panel['cases'])==10 and panel['all_certified']
    # Check every2x2 map under the full product action: no alignment of
    # independently selected Frobenius conjugacy representatives is assumed.
    mats=[matrix(GF(2),2,2,v) for v in itertools.product([0,1],repeat=4)]
    group=[M for M in mats if M.det()];assert len(group)==6
    orbits={}
    for M in mats:
        orbit={matrix_key(B*M*A.inverse()) for A in group for B in group}
        rank=int(M.rank());assert len(orbit)==[1,9,6][rank]
        if len(orbit)==1:assert M==0
        orbits.setdefault(rank,set()).add(len(orbit))
    # Explicit six-point torsion-isomorphism scheme. A solution sends a
    # target cubic root to u+v*theta+w*theta^2 in the source cubic algebra.
    S=PolynomialRing(QQ,names=('a','b','c','A','B','C','u','v','w'))
    a,b,c,A,B,C,u,v,w=S.gens()
    T=matrix(S,[[0,0,-c],[1,0,-b],[0,1,-a]])
    assert T**3+a*T**2+b*T+c*identity_matrix(S,3)==0
    N=u*identity_matrix(S,3)+v*T+w*T**2
    cp=N.charpoly();constraints=[cp[2]-A,cp[1]-B,cp[0]-C]
    assert [f.degree(u)+0 for f in constraints]==[1,2,3]
    scheme={'classification':'explicit degree-six torsion-isomorphism scheme, not a seed cover',
      'variables':list(S.variable_names()),
      'equations':[str(f) for f in constraints],
      'coefficient_substitution':'(a,b,c)=carrier_abc; (A,B,C)=original_abc from each frozen equation pair',
      'algebra_map':'original root -> u+v*theta+w*theta^2 in Q[theta]/(theta^3+a*theta^2+b*theta+c)',
      'proof':'For distinct source roots, evaluation of u+v*theta+w*theta^2 is an invertible Vandermonde map. The characteristic-polynomial equations give exactly the six bijections to the distinct target roots. They are reduced in characteristic zero. The full S3xS3 action is transitive, so the scheme is Spec of a degree-six field overQ. It has no rational point. This is a torsion-module isomorphism only, not an elliptic isogeny or a rational seed.',
      'map_orbit_sizes':{'rank0':1,'rank1':9,'rank2':6},
      'extension_degree_necessity':'A rank-one equivariant2-torsion map requires9 dividing[L:Q]; a rank-two map requires6 dividing[L:Q]. The six-point scheme supplies a degree-six torsion isomorphism field; geometric elliptic Hom remains zero over every extension.'}
    retain(DIR/'torsion-isomorphism-scheme.json',scheme)
    paths.append(DIR/'torsion-isomorphism-scheme.json')
    output={'status':'PASS_INDEPENDENT_TEN_GENUS1_GEOMETRIC_AND_TORSION_TRANSPORT_OBSTRUCTIONS',
      'classification':'new verified application and exact auxiliary-field obstruction',
      'cases':checks,'prime_exposures':640,'map_orbit_sizes':{'rank0':1,'rank1':9,'rank2':6},
      'boundary':'No homomorphism of curves or coefficient2-torsion modules transports the carrier class to the fixed original elliptic fibre overQ. These statements do not exclude nonlinear incidence or give a seed-production rule. No global Selmer group, Sha class or point search is computed.',
      'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths}}
    retain(DIR/'replay.json',output)
    print(output['status'],flush=True)
    for c in checks:print(c['label'],c['ordinary_pair'],flush=True)
if __name__=='__main__':
    signal.alarm(25);verify()

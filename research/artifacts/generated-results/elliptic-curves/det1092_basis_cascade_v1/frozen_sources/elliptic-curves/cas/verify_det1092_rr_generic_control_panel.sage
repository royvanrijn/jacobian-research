#!/usr/bin/env sage-python
"""Replay the entire frozen panel and verify the universal marked section.

V4 terminal labels are inspected only after every arithmetic replay succeeds.
No V4 file, process or protocol is changed.
"""
import hashlib,json,runpy,signal
from pathlib import Path
from collections import Counter
from sage.all import QQ,PolynomialRing,EllipticCurve
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
DIR=ART/'det1092_rr_generic_point_controls_v2'
CAS=ROOT/'elliptic-curves/cas'
OUT=DIR/'panel-replay.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify():
    checker=CAS/'verify_det1092_rr_generic_point_controls.sage'
    replay=runpy.run_path(str(checker))['verify']
    cases=[replay(i) for i in range(9)]
    assert all(r['status']=='PASS_GENERIC_ELLIPTIC_POINT_GIVES_NON_GENERIC_RATIONAL_JACOBIAN_CLASS' for r in cases)
    protocol_path=DIR/'protocol.json';protocol=json.loads(protocol_path.read_text())
    parents=[ART/n for n in ['curve302_recovered_mw17_parent_v1.json','det1092_first_centre_rr_net_v1.json',
             'det1092_universal_rr_descent_preflight_v1.json']]
    parent,net,uni=[json.loads(p.read_text()) for p in parents];uni=uni['universal_genus2']
    R=PolynomialRing(QQ,'t');t=R.gen();K=R.fraction_field()
    def rat(row):return K(R(row['numerator']))/R(row['denominator'])
    E=EllipticCurve(K,[rat(row) for row in parent['a_invariants']])
    basis=[E([rat(row) for row in P]) for P in parent['basis_weierstrass_coordinates']]
    centre=-sum((n*P for n,P in zip(net['trace_word'],basis)),E(0));h=R(uni['h'])
    assert centre[0].denominator()==h*h
    S0=basis[0];A,B=[[R(row) for row in net[key]] for key in ['A','B']]
    r=-(B[0]+B[1]*S0[0]+B[2]*S0[1])/(A[0]+A[1]*S0[0]+A[2]*S0[1])
    assert max(r.numerator().degree(),r.denominator().degree())==7
    assert list(map(str,r.numerator().list()))==protocol['selector']['u_numerator']
    assert list(map(str,r.denominator().list()))==protocol['selector']['u_denominator']
    # Here the RR curve coordinate T has been specialized to its marked t.
    f0,f1,f2=[B[j]+r*A[j] for j in range(3)]
    assert f0+f1*S0[0]+f2*S0[1]==0
    cx=centre[0]+E.b2()/12;cy=centre[1]+(E.a1()*centre[0]+E.a3())/2
    mnum=-f1+E.a1()*f2/2
    s=((2*(S0[0]+E.b2()/12)+cx)*f2**2-mnum**2)/h**3
    diagonal=sum(QQ(row['coefficient'])*t**row['T']*r**row['u'] for row in uni['sparse_q'] if row['v']==0)
    assert s*s==QQ(uni['scale'])*diagonal
    m=mnum/f2;omega=h**3*s/f2**2
    xshort=(m*m-cx+omega)/2;yshort=m*(xshort-cx)-cy
    xx=xshort-E.b2()/12;yy=yshort-(E.a1()*xx+E.a3())/2
    assert E([xx,yy])==S0
    def record(v):return {'numerator':list(map(str,v.numerator().list())),
                          'denominator':list(map(str,v.denominator().list()))}
    counts=Counter();table=[]
    for i,row in enumerate(cases):
        raw=json.loads((DIR/f'case-{i:02d}.json').read_text())
        counter=Counter(trial['status'] for trial in raw['trials']);counts.update(counter)
        table.append({'label':row['label'],'tau':row['tau'],
                      'local_blocks':row['passing_prime_blocks'],'local_skips':64-row['passing_prime_blocks'],
                      'fake_ranks':[16,17],'subgroup_ranks':[17,18],
                      'elliptic_class_mod_MW17':'ZERO','Jacobian_class_mod_generic_image':'NONZERO_AND_INDEPENDENT',
                      'known_rational_Kummer_subspace_dimension':18,'marked_Sha_image':'ZERO',
                      'theta_primes':row['theta_primes']})
    assert sum(counts.values())==576 and counts['PASS_LOCAL_KUMMER_BLOCK']==361
    # Outcome labels are read only now. Missing certificates are not nulls.
    v4=ROOT/'artifacts/local/elliptic-curves/det1092-v4-wide-bootstrap-v2'
    outcomes=[]
    for row in table[:8]:
        path=v4/row['label']/'v4-verified.json'
        if not path.exists():
            outcomes.append({'label':row['label'],'status':'NO_COMPLETED_V4_CERTIFICATE_AT_SNAPSHOT'})
        else:
            saved=json.loads(path.read_text())
            assert saved['status']=='PASS_INDEPENDENT_DET1092_V4_REPLAY'
            null=saved['rank_lower_bound']==17 and saved['gain']==0 and saved['stop_reason']=='COMPLETE_FRESH_BOOTSTRAP_NO_GAIN'
            outcomes.append({'label':row['label'],'status':'VERIFIED_BOUNDED_V4_NULL' if null else 'OTHER_VERIFIED_V4_OUTCOME',
                             'certificate':str(path.relative_to(ROOT)),'sha256':sha(path)})
    return {'classification':'verified application and new deduction; completed generic-input specificity panel',
            'status':'PASS_NINE_GENERIC_ELLIPTIC_POINTS_GIVE_NON_GENERIC_RATIONAL_JACOBIAN_CLASSES',
            'universal_family':{'equation':'C_t: s^2=c*q(T;r0(t),0), with q,c from the pinned generic universal RR certificate',
                                'r0':record(r),'degree_t_to_RR_pencil':7,
                                'marked_point':{'T':'t','s':record(s)},
                                'symbolic_point_and_map_identities':True,
                                'elliptic_image':'Exactly the original generic section0 over Q(t)',
                                'generic_Jacobian_rank_lower_bound':18,
                                'original_elliptic_generic_rank':17,
                                'rank_proof':'All18 subgroup generators extend over Q(t); independence at the verified first smooth specialization implies generic independence.',
                                'trace_relation':'Trace_{Q(t)/Q(u)}([P(t)-P0])=Phi(S0), where u=r0(t).'},
            'table':table,'prime_exposure':dict(counts),'post_arithmetic_V4_snapshot':outcomes,
            'class_distinction':{'inherited_basepoint':'Zero Abel-Jacobi class',
                                 'marked_Jacobian_classes':'Non-generic rational classes, independent of the full generic NS restriction image',
                                 'true_Kummer_classes':'In an explicitly certified18-dimensional rational subspace of Sel_2(J)',
                                 'Sha_images':'All zero; not merely classes surviving local tests',
                                 'Cassels_Tate_on_known_rational_subspace':'Zero, since the pairing factors through Sha',
                                 'full_Selmer_groups_and_other_Sha_classes':'NOT_COMPUTED'},
            'conclusion':'A non-generic rational Jacobian class relative to the original K3 Picard image does not imply that its marked elliptic image is an extra MW direction. The same property occurs uniformly in a generic-input construction with exact original generic elliptic rank17.',
            'scope':'This refutes that candidate discriminator. It does not prove absence of elliptic rank jumps on V4 controls or rule out a refined arithmetic discriminator.',
            'limits':protocol['limits'],
            'inputs':{str(p.relative_to(ROOT)):sha(p) for p in parents+[protocol_path,checker]+[DIR/f'case-{i:02d}.json' for i in range(9)]},
            'checker_sha256':sha(Path(__file__))}
if __name__=='__main__':
    signal.alarm(25);result=verify();payload=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if OUT.exists():
        # External V4 jobs may acquire later terminal certificates. Their old
        # snapshot is evidence, not a demand to mutate this frozen result.
        old=json.loads(OUT.read_text());result['post_arithmetic_V4_snapshot']=old['post_arithmetic_V4_snapshot']
        payload=json.dumps(result,indent=2,sort_keys=True)+'\n';assert OUT.read_text()==payload
    else:OUT.write_text(payload)
    print(result['status'],flush=True)

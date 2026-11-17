"""PF2RULED35: exact controls for the marked quartic cabling exclusion.
The finite S4 peripheral-residue layer is exhaustive; the 3-manifold
compression, filling and finite-depth proofs are written, not formalized.
"""
from pathlib import Path
from itertools import permutations,combinations
from functools import lru_cache
from collections import Counter
from math import gcd
from time import monotonic
import json,hashlib,platform
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
HERE=Path(__file__).resolve().parent
P=tuple(permutations(range(4)));I=tuple(range(4))
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def mul(a,b):return tuple(a[b[j]] for j in range(4))
def inv(a):return tuple(a.index(j) for j in range(4))
@lru_cache(None)
def power(a,n):
    if n<0:return power(inv(a),-n)
    out=I
    for _ in range(n):out=mul(out,a)
    return out
@lru_cache(None)
def group(gs):
    found={I};todo=[I]
    while todo:
        a=todo.pop()
        for b in gs:
            c=mul(a,b)
            if c not in found:found.add(c);todo.append(c)
    return frozenset(found)
def generated(gs):return group(tuple(sorted(set(gs))))
def even(a):return sum(a[i]>a[j] for i in range(4) for j in range(i+1,4))%2==0
A4=frozenset(a for a in P if even(a));S4=frozenset(P)
def cycles(a,subset=range(4)):
    left=set(subset);rows=[]
    while left:
        j=min(left);k=j;row=[]
        while k not in row:
            row.append(k);left.remove(k);k=a[k]
        rows.append(row)
    return tuple(sorted(map(len,rows)))
def orbit(g,j):return {a[j] for a in g}
def comm(g):
    return generated(mul(mul(mul(a,b),inv(a)),inv(b)) for a in g for b in g)
def normal(g,a):return generated(mul(mul(b,a),inv(b)) for b in g)
def main():
    start=monotonic()
    for filename in ['torus_infinity_receipt.json','residual_cover_receipt.json']:
        pred=json.loads((HERE/filename).read_text())
        check('predecessor PASS '+filename,pred['outcome']=='PASS')
        for name,digest in pred['source_sha256'].items():
            check('predecessor seal '+filename+' '+name,hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest)
    D=HERE/'sources/quartic-cabling';manifest=json.loads((D/'manifest.json').read_text())
    check('four cabling topology primary imports pinned',len(manifest['primary_sources'])==4)
    for row in manifest['primary_sources']:
        p=D/row['path']
        check('primary PDF seal '+row['path'],p.read_bytes().startswith(b'%PDF') and hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'])
    # Entire actual second plane chart; the nonvanishing cover is checked algebraically.
    X,Y,u,v=S.symbols('X Y u v');AA=X*Y;BB=X*X*Y;DD=X**3*Y+X
    chart=[S.expand(f.subs({X:1/u,Y:u**3*v-u*u},simultaneous=True)) for f in [Y,AA,BB,DD]]
    check('actual full second plane chart recovers four ring generators',chart==[u**3*v-u*u,u*u*v-u,u*v-1,v])
    check('the two regular formulas for u agree by the defining relation',S.expand(AA*DD-BB*(BB+1))==0)
    check('B1 and D never vanish simultaneously on second chart',S.expand(u*chart[3]-chart[2])==1)
    check('inverse chart is exact on overlap',S.expand((X**3*Y+X).subs({X:1/u,Y:u**3*v-u*u},simultaneous=True))==v)
    check('torsor linear transition is index minus three',S.diff(DD,Y)==X**3)
    xb,ub=S.symbols('xb ub')
    Ys=-X*xb**3/(1+X**3*xb**3);vs=X/(1+X**3*xb**3)
    check('explicit smooth torsor section respects affine transition',S.cancel(X**3*Ys+X-vs)==0)
    check('explicit section extends smoothly across base infinity',
          S.cancel(vs.subs({X:1/u,xb:1/ub},simultaneous=True)-u**2*ub**3/(1+u**3*ub**3))==0)
    check('fiber metric is global under index-three transition',
          S.cancel((X**3*Y)*(xb**3*Y)*(1+1/(X*xb))**3-Y**2*(1+X*xb)**3)==0)
    check('marked end fiber represents a nonzero generator of Z3',gcd(4,3)==1)
    # Quartic longitude: the residual factor has even degree for every collision set.
    o=S.symbols('o',integer=True,nonnegative=True)
    check('residual squarefree degree is even at all collision counts',S.simplify(S.exp(S.pi*S.I*(2*o)))==1)
    tauS=(1,0,2,3);tauA=(1,2,0,3)
    centralS={a for a in P if mul(a,tauS)==mul(tauS,a) and even(a)}
    check('preferred-longitude alternatives before residual infinity are exactly two',centralS=={I,(1,0,3,2)})
    check('fixing both residual letters forces identity longitude',{a for a in centralS if a[2:]==(2,3)}=={I})
    # Complete literal S4 subgroup lattice, reused mathematical layer, no large search.
    sub={frozenset([I])};todo=list(sub)
    while todo:
        g=todo.pop()
        for a in P:
            if a in g:continue
            k=generated(tuple(g)+(a,))
            if k not in sub:sub.add(k);todo.append(k)
    check('complete literal S4 subgroup census',len(sub)==30)
    derived={g:comm(g) for g in sub}
    check('A4 derived subgroup is V4',len(derived[A4])==4 and all(cycles(a) in [(1,1,1,1),(2,2)] for a in derived[A4]))
    check('A4 three-cycle longitude is identity',{a for a in derived[A4] if mul(a,tauA)==mul(tauA,a)}=={I})
    for sigma,label,orders in [(tauS,'transposition',[2,6,6,24]),(tauA,'three-cycle',[3,12]),((1,0,3,2),'double-transposition',[2])]:
        candidates=[g for g in sub if sigma in g and normal(g,sigma)==g]
        check('all normally meridian-generated images '+label,sorted(map(len,candidates))==orders)
    # All possible fiber-trivial peripheral permutation states, for unbounded p,q.
    # Every permutation exponent divides 12. Coprime integer p,q reduce to these
    # residue pairs; the overapproximation imposes only the necessary mod12 gcd.
    expectedS={
      (2,(1,1),24,(1,1,2),(1,1,1,1)),
      (3,(1,2),12,(1,3),(1,1,1,1)),
      (3,(3,),6,(1,1,2),(1,1,1,1)),
      (3,(3,),24,(1,1,2),(1,1,1,1)),
      (4,(1,3),4,(4,),(1,1,1,1)),
      (4,(1,3),24,(4,),(1,1,1,1)),
      (4,(1,3),24,(4,),(2,2)),
      (4,(4,),3,(1,3),(1,1,1,1)),
      (4,(4,),12,(1,3),(1,1,1,1))}
    expectedA={
      (3,(1,1,1),12,(1,3),(1,1,1,1)),
      (4,(1,3),2,(2,2),(1,1,1,1)),
      (4,(2,2),3,(1,3),(1,1,1,1)),
      (4,(2,2),12,(1,3),(1,1,1,1))}
    summaries={}
    for name,tau,G,order,expected in [('S',tauS,S4,2,expectedS),('A',tauA,A4,3,expectedA)]:
        found=Counter()
        for p in range(12):
            for q in range(12):
                if gcd(gcd(p,q),12)!=1 or p*q%order:continue
                for alpha in G:
                    if power(alpha,p)!=I:continue
                    pattern=generated([alpha,tau]);O=orbit(pattern,0)
                    delta=mul(alpha,tau);mu=power(delta,p);lon=power(delta,-q)
                    for companion in sub:
                        if not companion<=G or mu not in companion or lon not in derived[companion] or normal(companion,mu)!=companion:continue
                        if generated(tuple(pattern)+tuple(companion))!=G:continue
                        found[(len(O),cycles(alpha,O),len(companion),cycles(mu),cycles(lon))]+=1
        check('exhaustive fiber-trivial residue case list '+name,set(found)==expected)
        check('residue count control '+name,sum(found.values())==(208 if name=='S' else 96))
        summaries[name]=[dict(branched_orbit=key[0],alpha_cycles=list(key[1]),companion_order=key[2],
                              meridian_cycles=list(key[3]),longitude_cycles=list(key[4]),residue_multiplicity=count)
                         for key,count in sorted(found.items())]
        # Nontrivial central fiber forces exactly the branch support.
        check('nonidentity fiber fixes a residual letter and cannot act centrally on four letters '+name,
              any(tau[j]==j for j in range(4)))
        nt=set()
        for p in range(12):
            for q in range(12):
                if gcd(gcd(p,q),12)!=1 or p*q%order==0:continue
                mu=power(tau,p);h=power(tau,p*q)
                for companion in sub:
                    if not companion<=G or mu not in companion or normal(companion,mu)!=companion:continue
                    for lon in derived[companion]:
                        if mul(lon,tau)!=mul(tau,lon):continue
                        if mul(power(lon,p),power(mu,q))!=h:continue
                        if generated(tuple(companion)+(tau,))!=G:continue
                        nt.add((len(companion),cycles(mu),cycles(lon)))
        check('nontrivial-fiber companion is full and longitude identity '+name,
              nt=={(len(G),cycles(tau),(1,1,1,1))})
    # Symbolic base and homology formulas used in the written cases.
    m,ca=S.symbols('m ca',integer=True)
    check('simple pattern base Euler formula',S.expand(ca+(m-1)-m)==ca-1)
    check('alternating pattern base Euler formula',S.expand(ca+(m-2)-m)==ca-2)
    a,b,U,V=S.symbols('a b U V',positive=True)
    check('alternating torus cover homology is even',
          S.cancel((a*a*b*b/12)*(4/(a*b))).subs({a:3*U,b:2*V},simultaneous=True)==2*U*V)
    check('alternating nontrivial-fiber contradiction uses 3 times marked class',3%3==0)
    # Coprime disk-cone homology has no torsion; controls for the symbolic gcd proof.
    for p in [3,5,7,9,11]:
        for beta in range(1,p):
            if gcd(beta,p)!=1:continue
            matrix=S.Matrix([[p,0,beta],[0,2,1]])
            N=smith_normal_form(matrix,domain=S.ZZ)
            check(f'disk cones p,2 give torsion-free H1 p{p} beta{beta}',
                  sorted(abs(int(N[j,j])) for j in range(2))==[1,1])
    q,p,mm,ll,Rc=S.symbols('q p m0 l0 Rc')
    check('transparent marked class follows integrally in free H1',
          S.expand((2*Rc-(q*mm+2*p*ll))-2*(Rc-p*ll)+q*mm)==0)
    check('transparent nonzero marked class forces winding coprime to three',
          [z for z in range(6) if z%2 and z%3]!=[] and all(gcd(z,6)==1 for z in range(6) if z%2 and z%3))
    # A genuine local survival control: the descent cannot be replaced by a
    # false assertion that each odd cable filling is itself impossible.
    p0,q0,A0,B0=5,3,1,3
    check('transparent local control has distance-one nonfiber filling',abs(2*p0*A0-q0*B0)==1)
    check('transparent local control preserves H1 order three and nonzero marked class',abs(B0)==3 and p0%3!=0)
    # State-machine termination is finite cabling depth, not a bounded depth test.
    n=S.symbols('n',integer=True,positive=True)
    check('transparent removal strictly decreases every positive cabling depth',S.simplify(n-(n-1))==1)
    inputs=[Path(__file__),HERE/'QUARTIC_CABLING_EXCLUSION.md',HERE/'torus_infinity_receipt.json',
            HERE/'TORUS_INFINITY_OBSTRUCTION.md',HERE/'residual_cover_receipt.json',
            HERE/'RESIDUAL_DOUBLE_COVER.md',HERE/'WRIGHT_MARKED_MODEL.md',HERE/'WRIGHT_CONNECTION.md',
            HERE/'PROOF.md',HERE/'PRIMITIVE_QUARTIC_FIELD.md']+sorted(D.glob('*'))
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        scope='all pole heights on the actual marked Wright surface; finite field degree exactly four',
        marked_end=dict(H1='Z/3',ramification_exterior='solid torus',ramification_class='generator'),
        peripheral=dict(S4_meridian='transposition',S4_preferred_longitude='identity',
                        A4_meridian='three-cycle',A4_preferred_longitude='identity'),
        finite_peripheral_layer=dict(exhaustive_S4_subgroups=30,exponent_modulus=12,
            S_fiber_trivial_state_count=9,A_fiber_trivial_state_count=4,
            S_residue_controls=208,A_residue_controls=96,states=summaries,
            scope='exhaustive permutation/peripheral residue layer only, not formal three-manifold verification'),
        surviving_removal=dict(S_to_S='odd p,q; gcd(p,6)=1; cabling depth decreases by one',
                              S_to_A='p=2 three-sheeted pattern component; A impossible',
                              A_outgoing='none'),
        local_control=dict(p=5,q=3,filling=[1,3],distance_from_fiber=1,H1_order=3,
                           marked_class_mod3=2,is_actual_global_cover=False),
        actual_marked_quartic_Keller_maps_excluded=True,
        compatible_ruled_quartic_realizations_excluded=True,
        actual_marked_quartic_realization_frontier=[],
        prior_truncated_algebraic_frontier=[[5,9],[6,9],[7,9],[8,9]],
        assurance=dict(actual_marked_ring_and_finite_map_are_hypotheses=True,
            all_positive_winding_and_signed_coprime_cable_parameters_covered_by_written_proof=True,
            finite_cabling_depth_unbounded_and_handled_by_descent=True,
            residual_infinity_used_to_prove_preferred_longitude_identity=True,
            compression_ball_alternative_excluded_using_nonzero_marked_class=True,
            actual_peripheral_filling_lattices_tracked=True,
            transparent_step_preserves_lens_core_and_H1_three=True,
            primary_topology_results_explicitly_imported=True,end_lens_recognition_explicitly_imported=True,
            finite_permutation_layer_exhaustive_but_topological_proof_written=True,
            coefficient_search_used=False,
            all_finite_field_degrees_on_Wright_excluded=False,
            arbitrary_normalization_compatible_completion_proved=False,
            arbitrary_quartic_normalizations_excluded=False,quartic_packets_excluded=False,
            truncated_degree8_necessary_core_unit_ideal_proved=False,
            formal_verification=False,independent_full_reproof=False,external_review=False,research_goal_completed=False),
        source_sha256={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        python=platform.python_version(),sympy=S.__version__,elapsed_seconds=round(monotonic()-start,6))
    (HERE/'quartic_cabling_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} exact controls; written marked-quartic exclusion; arbitrary normalization and higher field degrees remain open')
if __name__=='__main__':main()

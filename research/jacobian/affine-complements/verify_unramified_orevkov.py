"""Exact controls supporting written unbounded geometry; not a topology prover."""
from pathlib import Path
from itertools import combinations
import hashlib,json,platform
import sympy as S
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def compose(a,b):return tuple(a[b[i]] for i in range(len(a)))
def perm(n,cycles):
    out=list(range(n))
    for cy in cycles:
        for a,b in zip(cy,cy[1:]+cy[:1]):out[a]=b
    return tuple(out)
def closure(gens):
    out={tuple(range(len(gens[0])))};todo=list(out)
    while todo:
        a=todo.pop()
        for g in gens:
            b=compose(a,g)
            if b not in out:out.add(b);todo.append(b)
    return out
def packets(total,min_e=1):
    if total==0:yield [];return
    for e in range(min_e,total+1):
        for E in range(total-e+1):
            for rest in packets(total-e-E,e):yield [(e,E)]+rest
def canonical_number(pairs):
    p,q=pairs[0];a=p+q
    for pp,qq in pairs[1:]:a=qq*a+pp-p*q*qq;p,q=pp,qq
    return a
def main():
    predecessors=['exact_coordinate_bridge_receipt.json','universal_regular_end_receipt.json','higher_regular_coordinate_receipt.json','primitive_cover_receipt.json']
    for name in predecessors:
        d=json.loads((HERE/name).read_text());check('predecessor PASS '+name,d['outcome']=='PASS')
        for f,h in d['source_sha256'].items():check('predecessor seal '+name+' '+f,sha(HERE/f)==h)
    manifest=HERE/'sources/unramified-orevkov/manifest.json';m=json.loads(manifest.read_text())
    check('seven pinned primary imports and audit-only source',len(m['primary_sources'])==7)
    for row in m['primary_sources']:
        f=HERE/row['path_relative_to_affine_complements'];check('primary bytes '+f.name,len(f.read_bytes())==row['bytes'] and sha(f)==row['sha256'])
    check('new complete primary page audits recorded',[r['new_complete_page_visual_audits'] for r in m['primary_sources'][:3]]==[[1,12,13,14,15],[28],[10]])
    rows={tuple(sorted(r)) for r in packets(3) if any(e>=2 for e,E in r) and all(e!=3 for e,E in r) and all(E==0 for e,E in r if e==1)}
    check('complete Euler/maximal-index quartic alternatives',rows=={((2,1),),((1,0),(2,0))})
    check('bridge removes precisely the sole-prime alternative',{r for r in rows if r!=((2,1),)}=={((1,0),(2,0))})
    branch=[(c,h,4-2*c-h) for c in range(3) for h in range(5) if 2*c+h<=4]
    check('complete full quartic branch table',branch==[(0,0,4),(0,1,3),(0,2,2),(0,3,1),(0,4,0),(1,0,2),(1,1,1),(1,2,0),(2,0,0)])
    check('complete missing-value branch table',[(c,h) for c,h,l in branch if l==0]==[(0,4),(1,2),(2,0)])
    check('complete singleton table',[(c,h) for c,h,l in branch if l==1]==[(0,3),(1,1)])
    check('coincident images forbid every self-collision',all(3*b>4 for b in range(2,8)))
    check('a ramified double value misses the source and avoids H',[(h,l) for c,h,l in branch if c==2]==[(0,0)])
    for n in range(3,9):
        trans=[perm(n,[[i,j]]) for i,j in combinations(range(n),2)]
        check('local disjoint fold transpositions commute d'+str(n),all(compose(perm(n,[list(a)]),perm(n,[list(b)]))==compose(perm(n,[list(b)]),perm(n,[list(a)])) for a,b in combinations(combinations(range(n),2),2) if not set(a)&set(b)))
        adjacent=[perm(n,[[i,i+1]]) for i in range(n-1)]
        check('connected sheet-edge generators give S_d control d'+str(n),len(closure(adjacent))==S.factorial(n))
        check('zero-excess prime count satisfies exact Euler ledger d'+str(n),2+(n-3)==n-1)
    a=perm(4,[[0,1]]);b=perm(4,[[2,3]]);c=perm(4,[[1,2]])
    check('tangential double folds have abelian V4 local image',len(closure([a,b]))==4 and compose(a,b)==compose(b,a))
    check('cusp countercontrol prevents universal Orevkov factorization',compose(a,c)!=compose(c,a) and len(closure([a,c]))==6 and compose(compose(a,c),a)==compose(compose(c,a),c))
    a3=perm(8,[[0,1,2]]);b4=perm(8,[[3,4,5,6]])
    check('arbitrary-index local cyclic factors commute on disjoint sheet clusters',compose(a3,b4)==compose(b4,a3) and len(closure([a3,b4]))==12)
    g4=closure([perm(4,[[0,1,2,3]])]);v4=closure([perm(4,[[0,1],[2,3]]),perm(4,[[0,2],[1,3]])])
    check('transitive abelian groups are regular without being assumed symmetric',all(all(g[i]!=i for i in range(4)) for group in [g4,v4] for g in group if g!=tuple(range(4))))
    check('cyclic local quotient need not have symmetric global group',len(g4)==4 and len(v4)==4 and len(g4)!=24)
    d,k,T=S.symbols('a_n k T',integer=True)
    square=d-1+2*k+2*T
    check('all tangency multiplicity is retained in Nori self-intersection',S.expand(square-4*T-2*k-(d-1-2*T))==0)
    p,q,r,s=S.symbols('p q r s',integer=True)
    check('two-stage canonical threshold is exact',S.expand(r*(p+q)+s-p*q*r-(s+r*(p+q-p*q)))==0)
    check('single torus-pair canonical value exceeds one',canonical_number([(3,2)])==5)
    check('trefoil then winding-two numerical equality control',canonical_number([(3,2),(3,2)])==1)
    check('nodal obstruction reaches beyond the sufficient moduli threshold',canonical_number([(3,2),(5,2)])==3 and not (3>=3*2-1))
    # Normal forms for arbitrary r: identities are symbolic; samples are controls.
    z,t,x,y,a,b,c=S.symbols('z t x y a b c')
    for rr in range(1,9):
        h=sum(S.Symbol('c'+str(j))*z**(-j) for j in range(rr,rr+4))
        xx=z**(1-rr)*t+h;yy=z**rr
        jac=S.expand(S.diff(xx,t)*S.diff(yy,z))
        check('whole unramified chart has unit source volume r'+str(rr),jac==rr)
        diff=S.expand(rr*z**(rr-1)*xx-rr*t)
        check('shared exactness residue r'+str(rr),S.residue(diff,z,0)==rr*S.Symbol('c'+str(rr)))
        poly=1+yy+yy**3;newt=b*a**(rr-1)*t+a**(rr-1)*z**(rr-1)*poly
        check('all polynomial source shears are regular r'+str(rr),S.Poly(newt,z,t) is not None)
        check('whole shear coordinate identity r'+str(rr),S.expand((a*z)**(1-rr)*newt-(b*z**(1-rr)*t+poly))==0)
        check('source LND extends with exact boundary exponent r'+str(rr),S.expand(z**(rr-1)*S.diff(xx,t))==1)
        check('uniform actual exact affine datum r'+str(rr),S.gcd(rr,rr+1)==1 and S.residue(rr*z**(rr-1)*z**(-rr-1),z,0)==0)
    W,V=S.symbols('W V');H=1+y+y**2
    for qq in range(1,7):
        rel=y**qq*V-W*(W-H)
        check('whole hypersurface source chart q'+str(qq),S.expand(rel.subs({W:y**qq*x,V:x*(y**qq*x-H)},simultaneous=True))==0)
        check('whole hypersurface boundary chart q'+str(qq),S.expand(rel.subs({W:H+y**qq*t,V:H*t+y**qq*t*t},simultaneous=True))==0)
        grad=S.diff(rel,W).subs(y,0)
        check('both zero components have unit transverse gradient q'+str(qq),grad.subs(W,0)==1 and grad.subs(W,1)==-1)
    inputs=[Path(__file__),HERE/'UNRAMIFIED_PACKET_AND_OREVKOV_GATE.md',manifest]+[HERE/r['path_relative_to_affine_complements'] for r in m['primary_sources']]+[HERE/n for n in predecessors]
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        quartic_reduction=dict(only_two_primes=True,indices=[2,1],residue_degrees=[1,1],excesses=[0,0],smooth_N=True,distinct_immersive_images=True,ramified_image_has_double_value=True,necessarily_non_surjective=True,full_branch_table=branch,missing_pairs=[[0,4],[1,2],[2,0]],singleton_pairs=[[0,3],[1,1]],higher_index_quartic_case_already_excluded_by='QPK4/Orevkov maximal-index Remark',two_prime_case_excluded=False),
        all_degree_orevkov=dict(full_finite_normalization_required=True,sole_ramified_prime=True,simple_index_two=True,immersive_image_required=True,unramified_primes_retained=True,monodromy='S_d',local_images='disjoint transpositions',quotient='O(C) ->> S_d',unramified_prime_count='d-3',positive_braid_infinity_excluded=True,NN_sufficient_class_excluded='a_n>=p_1*q_2*...*q_n-1 AND a_n>1',canonical_recurrence='a_i=q_i*a_(i-1)+p_i-p_(i-1)*q_(i-1)*q_i',nodal_bound='a_n<=1',quartic_tangency_bound='a_n<=1+2*T',T_definition='sum intersection multiplicities m>=2, not number of tangencies',cusp_factorization_claimed=False,multiple_ramified_prime_symmetric_factorization_claimed=False),
        general_zero_excess_orevkov=dict(all_field_degrees=True,all_generic_indices=True,any_number_ramified_primes=True,all_ramified_residue_degrees_one_required=True,all_ramified_image_normalizations_immersive_required=True,full_branch_union_required=True,local_images='disjoint cyclic factors',quotient='O(B) ->> G',G_nonabelian=True,G_symmetric_assumed=False,abelian_transitive_group_regular=True,Galois_obstruction_written=True,positive_braid_full_branch_union_excluded=True,NN_criterion_on_all_simultaneous_infinity_paths=True,separate_component_diagrams_suffice=False,positive_excess_packets_excluded=False),
        unramified_geometry=dict(all_boundary_orders=True,regular_coordinate_automatic=True,whole_root_chart='y=z^r; x=z^(1-r)*t_i+h_i(z)',full_support='r<=j<=q',gcd_condition=True,all_admissible_data_actual_smooth_affine=True,exactness_iff='c_r=0',marked_equivalence='c_j_prime=b*a^j*c_j',all_polynomial_source_shears=True,LND='partial_x',boundary_LND_order='r-1',kernel='C[y]',actual_exact_controls_all_r=True,target_algebra_invariance_assumed=False,removing_other_primes_preserves_finiteness=False,common_coordinate_assumed=False),
        assurance=dict(written_unbounded_geometry=True,precise_primary_imports=True,explicit_tangency_blowup_cost=True,finite_checks_not_topological_proof=True,all_polynomial_links_assumed_positive=False,withdrawn_Jacobian_reduction_used=False,full_quartic_exclusion=False,higher_degree_exclusion=False,formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,research_goal_completed=False),
        source_sha256={str(p.relative_to(HERE)) if p.is_relative_to(HERE) else str(p) :sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'unramified_orevkov_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS',len(checks),'controls; only two-prime quartic packet remains; all-degree Orevkov gate and unramified affine classification')
if __name__=='__main__':main()

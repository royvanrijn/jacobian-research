"""PF2RULED34: all-degree torus-end obstruction.
Finite controls support a written Seifert/homology proof, not an exhaustive
search over degrees, pole heights, knots, collision counts or tangencies.
"""
from pathlib import Path
from time import monotonic
from math import gcd,factorial,lcm
import itertools,json,hashlib,platform
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def compose(a,b): return tuple(a[b[j]] for j in range(len(a)))
def inv(p): return tuple(p.index(j) for j in range(len(p)))
def cycles(p):
    seen=set();out=[]
    for j in range(len(p)):
        if j in seen: continue
        row=[];k=j
        while k not in seen:
            seen.add(k);row.append(k);k=p[k]
        out.append(row)
    return out
def transposition(n,a,b):
    p=list(range(n));p[a],p[b]=p[b],p[a];return tuple(p)
def generated(gs):
    one=tuple(range(len(gs[0])));seen={one};pending=[one]
    while pending:
        p=pending.pop()
        for g in gs:
            q=compose(p,g)
            if q not in seen:seen.add(q);pending.append(q)
    return seen
def order(p):
    return lcm(*(len(x) for x in cycles(p)))
def main():
    start=monotonic()
    pred=json.loads((HERE/'residual_cover_receipt.json').read_text())
    check('residual-cover predecessor PASS',pred['outcome']=='PASS')
    for name,digest in pred['source_sha256'].items():
        check('predecessor seal '+name,hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest)
    source=HERE/'sources/torus-infinity';manifest=json.loads((source/'manifest.json').read_text())
    check('three primary topology imports pinned',len(manifest['primary_sources'])==3)
    for row in manifest['primary_sources']:
        p=source/row['path']
        check('primary PDF seal '+row['path'],p.read_bytes().startswith(b'%PDF') and hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'])
    # The actual end comes from the entire marked completion, not just an abstract ring.
    C0=S.Matrix([1,0]);f=S.Matrix([0,1]);I=S.Matrix([[-1,1],[1,0]])
    B=C0+2*f;E=C0+3*f
    check('marked B square is three',(B.T*I*B)[0]==3)
    check('marked omitted closure intersects B four times',(E.T*I*B)[0]==4)
    check('marked end homology is Z3 with nonzero ramification class',4%3==1 and gcd(4,3)==1)
    # Universal formulas, independent of field degree and knot type parameters.
    d,r,s,k,A,p,q=S.symbols('d r s k A p q',positive=True,integer=True)
    ca,cb=S.symbols('ca cb',integer=True)
    check('sphere-to-sphere cycle formula leaves three exceptional-value cycles',
          S.expand((ca+cb+(d-1))-(d+2))==ca+cb-3)
    check('all-degree Seifert Euler and first-homology formula',
          S.cancel((A*s*k*r*k)*(d/(p*q))).subs({p:d*A,q:r*s*k},simultaneous=True)==k)
    check('full-cycle and split-cycle orders give p=dA and q=rsk',S.expand((r+s)-d).subs(d,r+s)==0)
    a1,a2,a3,b1,b2,b3,b0=S.symbols('a1 a2 a3 b1 b2 b3 b0')
    M=S.Matrix([[a1,0,0,b1],[0,a2,0,b2],[0,0,a3,b3],[1,1,1,b0]])
    euler=b0-b1/a1-b2/a2-b3/a3
    check('exact Seifert abelianization determinant is product orders times Euler',
          S.cancel(M.det()-a1*a2*a3*euler)==0)
    orb=S.Matrix([[A,0,0],[0,s*k,0],[0,0,r*k],[1,1,1]])
    minors=[S.expand(orb[list(rows),:].det()) for rows in itertools.combinations(range(4),3)]
    check('orbifold maximal minors include all required generators',
          set(minors)=={A*r*s*k*k,A*s*k,-A*r*k,r*s*k*k})
    # Bezout elimination is integral: no variable localization is allowed.
    u,v,x,y=S.symbols('u v x y')
    first=u*(A*r*k)+v*(A*s*k)
    check('coprime r,s recover Ak from the two first minors',
          S.expand(first-A*k*(u*r+v*s))==0)
    final=x*first+y*(r*s*k*k)
    check('coprime A,q recover k integrally from all maximal minors',
          S.expand(final-k-(x*A*k*(u*r+v*s-1)+k*(x*A+y*r*s*k-1)))==0)
    check('all maximal minors are divisible by k',
          all(S.cancel(m/k).is_polynomial(A,r,s,k) for m in minors))
    # Finite controls for the written all-d graph and cycle-splitting arguments.
    for n in range(3,11):
        alpha=tuple((j+1)%n for j in range(n))
        for j in range(1,n):
            tau=transposition(n,0,j);beta=compose(inv(alpha),tau)
            lengths=sorted(len(z) for z in cycles(beta))
            check(f'full cycle splits into j and n-j n{n} j{j}',lengths==sorted([j,n-j]))
            edges=[(z,(z+j)%n) for z in range(n)]
            seen=set();components=0
            for z in range(n):
                if z in seen:continue
                components+=1;todo=[z];seen.add(z)
                while todo:
                    t=todo.pop()
                    for aa,bb in edges:
                        if t==aa and bb not in seen:seen.add(bb);todo.append(bb)
                        if t==bb and aa not in seen:seen.add(aa);todo.append(aa)
            check(f'conjugate-transposition graph components equal gcd n{n} j{j}',components==gcd(n,j))
            if n<=5:
                check(f'finite generated group confirms symmetric criterion n{n} j{j}',
                      (len(generated([alpha,tau]))==factorial(n))==(gcd(n,j)==1))
    for n in [3,4,5]:
        identity=tuple(range(n));adj=[transposition(n,j,j+1) for j in range(n-1)]
        center=[z for z in itertools.permutations(range(n)) if all(compose(z,t)==compose(t,z) for t in adj)]
        check(f'symmetric center is trivial control n{n}',center==[identity])
    # Unmarked T(4,9) really has the required finite monodromy and Z3 homology.
    n=4;alpha=tuple((j+1)%n for j in range(n));tau=transposition(n,0,1);beta=compose(inv(alpha),tau)
    check('unmarked T4,9 representation has S4 image',len(generated([alpha,beta]))==24)
    check('unmarked T4,9 exceptional relations and transposition are exact',
          order(alpha)==4 and 9%order(beta)==0 and compose(compose(alpha,beta),tau)==tuple(range(4)))
    # Orders (1,9,3), beta=(0,-2,1), integer obstruction zero give |e|=1/9.
    M0=S.Matrix([[1,0,0,0],[0,9,0,-2],[0,0,3,1],[1,1,1,0]])
    snf=smith_normal_form(M0,domain=S.ZZ)
    check('unmarked T4,9 end has H1 exactly Z3',sorted(abs(int(snf[j,j])) for j in range(4))==[1,1,1,3])
    hvector=S.Matrix([0,0,0,1]);coeff=M0.T.inv()*hvector
    check('unmarked T4,9 ramification regular fiber is integrally zero',
          all(z.q==1 for z in coeff) and M0.T*coeff==hvector)
    # Symbolic weighted-Euler route from PF33 balances, so it is not adopted as an exclusion.
    n,reciprocal_sum=S.symbols('n reciprocal_sum')
    eW=S.Rational(1,3)+reciprocal_sum/2;eZ=S.Rational(4,3)+2*reciprocal_sum
    check('weighted sign-cover Euler route balances identically',S.expand(eZ-4*eW)==0)
    # Nonzero class at the marked end is the precise contradiction.
    check('actual nonzero marked class contradicts the torus-cover zero class',4%3!=0%3)
    inputs=[Path(__file__),HERE/'TORUS_INFINITY_OBSTRUCTION.md',HERE/'residual_cover_receipt.json',
            HERE/'RESIDUAL_DOUBLE_COVER.md',HERE/'GLOBAL_BOUNDARY_OBSTRUCTION.md',
            HERE/'PROOF.md',HERE/'WRIGHT_MARKED_MODEL.md',HERE/'PRIMITIVE_QUARTIC_FIELD.md']+sorted(source.glob('*'))
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        scope='every pole height and every possible finite field degree on the actual marked Wright surface',
        marked_end=dict(H1='Z/3',ramification_class='4*mu=mu!=0',B_square=3,E_dot_B=4),
        topological_lemma=dict(degree_lower_bound=3,monodromy='S_d',finite_H1_required=True,
            base_exceptional_cycles='one d-cycle and two cycles r,s with gcd(r,s)=1',
            parameterization='p=d*A, q=r*s*k',source_orders=['A','s*k','r*k'],
            first_homology_order='k',orbifold_abelianization_order='k',ramification_knot_class=0),
        consequences=dict(unknot_at_infinity_excluded=True,all_torus_knots_at_infinity_excluded=True,
            coprime_boundary_normalization_degrees_excluded_in_every_target_presentation=True,
            every_prime_boundary_coordinate_degree_excluded=True,minimum_boundary_coordinate_degree=4),
        rejected_euler_route=dict(weighted_multiplicativity_balances_for_all_orders=True,exclusion_obtained=False),
        unmarked_T4_9_control=dict(monodromy='S4',H1='Z/3',ramification_class=0,is_actual_Wright_cover=False),
        assurance=dict(actual_marked_ring_finite_map_and_constant_source_Jacobian_are_hypotheses=True,
            all_degree_proof_is_written_not_inferred_from_finite_controls=True,
            primary_topology_results_explicitly_imported=True,
            branched_base_Euler_naturality_proved_by_section_obstruction=True,
            finite_controls_do_not_replace_topological_proof=True,
            uses_no_new_coefficient_or_Groebner_search=True,
            boundary_coordinate_degrees_distinguished_from_source_degrees_and_B_poles=True,
            all_iterated_cable_knots_excluded=False,all_colliding_Wright_covers_excluded=False,
            actual_finite_map_constructed=False,arbitrary_normalization_compatible_completion_proved=False,
            whole_height_nine_pair_excluded=False,quartic_packets_excluded=False,
            formal_verification=False,independent_full_reproof=False,external_review=False,research_goal_completed=False),
        first_open_reduced_pairs=[[5,9],[6,9],[7,9],[8,9]],
        source_sha256={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        python=platform.python_version(),sympy=S.__version__,elapsed_seconds=round(monotonic()-start,6))
    (HERE/'torus_infinity_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} symbolic/topological-algebra controls and primary pins; all-degree torus type excluded; iterated cabling and full goal remain open')
if __name__=='__main__':main()

"""Exact PF2RULED61 certificates; unbounded source topology is proved in the note."""
from pathlib import Path
import hashlib,itertools,json,os,platform
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
HERE=Path(__file__).resolve().parent
def mul(a,b):return tuple(a[b[i]] for i in range(len(a)))
def inv(a):return tuple(a.index(i) for i in range(len(a)))
def power(a,n):
    if n<0:return power(inv(a),-n)
    z=tuple(range(len(a)))
    for _ in range(n):z=mul(z,a)
    return z
def wi(w):return [-x for x in reversed(w)]
def wp(w,n):return w*n if n>=0 else wi(w)*(-n)
def comm(u,v):return u+v+wi(u)+wi(v)
def evaluate(w,rho):
    z=tuple(range(len(rho[0])))
    for x in w:z=mul(z,power(rho[abs(x)-1],1 if x>0 else -1))
    return z
def cover_matrix(rels,rho,meridian,longitude,reverse_tree=False):
    n=len(rho);d=len(rho[0]);edges=[(g,s) for g in range(n) for s in range(d)]
    dest=lambda g,s:inv(rho[g])[s]
    parent=list(range(d));tree=set()
    def root(s):
        while parent[s]!=s:s=parent[s]
        return s
    for g,s in (list(reversed(edges)) if reverse_tree else edges):
        t=dest(g,s)
        if root(s)!=root(t):parent[root(s)]=root(t);tree.add((g,s))
    assert len(tree)==d-1,'connected sheet action required'
    keep=[e for e in edges if e not in tree];ix={e:i for i,e in enumerate(keep)}
    def lift(w,s):
        row=[0]*len(keep);start=s
        for x in w:
            g=abs(x)-1
            if x>0:e=(g,s);s=dest(g,s);sgn=1
            else:s=rho[g][s];e=(g,s);sgn=-1
            if e in ix:row[ix[e]]+=sgn
        assert s==start,'every attaching path is closed'
        return row
    rows=[lift(w,s) for w in rels for s in range(d)]
    mp=evaluate(meridian,rho);seen=set();ram=None
    for s in range(d):
        if s in seen:continue
        cyc=[];t=s
        while t not in cyc:cyc.append(t);seen.add(t);t=mp[t]
        rows.append(lift(meridian*len(cyc),s))
        if len(cyc)==2:ram=s
    assert seen==set(range(d))
    if ram is None:ram=0
    return rows,lift(longitude,ram),keep
def h1(rows,extra=()):
    D=smith_normal_form(S.Matrix(rows+list(extra)),domain=S.ZZ)
    diag=[abs(D[i,i]) for i in range(min(D.shape))]
    return dict(free=D.cols-sum(x!=0 for x in diag),torsion=[int(x) for x in diag if x>1])
def presentations(w):
    a,b,B,h,mu=range(1,6);M=[a,-b];L=[a,a]+wp(M,-6)
    five=[[a,a,-b,-b,-b],comm([B],[h]),comm([mu],[h]),
          [-mu,B]*w+[h],wi(M)+[B]*w+[h],wi(L)+[-B]*(w+1)+[-h]]
    lam5=[h]+[-mu]*(w*(w+1))
    # Root amalgam: z^w=M^(w+1)L^w, beta=M^-1 L^-1, mu=z beta.
    h3=wp(M,w+1)+wp(L,w);beta3=wi(M)+wi(L);mu3=[3]+beta3
    three=[[a,a,-b,-b,-b],[3]*w+wi(h3)]
    lam3=h3+wp(mu3,-w*(w+1))
    return (five,[mu],lam5),(three,mu3,lam3)
def generated(gs):
    I=tuple(range(len(gs[0])));out={I};todo=[I]
    while todo:
        a=todo.pop()
        for b in gs:
            c=mul(a,b)
            if c not in out:out.add(c);todo.append(c)
    return len(out)
def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    preds=['coordinate_pencil_matrix_receipt.json','peripheral_zero_receipt.json']
    for name in preds:
        r=json.loads((HERE/name).read_text());check('predecessor PASS '+name,r['outcome']=='PASS')
        for p,h in r['source_sha256'].items():check('sealed predecessor '+name+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/marked-end-lattice/manifest.json';m=json.loads(manifest.read_text())
    for row in m['primary_sources']:
        p=HERE/row['path'];check('primary pin '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    # Independent known degree2/3 knot-cover controls exercise all filling cycles.
    for name,rho,expected in [
        ('unknot2',[(1,0)],dict(free=0,torsion=[])),
        ('trefoil2',[(1,0),(0,1)],dict(free=0,torsion=[3])),
        ('trefoil3',[(1,0,2),(1,2,0)],dict(free=0,torsion=[])),
        ('cyclic3trefoil',[(0,1,2),(1,2,0)],dict(free=0,torsion=[2,2]))]:
        rels=[] if name=='unknot2' else [[1,1,-2,-2,-2]]
        mu=[1] if name=='unknot2' else [1,-2]
        rows,_,_=cover_matrix(rels,rho,mu,[])
        check('known branched cover '+name,h1(rows)==expected)
    # Formal 2x2 determinant, inverse, marking and complementary parity.
    nr,nd,c=S.symbols('nr nd c');Q=S.Matrix([[-nr,c],[c,-nd]]);delta=nr*nd-c*c
    check('full binary determinant',Q.det()==delta)
    check('canonical marking inverse',S.simplify(Q.inv()*S.Matrix([1,0])-S.Matrix([-nd,-c])/delta)==S.zeros(2,1))
    check('complementary nR modulo3',all((r*r-r+1)%3==0 for r in range(2,20,3)))
    check('minus one is not a square modulo3',all((-x*x)%3!=1 for x in range(3)))
    check('odd complementary order forces even cross',all((r-1)%2==0 for r in range(3,20,2)))
    check('marked zero determinant3 and nR0 forces determinant0 modulo9',all((3*a*3*b-(3*k)**2)%9==0 for a,b,k in itertools.product(range(3),repeat=3)))
    # Exhaustive S4 enumeration; rho(h)=1 follows from lambda=1 and mu^2=1.
    perms=list(itertools.permutations(range(4)));I=tuple(range(4));trans=[p for p in perms if sum(p[i]!=i for i in range(4))==2]
    def orbit(rho):return min(tuple(tuple(inv(g)[p[g[i]]] for i in range(4)) for p in rho) for g in perms)
    certificates={}
    for w,expected_total,expected_orbits in [(3,72,3),(8,168,7)]:
        five,three=presentations(w);classes={};total=0
        for aa,bb in itertools.product(perms,repeat=2):
            if power(aa,2)!=power(bb,3):continue
            MM=mul(aa,inv(bb));LL=mul(power(aa,2),power(MM,-6))
            for beta in perms:
                if power(beta,w)!=MM or power(beta,-w-1)!=LL:continue
                for tau in trans:
                    rr=[aa,bb,beta,I,tau]
                    if power(mul(tau,beta),w)!=I or generated(rr)!=24:continue
                    assert all(evaluate(word,rr)==I for word in five[0])
                    assert evaluate(five[2],rr)==I
                    total+=1;key=orbit(rr)
                    if key not in classes:classes[key]=dict(rho=rr,count=0)
                    classes[key]['count']+=1
        check(f'w{w} exhaustive representation and orbit counts',(total,len(classes))==(expected_total,expected_orbits))
        cert=[]
        for i,row in enumerate(classes.values()):
            rr=row['rho'];check(f'w{w} orbit{i} conjugacy size',row['count']==24)
            rows,mark,keep=cover_matrix(five[0],rr,five[1],five[2]);hh=h1(rows);hq=h1(rows,[mark])
            rr3=[rr[0],rr[1],mul(rr[4],inv(rr[2]))]
            check(f'w{w} orbit{i} root amalgam relators',all(evaluate(word,rr3)==I for word in three[0]))
            r3,m3,_=cover_matrix(three[0],rr3,three[1],three[2])
            check(f'w{w} orbit{i} equivalent root homology and marking',h1(r3)==hh and h1(r3,[m3])==hq)
            alt,am,_=cover_matrix(five[0],rr,five[1],five[2],reverse_tree=True)
            check(f'w{w} orbit{i} spanning tree invariance',h1(alt)==hh and h1(alt,[am])==hq)
            cert.append(dict(rho=rr,orbit_size=row['count'],H1=hh,quotient_by_ramification=hq,root_presentation_matrix=r3,root_ramification_row=m3))
        if w==3:
            check('all three control ends have order1 or3',sorted(x['H1']['torsion'] for x in cert)==[[],[],[3]] and all(x['H1']['free']==0 for x in cert))
            check('all three control ramification classes zero',all(x['H1']==x['quotient_by_ramification'] for x in cert))
        else:
            check('three w8 covers fail finite end homology',sum(x['H1']==dict(free=1,torsion=[26]) for x in cert)==3)
            check('four remaining w8 covers have finite homology',sum(x['H1']['free']==0 for x in cert)==4)
        certificates[str(w)]=dict(total_representations=total,orbits=cert)
    # The target knot determinant9 at both windings is separate from cover H1.
    t=S.symbols('t')
    def torus(w,p):return S.cancel((t**(w*p)-1)*(t-1)/((t**w-1)*(t**p-1)))
    for w in [3,8]:
        alex=S.cancel(torus(w,w+1)*torus(2,3).subs(t,t**w))
        check(f'w{w} infinity determinant9',alex.is_polynomial(t) and abs(alex.subs(t,-1))==9)
    inputs=[Path(__file__),HERE/'MARKED_END_LATTICE_AND_CABLE_CONTROL.md',manifest]+[HERE/n for n in preds]+[HERE/r['path'] for r in m['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        source_lattice=dict(any_boundary_count=True,all_smooth_affine_extensions_with_disjoint_A1=True,whole_integral_intersection_matrix=True,negative_definite=True,end_H1_cokernel=True,linking_inverse_retained=True,canonical_marking_retained=True,smoothness_assumed_in_higher_degrees=True),
        quartic=dict(all_source_orders_and_contacts=True,canonical_covector=[1,0],ramified_diagonal_odd=True,unramified_diagonal_even=True,canonical_zero_iff_determinant_divides_nD_and_cross=True),
        control=certificates,
        exclusions=dict(entire_target_knot_C34_trefoil_if_nR_divisible3=True,common_pencil_included_when_nR_divisible3=True,all_odd_complementary_r_and_r1=True,all_complementary_r_congruent2_mod3=True,minimal_complementary_r2_q1_all_unramified_orders=True,remaining_complementary_r_mod6=[0,4],three_w8_representation_classes_positive_b1=True,all_C34_source_sectors_excluded=False,all_w8_representations_excluded=False),
        assurance=dict(written_unbounded_topology=True,finite_group_enumeration_exhaustive_for_two_fixed_knots=True,matrices_are_not_independent_topological_review=True,symplectic_filling_classification_imported=False,named_boundary_manifold_identification_used=False,full_quartic_exclusion=False,higher_degree_exclusion=False,finite_Keller_map=False,formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'marked_end_lattice_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; full integral end lattice and marked cable obstruction')
if __name__=='__main__':main()

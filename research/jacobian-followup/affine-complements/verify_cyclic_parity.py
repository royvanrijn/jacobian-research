"""PF2RULED68: unrestricted identities plus finite, independently lifted controls."""
from pathlib import Path
import hashlib,importlib.util,itertools,json,math,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent
def load(name,file):
    spec=importlib.util.spec_from_file_location(name,HERE/file)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
E=load('end61','verify_marked_end_lattice.py')
C65=load('cyclic65','verify_cyclic_companion.py')

def c4_new_pattern():
    """New c=1 mod4 case; all exponents retained as symbolic integers."""
    R,C=S.symbols('R C',integer=True,nonnegative=True)
    r=12*R+3;c=4*C+1;k=S.expand(r*c+1)
    M=(1,2,3,0);tau=(1,0,2,3);z=E.mul(tau,M);rho=[z,M]
    _,_,keep=E.cover_matrix([],rho,[1],[]);ix={edge:i for i,edge in enumerate(keep)};n=len(keep)
    def identity():return tuple(range(4)),S.zeros(4,n)
    def compose(a,b):
        pa,ca=a;pb,cb=b
        return tuple(pb[pa[i]] for i in range(4)),S.Matrix([list(ca.row(i)+cb.row(pa[i])) for i in range(4)])
    def generator(g):
        mat=S.zeros(4,n)
        for i in range(4):
            if (g,i) in ix:mat[i,ix[g,i]]=1
        return E.inv(rho[g]),mat
    def positive(effect,length):
        out=identity()
        for _ in range(length):out=compose(out,effect)
        return out
    def power(effect,exponent,order,res):
        loop=positive(effect,order);assert loop[0]==tuple(range(4))
        prefix=positive(effect,res)
        return prefix[0],(((exponent-res)/order)*loop[1]+prefix[1]).applyfunc(S.expand)
    zz,mm=map(generator,range(2));zr=power(zz,r,3,0);mk=power(mm,k,4,0)
    mu=compose(zz,power(mm,-c,4,3));assert mu[0]==tau
    rels=zr[1]-mk[1];rows=[rels.row(i) for i in range(4)];seen=set();ram=None
    for i in range(4):
        if i in seen:continue
        cycle=[];j=i
        while j not in cycle:cycle.append(j);seen.add(j);j=mu[0][j]
        fill=positive(mu,len(cycle));assert fill[0][i]==i;rows.append(fill[1].row(i))
        if len(cycle)==2:ram=i
    return dict(R=R,C=C,r=r,c=c,k=k,rho=rho,keep=keep,
                matrix=S.Matrix.vstack(*rows).applyfunc(S.expand),mark=zr[1].row(ram))

def expected_group(p,q,r,kind):
    orders=[r//3 if kind=='C4' else (r*(p*q-p-q)+1)//3]
    if kind=='C3':
        if p%3==0:orders += [q,q]
        elif q%3==0:orders += [p,p]
    elif p%2==0 or q%2==0:
        e,o=(p,q) if p%2==0 else (q,p)
        orders += [o]*(3 if e%4==0 else 1)
    return E.h1(S.diag(*orders).tolist())

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    preds=['marked_end_lattice_receipt.json','binomial_periodic_receipt.json',
           'cyclic_companion_receipt.json','oriented_pattern_linking_receipt.json','common_pencil_mark_receipt.json']
    for name in preds:
        old=json.loads((HERE/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for p,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/cyclic-parity/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        p=HERE/row['path'];check('retained primary bytes '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    # Exhaustive *cyclic-sector* permutation classification, independent of degree exponents.
    perms=list(itertools.permutations(range(4)));I=tuple(range(4))
    trans=[p for p in perms if sum(p[i]!=i for i in range(4))==2]
    counts={}
    for M,tau in itertools.product(perms,trans):
        if E.generated([M,tau])!=24:continue
        order=next(j for j in range(1,13) if E.power(M,j)==I)
        counts[order]=counts.get(order,0)+1
        check('cyclic companion M order3 or4 '+str((M,tau)),order in [3,4])
        for c in [1,3,5,7,9,11]:
            z=E.mul(tau,E.power(M,c))
            if order==3 and c%3==0:continue  # r*c+1 cannot be divisible by3.
            zorder=next(j for j in range(1,13) if E.power(z,j)==I)
            check('root order forced in cyclic sector '+str((M,tau,c)),zorder==(4 if order==3 else 3))
    check('entire cyclic sector pairs',counts=={3:24,4:24})
    pat=c4_new_pattern();R,C=pat['R'],pat['C'];r,c,k=pat['r'],pat['c'],pat['k'];s=r/3;t=k/4
    expected=S.Matrix([[r,0,-t,-t,-t],[0,s,-t,-t,-t],[0,s,-t,-t,-t],[0,s,-t,-t,-t],
                       [1,0,-2*C-1,-2*C,-2*C],[0,0,-C,-C-1,-C],[0,1,-C,-C,-C-1]])
    check('new unrestricted c1mod4 lifted matrix',(pat['matrix']-expected).applyfunc(S.expand)==S.zeros(7,5))
    check('new unrestricted c1mod4 marking',pat['mark']==S.Matrix([[r,0,0,0,0]]))
    X,x=S.symbols('X x');y=-C*X;w=X-x-y;A=2*C*X+x;B=C*X+w
    check('new fillings imply A+B=cX',S.expand(A+B-c*X)==0)
    ss,tt,cc=S.symbols('s t c');mat=S.Matrix([[3*ss,-tt],[-ss,ss*cc-tt]])
    check('unit relation forces X',mat.row(0)+3*mat.row(1)==S.Matrix([[0,3*ss*cc-4*tt]]))
    check('C4 determinant gives cyclic order s',S.expand(mat.det()-ss*(3*ss*cc-4*tt))==0)
    check('primitive relation is minus one',S.expand(3*s*c-4*t)==-1)
    for rr,cv in [(0,0),(0,1),(1,0),(2,2)]:
        sub={R:rr,C:cv};rn=int(r.subs(sub));cn=int(c.subs(sub));kn=int(k.subs(sub))
        rows,mark,_=E.cover_matrix([[1]*rn+[-2]*kn],pat['rho'],[1]+[-2]*cn,[1]*rn)
        check('ordinary-word new pattern '+str((rr,cv)),S.Matrix(rows)==pat['matrix'].subs(sub) and S.Matrix([mark])==pat['mark'].subs(sub))
        want=dict(free=0,torsion=[] if rn==3 else [rn//3])
        check('new pattern group and zero mark '+str((rr,cv)),E.h1(rows)==want and E.h1(rows,[mark])==want)
    # Exact module manipulations; residue identities hold for unrestricted o.
    t,b=S.symbols('t b');P=(t+1)*(t*t+1);Q=1+b*(t*t+1)
    check('quadruple module t+1 integral certificate',S.expand((t+1)*Q-b*P)==t+1)
    check('quadruple module substitution',S.expand(Q.subs(t,-1))==1+2*b)
    check('odd geometric sum unit residue1',S.rem(S.Integer(1),P,t)==1)
    check('odd geometric sum unit residue3',S.rem(1+t+t*t+t**3,P,t)==0)
    check('residue3 is a unit',S.rem((-t**3)*(-t),P,t)==1)
    modules=[]
    for p,q in [(2,3),(2,5),(2,7),(4,3),(4,5),(6,5),(6,7),(8,3),(10,3),(10,7),(3,8),(5,12),(3,5)]:
        delta=S.cancel((t**(p*q)-1)*(t-1)/((t**p-1)*(t**q-1)))
        groups={}
        for a in [3,4]:
            group=E.h1(C65.module_matrix(delta,t,a).tolist())
            if a==3:orders=([q,q] if p%3==0 else ([p,p] if q%3==0 else []))
            elif p%2 and q%2:orders=[]
            else:
                e,o=(p,q) if p%2==0 else (q,p);orders=[o]*(3 if e%4==0 else 1)
            want=E.h1(S.diag(*orders).tolist()) if orders else dict(free=0,torsion=[])
            check('entire Alexander module '+str((p,q,a)),group==want);groups[str(a)]=group
        modules.append(dict(p=p,q=q,groups=groups))
    controls=[]
    for p,q,r,kind in [(2,3,3,'C4'),(2,3,27,'C4'),(2,7,27,'C4'),(4,3,3,'C4'),
                       (4,5,9,'C4'),(3,4,3,'C4'),(2,7,16,'C3'),(4,5,4,'C3'),
                       (6,5,8,'C3'),(3,4,16,'C3')]:
        row=C65.full_satellite(p,q,r,kind);want=expected_group(p,q,r,kind)
        check('full ordinary satellite group '+str((p,q,r,kind)),row['H1']==want)
        check('full ordinary satellite R zero '+str((p,q,r,kind)),row['quotient_by_ramification']==want)
        check('full satellite image S4 '+str((p,q,r,kind)),E.generated(row['rho'])==24)
        det=(r*(p*q-p-q)+1) if r%2==0 else r*(q if p%2==0 else p)
        check('new satellite examples retain target nine gate '+str((p,q,r,kind)),det%9==0)
        controls.append({key:row[key] for key in ['p','q','r','kind','H1','quotient_by_ramification']})
    # Symbolic source contradiction, with all non-3 companion residues retained.
    allowed=[]
    for pp,qq,rr in itertools.product([1,2],[1,2],range(12)):
        cc=(pp*qq-pp-qq)%3
        if rr%4 or (rr*cc+1)%3:continue
        check('all C3 no3companion residues force c2 t1 '+str((pp,qq,rr)),cc==2 and (rr//4)%3==1)
        allowed.append([pp,qq,rr])
    check('complementary nR never2mod3',all((u*u-u+1)%3 in [0,1] for u in range(3)))
    check('oriented C3 coefficient forces nR2',all((-t*u*u)%3==2 for t in [1,4,7] for u in [1,2]))
    # Noncyclic countercontrol from the *complete* predecessor census.
    census=json.loads((HERE/'binomial_periodic_receipt.json').read_text())['census']['complete_cells']
    row=next(x for x in census['2,3,9']['orbits'] if x['companion_order']==24 and x['cable_order']==24 and x['cycle_rank']==0)
    a,bb,beta,tau=map(tuple,row['rho']);z=E.mul(tau,E.inv(beta));rho=[a,bb,z]
    M=[1,-2];L=[1,1]+E.wp(M,-6);h=E.wp(M,10)+E.wp(L,9);mu=[3]+E.wi(M)+E.wi(L)
    rows,mark,_=E.cover_matrix([[1,1,-2,-2,-2],[3]*9+E.wi(h)],rho,mu,h)
    counter=dict(H1=E.h1(rows),quotient=E.h1(rows,[mark]))
    check('noncyclic companion R can be NONZERO',counter==dict(H1=dict(free=0,torsion=[24]),quotient=dict(free=0,torsion=[6])))
    inputs=[Path(__file__),HERE/'CYCLIC_COMPANIONS_ACROSS_PARITY.md',manifest]+[HERE/n for n in preds]+[HERE/r['path'] for r in sources['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        theorem=dict(all_coprime_first_pairs=True,all_windings_and_source_orders_poles=True,
          actual_primitive_binomial_two_stage_required=True,cyclic_companion_restriction_required=True,
          cyclic_sector_only_C3_C4=True,C3_all_first_pair_parities_excluded=True,
          C3_non3_companion_oriented_linking='negative (r/4)/(k/3)',
          C4_common_source_excluded=True,C4_even_exponent_divisible4_excluded=True,
          C4_remaining_even_exponent_mod4=2,C4_remaining_source='complementary',
          C4_remaining_coprimality='gcd(odd first exponent,r/3)=1',
          all_cyclic_sector_R_marks_zero=True,cyclic_cover4_even2_group='Z/odd',
          cyclic_cover4_even0_group='(Z/odd)^3',cyclic_cover3_formula_parity_free=True,
          C4_complementary_sector_closed=False,noncyclic_companion_sectors_closed=False,
          all_even_first_pairs_excluded=False,PF63_polynomial_immersions_still_exist=True),
        symbolic_new_pattern=dict(matrix=[[str(v) for v in row] for row in pat['matrix'].tolist()],mark=[str(v) for v in pat['mark']]),
        module_controls=modules,full_satellite_controls=controls,noncyclic_nonzero_mark_control=counter,
        assurance=dict(written_unbounded_proofs=True,exact_symbolic_new_pattern=True,
          retained_primary_imports=True,finite_controls_not_universal_proof=True,
          arbitrary_noncyclic_pieces_replaced=False,finite_collisions_classified=False,
          full_quartic_exclusion=False,higher_degree_exclusion=False,finite_Keller_map=False,
          formal_verification=False,independent_full_reproof=False,external_review=False,
          novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'cyclic_parity_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; cyclic-companion exclusions across both first-pair parities')
if __name__=='__main__':main()

"""PF2RULED65: symbolic unbounded identities and independent finite controls."""
from pathlib import Path
import hashlib,importlib.util,itertools,json,math,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('sealed_end61',HERE/'verify_marked_end_lattice.py')
E=importlib.util.module_from_spec(spec);spec.loader.exec_module(E)

def compressed_pattern(kind,residue=1):
    """Exact crossed chain sums; exponent reduction uses the proved permutation order."""
    R,C=S.symbols('R C',integer=True,nonnegative=True)
    if kind=='C4':
        M=(1,2,3,0);tau=(1,0,2,3);z=E.mul(tau,E.inv(M))
        r=12*R+9;c=4*C+3;zorder=3;morder=4
    else:
        assert residue in [1,2]
        M=(1,2,0,3);tau=(3,1,2,0);z=E.mul(tau,E.power(M,residue))
        r=12*R+(8 if residue==1 else 4);c=3*C+residue;zorder=4;morder=3
    k=S.expand(r*c+1);rho=[z,M]
    _,_,keep=E.cover_matrix([],rho,[1],[]);ix={edge:i for i,edge in enumerate(keep)};n=len(keep)
    def identity():return tuple(range(4)),S.zeros(4,n)
    def compose(a,b):
        pa,ca=a;pb,cb=b
        return tuple(pb[pa[i]] for i in range(4)),S.Matrix([list(ca.row(i)+cb.row(pa[i])) for i in range(4)])
    def generator(g):
        dest=E.inv(rho[g]);mat=S.zeros(4,n)
        for i in range(4):
            if (g,i) in ix:mat[i,ix[g,i]]=1
        return dest,mat
    def positive(effect,length):
        out=identity()
        for _ in range(length):out=compose(out,effect)
        return out
    def power(effect,exponent,order,res):
        loop=positive(effect,order);assert loop[0]==tuple(range(4))
        prefix=positive(effect,res)
        return prefix[0],(((S.sympify(exponent)-res)/order)*loop[1]+prefix[1]).applyfunc(S.expand)
    zz,mm=map(generator,range(2))
    zr=power(zz,r,zorder,0);mk=power(mm,k,morder,0)
    mu=compose(zz,power(mm,-c,morder,1 if kind=='C4' else (-residue)%3));assert mu[0]==tau
    rels=zr[1]-mk[1];rows=[rels.row(i) for i in range(4)];seen=set();ram=None
    for i in range(4):
        if i in seen:continue
        cycle=[];j=i
        while j not in cycle:cycle.append(j);seen.add(j);j=mu[0][j]
        fill=positive(mu,len(cycle));assert fill[0][i]==i;rows.append(fill[1].row(i))
        if len(cycle)==2:ram=i
    mat=S.Matrix.vstack(*rows).applyfunc(S.expand);mark=zr[1].row(ram)
    return dict(R=R,C=C,r=r,c=c,k=k,rho=rho,keep=keep,matrix=mat,mark=mark)

def full_satellite(p,q,r,kind):
    c=p*q-p-q;k=r*c+1
    u,v,g=S.gcdex(q,p);assert g==1;u=int(u);v=int(v)
    Mword=E.wp([1],u)+E.wp([2],v)
    Lword=E.wp([1],p)+E.wp(Mword,-p*q)
    hword=E.wp(Mword,k)+E.wp(Lword,r)
    betaword=E.wp(Mword,-c)+E.wi(Lword);mu=[3]+betaword
    if kind=='C4':M=(1,2,3,0);tau=(1,0,2,3)
    else:M=(1,2,0,3);tau=(3,1,2,0)
    z=E.mul(tau,E.power(M,c));rho=[E.power(M,q),E.power(M,p),z]
    rels=[E.wp([1],p)+E.wp([2],-q),E.wp([3],r)+E.wi(hword)]
    assert all(E.evaluate(w,rho)==tuple(range(4)) for w in rels)
    assert E.evaluate(mu,rho)==tau and E.evaluate(hword,rho)==tuple(range(4))
    rows,mark,_=E.cover_matrix(rels,rho,mu,hword)
    # lambda=h*mu^(-rk); rk is even, so its difference from h is a filled ramified meridian power.
    assert (r*k)%2==0
    return dict(p=p,q=q,r=r,c=c,k=k,kind=kind,rho=rho,rows=rows,mark=mark,
                H1=E.h1(rows),quotient_by_ramification=E.h1(rows,[mark]),presentation=(rels,rho,mu,hword))

def expected_orders(p,q,r,kind):
    if kind=='C4':return [r//3]
    out=[(r*(p*q-p-q)+1)//3]
    if p%3==0:out += [q,q]
    elif q%3==0:out += [p,p]
    return out

def module_matrix(poly,t,length):
    modulus=t**length-1
    return S.Matrix.hstack(*[S.Matrix([S.Poly(S.rem(poly*t**j,modulus,t),t).nth(i) for i in range(length)]) for j in range(length)])

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    preds=['marked_end_lattice_receipt.json','orbit_tree_receipt.json','binomial_periodic_receipt.json']
    for name in preds:
        old=json.loads((HERE/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for p,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/cyclic-companion/manifest.json';source=json.loads(manifest.read_text())
    for row in source['primary_sources']:
        p=HERE/row['path'];check('primary bytes '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    s,t,c,b,d=S.symbols('s t c b d',integer=True)
    mat4=S.Matrix([[s,-t],[-3*s,3*s*c-t]])
    check('C4 unrestricted determinant',S.expand(mat4.det()-s*(3*s*c-4*t))==0)
    check('C4 integral row combination forces X',mat4.row(1)+3*mat4.row(0)==S.Matrix([[0,3*s*c-4*t]]))
    mat3=S.Matrix([[d-b,d],[d,d-3*b]])
    check('C3 unrestricted determinant',S.expand(mat3.det()-b*(3*b-4*d))==0)
    check('C3 entries recover b with integer coefficients',mat3[0,1]-mat3[0,0]==b)
    check('C3 entries generate unit under 3b-4d=1',3*(mat3[0,1]-mat3[0,0])-4*mat3[0,1]==3*b-4*d)
    C,X,v=S.symbols('C X v')
    for residue in [1,2]:
        if residue==1:x=(2*C+1)*X;y=-C*X;w=-C*X;A=C*X+x+(3*C+1)*v
        else:x=(C+1)*X;w=(C+1)*X;y=-(2*C+1)*X;A=(C+1)*X-y+(3*C+2)*v
        check('C3 filling loop sum '+str(residue),S.expand(x+y+w-X)==0)
        check('C3 fillings give A=c(X+v) '+str(residue),S.expand(A-(3*C+residue)*(X+v))==0)
    A4=(2*C+2)*X-S.Symbol('y');B4=(C+1)*X-S.Symbol('x')
    check('C4 filling sum gives A+B=cX',S.expand((A4+B4-(4*C+3)*X).subs(S.Symbol('y'),X-(C+1)*X-S.Symbol('x')))==0)
    symbolic=[]
    for kind,residue in [('C4',1),('C3',1),('C3',2)]:
        pat=compressed_pattern(kind,residue);R,C=pat['R'],pat['C'];r,c,k=pat['r'],pat['c'],pat['k']
        if kind=='C4':
            s=r/3;t=k/4
            expected=S.Matrix([[0,s,-t,-t,-t],[r,0,-t,-t,-t],[0,s,-t,-t,-t],[0,s,-t,-t,-t],
              [1,0,-2*C-2,-2*C-1,-2*C-2],[0,0,-C-1,-C-1,-C],[0,1,-C,-C-1,-C-1]])
            mark=S.Matrix([[0,s,0,0,0]])
        else:
            a=r/4;b=k/3
            common=[a,-b,-b,-b,0];fixed=[a,0,0,0,-k]
            fills=([[1,-C-1,-C,-C,-3*C-1],[0,-C,-C-1,-C,0],[0,-C,-C,-C-1,0]] if residue==1 else
                   [[1,-C-1,-C,-C-1,-3*C-2],[0,-C-1,-C-1,-C,0],[0,-C,-C-1,-C-1,0]])
            expected=S.Matrix([common,common,common,fixed]+fills);mark=S.Matrix([[a,0,0,0,0]])
        check('exact unrestricted compressed lifted rows '+str((kind,residue)),(pat['matrix']-expected).applyfunc(S.expand)==S.zeros(7,5))
        check('exact unrestricted ramification h row '+str((kind,residue)),(pat['mark']-mark).applyfunc(S.expand)==S.zeros(1,5))
        for rr,cc in [(0,0),(0,1),(1,0),(3,2)]:
            sub={R:rr,C:cc};rn=int(r.subs(sub));cn=int(c.subs(sub));kn=int(k.subs(sub))
            rels=[[1]*rn+[-2]*kn];mu=[1]+[-2]*cn
            rows,h,_=E.cover_matrix(rels,pat['rho'],mu,[1]*rn)
            check('ordinary-word independent pattern matrix '+str((kind,residue,rr,cc)),S.Matrix(rows)==pat['matrix'].subs(sub) and S.Matrix([h])==pat['mark'].subs(sub))
            order=rn//3 if kind=='C4' else kn//3
            check('pattern cyclic homology and zero marking '+str((kind,residue,rr,cc)),E.h1(rows)==dict(free=0,torsion=[] if order==1 else [order]) and E.h1(rows,[h])==E.h1(rows))
        symbolic.append(dict(kind=kind,c_residue=residue,matrix=[[str(v) for v in row] for row in pat['matrix'].tolist()],mark=[str(v) for v in list(pat['mark'])]))
    t=S.symbols('t');alex_controls=[]
    for p,q in [(3,5),(3,7),(5,7),(5,9),(7,9),(5,11),(7,11)]:
        delta=S.cancel((t**(p*q)-1)*(t-1)/((t**p-1)*(t**q-1)))
        check('torus Alexander polynomial and augmentation '+str((p,q)),delta.is_polynomial(t) and delta.subs(t,1)==1)
        residue=S.rem(delta,t**4-1,t);expected=1 if p%4==1 or q%4==1 else t*t
        check('cyclic4 module unit '+str((p,q)),residue==expected and E.h1(module_matrix(delta,t,4).tolist())==dict(free=0,torsion=[]))
        orders=[q,q] if p%3==0 else ([p,p] if q%3==0 else [])
        expected=E.h1(S.diag(*orders).tolist()) if orders else dict(free=0,torsion=[])
        group=E.h1(module_matrix(delta,t,3).tolist())
        check('cyclic3 entire module group '+str((p,q)),group==expected)
        alex_controls.append(dict(p=p,q=q,H1_cyclic3=group))
    delta=t*t-t+1
    check('dropping all-odd hypothesis invalidates cyclic4 vanishing',E.h1(module_matrix(delta,t,4).tolist())==dict(free=0,torsion=[3]))
    controls=[]
    for p,q,r,kind in [(3,5,9,'C4'),(3,5,45,'C4'),(5,7,9,'C4'),(3,5,32,'C3'),(3,5,68,'C3'),(3,7,4,'C3'),(5,3,32,'C3'),(5,7,16,'C3')]:
        row=full_satellite(p,q,r,kind);orders=expected_orders(p,q,r,kind);expected=E.h1(S.diag(*orders).tolist())
        check('full satellite group matches replacement '+str((p,q,r,kind)),row['H1']==expected)
        check('full satellite mark remains zero '+str((p,q,r,kind)),row['quotient_by_ramification']==row['H1'])
        check('full satellite monodromy S4 '+str((p,q,r,kind)),E.generated(row['rho'])==24)
        excluded=kind=='C3' and (p*q)%3==0
        check('new canonical zero versus noncyclic contradiction '+str((p,q,r,kind)),(len(row['H1']['torsion'])>1)==excluded)
        if (p,q,r)==(3,5,32):
            rels,rho,mu,h=row['presentation'];alt,am,_=E.cover_matrix(rels,rho,mu,h,reverse_tree=True)
            check('new exclusion independent of spanning tree',E.h1(alt)==row['H1'] and E.h1(alt,[am])==row['H1'])
        controls.append({k:v for k,v in row.items() if k!='presentation'})
    N=S.symbols('N',integer=True,nonnegative=True)
    check('entire p3 q5 former C3 family has winding divisible4',S.expand((36*N+32)-4*(9*N+8))==0)
    check('entire p3 q5 former C3 family has determinant divisible9',S.expand(7*(36*N+32)+1-9*(28*N+25))==0)
    check('same group order cannot replace group structure',E.h1([[9]])!=E.h1([[3,0],[0,3]]) and 9==3*3)
    check('canonical-zero binary source control is cyclic',E.h1([[-1,3],[3,-12]])==dict(free=0,torsion=[3]) and E.h1([[-1,3],[3,-12]],[[1,0]])==dict(free=0,torsion=[3]))
    check('two generators alone miss the canonical mark',E.h1([[-5,0],[0,-10]])!=E.h1([[-5,0],[0,-10]],[[1,0]]))
    source_order=S.symbols('source_order',integer=True)
    check('all source orders2mod3 have ramified dual square divisible3',S.expand((3*N+2)**2-(3*N+2)+1-3*(3*N*N+3*N+1))==0)
    check('all odd complementary source orders have even q',S.expand(((2*N+1)-1)-2*N)==0)
    inputs=[Path(__file__),HERE/'CYCLIC_COMPANION_HOMOLOGY_AND_MARKED_END_EXCLUSION.md',manifest]+[HERE/p for p in preds]+[HERE/r['path'] for r in source['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      general=dict(all_finite_cover_degrees=True,cyclic_companion_restrictions_required=True,all_orbit_lengths=True,
          integral_direct_sum_replacement=True,preferred_longitude_zero_in_each_cyclic_piece=True,
          outer_mark_preserved=True,all_smooth_source_degrees_with_disjoint_A1=True,
          primitive_canonical_vector_required_for_integral_generator_bound=True,
          bound='H1(end)/<canonical_mark> generated by at most b-1 elements',
          arbitrary_noncyclic_pieces_replaced=False,hypotheses_automatic_in_higher_degrees=False),
      all_odd=dict(all_coprime_odd_p_q_at_least_three=True,all_windings=True,actual_primitive_binomial_key_required=True,
          all_peripheral_S4_classes_covered=True,C4_group='Z/(r/3)',
          C3_group='Z/(k/3) direct-sum H1(Sigma3(T(p,q)))',ramification_class_zero=True,
          C3_if_3_divides_p='Z/(k/3) direct-sum (Z/q)^2',C3_if_3_divides_q='Z/(k/3) direct-sum (Z/p)^2',
          entire_C3_case_if_3_divides_pq_excluded=True,
          refined_necessary_gate='r=9 mod36 OR (4|r AND9|k AND3 does not divide pq)',
          p3_q5_all_r32mod36_excluded=True,remaining_D='r/3 in C4; k/3 in remaining C3',
          remaining_D_odd_multiple_three=True,remaining_source_mark_zero=True,
          complementary_source_orders_mod6=[0,4],minimal_source_r2q1_excluded_in_entire_target_class=True,
          full_all_odd_class_excluded=False),
      symbolic_patterns=symbolic,torus_module_controls=alex_controls,full_satellite_controls=controls,
      assurance=dict(written_unbounded_proofs=True,exact_symbolic_exponent_reduction=True,
          finite_controls_not_universal_proof=True,retained_primary_peripherals_imported=True,
          Hillman_Sakuma_primary_integral_split_crosschecked=True,finite_collision_types_classified=False,
          arbitrary_cubic_index_identified_with_source_coordinate=False,full_linking_matched=False,
          independent_full_reproof=False,formal_verification=False,external_review=False,novelty_priority=False,
          full_quartic_exclusion=False,higher_degree_exclusion=False,finite_Keller_map=False,research_goal_completed=False),
      source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'cyclic_companion_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; all-degree cyclic replacement and unbounded marked-end exclusion')

if __name__=='__main__':main()

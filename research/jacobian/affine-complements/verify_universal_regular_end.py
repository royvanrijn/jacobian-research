"""PF2RULED43 exact supporting controls; universal geometric proofs are in the note.
No bounded support census proves the affine or end theorems.
"""
from pathlib import Path
from math import gcd
from functools import reduce
from heapq import heappop,heappush
import hashlib,json,platform
import sympy as S
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def data(r,support):
    support=sorted(set(support),reverse=True);q=support[0];g=r;chain=[];partial=0;generators=[r];pairs=[]
    for j in support:
        d=gcd(g,j)
        if d==g:continue
        beta=2*q-j;partial+=(g-d)*beta
        wind=g//d
        alpha=beta//d if not chain else wind*pairs[-1][0]*pairs[-1][1]+(beta-chain[-1]['beta'])//d
        if not chain:generators.append(beta)
        else:generators.append((chain[-1]['before']//chain[-1]['after'])*generators[-1]+beta-chain[-1]['beta'])
        pairs.append([wind,alpha]);chain.append(dict(j=j,before=g,after=d,beta=beta))
        g=d
        if g==1:break
    assert g==1
    c=partial-r+1;F=pairs[-1][0]*pairs[-1][1];k=chain[-1]['j']-r+2
    return dict(r=r,q=q,support=support,chain=chain,S=partial,conductor=c,F=F,kappa=k,n=F+k,pairs=pairs,generators=generators)
def apery(generators):
    m=min(generators);dist=[10**12]*m;dist[0]=0;heap=[(0,0)]
    while heap:
        d,a=heappop(heap)
        if d!=dist[a]:continue
        for g in generators:
            b=(a+g)%m;v=d+g
            if v<dist[b]:dist[b]=v;heappush(heap,(v,b))
    assert all(d<10**12 for d in dist)
    conductor=max(dist)-m+1
    genus=sum(dist)//m-(m-1)//2 if m%2 else S.Rational(sum(dist),m)-S.Rational(m-1,2)
    return conductor,int(genus)
z=S.symbols('z')
def order(expr):
    a,b=S.fraction(S.cancel(expr))
    def val(poly):return min(mon[0] for mon,c in S.Poly(poly,z).terms() if c)
    return val(a)-val(b)
def leading(expr):
    a,b=S.fraction(S.cancel(expr));pa=S.Poly(a,z);pb=S.Poly(b,z)
    va=min(m[0] for m,c in pa.terms() if c);vb=min(m[0] for m,c in pb.terms() if c)
    return pa.nth(va)/pb.nth(vb)
def resolution(r,support):
    q=max(support);H=sum(z**(q-j) for j in support)
    aa=z**q/H;bb=z**r;la='Dx';lb='L';weights={'Dx':0,'L':0};edges={frozenset(['Dx','L'])}
    vw={'Dx':1,'L':0};vy={'Dx':0,'L':1};K={'Dx':0,'L':0};exceptional=[]
    for it in range(1,160):
        labs=[v for v in [la,lb] if v]
        if len(labs)==1 and labs[0] in exceptional:
            if order(aa if la else bb)==1:return weights,edges,exceptional,labs[0],vw,vy,K
        va,vb=order(aa),order(bb)
        assert va>0 and vb>0
        new='E'+str(it)
        for lab in labs:weights[lab]-=1
        if len(labs)==2:edges.remove(frozenset(labs))
        weights[new]=-1
        for lab in labs:edges.add(frozenset([lab,new]))
        vw[new]=sum(vw[k] for k in labs);vy[new]=sum(vy[k] for k in labs);K[new]=1+sum(K[k] for k in labs)
        exceptional.append(new)
        if va<=vb:
            ratio=S.cancel(bb/aa)
            if va==vb:bb=S.cancel(ratio-leading(ratio));lb=None
            else:bb=ratio
            la=new
        else:
            aa=S.cancel(aa/bb);lb=new
    raise AssertionError('resolution did not terminate')
def matrix(names,weights,edges):
    return S.Matrix([[-weights[a] if a==b else -int(frozenset([a,b]) in edges) for b in names] for a in names])
def cofactor(C,i):return C.minor_submatrix(i,i).det() if C.rows>1 else 1
def main():
    pred_names=['coprime_regular_end_receipt.json','higher_regular_coordinate_receipt.json','torus_infinity_receipt.json','regular_sector_end_receipt.json']
    for name in pred_names:
        report=json.loads((HERE/name).read_text());check('predecessor PASS '+name,report['outcome']=='PASS')
        for f,h in report['source_sha256'].items():check('predecessor seal '+name+' '+f,digest(HERE/f)==h)
    manifest=HERE/'sources/universal-regular-end/manifest.json';pins=json.loads(manifest.read_text())
    check('five primary sources with explicit restricted imports',len(pins['primary_sources'])==5)
    for row in pins['primary_sources']:
        p=HERE/row['path_relative_to_affine_complements']
        check('primary bytes '+p.name,p.read_bytes().startswith(b'%PDF') and p.stat().st_size==row['bytes'] and digest(p)==row['sha256'])
    examples=[(2,[1]),(3,[2]),(4,[3]),(3,[4]),(4,[6,5]),(4,[6,3]),(6,[10,8,5]),(6,[12,10,5]),(8,[12,10,9]),(8,[12,10,7]),(12,[18,16,13]),(12,[18,16,11]),(5,[10,9]),(6,[18,16,15,11])]
    controls=[]
    for r,support in examples:
        d=data(r,support);q=d['q'];tag=f' r{r} support{support}'
        check('admissible full Laurent data'+tag,min(support)>=r-1 and reduce(gcd,[r]+support)==1)
        contacts=[2*q-max(j for j in support if a*j%r) for a in range(1,r)]
        check('individual cyclic contacts equal gcd-block contacts'+tag,sum(contacts)==d['S'] and all(contacts.count(c['beta'])==c['before']-c['after'] for c in d['chain']))
        for a in range(1,r):
            j=max(j for j in support if a*j%r)
            check('both boundary directional translations are polar'+tag+f' a{a}',j-(r-2)>=1)
        check('both source-boundary directional translations are polar'+tag,q>=1 and q-r+2>=1)
        conductor,genus=apery(d['generators'])
        check('independent Apery conductor and symmetry'+tag,conductor==d['conductor'] and 2*genus==conductor)
        check('last cabling framing equals contact sum'+tag,d['F']==d['S']+d['chain'][-1]['beta'])
        check('odd universal surgery and terminal-distance gap'+tag,d['n']==conductor+2*q+1 and d['n']%2==1 and d['n']>=3 and d['kappa']>=1 and d['kappa']!=2)
        higher=[j for j in support if j>r-1]
        check('distance one iff essential lowest symmetry breaker'+tag,(d['kappa']==1)==((r-1 in support) and reduce(gcd,[r]+higher)>1))
        if len(d['chain'])==1:
            beta=d['chain'][0]['beta'];check('single-pair formula agrees with PF2RULED42'+tag,d['n']==r*beta+d['kappa'])
        if len(d['chain'])>1 and d['kappa']==1:
            wind=d['pairs'][-1][0];prev=d['pairs'][-2][0]*d['pairs'][-2][1]
            gap=d['n']-wind**2*prev
            check('distance-one collapse has reduced denominator and next distance at least seven'+tag,gcd(d['n'],wind**2)==1 and gap>=wind*(wind+1)+1>=7)
            d['collapsed_companion_slope']=[d['n'],wind**2];d['next_fiber_distance']=gap
        controls.append(d)
    for r,support in [(4,[6,5]),(4,[6,3]),(6,[12,10,5])]:
        d=data(r,support);q=max(support);w,y=S.symbols('w y');H=sum(z**(q-j) for j in support)
        f=S.resultant(y-z**r,w*H-z**q,z)
        polar=S.cancel(S.diff(f,w).subs({y:z**r,w:z**q/H}))
        check('actual reduced reference resultant polar '+str((r,support)),order(polar)==d['S'])
    # Actual successive point blowups, independent of the support recurrence.
    graphs=[]
    for r,support in [(2,[1]),(3,[4]),(4,[6,5]),(4,[6,3]),(6,[12,10,5]),(8,[12,10,9])]:
        d=data(r,support);tag=str((r,support));weights,edges,names,E,vw,vy,K=resolution(r,support)
        C0=matrix(names,weights,edges);idx=names.index(E)
        check('actual embedded-resolution exceptional tree is unimodular '+tag,C0.det()==1 and weights[E]==-1)
        check('last curvette has correct source orders and discrepancy '+tag,vw[E]==d['q'] and vy[E]==r and K[E]==r+d['chain'][-1]['beta']-1)
        check('actual last curvette framing equals inverse marked entry '+tag,cofactor(C0,idx)==d['F'])
        # Include the zero leaf at the original infinity curve: one remaining
        # exceptional attachment is precisely the retained R6 situation.
        check('original deleted infinity curve has one exceptional attachment '+tag,sum('Dx' in e and bool(set(e)&set(names)) for e in edges)==1)
        prev=E
        for k in range(1,d['kappa']+1):
            new='F'+str(k);weights[prev]-=1;weights[new]=-1;edges.add(frozenset([prev,new]));prev=new
        mark=E if d['kappa']==1 else 'F'+str(d['kappa']-1)
        deleted=names+['F'+str(k) for k in range(1,d['kappa'])];C=matrix(deleted,weights,edges)
        check('actual whole end determinant and primitive marked cofactor '+tag,C.det()==d['n'] and cofactor(C,deleted.index(mark))==d['n']-1)
        graphs.append(dict(r=r,support=support,exceptional_vertices=len(names),F=int(cofactor(C0,idx)),n=int(C.det()),marked_cofactor=int(cofactor(C,deleted.index(mark)))))
    # Symbolic chain algebra proves the determinant identity for arbitrary
    # base inverse entry F after eliminating the unimodular reference tree.
    F=S.symbols('F',positive=True)
    for k in range(1,10):
        if k==1:det=1+F;co=F
        else:
            T=S.diag(*([2]*(k-1)))
            for i in range(k-2):T[i,i+1]=T[i+1,i]=-1
            T[0,0]-=F/(1+F)
            det=S.cancel((1+F)*T.det());co=S.cancel((1+F)*(T[:-1,:-1].det() if k>2 else 1))
        check('symbolic Schur determinant and marked cofactor k'+str(k),S.expand(det-(F+k))==0 and S.expand(co-(F+k-1))==0)
    a,b,t,x=S.symbols('a b t x',nonzero=True)
    for r in [2,3,4,6]:
        q=r+2;h=z**(-q)+2*z**(-(r-1));poly=1+z**r+3*z**(2*r)
        lhs=b*(z**(2-r)*t+h)+poly
        rhs=(a*z)**(2-r)*(b*a**(r-2)*t+a**(r-2)*z**(r-2)*poly)+b*h
        check('marked equivalence is an actual polynomial boundary chart substitution r'+str(r),S.cancel(lhs-rhs)==0 and r-2>=0)
    inputs=[Path(__file__),HERE/'UNIVERSAL_REGULAR_CLASSIFICATION_AND_END.md',manifest]+[HERE/r['path_relative_to_affine_complements'] for r in pins['primary_sources']]+[HERE/n for n in pred_names]
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        affine_classification=dict(necessary_and_sufficient=True,all_Laurent_supports=True,r='r>=2',range='q>=r-1',support_gcd=1,Dubouloz_section='2.8',both_directions_all_pairs_polar=True,smooth_affine_quotient=True,marked_equivalence="c'_j=b*a^j*c_j",all_polynomial_source_shears_extend=True,volume_exact_iff='c_r=0',target_algebra_invariance_assumed=False),
        universal_end=dict(source_reference='Gamma: w=1/h(z), y=z^r',characteristic_beta='2q-j_i at each strict support gcd drop',conductor='2delta=sum_i(g_(i-1)-g_i)*beta_i-r+1',F='sum_i(g_(i-1)-g_i)*beta_i+beta_last',kappa='j_last-r+2',surgery='n=2delta+2q+1',H1='Z/n',n_odd=True,marked_class_primitive=True,plumbing_pair_identified=True,stable_peripheral_data_transfer_to_radial_section=True,arbitrary_multicharacteristic_radial_homeomorphism_claimed=False),
        terminal_structure=dict(distance_two_impossible=True,distance_one_iff='essential c_(r-1) and higher support gcd>1',multi_distance_one_companion_filling='n/p_last^2',next_fiber_distance_at_least=7,no_consecutive_distance_one_collapse=True,distance_at_least_three_model_piece='Seifert D2(p_last,kappa)',actual_radial_JSJ_recognition_claimed=False),
        all_degree_consequences=dict(entire_regular_coordinate_sector=True,target_torus_or_unknot_infinity_excluded=True,target_coprime_boundary_degrees_excluded=True,target_prime_boundary_degrees_excluded=True,all_Keller_maps_excluded=False),
        controls=controls,resolution_graph_controls=graphs,
        assurance=dict(written_unbounded_proofs=True,primary_imports_explicit=True,finite_checks_not_complete_geometric_proof=True,general_affine_converse_closed=True,multi_characteristic_quartic_sector_excluded=False,no_regular_coordinate_case_excluded=False,extra_omitted_primes_excluded=False,finite_Keller_map_constructed=False,formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,research_goal_completed=False),
        source_sha256={str(p.relative_to(HERE)):digest(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'universal_regular_end_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS '+str(len(checks))+' universal regular-coordinate controls; affine classification closed; finite-map goal active')
if __name__=='__main__':main()

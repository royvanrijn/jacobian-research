#!/usr/bin/env sage-python
"""Independent full-frame / four-chart replay and smooth-reduction inputs.

Manual short group law and polynomial division; no producer import,
EllipticCurve, nullspace solver, point search or new orbit/prime exposure.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,vector,block_diagonal_matrix
signal.alarm(25)
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_uniform_modular_rr_v2';OLD=ART/'det1092_modular_bisection_direct_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
FROB=ART/'curve302_parent_geometric_picard19_v1.json'
BRAUER=ART/'det1092_surface_brauer_triviality_v1/arithmetic-lattice.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bind(d):
    for name,h in d['inputs'].items():assert sha(ROOT/name)==h,name
protocol=read(OUT/'protocol.json');summary=read(OUT/'summary.json');bind(protocol);bind(summary)
bind(read(OLD/'summary.json'))
assert protocol['chart_order']==['infinity',1,2,3]
assert protocol['limits']['old_target_trials']==30 and protocol['limits']['new_orbits']==0
assert protocol['limits']['good_primes']==[149,151]
parent=read(PARENT);geometry=read(OUT/'geometry.json')['primes']
frob=read(FROB);assert frob['status']=='PASS'
for name,h in frob['sources'].items():assert sha(ROOT/name)==h
G=matrix(ZZ,parent['generic_height_gram']);assert G.det()==1092
NS=block_diagonal_matrix(matrix(ZZ,[[-2,1],[1,0]]),-G)
lat=read(BRAUER);v=vector(ZZ,lat['unique_order_two_dual_word'])
assert NS.det()==1092 and all(c%2==0 for c in NS*v)
assert NS.change_ring(GF(2)).right_kernel().dimension()==1 and v*NS*v==-180
assert QQ(v*NS*v/4)==-45  # The unique index-two candidate is not even.
assert 1092==4*3*7*13  # Only index2 could enlarge an equal-rank integral lattice.
masks=sorted(summary['new_masks']+summary['regression_masks']);assert len(masks)==15
assert summary['new_masks']==read(OLD/'summary.json')['new_masks']
RX=PolynomialRing(QQ,'z');z=RX.gen();checks=[];certificates=[];recovered=[]
for geom in geometry:
    p=geom['prime'];field=GF(p);R=PolynomialRing(field,'t');t=R.gen();K=R.fraction_field()
    def dec(record):
        numerator=R([field(QQ(c)) for c in record['numerator']])
        denominator=R([field(QQ(c)) for c in record['denominator']]);assert denominator
        return K(numerator)/denominator
    a1,a2,a3,a4,a6=[dec(r) for r in parent['a_invariants']]
    b2=a1*a1+4*a2;b4=a1*a3+2*a4;b6=a3*a3+4*a6
    A=R(b4/2-b2*b2/48);B=R(b6/4-b2*b4/24+b2**3/864)
    D=-16*(4*A**3+27*B**2)
    assert A==R(geom['A']) and B==R(geom['B'])
    assert D.degree()==geom['discriminant_degree']==24
    g=D.gcd(D.derivative());assert g==R(geom['discriminant_derivative_gcd'])
    if p==157:
        assert g.degree()==3 and not geom['rootless_24I1_gate'];continue
    assert p in [149,151] and g.degree()==0 and A.gcd(D).degree()==0
    basis=[]
    for pair in parent['basis_weierstrass_coordinates']:
        x,y=map(dec,pair);X=x+b2/12;Y=y+(a1*x+a3)/2
        assert Y*Y==X**3+A*X+B;basis.append((X,Y))
    def add(P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;u,v=Q
        if x==u:
            if y==-v:return None
            assert y==v and y
            slope=(3*x*x+A)/(2*y)
        else:slope=(v-y)/(u-x)
        r=slope*slope-x-u
        return (r,-y+slope*(x-r))
    def mul(n,P):
        if n<0:return mul(-n,(P[0],-P[1]))
        Q=None
        while n:
            if n&1:Q=add(Q,P)
            n//=2;P=add(P,P)
        return Q
    def height(P):
        if P is None:return 0
        X,Y=P;dx=int(X.denominator().degree());dy=int(Y.denominator().degree())
        assert dx%2==0 and dy*2==dx*3
        infinity_x=max(0,int(X.numerator().degree())-dx-4)
        infinity_y=max(0,int(Y.numerator().degree())-dy-6)
        assert infinity_x%2==0 and infinity_y*2==infinity_x*3
        return 4+dx+infinity_x
    heights=[height(P) for P in basis];assert heights==geom['basis_heights']
    gram=matrix(QQ,17,17)
    for i in range(17):gram[i,i]=heights[i]
    pairs=[]
    for i in range(17):
        for j in range(i):
            h=height(add(basis[j],basis[i]));pairs.append([i,j,h])
            gram[i,j]=gram[j,i]=QQ(h-heights[i]-heights[j])/2
    assert pairs==geom['pair_sum_heights'] and gram==G==matrix(ZZ,geom['Gram'])
    # Old certified H2 polynomial, not a new point count. Tate/Artin--Tate
    # and the effective-divisor implication are supplied in the written proof.
    old=next(r for r in frob['reductions'] if r['prime']==p)
    characteristic=RX(old['full_H2_polynomial'])
    tail,rem=characteristic.quo_rem((z-p)**19*(z+p));assert not rem
    ap=-tail[1];assert tail==z*z-ap*z+p*p
    assert (ap/QQ(p)).denominator()!=1
    artin_tate_value=4*p*p-ap*ap
    anti_squareclass=QQ(artin_tate_value)/G.det()
    assert anti_squareclass==({149:QQ(25),151:QQ(29)}[p])
    assert not (QQ(6)/anti_squareclass).is_square()
    checks.append({'prime':p,'pair_trace':int(ap),'arithmetic_Picard_rank':19,
        'geometric_Picard_rank':20,'arithmetic_Picard_index':1,
        'Artin_Tate_Fp_squared_value':int(artin_tate_value),
        'anti_height_squareclass_representative':str(anti_squareclass),
        'anti_norm_six_possible':False,'basis_Gram_replayed':True})
    for mask in masks:
        row=read(OUT/('orbit-%d-p%d.json'%(mask,p)))
        old=read(OLD/('orbit-%d.json'%mask));assert row['word']==old['word']
        w=vector(ZZ,row['word']);assert w*G*w==10
        P=None
        for n,Q in reversed(list(zip(w,basis))):P=add(P,mul(-int(n),Q))
        assert height(P)==10
        X,Y=P;attempts=[]
        for chart in protocol['chart_order']:
            if chart=='infinity':x,y,aa,bb,target=X,Y,A,B,field(0)
            else:
                subst=K(chart+1/t)
                x=t**4*X(subst);y=t**6*Y(subst)
                aa=R(t**8*A(subst));bb=R(t**12*B(subst));target=-field(1)/chart
                assert chart+1/target==0
                if X.denominator()(0) and x.denominator()(target):
                    assert x(target)==X(0)*target**4
            degree=int(x.denominator().degree());attempts.append({'chart':chart,'denominator_degree':degree})
            if degree==6:break
        assert degree==6 and attempts==row['charts'] and chart==row['selected_chart']
        assert int(target)==row['target_coordinate']
        h,nx,ny=[R(row[k]) for k in ['pole','nx','ny']]
        assert h.is_monic() and h.degree()==3 and h*h==x.denominator().monic()
        assert x==nx/h**2 and y==ny/h**3 and nx.degree()<=10 and ny.degree()<=15
        columns=[h**3*t**j for j in range(10)]+[nx*h*t**j for j in range(6)]+[ny*t**j for j in range(4)]
        mat=matrix(field,[[f[j] for f in columns] for j in range(19)])
        kernel=vector(field,row['kernel']);assert kernel and mat*kernel==0
        cols=row['pivot_columns'];assert len(cols)==len(set(cols))==19
        minor=mat.matrix_from_columns(cols).det();assert minor and int(minor)==row['pivot_minor_determinant']
        f0,f1,f2=[R(list(a)) for a in [kernel[:10],kernel[10:16],kernel[16:]]]
        assert f0*h**3+f1*nx*h+f2*ny==0
        assert -16*(4*aa(target)**3+27*bb(target)**2)
        S=PolynomialRing(field,'X');u=S.gen()
        if h(target):
            assert row['mode']=='finite-trace-quadratic' and f2(target)
            cubic=(f0(target)+f1(target)*u)**2-f2(target)**2*(u**3+aa(target)*u+bb(target))
            residual,remainder=cubic.quo_rem(u-nx(target)/h(target)**2);assert not remainder
            value=residual.discriminant()
        else:
            assert row['mode']=='zero-trace-vertical-line' and not f2(target) and f1(target)
            xc=-f0(target)/f1(target);value=xc**3+aa(target)*xc+bb(target)
            residual=u*u-value
        assert list(map(int,residual.list()))==row['residual_coefficients'] and int(value)==row['square_test_value']
        roots=[int(a) for a in field if residual(a)==0]
        status=('EXCLUDED_NONSQUARE' if not roots else
                'UNKNOWN_REPEATED_REDUCTION' if len(roots)==1 else 'LOCALLY_SPLIT_SIMPLE_UNKNOWN_OVER_Q')
        assert status==row['status']
        if status=='EXCLUDED_NONSQUARE':certificates.append([mask,p])
        if chart!='infinity':recovered.append({'mask':mask,'prime':p,'chart':chart,'status':status})
        prior=next(r for r in old['trials'] if r['prime']==p)
        if prior['status']=='EXCLUDED_NONSQUARE':assert status=='EXCLUDED_NONSQUARE'
        elif prior['status']=='UNKNOWN_SPLIT_OR_REPEATED_REDUCTION':assert status!='EXCLUDED_NONSQUARE'
assert len(checks)==2 and recovered==summary['pole_gate_recovered']
assert summary['new_excluded']==[m for m in summary['new_masks'] if any(c[0]==m for c in certificates)]
# Symbolic vertical-line intersection: no new arithmetic parameter evaluated.
V=PolynomialRing(QQ,['aa','bb','l0','l1','yy']);aa,bb,l0,l1,yy=V.gens();L=V.fraction_field()
xc=-L(l0)/l1
assert l1**3*(yy**2-(xc**3+aa*xc+bb))==l1**3*yy**2+l0**3+aa*l0*l1**2-bb*l1**3
report={'status':'PASS_INDEPENDENT_UNIFORM_MODULAR_RR_AND_SMOOTH_REDUCTION_INPUTS',
    'classification':'verified exact applications; uniform divisor, four-chart and Artin--Tate implications are written proofs',
    'good_prime_checks':checks,'geometry_regression':{'prime':157,'discriminant_derivative_gcd_degree':3},
    'old_orbit_trials':30,'pole_gate_recovered':recovered,'new_excluded':summary['new_excluded'],
    'new_unresolved':summary['new_unresolved'],'vertical_chart_symbolic_identity':True,
    'boundary':'No full atlas execution, new prime, target, point search or seed. The residual-rank and smooth-reduction theorem covers all minimum norm10 classes by geometry, not by an equation census.',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [OUT/'protocol.json',OUT/'summary.json',FROB,BRAUER,OLD/'summary.json',Path(__file__)]}}
dest=OUT/'independent-replay.json';encoded=json.dumps(report,indent=2,sort_keys=True)+'\n'
if dest.exists():assert read(dest)==report
else:
    with dest.open('x') as stream:stream.write(encoded)
print(report['status'],checks,'recovered',recovered,flush=True)

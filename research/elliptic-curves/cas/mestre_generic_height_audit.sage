#!/usr/bin/env sage-python
"""Exact section heights on six existing parents; no fibration enumeration."""
import argparse,sys,json
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector,prod
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
from mestre_parent_adapter import rational_function
from probe_mestre_fermigier_two_section_local_continuation import reconstructed_second_line,normalized_data
from research_runtime.store import checkpoint,digest
ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/mestre-generic-height-v1'

def construct(outer_u):
    data=cert.read(ART/'mestre_component_coherent_sections_v1.json');u=cert.F(outer_u)
    R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field();Sring=PolynomialRing(R,'x');x=Sring.gen()
    roots=[QQ(str(rational_function(r,u))) for r in data['roots']]
    F=prod((x-r-T)*(x-r+T) for r in roots);S=x**6
    for j in range(5,-1,-1):S+=R((F[6+j]-(S*S)[6+j])/2)*x**j
    rem=S*S-F
    if rem.degree()>4:raise ArithmeticError('square approximation failed')
    quartic=[]
    for j in range(5):
        q,res=rem[j].quo_rem(T*T)
        if res:raise ArithmeticError('quartic normalization failed')
        quartic.append(K(q))
    qp=[]
    for j in data['fixed_source_root_order']:
        for sign in (-1,1):
            xx=K(roots[j]+sign*T);qp.append((xx,K(S(xx))/T))
    v,c2,m2=reconstructed_second_line(u);c1,m1=normalized_data(u,v)[4:6]
    for i,(c,m) in enumerate(((c1,m1),(c2,m2))):
        yy=R([QQ(str(rational_function(a,u))) for a in data['extra_ordinate_coefficients'][i]])
        qp.append((K(QQ(str(c))+QQ(str(m))*T),K(yy)))
    e,d,c,b,a=quartic;I=12*a*e-3*b*d+c*c;J=72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c*c*c
    A=R(-27*I);B=R(-27*J);E=EllipticCurve(K,[A,B]);points=[]
    for xx,yy in qp:
        if yy*yy!=sum(q*xx**j for j,q in enumerate(quartic)):raise ArithmeticError('labelled rational section missed quartic')
        g0=b*b/16-a*c/6;g1=b*c/12-a*d/2;g2=c*c/12-b*d/8-a*e;g3=c*d/12-b*e/2;g4=d*d/16-c*e/6
        g=g0*xx**4+g1*xx**3+g2*xx*xx+g3*xx+g4
        gx=4*g0*xx**3+3*g1*xx*xx+2*g2*xx+g3;gy=g1*xx**3+2*g2*xx*xx+3*g3*xx+4*g4
        fx=4*a*xx**3+3*b*xx*xx+2*c*xx+d;fy=b*xx**3+2*c*xx*xx+3*d*xx+4*e
        h=(fx*gy-fy*gx)/8;points.append(E([36*g/(yy*yy),108*h/(yy**3)]))
    parent=next(r for r in cert.read(ART/'mestre_parent_portfolio_intake_v1.json')['rows'] if r['outer_u']==str(u))
    if A!=R(parent['A_coefficients']) or B!=R(parent['B_coefficients']):raise ArithmeticError('exact parent model binding differs')
    return E,A,B,points

def setup_height(A,B):
    R=A.parent();T=R.gen();delta=-16*(4*A**3+27*B**2)
    if A.degree()!=8 or B.degree()!=12 or delta.degree()!=20 or delta.gcd(A).degree()!=0:raise ArithmeticError('minimal semistable K3 model required')
    factors=delta.squarefree_decomposition()
    if sorted((int(f.degree()),int(m)) for f,m in factors)!=[(2,2),(16,1)]:raise ArithmeticError('fixed fibre configuration differs')
    double=next(f for f,m in factors if m==2);bases=double.roots(QQ,multiplicities=False)
    if len(bases)!=2:raise ArithmeticError('two rational I2 values required')
    Xring=PolynomialRing(QQ,'x');x=Xring.gen();nodes=[]
    for t in sorted(bases):
        f=x**3+A(t)*x+B(t);g=f.gcd(f.derivative()).monic()
        if g.degree()!=1:raise ArithmeticError('unique finite node required')
        nodes.append((t,-g[0]))
    f=x**3+A[8]*x+B[12];g=f.gcd(f.derivative()).monic()
    if g.degree()!=1:raise ArithmeticError('unique infinity node required')
    node=-g[0];c=(3*node).sqrt()
    if c not in QQ or not c:raise ArithmeticError('split infinity required')
    if A[7] or B[11]:raise ArithmeticError('infinity node has no linear drift in this chart')
    # The first blown-up node has no constant t^2 term. Its intersection
    # is the component2 point; the two tangent directions label1 and3.
    if A[6]*node+B[10]:raise ArithmeticError('I4 first blowup equation differs')
    return nodes,node,QQ(c)

def infinity_value(f,weight):
    if not f:return QQ(0)
    n,d=f.numerator(),f.denominator();order=d.degree()-n.degree()+weight
    if order<0:return None
    return QQ(0) if order>0 else n.leading_coefficient()/d.leading_coefficient()

def pole_intersection(P):
    if P.is_zero():raise ValueError('zero section self-intersection is separate')
    x=P[0];n,d=x.numerator(),x.denominator()
    twice=max(int(d.degree()),int(n.degree())-4)
    if twice<0 or twice%2:raise ArithmeticError('integral pole intersection required')
    return twice//2

def height(P,A,B,geometry):
    if P.is_zero():return QQ(0),{'zero':True,'components':[0,0,0]}
    nodes,node,tangent=geometry;T=A.parent().gen();x,y=P.xy();components=[]
    for base,nodex in nodes:
        if x.denominator()(base)==0:components.append(0)
        else:components.append(int(x(base)==nodex and y(base)==0))
    X,Y=infinity_value(x,4),infinity_value(y,6);component=0
    if X==node and Y==0:
        f=x-node*T**4
        vx=10**9 if not f else 4+f.denominator().degree()-f.numerator().degree()
        vy=10**9 if not y else 6+y.denominator().degree()-y.numerator().degree()
        if min(vx,vy)<1:raise ArithmeticError('node valuation failed')
        if min(vx,vy)>=2:component=2
        else:
            ratio=infinity_value(y/f,2)
            if ratio==tangent:component=1
            elif ratio==-tangent:component=3
            else:raise ArithmeticError('first blowup tangent differs')
    components.append(component);intersection=pole_intersection(P)
    correction=QQ(components[0]+components[1])/2+QQ(component*(4-component))/4
    h=4+2*intersection-correction
    if h<=0:raise ArithmeticError('nonzero supplied combination has nonpositive height')
    return h,{'zero':False,'zero_section_intersection':intersection,'components':components,'correction':str(correction)}

def worker(index):
    protocol=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in protocol['sources'].items()):raise ArithmeticError('frozen sources changed')
    u=protocol['outer_parameters'][index];out=D/('u'+u)/'height.json'
    if out.exists():raise FileExistsError('preserve height computation')
    data={'status':'RUNNING','stage':'construct','protocol_hash':digest(protocol),'outer_u':u};checkpoint(out,data)
    E,A,B,P=construct(u);geometry=setup_height(A,B)
    heights=[];profiles=[]
    for point in P:
        h,profile=height(point,A,B,geometry);heights.append(h);profiles.append(profile)
        h2,profile2=height(2*point,A,B,geometry)
        if h2!=4*h or profile2['components']!=[(2*profile['components'][j])%n for j,n in enumerate((2,2,4))]:raise ArithmeticError('doubling height/component check failed')
    G=matrix(QQ,14,14);pairs=[]
    data['stage']='height_pairs';checkpoint(out,data)
    for i in range(14):
        G[i,i]=heights[i]
        for j in range(i):
            hplus,pplus=height(P[i]+P[j],A,B,geometry);hminus,pminus=height(P[i]-P[j],A,B,geometry)
            if hplus+hminus!=2*(heights[i]+heights[j]):raise ArithmeticError('exact parallelogram identity failed')
            for n,k in ((2,0),(2,1),(4,2)):
                if pplus['components'][k]!=(profiles[i]['components'][k]+profiles[j]['components'][k])%n:raise ArithmeticError('component addition check failed')
            G[i,j]=G[j,i]=(hplus-heights[i]-heights[j])/2
            pairs.append({'indices':[j,i],'sum_height':str(hplus),'difference_height':str(hminus),'sum_profile':pplus,'difference_profile':pminus})
        print('u'+u,'HEIGHT ROW',i+1,'of14',flush=True)
    if G.rank()!=11:raise ArithmeticError('rank11 covariant Gram required')
    # Replace covariant images by the same divisor cloud as the pilot.
    C=matrix(QQ,14,14);C[0,0]=1
    for j in range(1,14):C[j,0]=-QQ(1)/2;C[j,j]=QQ(1)/2
    cloud_gram=C*G*C.transpose()
    seed_path=ROOT/'artifacts/local/elliptic-curves/mestre-parent-calibration-v1'/('u'+u+'-unit')/'seed.json';seed=cert.read(seed_path)
    indices=seed['independent_column_indices'];basis_gram=cloud_gram.matrix_from_rows_and_columns(indices,indices)
    if not basis_gram.is_positive_definite():raise ArithmeticError('independent divisor Gram must be positive definite')
    span=next(r for r in cert.read(ART/'mestre_parent_calibration_input_span_v1.json')['rows'] if r['id']=='u'+u+'-unit');relations=[]
    for r in span['relations']:
        w=vector(QQ,14);w[r['target_index']]=r['target_multiplier']
        for i,a in zip(indices,r['basis_coefficients']):w[i]-=a
        raw=2*w*C
        if any(a.denominator()!=1 for a in raw):raise ArithmeticError('integer covariant relation required')
        if sum((ZZ(a)*point for a,point in zip(raw,P)),E(0))!=E(0):raise ArithmeticError('specialized relation failed over Q(T)')
        if G*raw!=0:raise ArithmeticError('relation not in height kernel')
        relations.append(list(map(int,raw)))
    if matrix(QQ,relations).rank()!=3:raise ArithmeticError('three independent generic relations required')
    data.update(status='PASS',stage='complete',curve_A_coefficients=list(map(str,A.list())),curve_B_coefficients=list(map(str,B.list())),
      covariant_points=[[str(c) for c in point.xy()] for point in P],covariant_height_gram=[[str(c) for c in row] for row in G.rows()],
      covariant_profiles=profiles,pair_checks=pairs,generic_covariant_relations=relations,
      divisor_change_matrix=[[str(c) for c in row] for row in C.rows()],divisor_height_gram=[[str(c) for c in row] for row in cloud_gram.rows()],
      seed_indices=indices,seed_height_gram=[[str(c) for c in row] for row in basis_gram.rows()],
      seed_height_determinant=str(basis_gram.det()),known_rank18_divisor_lattice_absolute_determinant=str(16*basis_gram.det()),
      supplied_generic_section_subgroup_rank=11,
      scope='Exact Q(T) section equations, local intersection heights, doubling/parallelogram/component checks and three rational-function group relations. The supplied section subgroup has rank exactly11; the ambient generic MW rank and full NS rank remain unknown. The rank18 divisor lattice determinant includes U and all five reducible-fibre components, but no saturation or full NS identification is claimed. No new fibration, parent or point search.')
    checkpoint(out,data);print('PASS u'+u,'MW DET',basis_gram.det(),'NS SUBLATTICE DET',16*basis_gram.det(),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--index',type=int,required=True);args=parser.parse_args();worker(args.index)

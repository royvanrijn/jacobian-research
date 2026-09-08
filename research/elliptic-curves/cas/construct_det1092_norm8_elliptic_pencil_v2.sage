#!/usr/bin/env sage-python
"""One generic-only norm8 elliptic pencil, not a point search.

Select the minimum(l1,linfinity,orbit) norm8 row of the existing census.
One exact parity check; one quartet; two known generic section restrictions;
one predeclared auxiliary fibre z=0, with small finite-group certificates.
"""
import csv,hashlib,json,runpy,sys
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector,gcd
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';CAS=ROOT/'elliptic-curves/cas'
sys.path.insert(0,str(CAS))
from visibility_lattice_v2 import ExactParity
OUT=ART/'det1092_norm8_elliptic_pencil_v2.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def build():
    pp=ART/'curve302_recovered_mw17_parent_v1.json';op=ART/'curve302_parent_degree2_multisection_orbits_v1.tsv';loader=CAS/'load_curve302_recovered_parent.sage'
    with op.open() as stream:rows=[r for r in csv.DictReader(stream,delimiter='\t') if r['minimum_norm']=='8']
    def key(row):
        w=list(map(int,row['parent_MW17_w'].split()));return sum(map(abs,w)),max(map(abs,w)),int(row['orbit_mask'])
    chosen=min(rows,key=key);w=vector(ZZ,list(map(int,chosen['parent_MW17_w'].split())))
    E,base,_=runpy.run_path(str(loader))['load_curve302_recovered_parent'](pp)
    parent=json.loads(pp.read_text());G=matrix(ZZ,parent['generic_height_gram']);assert w*G*w==8
    U=G.LLL_gram();wr=vector(ZZ,U.inverse()*w);cvp=ExactParity((U.transpose()*G*U).rows()).solve(wr,wr,node_limit=200000);assert cvp['norm']==8
    R=E.base_ring().ring();t=R.gen();C=sum((n*P for n,P in zip(w,base)),E(0))
    a=-E.c4()/48;b=-E.c6()/864;cx=C[0]+E.b2()/12;cy=C[1]+(E.a1()*C[0]+E.a3())/2
    h=R(C[0].denominator().sqrt());nx=R(cx*h*h);ny=R(cy*h**3)
    assert h.degree()==2 and nx.degree()==8 and ny.degree()==12
    shift=(ny*nx.inverse_mod(h*h))%(h*h)
    Z=PolynomialRing(QQ,'z');z=Z.gen();S=PolynomialRing(Z,'t')
    N=S(h*h)*z-S(shift)
    raw=N**4-6*S(nx)*N*N-8*S(ny)*N-3*S(nx)**2-4*S(R(a))*S(h)**4
    f,rem=raw.quo_rem(S(h)**6);assert not rem;f=S(f);assert f.degree()==4
    e,d,c,bb,aa=f.list();I=12*aa*e-3*bb*d+c*c;J=72*aa*c*e+9*bb*c*d-27*aa*d*d-27*bb*bb*e-2*c**3
    ja=-27*I;jb=-27*J;jdelta=-16*(4*ja**3+27*jb**2)
    assert jdelta and ja.degree()<=8 and jb.degree()<=12
    F=Z.fraction_field();jf=F(6912*ja**3)/(4*ja**3+27*jb**2);assert jf.derivative()
    # The chosen short word is checked, not silently assumed from selection.
    assert list(w)==[0]*14+[1,-1,0]
    sections=[]
    for i,sgn in [(14,1),(15,-1)]:
        word=vector(ZZ,[sgn if j==i else 0 for j in range(17)]);P=sgn*base[i]
        assert word*G*word-word*G*w==1
        xx=P[0]+E.b2()/12;yy=P[1]+(E.a1()*P[0]+E.a3())/2
        old_slope=(yy+cy)/(xx-cx);param=(h*old_slope+shift)/(h*h)
        assert max(param.numerator().degree(),param.denominator().degree())==1
        n0,n1=[param.numerator()[j] for j in range(2)];d0,d1=[param.denominator()[j] for j in range(2)]
        T=F(n0-z*d0)/(z*d1-n1)
        def at(v):return F(v.numerator()(T))/v.denominator()(T)
        W=(2*at(xx)+at(cx)-at(old_slope)**2)/h(T)
        assert W*W==sum(F(v)*T**j for j,v in enumerate(f.list()))
        sections.append({'word':list(map(int,word)),'t_of_z':T,'W_of_z':W})
    assert sections[0]['t_of_z']==sections[1]['t_of_z'] and sections[0]['W_of_z']==-sections[1]['W_of_z']
    # Intersection on the original K3, hence on the alternate elliptic K3.
    intersection=QQ(G[14,14]+G[15,15])/2-2+G[14,15];assert intersection==4
    # Fixed fibre z=0, selected before any solubility/rank outcome.
    f0=R([v(0) for v in f.list()]);assert f0.degree()==4 and f0.discriminant()
    t0=sections[0]['t_of_z'](0);W0=sections[0]['W_of_z'](0);assert W0 and W0*W0==f0(t0)
    jac=EllipticCurve(QQ,[ja(0),jb(0)])
    ee,dd,cc,bbb,aaa=f0.list()
    gs=[bbb*bbb/16-aaa*cc/6,bbb*cc/12-aaa*dd/2,cc*cc/12-bbb*dd/8-aaa*ee,cc*dd/12-bbb*ee/2,dd*dd/16-cc*ee/6]
    g0,g1,g2,g3,g4=gs;gv=g0*t0**4+g1*t0**3+g2*t0*t0+g3*t0+g4
    gx=4*g0*t0**3+3*g1*t0*t0+2*g2*t0+g3;gy=g1*t0**3+2*g2*t0*t0+3*g3*t0+4*g4
    ux=4*aaa*t0**3+3*bbb*t0*t0+2*cc*t0+dd;uy=bbb*t0**3+2*cc*t0*t0+3*dd*t0+4*ee
    hv=(ux*gy-uy*gx)/8;P=jac([36*gv/W0**2,108*hv/W0**3])
    finite=[];bound=ZZ(0)
    for p in [17,47,53,61,67,71,79,83,89,101,107,113,127,137,149,179,191,197]:
        if any(v.denominator()%p==0 for v in [*jac.a_invariants(),*P[:2]]):continue
        if jac.discriminant().valuation(p)!=0:continue
        field=GF(p);ep=EllipticCurve(field,jac.a_invariants())
        if not ep.discriminant():continue
        order=ZZ(ep.cardinality());bound=gcd(bound,order);finite.append({'prime':p,'order':int(order)})
        if len(finite)>=2 and bound*ep(list(P[:2]))!=ep(0):
            witness={'prime':p,'torsion_order_bound':int(bound),'multiple_nonzero':True};break
    else:raise ArithmeticError('fixed finite pool did not certify the point nontorsion')
    def rec(v):return {'numerator':list(map(str,v.numerator().list())),'denominator':list(map(str,v.denominator().list()))}
    return {'classification':'verified application and new deduction','status':'PASS_GENERIC_NORM8_ELLIPTIC_PENCIL_WITH_POSITIVE_RANK_BASE',
        'selection':{'orbit':int(chosen['orbit_mask']),'rule':'minimum(l1,linfinity,orbit) of stored norm8 rows','word':list(map(int,w))},
        'exact_parity_minimum':8,'exact_nodes':cvp['nodes'],'pole_h':list(map(str,h.list())),
        'nx':list(map(str,nx.list())),'ny':list(map(str,ny.list())),'shift':list(map(str,shift.list())),
        'quartic_t_coefficients_in_z':[list(map(str,v.list())) for v in f.list()],
        'Jacobian_A':list(map(str,ja.list())),'Jacobian_B':list(map(str,jb.list())),
        'sections':[{'word':row['word'],'t_of_z':rec(row['t_of_z']),'W_of_z':rec(row['W_of_z'])} for row in sections],
        'section_intersection':4,'generic_section_difference_nontorsion':'Distinct torsion sections cannot meet the zero section over a characteristic-zero base; the two displayed sections meet four times.',
        'forward_map':'m=(Y+cy)/(X-cx); z=(h*m+shift)/h^2; W=(2*X+cx-m^2)/h',
        'inverse_map':'m=h*z-shift/h; X=(h*W-cx+m^2)/2; Y=m*(X-cx)-cy; x=X-b2/12; y=Y-(a1*x+a3)/2',
        'fixed_base':{'z':0,'quartic_coefficients':list(map(str,f0.list())),'rational_point':[str(t0),str(W0)],
            'Jacobian_model':list(map(str,jac.a_invariants())),'nontorsion_Jacobian_point':list(map(str,P[:2])),
            'finite_group_orders':finite,'nontorsion_witness':witness},
        'limits':{'new_pencils':1,'parity_cosets':1,'node_limit':200000,'fixed_auxiliary_fibres':1,'point_searches':0,'exceptional_point_inputs':0},
        'boundary':'A source-only alternate elliptic fibration and one positive-rank elliptic base. This does not yet identify the first302 point as a generically generated section on the alternate fibration.',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [pp,op,loader,CAS/'visibility_lattice_v2.py',Path(__file__)]}}
if __name__=='__main__':
    if OUT.exists():raise FileExistsError(OUT)
    d=build();OUT.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(d['status'],'orbit',d['selection']['orbit'],'fixed base orders',d['fixed_base']['finite_group_orders'],flush=True)


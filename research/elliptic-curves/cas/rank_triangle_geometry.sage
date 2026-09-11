#!/usr/bin/env sage-python
"""Exact bounded multisection atlas through the fixed14/8 representatives.

Vertical x and constant-slope chords through all signed generic sections;
add the existing class1 genus-one pencil on302. No claim of minimal genus
outside this atlas. Covers are compared over the fixed original base.
"""
import gzip,json,hashlib
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';OUT=ART/'rank_triangle_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def save(p,d):
    p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():assert read(p)==d
    else:p.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
R=PolynomialRing(QQ,'t');F=R.fraction_field();t=R.gen()
def dec(d):
    if 'numerator' in d:return F(R(d['numerator'])/R(d['denominator']))
    return F(R(d['numerator_coefficients_low_to_high'])/R(d['denominator_coefficients_low_to_high']))
def at(f,t0):return F(f).numerator()(t0)/F(f).denominator()(t0)
def enc(f):
    f=F(f);return {'numerator':list(map(str,f.numerator().list())),'denominator':list(map(str,f.denominator().list()))}
def squarefree_cover(raw):
    raw=F(raw);n,d=raw.numerator(),raw.denominator();assert raw
    c=n.leading_coefficient()/d.leading_coefficient();q=R(1);s=F(1)
    for poly,sign in [(n/n.leading_coefficient(),1),(d/d.leading_coefficient(),-1)]:
        for factor,e in poly.squarefree_decomposition():
            factor=R(factor);q*=factor**(e%2)
            s*=F(factor)**(e//2 if sign==1 else -((e+1)//2))
    assert raw==c*q*s*s and q.gcd(q.derivative())==1
    return c,q,s
def parent(cid):
    if cid=='302':
        p=read(ART/'curve302_recovered_mw17_parent_v1.json');a1,a2,a3,a4,a6=map(dec,p['a_invariants'])
        b2=a1*a1+4*a2;b4=a1*a3+2*a4;b6=a3*a3+4*a6
        A=-(b2*b2-24*b4)/48;B=-(-b2**3+36*b2*b4-216*b6)/864
        points=[]
        for P in p['basis_weierstrass_coordinates']:
            x,y=map(dec,P);points.append((x+b2/12,y+(a1*x+a3)/2))
        return A,B,points,QQ(0)
    p=read(ROOT/'artifacts/local/elliptic-curves/broad-rank-v1/runtime/research/broad-inputs/parents/11952.json')
    return F(R(p['A_coefficients_low_to_high'])),F(R(p['B_coefficients_low_to_high'])),[(dec(P['X']),dec(P['Y'])) for P in p['sections']],QQ(921)/653
def class1(A,B,P,t0):
    tr,rr,eq=[read(ART/('x1092_class1_realization_'+s+'_v1.json')) for s in ['trace','rr','equation']]
    assert A==R(tr['A']) and B==R(tr['B'])
    h,nx,ny=[R(tr[k]) for k in ['h','Nx','Ny']];a0,b0,a1,b1=[R(rr[k]) for k in ['a0','b0','a1','b1']]
    x,y=P;u=(a1(t0)*(x*h(t0)**2-nx(t0))+b1(t0)*(y*h(t0)**3+ny(t0)))/(a0(t0)*(x*h(t0)**2-nx(t0))+b0(t0)*(y*h(t0)**3+ny(t0)))
    ev=lambda d:at(dec(d),u)
    q=R([ev(c) for c in eq['quartic_coefficients']]);sf=R([ev(c) for c in eq['radical_square_factor']]);ds=R([ev(c) for c in eq['radical_denominator_sqrt']])
    m=F((a1-u*a0)/((u*b0-b1)*h));xp=F(nx/h**2);yp=F(ny/h**3)
    x0=(m*m-xp)/2;x1=F(sf/(2*ds));y0=m*(x0-xp)-yp;y1=m*x1
    v=(x-at(x0,t0))/at(x1,t0)
    return q,(x0,x1,y0,y1),v,{'family':'class1_fibre','parameter':str(u)}
def main():
    paths=[Path(__file__),ART/'rank_accessibility_subsets_v1/inputs.json.gz',ART/'curve302_recovered_mw17_parent_v1.json',
           ROOT/'artifacts/local/elliptic-curves/broad-rank-v1/runtime/research/broad-inputs/parents/11952.json']
    paths += [ART/('x1092_class1_realization_'+s+'_v1.json') for s in ['trace','rr','equation']]
    save(OUT/'geometry-protocol-v2.json',{'sha256':{str(p):sha(p) for p in paths},'atlas':'vertical x; constant slope through +/- each of17 generic sections; plus302 class1 pencil',
          'selection':'Fixed representatives from the prior accessibility input, not the different historical seed-universality basis.',
          'global_minimum_genus':'UNKNOWN','base_reparametrizations_quotiented_out':False})
    allrows=[]
    for data in read(paths[1])['curves']:
        cid=data['id'];A,B,sections,t0=parent(cid);Es=EllipticCurve(QQ,list(map(QQ,data['model'])));E0=EllipticCurve(QQ,[at(A,t0),at(B,t0)])
        images=[Es(list(map(QQ,P))) for P in data['basis_points']]
        def signkey(P):return (P[0],min(P[1],-P[1]))
        desired={signkey(E0(at(x,t0),at(y,t0))) for x,y in sections}
        iso=next(i for i in Es.isomorphisms(E0) if {signkey(i(P)) for P in images[:17]}==desired)
        for x,y in sections:assert y*y==x**3+A*x+B
        for j,original in enumerate(images[17:]):
            P=iso(original);px,py=P.xy();candidates=[]
            candidates.append((F(px**3+A*px+B),(F(px),F(0),F(0),F(1)),py,{'family':'vertical_x'}))
            for i,(tx,ty) in enumerate(sections):
                for sign in [-1,1]:
                    sy=sign*ty;m=(py+at(sy,t0))/(px-at(tx,t0));D=m**4-6*tx*m*m-8*sy*m-3*tx*tx-4*A
                    x0=(m*m-tx)/2;y0=m*(x0-tx)-sy
                    candidates.append((D,(x0,F(QQ(1)/2),y0,F(m/2)),2*px-(m*m-at(tx,t0)),{'family':'constant_slope','generic_index':i+1,'sign':sign,'slope':str(m)}))
            if cid=='302':candidates.append(class1(A,B,(px,py),t0))
            covers=[]
            for raw,(x0,x1,y0,y1),v0,tag in candidates:
                raw=F(raw);assert y0*y0+y1*y1*raw==x0**3+3*x0*x1*x1*raw+A*x0+B
                assert 2*y0*y1==3*x0*x0*x1+x1**3*raw+A*x1
                assert at(raw,t0)==v0*v0 and at(x0,t0)+at(x1,t0)*v0==px and at(y0,t0)+at(y1,t0)*v0==py
                c,q,s=squarefree_cover(raw);assert q.degree()>0 and at(s,t0)!=0
                v=v0/at(s,t0);assert v*v==c*q(t0)
                covers.append({'tag':tag,'degree_over_base':2,'normalization_genus':int((q.degree()-1)//2),
                               'branch_monic_polynomial':list(map(str,q.list())),'branch_at_infinity':bool(q.degree()%2),
                               'constant_twist_representative':str(c),'square_multiplier':enc(s),
                               'lift':[str(t0),str(v)],'raw_radical':enc(raw),
                               'maps':{k:enc(f) for k,f in zip(['x0','x1','y0','y1'],[x0,x1*s,y0,y1*s])},
                               'rational_bisection_orbit_class':'NOT_APPLICABLE_POSITIVE_GENUS' if q.degree()>2 else 'UNKNOWN_ORBIT_LABEL'})
            row={'fibre':cid,'target':j+1,'source_point':list(map(str,original.xy())),
                 'parent_point':list(map(str,P.xy())),'fibre_isomorphism':list(map(str,iso.tuple())),
                 'covers':covers,'best_genus_in_atlas':min(c['normalization_genus'] for c in covers),
                 'global_minimum_genus':'UNKNOWN','degree_upper_bound':2}
            save(OUT/'geometry'/f'{cid}-E{j+1}.json',row);allrows.append(row)
            print('ANCESTRY',cid,j+1,'best genus',row['best_genus_in_atlas'],flush=True)
    groups=[]
    for row in allrows:
        for k,c in enumerate(row['covers']):
            match=next((g for g in groups if g['fibre']==row['fibre'] and g['branch']==c['branch_monic_polynomial'] and (QQ(g['constant'])/QQ(c['constant_twist_representative'])).is_square()),None)
            if match is None:
                match={'fibre':row['fibre'],'branch':c['branch_monic_polynomial'],'constant':c['constant_twist_representative'],'members':[]};groups.append(match)
            match['members'].append([row['target'],k])
    for g in groups:g['targets']=sorted({r[0] for r in g['members']})
    save(OUT/'geometry.json',{'status':'PASS_FIXED_ATLAS','rows':allrows,'cover_equivalence_classes':groups,
         'equivalence':'Same quadratic extension over the fixed base iff identical monic branch polynomial and constant ratio a rational square. No full coefficient factorization required.',
         'boundary':'No arbitrary-base PGL2 canonicalization or genus-minimality proof. One pencil is not one shared cover. Fixed representatives, hence not invariant under arbitrary exceptional rebasing.'})
if __name__=='__main__':main()

#!/usr/bin/env sage-python
"""Independent full-Picard specificity replay on all nine fixed controls.

Exact divisor identities and finite multiplication determinants, no
constructor import or resultant. No exceptional coordinate, point search,
number field or Selmer group.25-second cap.
"""
import hashlib,json,runpy,signal
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
FRAME=ART/'det1092_genus1_picard_image_v1';SOURCE=ART/'det1092_norm8_seed_cover_v2'
DIR=ART/'det1092_genus1_picard_controls_v1'
HELPER=ROOT/'elliptic-curves/cas/verify_det1092_genus1_picard_image.sage'
finite_norm=runpy.run_path(str(HELPER))['finite_norm']
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def provenance(d):
    for path,digest in d['inputs'].items():assert sha(ROOT/path)==digest,(path,'HASH_MISMATCH')
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def residue(v,p):
    v=QQ(v);assert v.denominator()%p
    return int(v.numerator()*v.denominator().inverse_mod(p)%p)
def verify():
    paths=[DIR/'protocol.json',DIR/'selection.json',DIR/'panel.json',
      ART/'curve302_recovered_mw17_parent_v1.json',SOURCE/'generic.json',
      FRAME/'frame.json',FRAME/'generic-restrictions.json',HELPER,Path(__file__)]
    protocol,selection,panel,parent,pencil,frame,functions=map(read,paths[:7])
    for d in [protocol,selection,panel,frame,functions]:provenance(d)
    assert protocol['limits']['old_addresses']==9 and protocol['limits']['exceptional_point_inputs']==0
    G=matrix(QQ,parent['generic_height_gram']);NS=matrix(QQ,19)
    NS[0,0]=-2;NS[0,1]=NS[1,0]=1
    for i in range(17):
        for j in range(17):NS[i+2,j+2]=-G[i,j]
    D=vector(QQ,frame['D']);zero=vector(QQ,frame['new_zero_section']);kernel=[D,zero]
    assert D*NS*D==0 and D*NS*zero==1
    for row in frame['vertical_old_sections']:
        x=vector(QQ,row['word']);assert all(v in ZZ for v in x)
        S=vector(QQ,[1,x*G*x/2,*x]);assert D*NS*S==0 and S*NS*S==-2;kernel.append(S)
    assert matrix(QQ,kernel).rank()==7 and 19-7==12
    R=PolynomialRing(QQ,'t');K=R.fraction_field();t=R.gen();Z=PolynomialRing(QQ,'z')
    def dec(d):return K(R(d['numerator']))/R(d['denominator'])
    ai=list(map(dec,parent['a_invariants']));E=EllipticCurve(K,ai)
    basis=[E(list(map(dec,P))) for P in parent['basis_weierstrass_coordinates']]
    C=basis[14]-basis[15];cx=C[0]+E.b2()/12;cy=C[1]+(E.a1()*C[0]+E.a3())/2
    h=R(pencil['pole_h']);shift=R(pencil['shift'])
    assert h*h==C[0].denominator() and h.degree()==2
    section_functions=[]
    for i,P in enumerate(basis):
        xx=P[0]+E.b2()/12;yy=P[1]+(E.a1()*P[0]+E.a3())/2
        slope=(yy+cy)/(xx-cx);z=(h*slope+shift)/(h*h);ww=(2*xx+cx-slope*slope)/h
        assert z==dec(functions['sections'][i]['z']) and ww==dec(functions['sections'][i]['W'])
        section_functions.append((z,ww))
    selector=section_functions[0][0]
    assert max(selector.numerator().degree(),selector.denominator().degree())==selection['selector_degree']==4
    assert len(selection['cases'])==len(protocol['cases'])==len(panel['cases'])==9
    coefficients=list(map(Z,pencil['quartic_t_coefficients_in_z']));checks=[];casepaths=[]
    for index,(sel,case) in enumerate(zip(selection['cases'],protocol['cases'])):
        dp=DIR/f'case-{index:02d}-generic.json';cp=DIR/f'case-{index:02d}-codes.json';mp=DIR/f'case-{index:02d}-marked.json'
        data,codes,marked=map(read,[dp,cp,mp]);provenance(marked);casepaths +=[dp,cp,mp]
        assert sel['label']==case['label']==data['label']==marked['label']
        assert sel['tau']==case['parameter']==data['tau']==marked['tau']
        tau=QQ(sel['tau']);zvalue=selector(tau);assert str(zvalue)==sel['z']==data['z']==marked['z']
        q=R(data['quartic']);assert q==R([v(zvalue) for v in coefficients])
        assert q.degree()==4 and q.gcd(q.derivative())==1 and q.gcd(R(E.discriminant()))==1
        M=h*h*zvalue-shift
        nx=cx*h*h;ny=cy*h**3
        assert M**4-6*nx*M*M-8*ny*M-3*nx*nx-4*(-E.c4()/48)*h**4==h**6*q
        t0,s=map(QQ,data['basepoint']);assert s and q(t0)==s*s
        assert section_functions[14][0](t0)==zvalue and section_functions[14][1](t0)==s
        q0,q1,q2,q3,q4=q(t+t0).list();assert q0==s*s
        abc=list(map(QQ,data['cubic']));assert abc==[q2,q1*q3-4*s*s*q4,s*s*q3*q3+q1*q1*q4-4*s*s*q2*q4]
        J=EllipticCurve(QQ,[0,abc[0],0,abc[1],abc[2]])
        HF=list(map(QQ,data['F_point']));xb=q1*q1/(4*s*s)-q2;yb=-(q1*xb+2*s*s*q3)/(2*s)
        assert HF==[xb,yb] and J(HF)
        pairs=[];assert len(data['divisors'])==17
        for i,row in enumerate(data['divisors']):
            assert row['index']==i and row['degree']==functions['sections'][i]['intersection_degree']
            z,wi=section_functions[i]
            if i==14:
                assert row['zero'] and row['degree']==1;pairs.append(None);continue
            if row['degree']==0:
                assert row['zero'] and z!=zvalue;pairs.append(None);continue
            assert not row['zero']
            g=R(row['g']);ww=R(row['W']);xx=R(row['X']);yy=R(row['Y'])
            assert g==(z.numerator()-zvalue*z.denominator()).monic()
            assert g.degree()==row['degree'] and g.gcd(g.derivative())==1 and g(t0)
            assert wi.denominator().gcd(g)==1 and (ww*wi.denominator()-wi.numerator())%g==0
            assert (ww*ww-q)%g==0
            du=t-t0
            assert (xx*du*du-2*s*(ww+s)-q1*du)%g==0
            assert (2*s*yy-(xx*xx-4*s*s*q4)*du+q1*xx+2*s*s*q3)%g==0
            assert (yy*yy-xx**3-abc[0]*xx*xx-abc[1]*xx-abc[2])%g==0
            pairs.append((g,xx))
        assert [r['p'] for r in codes['trials']]==protocol['primes']
        rows=[];owners=[];no2=None;local=[]
        for trial in codes['trials']:
            p=trial['p'];status=trial['status']
            if any(v.denominator()%p==0 for v in abc):
                assert status=='SKIP_JACOBIAN_DENOMINATOR';continue
            a,b,c=[residue(v,p) for v in abc]
            if (a*a*b*b-4*b**3-4*a**3*c-27*c*c+18*a*b*c)%p==0:
                assert status=='SKIP_BAD_JACOBIAN_REDUCTION';continue
            roots=[x for x in range(p) if (x**3+a*x*x+b*x+c)%p==0]
            assert roots==trial['roots'] and trial['cubic']==[c,b,a,1]
            if not roots:
                assert status=='COMPLETE_ZERO_QUOTIENT' and trial['rows']==[]
                no2=p if no2 is None else no2;local.append({'p':p,'generic_rank':0,'full_dimension':0});continue
            if any(v.denominator()%p==0 for pair in pairs if pair for f in pair for v in f):
                assert status=='SKIP_FULL_DIVISOR_DENOMINATOR';continue
            S=PolynomialRing(GF(p),'t');norms=[];raw=[];failed=False
            for root in roots:
                v=1 if HF[0].valuation(p)<0 else (residue(HF[0],p)-root)%p
                if not v:v=(3*root*root+2*a*root+b)%p
                vals=[v]+[1 if pair is None else finite_norm(S(pair[1])-root,S(pair[0])) for pair in pairs]
                norms.append(vals)
                if any(v==0 for v in vals):failed=True;break
                raw.append([int(pow(v,(p-1)//2,p)==p-1) for v in vals])
            if failed:
                assert status=='SKIP_NONUNIT_DIVISOR_CLASS' and trial['partial_norms']==norms;continue
            assert status=='COMPLETE_FULL_PICARD_REDUCTION' and trial['rows']==raw and trial['norms']==norms
            rows+=raw;owners +=[p]*len(raw)
            local.append({'p':p,'generic_rank':int(matrix(GF(2),raw).rank()),'full_dimension':1 if len(roots)==1 else 2})
        A=matrix(GF(2),rows)
        assert rows==codes['rows'] and owners==codes['row_primes'] and A.rank()==codes['rank']==12
        assert no2==codes['no_two_torsion_prime'] and no2 is not None
        wm=section_functions[0][1](tau);assert wm*wm==q(tau)
        assert marked['marked_quartic_point']==list(map(str,[tau,wm]))
        du=tau-t0;assert du
        xm=(2*s*(wm+s)+q1*du)/(du*du);ym=((xm*xm-4*s*s*q4)*du-q1*xm-2*s*s*q3)/(2*s)
        assert marked['marked_Jacobian_point']==list(map(str,[xm,ym])) and J([xm,ym])
        m=h(tau)*zvalue-shift(tau)/h(tau);oldX=(h(tau)*wm-cx(tau)+m*m)/2
        oldx=oldX-E.b2()(tau)/12;oldY=m*(oldX-cx(tau))-cy(tau)
        oldy=oldY-(E.a1()(tau)*oldx+E.a3()(tau))/2
        assert oldx==basis[0][0](tau) and oldy==basis[0][1](tau)
        assert marked['old_elliptic_point']==list(map(str,[oldx,oldy]))
        column=[]
        for trial in codes['trials']:
            if not trial['status'].startswith('COMPLETE_'):continue
            p=trial['p'];c,b,a,one=trial['cubic'];assert one==1
            for root in trial['roots']:
                v=1 if xm.valuation(p)<0 else (residue(xm,p)-root)%p
                if not v:v=(3*root*root+2*a*root+b)%p
                assert v;column.append(int(pow(v,(p-1)//2,p)==p-1))
        assert column==marked['marked_column']
        aug=A.augment(matrix(GF(2),len(column),1,column));sep=vector(GF(2),marked['separating_row'])
        assert aug.rank()==13 and sep*A==0 and sep*vector(GF(2),column)==1
        assert marked['status']=='GENERIC_OLD_POINT_BUT_NON_GENERIC_CARRIER_CLASS'
        checks.append({'index':index,'label':case['label'],'generic_Picard_image_rank':12,
          'rank_with_marked_class':13,'no_rational_2_torsion_prime':no2,
          'old_elliptic_class':'GENERIC_SECTION0; ZERO_RELATIVE_CLASS',
          'carrier_class':'RATIONAL_OUTSIDE_FULL_PICARD_IMAGE; SELMER_WITH_ZERO_SHA_IMAGE',
          'local_comparisons':local})
    result={'status':'PASS_INDEPENDENT_NINE_GENUS1_PICARD_SPECIFICITY_COUNTEREXAMPLES',
      'classification':'new exact specificity obstruction and verified application',
      'cases':checks,'selector_degree':4,
      'generic_explanation':'The source-only degree4 base change z=z_section0(tau) gives a rational marked section of the carrier Jacobian. Its independence from the full inherited image at a certified specialization implies the same independence overQ(tau), while its old-elliptic image is identically the generic section0.',
      'boundary':'The first seed also gives rank12->13 on its calibrated carrier, but this property holds for generic old points on all eight controls and the source-only302 member. It is not an original-fibre seed discriminator. No full Selmer dimension or nonzero Sha class is inferred.',
      'limits':protocol['limits'],
      'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+casepaths}}
    retain(DIR/'replay.json',result)
    print(result['status'],'all nine: carrier12->13, original point generic',flush=True)
if __name__=='__main__':
    signal.alarm(25);verify()

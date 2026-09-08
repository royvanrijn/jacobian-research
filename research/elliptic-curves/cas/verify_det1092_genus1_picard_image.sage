#!/usr/bin/env sage-python
"""Independent full-NS restriction and marked-class finite-norm replay.

No CVP, constructor import, polynomial resultant, class group or point search.
The finite norm is an explicit multiplication determinant.25-second cap.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
DIR=ART/'det1092_genus1_picard_image_v1';SOURCE=ART/'det1092_norm8_seed_cover_v2'
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
def finite_norm(poly,g):
    S=g.parent();t=S.gen();degree=g.degree();assert g.is_monic()
    columns=[list((poly*t**j%g).list()) for j in range(degree)]
    columns=[c+[S.base_ring()(0)]*(degree-len(c)) for c in columns]
    return int(matrix(S.base_ring(),columns).det())
def verify():
    paths=[DIR/f for f in ['protocol.json','frame.json','generic-restrictions.json','codes-frozen.json','first-class.json']]
    protocol,frame,functions,frozen,marked=map(read,paths)
    for d in [protocol,frame,functions,frozen,marked]:provenance(d)
    assert protocol['limits']['NS_generators']==19 and protocol['limits']['fixed_members']==2
    parent=read(ART/'curve302_recovered_mw17_parent_v1.json');generic=read(SOURCE/'generic.json')
    G=matrix(QQ,parent['generic_height_gram']);assert G.det()==1092 and G.nrows()==17
    NS=matrix(QQ,19);NS[0,0]=-2;NS[0,1]=NS[1,0]=1
    for i in range(17):
        for j in range(17):NS[i+2,j+2]=-G[i,j]
    D=vector(QQ,frame['D']);zero=vector(QQ,frame['new_zero_section'])
    w=vector(QQ,[0]*14+[1,-1,0]);assert D==vector(QQ,[2,4,*w])
    assert D*NS*D==0 and D*NS*zero==1
    assert zero==vector(QQ,[1,G[14,14]/2,*[int(i==14) for i in range(17)]])
    kernel=[D,zero];seen=set()
    for row in frame['vertical_old_sections']:
        x=vector(QQ,row['word']);assert all(v in ZZ for v in x)
        assert tuple(x) not in seen;seen.add(tuple(x))
        S=vector(QQ,[1,x*G*x/2,*x])
        assert list(map(str,S))==row['NS_coordinates'] and S*NS*S==-2
        assert S*NS*D==0 and int(x*G*x)==row['height']
        kernel.append(S)
    assert tuple([0]*17) in seen
    assert matrix(QQ,kernel).rank()==frame['kernel_rank']==7
    assert 19-frame['kernel_rank']==frame['Picard_image_rank_upper_bound']==12
    # Only these exhibited kernel classes are needed for the upper bound;
    # no completeness of the CVP list is assumed by the independent proof.
    R=PolynomialRing(QQ,'t');K=R.fraction_field();t=R.gen()
    def dec(d):return K(R(d['numerator']))/R(d['denominator'])
    ai=list(map(dec,parent['a_invariants']));E=EllipticCurve(K,ai)
    basis=[E(list(map(dec,P))) for P in parent['basis_weierstrass_coordinates']]
    C=basis[14]-basis[15];cx=C[0]+E.b2()/12;cy=C[1]+(E.a1()*C[0]+E.a3())/2
    h=R(generic['pole_h']);shift=R(generic['shift']);assert h*h==C[0].denominator()
    assert len(functions['sections'])==17
    for i,(row,P) in enumerate(zip(functions['sections'],basis)):
        xx=P[0]+E.b2()/12;yy=P[1]+(E.a1()*P[0]+E.a3())/2
        slope=(yy+cy)/(xx-cx);z=(h*slope+shift)/(h*h);W=(2*xx+cx-slope*slope)/h
        assert row['index']==i and dec(row['z'])==z and dec(row['W'])==W
        degree=G[i,i]-(G*w)[i]
        assert degree==row['intersection_degree']==max(z.numerator().degree(),z.denominator().degree())
    checks=[];case1=None
    for index in range(2):
        dp=DIR/f'case-{index:02d}-generic.json';cp=DIR/f'case-{index:02d}-codes.json'
        data,codes=map(read,[dp,cp]);provenance(data);provenance(codes)
        payload=read(SOURCE/f'equation-only-cover-{index:02d}.json')
        q=R(data['quartic']);assert q==R(payload['quartic_coefficients'])
        zvalue=QQ(data['z']);assert str(zvalue)==payload['cover_parameter_z']
        t0,s=map(QQ,data['basepoint']);assert s*s==q(t0) and s
        assert t0==QQ(payload['infinite_base_point_source']['t0'])
        q0,q1,q2,q3,q4=q(t+t0).list();assert q0==s*s
        abc=list(map(QQ,data['Jacobian_cubic_coefficients']))
        assert abc==[q2,q1*q3-4*s*s*q4,s*s*q3*q3+q1*q1*q4-4*s*s*q2*q4]
        J=EllipticCurve(QQ,[0,abc[0],0,abc[1],abc[2]])
        HF=list(map(QQ,data['old_F_restriction']['point']))
        xbar=q1*q1/(4*s*s)-q2;ybar=-(q1*xbar+2*s*s*q3)/(2*s)
        assert HF==[xbar,ybar] and J(HF)
        assert data['old_O_restriction']=='ZERO_D_DOT_O_IS_ZERO'
        assert len(data['old_section_restrictions'])==17
        pairs=[]
        for i,(row,func) in enumerate(zip(data['old_section_restrictions'],functions['sections'])):
            assert row['index']==i and row['degree']==func['intersection_degree']
            z=dec(func['z']);W=dec(func['W'])
            if i==14:
                assert row['status']=='ZERO_NORMALIZED_BASEPOINT' and z(t0)==zvalue and W(t0)==s
                pairs.append(None);continue
            if row['degree']==0:
                assert row['status']=='ZERO_VERTICAL_SECTION' and z!=zvalue
                pairs.append(None);continue
            g=R(row['g']);ww=R(row['Wmod']);xx=R(row['Jacobian_Xmod']);yy=R(row['Jacobian_Ymod'])
            assert g==(z.numerator()-zvalue*z.denominator()).monic()
            assert g.degree()==row['degree'] and g.gcd(g.derivative())==1
            assert W.denominator().gcd(g)==1 and g(t0)
            assert (ww*W.denominator()-W.numerator())%g==0 and (ww*ww-q)%g==0
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
            discriminant=a*a*b*b-4*b**3-4*a**3*c-27*c*c+18*a*b*c
            if discriminant%p==0:
                assert status=='SKIP_BAD_JACOBIAN_REDUCTION';continue
            roots=[x for x in range(p) if (x**3+a*x*x+b*x+c)%p==0]
            assert roots==trial['rational_roots'] and trial['cubic']==[c,b,a,1]
            if not roots:
                assert status=='COMPLETE_ZERO_LOCAL_QUOTIENT' and not trial['rows']
                no2=p if no2 is None else no2
                local.append({'p':p,'generic_local_rank':0,'full_local_dimension':0});continue
            if any(v.denominator()%p==0 for pair in pairs if pair for f in pair for v in f):
                assert status=='SKIP_FULL_DIVISOR_DENOMINATOR';continue
            S=PolynomialRing(GF(p),'x');norms=[];raw=[];failed=False
            for root in roots:
                if HF[0].valuation(p)<0:values=[1]
                else:
                    v=(residue(HF[0],p)-root)%p
                    if not v:v=(3*root*root+2*a*root+b)%p
                    assert v;values=[v]
                for pair in pairs:
                    values.append(1 if pair is None else finite_norm(S(pair[1])-root,S(pair[0])))
                norms.append(values)
                if any(v==0 for v in values):failed=True;break
                raw.append([int(pow(v,(p-1)//2,p)==p-1) for v in values])
            if failed:
                assert status=='SKIP_NONUNIT_DIVISOR_CLASS' and trial['partial_norms']==norms;continue
            assert status=='COMPLETE_FULL_PICARD_REDUCTION' and trial['norms']==norms and trial['rows']==raw
            rows+=raw;owners +=[p]*len(raw)
            dimension=1 if len(roots)==1 else 2
            rank=int(matrix(GF(2),raw).rank());assert rank<=dimension
            local.append({'p':p,'generic_local_rank':rank,'full_local_dimension':dimension})
        M=matrix(GF(2),rows)
        assert rows==codes['matrix_rows'] and owners==codes['row_primes']
        assert M.rank()==codes['rank']==data['rank_upper_bound']==12 and no2 is not None
        checks.append({'index':index,'Picard_image_rank':12,'kernel_rank':7,
          'no_rational_2_torsion_prime':no2,'complete_local_comparisons':local,
          'integral_generators':'O,F and all17 old sections; none omitted from an accepted place',
          'odd_rational_saturation_index':True})
        if index==1:case1=(data,codes,M)
    data,codes,M=case1;q=R(data['quartic']);t0,s=map(QQ,data['basepoint'])
    W=q(0).sqrt();q0,q1,q2,q3,q4=q(t+t0).list();u=-t0
    X=(2*s*(W+s)+q1*u)/(u*u);Y=((X*X-4*s*s*q4)*u-q1*X-2*s*s*q3)/(2*s)
    assert marked['quartic_point']==['0',str(W)] and marked['Jacobian_point']==list(map(str,[X,Y]))
    abc=list(map(QQ,data['Jacobian_cubic_coefficients']))
    J=EllipticCurve(QQ,[0,abc[0],0,abc[1],abc[2]]);assert J([X,Y])
    column=[];evaluations=[]
    for trial in codes['trials']:
        if not trial['status'].startswith('COMPLETE_'):continue
        p=trial['p'];c,b,a,one=trial['cubic'];assert one==1
        bits=[]
        for root in trial['rational_roots']:
            if X.valuation(p)<0:v=1
            else:
                v=(residue(X,p)-root)%p
                if not v:v=(3*root*root+2*a*root+b)%p
            assert v;bits.append(int(pow(v,(p-1)//2,p)==p-1))
        evaluations.append({'p':p,'bits':bits});column+=bits
    assert evaluations==marked['prime_evaluations']
    A=M.augment(matrix(GF(2),len(column),1,column));sep=vector(GF(2),marked['separating_row'])
    assert M.rank()==12 and A.rank()==13 and sep*M==0 and sep*vector(GF(2),column)==1
    assert marked['status']=='FIRST_CLASS_OUTSIDE_FULL_PICARD_IMAGE'
    support=[{'p':p,'root_index':j} for trial in codes['trials'] if trial['status'].startswith('COMPLETE_')
             for j,p in enumerate([trial['p']]*len(trial['rational_roots']))]
    support=[s for s,b in zip(support,sep) if b]
    result={'status':'PASS_INDEPENDENT_FULL_PICARD_RANK12_AND_FIRST_CLASS_RANK13',
      'classification':'new exact arithmetic obstruction on the lower-genus carrier',
      'generic_cases':checks,'first_augmented_rank':13,'separating_support':support,
      'first_class':'RATIONAL_NON_GENERIC_KUMMER_CLASS; SELMER WITH ZERO SHA IMAGE',
      'obstruction':'Every generic section of this alternate fibration specializes inside the full Picard image. The first point and its hyperelliptic conjugate lie outside its rational span. Therefore the infinite rational-point source generated by inherited classes cannot reach either point above302.',
      'boundary':'This is one calibrated carrier. No prospective curve selector, full Selmer dimension, Sha obstruction or matched generic-point control comparison on this new carrier family is claimed.',
      'limits':protocol['limits'],
      'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+[DIR/'case-00-generic.json',DIR/'case-01-generic.json',DIR/'case-00-codes.json',DIR/'case-01-codes.json',Path(__file__)]}}
    retain(DIR/'replay.json',result)
    print(result['status'],'all19 Picard generators accounted for',flush=True)
if __name__=='__main__':
    signal.alarm(25);verify()

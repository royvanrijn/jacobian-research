#!/usr/bin/env sage-python
"""Same full-Picard test on nine generic-point controls, no exceptional input.

Freeze z(tau)=the old section0 pencil parameter and the old64 primes.
Every case is checkpointed; total25-second cap. No new original address,
point search, later point, V3 artifact, class group or Selmer dimension.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
FRAME=ART/'det1092_genus1_picard_image_v1';SOURCE=ART/'det1092_norm8_seed_cover_v2'
DIR=ART/'det1092_genus1_picard_controls_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def residue(v,p):
    v=QQ(v);assert v.denominator()%p
    return int(v.numerator()*v.denominator().inverse_mod(p)%p)
def codes_for(data,primes):
    abc=list(map(QQ,data['cubic']));J=EllipticCurve(QQ,[0,abc[0],0,abc[1],abc[2]])
    HF=list(map(QQ,data['F_point']));R=PolynomialRing(QQ,'t')
    pairs=[None if d['zero'] else (R(d['g']),R(d['X'])) for d in data['divisors']]
    trials=[];rows=[];owners=[];no2=None
    for p in primes:
        trial={'p':p}
        if any(v.denominator()%p==0 for v in abc):
            trial['status']='SKIP_JACOBIAN_DENOMINATOR';trials.append(trial);continue
        if J.discriminant().valuation(p)!=0:
            trial['status']='SKIP_BAD_JACOBIAN_REDUCTION';trials.append(trial);continue
        a,b,c=[residue(v,p) for v in abc];roots=[x for x in range(p) if (x**3+a*x*x+b*x+c)%p==0]
        trial.update(cubic=[c,b,a,1],roots=roots)
        if not roots:
            no2=p if no2 is None else no2
            trial.update(status='COMPLETE_ZERO_QUOTIENT',rows=[],norms=[]);trials.append(trial);continue
        if any(v.denominator()%p==0 for pair in pairs if pair for f in pair for v in f):
            trial['status']='SKIP_FULL_DIVISOR_DENOMINATOR';trials.append(trial);continue
        S=PolynomialRing(GF(p),'t');norms=[];raw=[];failed=False
        for root in roots:
            value=1 if HF[0].valuation(p)<0 else (residue(HF[0],p)-root)%p
            if not value:value=(3*root*root+2*a*root+b)%p
            vals=[value]
            for pair in pairs:
                vals.append(1 if pair is None else int(S(pair[0]).resultant(S(pair[1])-root)))
            norms.append(vals)
            if any(v==0 for v in vals):failed=True;break
            raw.append([int(pow(v,(p-1)//2,p)==p-1) for v in vals])
        if failed:
            trial.update(status='SKIP_NONUNIT_DIVISOR_CLASS',partial_norms=norms);trials.append(trial);continue
        if len(roots)==3:assert all(sum(r[j] for r in raw)%2==0 for j in range(18))
        trial.update(status='COMPLETE_FULL_PICARD_REDUCTION',norms=norms,rows=raw)
        rows+=raw;owners +=[p]*len(raw);trials.append(trial)
    M=matrix(GF(2),len(rows),18,[v for r in rows for v in r])
    assert M.rank()<=12
    return {'trials':trials,'rows':rows,'row_primes':owners,'rank':int(M.rank()),'no_two_torsion_prime':no2}
def main():
    paths=[ART/'curve302_recovered_mw17_parent_v1.json',SOURCE/'generic.json',
       FRAME/'frame.json',FRAME/'generic-restrictions.json',FRAME/'replay.json',
       ART/'det1092_rr_generic_point_controls_v2/protocol.json',Path(__file__)]
    parent,pencil,frame,functions,proof,roster=map(read,paths[:6])
    protocol={'classification':'frozen source-only specificity controls',
       'rule':'Use the same orbit20124 pencil. For each of the nine unchanged original addresses tau, take exactly z=z_section0(tau) from the generic restriction function already constructed, and mark that generic section0 point. Compare its Jacobian class with the full inherited Picard image using exactly the old64-prime pool. Do not replace singular, degenerate or unresolved cases.',
       'limits':{'wall_seconds':25,'old_addresses':9,'new_addresses':0,'generic_sections_used_for_selection':[0],
          'NS_generators':19,'primes':64,'exceptional_point_inputs':0,'point_searches':0,
          'V3_inputs':0,'pilot_changes':0,'class_groups':0,'Selmer_dimensions':0},
       'primes':roster['primes'],'cases':[{k:c[k] for k in ['label','parameter']} for c in roster['cases']],
       'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths}}
    retain(DIR/'protocol.json',protocol)
    for path,digest in proof['inputs'].items():assert sha(ROOT/path)==digest
    assert frame['Picard_image_rank_upper_bound']==12
    R=PolynomialRing(QQ,'t');K=R.fraction_field();t=R.gen();Z=PolynomialRing(QQ,'z')
    def dec(d):return K(R(d['numerator']))/R(d['denominator'])
    coefficients=list(map(Z,pencil['quartic_t_coefficients_in_z']))
    ai=list(map(dec,parent['a_invariants']));E=EllipticCurve(K,ai)
    S0=E(list(map(dec,parent['basis_weierstrass_coordinates'][0])))
    C=E(list(map(dec,parent['basis_weierstrass_coordinates'][14])))-E(list(map(dec,parent['basis_weierstrass_coordinates'][15])))
    cx=C[0]+E.b2()/12;cy=C[1]+(E.a1()*C[0]+E.a3())/2
    h=R(pencil['pole_h']);shift=R(pencil['shift'])
    selector=dec(functions['sections'][0]['z']);W0=dec(functions['sections'][0]['W'])
    assert max(selector.numerator().degree(),selector.denominator().degree())==4
    selection=[{'label':case['label'],'tau':case['parameter'],'z':str(selector(QQ(case['parameter'])))} for case in roster['cases']]
    retain(DIR/'selection.json',{'classification':'equations selected before any marked-class arithmetic',
       'selector_degree':4,'cases':selection,'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [DIR/'protocol.json',paths[3]]}})
    results=[]
    for index,selected in enumerate(selection):
        tau=QQ(selected['tau']);zvalue=QQ(selected['z']);q=R([v(zvalue) for v in coefficients])
        assert q.degree()==4 and q.gcd(q.derivative())==1 and q.gcd(R(E.discriminant()))==1
        sec=pencil['sections'][0]
        def at(d):return Z(d['numerator'])(zvalue)/Z(d['denominator'])(zvalue)
        t0,s=at(sec['t_of_z']),at(sec['W_of_z']);assert s and s*s==q(t0)
        q0,q1,q2,q3,q4=q(t+t0).list();assert q0==s*s
        cubic=[q2,q1*q3-4*s*s*q4,s*s*q3*q3+q1*q1*q4-4*s*s*q2*q4]
        xf=q1*q1/(4*s*s)-q2;yf=-(q1*xf+2*s*s*q3)/(2*s)
        J=EllipticCurve(QQ,[0,cubic[0],0,cubic[1],cubic[2]]);J([xf,yf])
        divs=[]
        for i,func in enumerate(functions['sections']):
            degree=func['intersection_degree']
            if i==14 or degree==0:
                divs.append({'index':i,'zero':True,'degree':degree});continue
            z=dec(func['z']);W=dec(func['W']);g=R(z.numerator()-zvalue*z.denominator()).monic()
            assert g.degree()==degree and g.gcd(g.derivative())==1 and g(t0)
            assert g.gcd(W.denominator())==1
            ww=(W.numerator()*W.denominator().inverse_mod(g))%g
            du=t-t0;xx=((2*s*(ww+s)+q1*du)*(du*du).inverse_mod(g))%g
            yy=(((xx*xx-4*s*s*q4)*du-q1*xx-2*s*s*q3)/(2*s))%g
            assert (ww*ww-q)%g==0 and (yy*yy-xx**3-cubic[0]*xx*xx-cubic[1]*xx-cubic[2])%g==0
            divs.append({'index':i,'zero':False,'degree':degree,'g':list(map(str,g.list())),
                         'W':list(map(str,ww.list())),'X':list(map(str,xx.list())),'Y':list(map(str,yy.list()))})
        data={'index':index,**selected,'quartic':list(map(str,q.list())),
           'basepoint':[str(t0),str(s)],'cubic':list(map(str,cubic)),
           'F_point':list(map(str,[xf,yf])),'divisors':divs,'rank_upper_bound':12}
        retain(DIR/f'case-{index:02d}-generic.json',data)
        codes=codes_for(data,protocol['primes']);retain(DIR/f'case-{index:02d}-codes.json',codes)
        # Only after freezing this case's complete inherited code, evaluate
        # the marked rational point from the already selected generic S0.
        wm=W0(tau);assert wm*wm==q(tau);du=tau-t0;assert du
        xm=(2*s*(wm+s)+q1*du)/(du*du)
        ym=((xm*xm-4*s*s*q4)*du-q1*xm-2*s*s*q3)/(2*s);J([xm,ym])
        m=h(tau)*zvalue-shift(tau)/h(tau)
        oldX=(h(tau)*wm-cx(tau)+m*m)/2
        oldx=oldX-E.b2()(tau)/12;oldY=m*(oldX-cx(tau))-cy(tau)
        oldy=oldY-(E.a1()(tau)*oldx+E.a3()(tau))/2
        assert oldx==S0[0](tau) and oldy==S0[1](tau)
        column=[]
        for trial in codes['trials']:
            if not trial['status'].startswith('COMPLETE_'):continue
            p=trial['p'];c,b,a,one=trial['cubic'];assert one==1
            for root in trial['roots']:
                v=1 if xm.valuation(p)<0 else (residue(xm,p)-root)%p
                if not v:v=(3*root*root+2*a*root+b)%p
                assert v;column.append(int(pow(v,(p-1)//2,p)==p-1))
        M=matrix(GF(2),codes['rows']);A=M.augment(matrix(GF(2),len(column),1,column))
        separator=None
        if A.rank()>M.rank():separator=next(v for v in M.left_kernel().basis() if v*vector(GF(2),column))
        result={'index':index,'label':selected['label'],'tau':str(tau),'z':str(zvalue),
           'status':'GENERIC_OLD_POINT_BUT_NON_GENERIC_CARRIER_CLASS' if M.rank()==12 and A.rank()==13 and codes['no_two_torsion_prime'] else 'UNRESOLVED_FIXED_CONTROL',
           'generic_rank':int(M.rank()),'augmented_rank':int(A.rank()),
           'old_elliptic_point':list(map(str,[oldx,oldy])),
           'marked_quartic_point':[str(tau),str(wm)],'marked_Jacobian_point':list(map(str,[xm,ym])),
           'marked_column':column,'separating_row':list(map(int,separator)) if separator is not None else None,
           'inputs':{str(p.relative_to(ROOT)):sha(p) for p in
              [DIR/f'case-{index:02d}-generic.json',DIR/f'case-{index:02d}-codes.json',DIR/'protocol.json']}}
        retain(DIR/f'case-{index:02d}-marked.json',result)
        results.append({'index':index,'label':selected['label'],'status':result['status'],
                        'generic_rank':int(M.rank()),'augmented_rank':int(A.rank())})
        print('CHECKPOINT_GENERIC_POINT_CONTROL',index,M.rank(),A.rank(),flush=True)
    retain(DIR/'panel.json',{'status':'COMPLETED_FROZEN_GENUS1_GENERIC_POINT_CONTROLS',
       'classification':'source-only specificity comparison; independent replay required',
       'cases':results,'selector_degree':4,'limits':protocol['limits'],
       'boundary':'No exceptional old-elliptic point was selected. A non-generic carrier Jacobian class is not by itself a non-generic original-fibre point. No full Selmer dimension or Sha event is claimed.',
       'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+[DIR/'protocol.json',DIR/'selection.json']}})
if __name__=='__main__':
    signal.alarm(25);DIR.mkdir(parents=True,exist_ok=True);main()

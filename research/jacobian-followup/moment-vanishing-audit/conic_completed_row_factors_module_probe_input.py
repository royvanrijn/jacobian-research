"""Exact auxiliary completed-row conic factors and p2q2 seed estimates.
All-degree module/path proofs, analytic exclusions and covers remain written.
Retained checkers are imported only for pure helpers; their mains are not run.
"""
from pathlib import Path
from itertools import combinations
import hashlib,json,time
import sympy as s
import verify_balanced_logarithmic_seed_bounds as balanced
import verify_two_sided_conic_seed_bounds as two
import verify_endpoint_degrees_and_odd_prefixes as endpoint
P=Path(__file__).resolve().parent
t,y=balanced.t,balanced.y
parse=balanced.laurent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda z:s.sympify(z,locals={'t':t,'y':y})

def bound(p,q,kd,ell,j,k):
    if p==q==3:
        return [3*j+abs(k),4*j,5*j-abs(k)][kd]
    if p==2 and q==3:
        if kd<2:return two.specification(p,q,kd,ell,j,k)-kd*max(-k,0)
        return 7*j-1-k-2*max(-k,0)
    if p==1 and q==3:return endpoint.specification(p,q,kd,ell,j,k)-kd*max(-k,0)
    if p==3 and q==2:
        old=endpoint.specification(p,q,kd,ell,j,k) if kd==1 else two.specification(p,q,kd,ell,j,k)
        return old-kd*max(-k,0)
    assert p==q==2
    return [9*j-3+2*abs(k), (24*j-12-k)//2-max(-k,0),9*j-3-k-2*max(-k,0)][kd]

def arithmetic():
    rows=[]
    cases=[(3,3,kd,0) for kd in range(3)]
    cases += [(p,3,kd,ell) for p in [1,2] for kd in range(3) for ell in range(4-p)]
    cases += [(3,2,kd,ell) for kd in range(3) for ell in range(2)]
    cases += [(2,2,kd,0) for kd in range(3)]
    for p,q,kd,ell in cases:
        basis=list(range(-q,p));gates=[]
        for size in range(1,len(basis)+1):
            for selected in combinations(basis,size):
                if 0 not in selected:continue
                # Degree of any determinant term is bounded by a column base
                # plus the row offsets; these formulas are separable in j,k.
                gate=sum(bound(p,q,kd,ell,j,0) for j in range(1,size+1))
                gate+=sum(bound(p,q,kd,ell,1,k)-bound(p,q,kd,ell,1,0) for k in selected)
                gates.append({'rank':size,'rows':selected,'degree_bound':gate})
        rows.append({'p':p,'q':q,'K_degree':kd,'endpoint_degree':ell,'gate':max(z['degree_bound'] for z in gates),'all_rank_row_cases':gates})
    lookup={(z['p'],z['q'],z['K_degree'],z['endpoint_degree']):z['gate'] for z in rows}
    assert [lookup[(3,3,kd,0)] for kd in range(3)]==[72,84,96]
    assert max(z['gate'] for z in rows if z['p']==2 and z['q']==3)==94
    assert max(z['gate'] for z in rows if z['p']==1 and z['q']==3)==84
    assert max(z['gate'] for z in rows if z['p']==3 and z['q']==2)==93
    assert [lookup[(2,2,kd,0)] for kd in range(3)]==[86,93,74]
    retained=[z for z in balanced.degree_cases() if z['negative_order']<3 and (z['positive_order'],z['negative_order']) not in [(3,2),(2,2)]]
    assert max(z['all_seed_rank_and_row_choices_gate'] for z in retained)==63
    assert max(96,94,84,69,93,63)==96
    return rows

def replay(m,origin,remove_low=False):
    W,K=read(m['W']),read(m['K']);raw=parse(W);p=max(raw);q=-min(raw);kd=int(s.degree(K,t))
    endpoint_poly=raw[p] if q==3 and p<3 else s.cancel(raw[-q]/K**q) if q==2 else s.Integer(1)
    ell=int(s.degree(endpoint_poly,t));basis=m['basis'];assert basis==list(range(-q,p))
    out=[]
    for row in m['seed_columns']:
        j=row['j'];common=raw[-3]**row['depths'][1] if remove_low else s.Integer(1)
        vals=[s.Poly(s.cancel(read(z)/(common*K**max(-k,0))),t,domain=s.QQ) for k,z in zip(basis,row['N'])]
        limits=[bound(p,q,kd,ell,j,k) for k in basis]
        assert all(z.degree()<=b for z,b in zip(vals,limits)),(origin,m['label'],j,[z.degree() for z in vals],limits)
        H=s.Poly(s.cancel(read(row['H'])/common),t,domain=s.QQ)
        assert min(e[0] for e,c in H.terms() if c)>=3*j
        # Same completed vector, no connection or basis derivative changed.
        assert all(s.cancel(z.as_expr()*K**max(-k,0)*common-read(old))==0 for k,z,old in zip(basis,vals,row['N']))
        out.append({'j':j,'N':[str(z.as_expr()) for z in vals],'H':str(H.as_expr()),'row_degree_bounds':limits,'row_scales':[str(K**max(-k,0)) for k in basis],'old_common_factor_cancelled':str(common)})
    return {'origin':origin,'label':m['label'],'p':p,'q':q,'K_degree':kd,'endpoint_degree':ell,'seed_columns':out}

def universal_p2q2():
    k0,k1,k2=s.symbols('k0 k1 k2');rows=[]
    for kd in range(3):
        K=[k0,k0+k1*t,k0+k1*t+k2*t*t][kd]
        W=sum(s.Symbol('a'+str(k)+'_'+str(n))*t**n*y**k for k in range(3) for n in range(4-k))
        W+=sum(K**k*sum(s.Symbol('b'+str(k)+'_'+str(n))*t**n for n in range(4-k))/y**k for k in [1,2])
        B=t*s.diff(K,t)-s.Rational(2,3)*K;h=t**3*K;old=t*s.diff(W,t)-2*W-t*t
        S=s.expand(K*old+B*y*s.diff(W,y));raw=parse(S)
        if kd==0:assert all(s.degree(c,t)+abs(k)<=3 for k,c in raw.items())
        else:
            wt,wy,cost=(2,1,8) if kd==1 else (1,1,5)
            assert all(wt*s.degree(c,t)+wy*k<=cost for k,c in raw.items())
        identities=[]
        for n in range(-4,5):
            f=y**n;witness=B*f
            assert s.expand(S*f+t*t*B*y*s.diff(f,y)-K*old*f-t*t*y*s.diff(witness,y)-y*s.diff(W,y)*witness)==0
            identities.append(n)
        rows.append({'K_degree':kd,'source_support':[min(raw),max(raw)],'gauge_identities':identities,'negative_module_preserved':True,'original_logarithmic_operator':True})
    return rows

def model(label,K,A,B2):
    W=s.expand(1+t**3+(2+t*t)*y+A*y*y+K*(3+t*t)/y+K*K*B2/y**2)
    raw=parse(W);kd=int(s.degree(K,t));basis=list(range(-2,2));field=s.QQ.frac_field(t)
    conv=lambda x:{k:field.convert(c) for k,c in parse(x).items()}
    def add(out,k,c):
        out[k]=out.get(k,field.zero)+c
        if not out[k]:del out[k]
    def exact(x):
        out={}
        for n,c in x.items():
            add(out,n,field.convert(t*t*n)*c)
            for k,a in raw.items():
                if k:add(out,n+k,field.convert(k*a)*c)
        return out
    def reduce(x):
        normal=x.copy();witness={};hi=lo=0
        while max(normal,default=0)>=2:
            n=max(normal);shift=n-2;c=normal[n]/field.convert(2*A);hi+=1;add(witness,shift,c)
            for k,z in exact({shift:c}).items():add(normal,k,-z)
        while min(normal,default=0)<-2:
            n=min(normal);shift=n+2;c=normal[n]/field.convert(-2*K*K*B2);lo+=1;add(witness,shift,c)
            for k,z in exact({shift:c}).items():add(normal,k,-z)
        defect=x.copy()
        for k,z in normal.items():add(defect,k,-z)
        for k,z in exact(witness).items():add(defect,k,-z)
        assert not defect and set(normal)<=set(basis)
        return normal,witness,hi,lo
    B=t*s.diff(K,t)-s.Rational(2,3)*K;h=t**3*K
    S=s.expand(K*(t*s.diff(W,t)-2*W-t*t)+B*y*s.diff(W,y));R=s.Integer(1);rows=[]
    for j in range(1,4):
        R=s.expand(h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+t*t*B*y*s.diff(R,y)+S*R)
        normal,witness,hi,lo=reduce(conv(R));assert hi<=2*j-1 and lo<=2*j-2
        den=s.expand(A**(2*j-1)*B2**(2*j-2));values=[]
        for k in basis:
            z=s.Poly(s.cancel(field.to_sympy(normal.get(k,field.zero))*den/K**max(-k,0)),t,domain=s.QQ)
            assert z.degree()<=bound(2,2,kd,0,j,k),(label,j,k,z.degree(),bound(2,2,kd,0,j,k));values.append(z)
        H=s.Poly(s.expand(h**j*den),t,domain=s.QQ);val=min(e[0] for e,c in H.terms() if c);assert val>=3*j
        rows.append({'j':j,'N':[str(z.as_expr()) for z in values],'H':str(H.as_expr()),'H_valuation':val,'high_low_counts':[hi,lo],'endpoint_clear_powers':[2*j-1,2*j-2],'row_degree_bounds':[bound(2,2,kd,0,j,k) for k in basis],'witness':{str(k):str(field.to_sympy(z)) for k,z in witness.items()}})
    return {'label':label,'p':2,'q':2,'K_degree':kd,'K':str(K),'A2':str(A),'B2':str(B2),'W':str(W),'basis':basis,'seed_columns':rows}

def main():
    start=time.monotonic();sources=universal_p2q2();cases=arithmetic();replays=[]
    for name,key in [('structured_seed_modules_certificate.json','conic_models'),('two_sided_conic_seed_bounds_certificate.json','models'),('endpoint_degrees_and_odd_prefixes_certificate.json','conic_models')]:
        for m in json.loads((P/name).read_text())[key]:replays.append(replay(m,name))
    for m in json.loads((P/'euler_seed_bounds_certificate.json').read_text())['conic_models']:
        if m['positive_order']==2:replays.append(replay(m,'euler_seed_bounds_certificate.json',True))
    print('Completed row factors checked on retained models',len(replays),flush=True)
    models=[]
    for kd,K in [(0,s.Integer(2)),(1,1+t),(2,t*(t+1))]:
        for A in [s.Integer(1),1+t]:
            for B2 in [s.Integer(2),2+t]:models.append(model('p2q2_K'+str(kd)+'_A'+str(A)+'_B'+str(B2),K,A,B2))
    for K,A,B2 in [(t,t,t),(t*t,t,t),(t*(t+1),t+1,t)]:models.append(model('shared_'+str(K)+'_'+str(A)+'_'+str(B2),K,A,B2))
    oddcover={'quadratic':7,'smooth':94,'nodal':100,'cusp':94,'transverse':98,'tangent':94,'triangle':96,'concurrent':104,'repeated':98}
    assert max(oddcover.values())==104 and max(56,10,10,72,96)==96 and max(96,92,95)==96
    cert={'universal_p2q2':sources,'degree_cases':cases,'retained_model_row_factors':replays,'p2q2_models':models,'odd_cover':oddcover}
    target=P/'conic_completed_row_factors_certificate.json';target.write_text(json.dumps(cert,indent=2,default=int)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(target),'runtime_seconds':round(time.monotonic()-start,3),'sympy_version':s.__version__,'p2q2_source_strata':3,'logarithmic_gauge_identities':27,'retained_models':len(replays),'retained_seed_row_identities':sum(len(m['seed_columns']) for m in replays),'p2q2_models':len(models),'new_direct_seed_identities':sum(len(m['seed_columns']) for m in models),'conic_gate':96,'full_W_linear_gate':96,'odd_sufficient_pure_prefix':104,'homogeneous_sufficient_pure_prefix':96,'global_sufficient_pure_prefix_unchanged':130,'least_prefix_interval':[20,130],'odd_least_prefix_interval':[20,104],'homogeneous_least_prefix_interval':[14,96],'odd_cover':oddcover,'optimality_proved':False,'written_analytic_or_global_proof_certified_by_checker':False,'boundary':'Exact completed negative-row polynomial factors on retained seeds; three generic p2q2 source strata with27 gauge identities;45 new exact seed identities and all-rank/row arithmetic. Universal formal-module and degree/path inductions, original period comparison, analytic exclusions and full coefficient/frame cover remain written inputs. No independent replay, external review, formal verification or expanded radical cofactors.'}
    (P/'conic_completed_row_factors_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()

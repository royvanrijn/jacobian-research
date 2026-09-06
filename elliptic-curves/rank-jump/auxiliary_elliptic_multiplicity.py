#!/usr/bin/env python3
"""Small-field exclusion of hidden elliptic factors in retained auxiliaries."""
import argparse
from pathlib import Path
import subprocess
import retrospective as r
import residual_alignment as source

PROTOCOL=Path(__file__).with_name('AUXILIARY_ELLIPTIC_MULTIPLICITY_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_auxiliary_elliptic_multiplicity_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-auxiliary-elliptic-multiplicity-v1'
MODELS={'rational':[0,11,0,-14,1],'Sha':[0,-11,0,-14,-1]}


def count(coeff,p,n):
    from sage.all import GF,PolynomialRing
    F=GF(p**n,'a');R=PolynomialRing(F,'x');f=R(list(map(F,coeff)))
    total=0;zeros=0;residues=0
    for x in F:
        y=f(x)
        if y==0:total+=1;zeros+=1
        elif y.is_square():total+=2;residues+=1
    infinity=2 if f.leading_coefficient().is_square() else 0
    return {'extension_degree':n,'field_size':p**n,'affine_branch_points':zeros,
        'nonzero_square_values':residues,'points_at_infinity':infinity,'points':total+infinity}


def weil_from_counts(p,counts):
    from sage.all import QQ,PolynomialRing
    R=PolynomialRing(QQ,'T');T=R.gen()
    s1,s2,s3=[QQ(p**i+1-row['points']) for i,row in enumerate(counts,1)]
    e1=s1;e2=(s1*s1-s2)/2;e3=(s1**3-3*s1*s2+2*s3)/6
    f=T**6-e1*T**5+e2*T**4-e3*T**3+p*e2*T*T-p*p*e1*T+p**3
    assert all(c.denominator()==1 for c in f)
    return f


def compute(index):
    from sage.all import QQ,GF,PolynomialRing,EllipticCurve
    from sage.version import version
    row=r.read(source.OUTPUT)['rows'][index]
    R=PolynomialRing(QQ,'x');x=R.gen()
    f,g=[sum(QQ(c)*x**(4-i) for i,c in enumerate(h)) for h in row['quartic_coefficients']]
    h=f*g;assert h.degree()==8 and h.gcd(h.derivative())==1
    records=[];conclusion='UNKNOWN'
    for p in r.read(PROTOCOL)['primes']:
        fp=h.change_ring(GF(p))
        if fp.degree()!=8 or fp.gcd(fp.derivative())!=1:
            records.append({'prime':p,'status':'SKIP_BAD_DISPLAYED_REDUCTION'});continue
        assert all(EllipticCurve(QQ,v).discriminant()%p for v in MODELS.values())
        counts=[count(h.list(),p,n) for n in [1,2,3]]
        W=weil_from_counts(p,counts);T=W.parent().gen();elliptic=[]
        for name,model in MODELS.items():
            E=EllipticCurve(GF(p),model);a=p+1-E.cardinality();P=T*T-a*T+p
            elliptic.append({'name':name,'point_count':int(E.cardinality()),'trace':int(a),
                'Frobenius_coefficients_ascending':list(map(str,P)),
                'remainder_coefficients_ascending':list(map(str,W%P)),'factor_excluded':bool(W%P)})
        records.append({'prime':p,'status':'PASS','counts':counts,
            'Jacobian_Frobenius_coefficients_ascending':list(map(str,W)),'elliptic_checks':elliptic})
        if all(x['factor_excluded'] for x in elliptic):conclusion='BOTH_Q_HOM_SPACES_ZERO';break
    return {'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),PROTOCOL,source.OUTPUT)},'alignment':row['alignment'],'status':conclusion,
        'software':{'sage':version},'product_polynomial_ascending':list(map(str,h)),'records':records,
        'boundary':'No geometric isogeny exclusion or rational point/rank claim. A nonzero remainder is used only as a necessary-condition obstruction to a Q-isogeny factor.'}


def capture(check=False):
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for i in range(2):
        if check:assert compute(i)==r.read(OUTPUT)['rows'][i];print('PASS elliptic multiplicity',i);continue
        path=WORK/f'case-{i}.json'
        if not path.exists():
            with (WORK/f'case-{i}.log').open('x') as log:
                try:
                    p=subprocess.run(['sage','-python',str(Path(__file__).resolve()),'worker','--index',str(i),
                        '--destination',str(path)],cwd=r.ROOT,stdout=log,stderr=log,timeout=30)
                    reason=None if p.returncode==0 else 'worker failure'
                except subprocess.TimeoutExpired:reason='30-second timeout'
                if reason and not path.exists():r.write_new(path,{'status':'UNKNOWN','reason':reason})
        row=r.read(path);rows.append(row);print(row['alignment'],row['status'],flush=True)
    if not check:r.write_new(OUTPUT,{'schema':'rank-jump.auxiliary-elliptic-multiplicity.v1','rows':rows})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker','check'])
    p.add_argument('--index',type=int);p.add_argument('--destination',type=Path);a=p.parse_args()
    if a.mode=='worker':r.write_new(a.destination,compute(a.index))
    else:capture(a.mode=='check')

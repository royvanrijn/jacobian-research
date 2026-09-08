#!/usr/bin/env sage-python
"""Complete real Kummer images from inherited divisors and trace signatures."""
import argparse,hashlib,json,signal
from pathlib import Path
from sage.all import QQ,PolynomialRing,matrix,identity_matrix,pari
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_rr_real_images_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def congruence(H):
    B=H.copy();V=identity_matrix(QQ,H.nrows());n=H.nrows()
    for j in range(n):
        if not B[j,j]:
            options=[k for k in range(j+1,n) if B[k,k]]
            C=identity_matrix(QQ,n)
            if options:
                k=options[0];C.swap_columns(j,k)
            else:
                options=[(k,l) for k in range(j,n) for l in range(k+1,n) if B[k,l]]
                assert options
                k,l=options[0];C.swap_columns(j,k);C[l,j]+=1
            B=C.transpose()*B*C;V=V*C
        assert B[j,j]
        C=identity_matrix(QQ,n)
        for k in range(j+1,n):C[j,k]=-B[j,k]/B[j,j]
        B=C.transpose()*B*C;V=V*C
    assert B.is_diagonal() and V.det() and V.transpose()*H*V==B
    return V,B
def construct(i):
    paths=[ART/'det1092_rr_good_local_images_v1'/('case-%02d.json'%j) for j in range(10)]
    protocol={'classification':'frozen real-place local descent',
              'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths},
              'limits':{'cases':10,'seconds_per_case':25,'generic_divisors_per_case':17,
                        'point_searches':0,'class_groups':0,'remote_arithmetic':0,'pilot_changes':0},
              'selector':'Use inherited divisors in existing order until the required real dimension is reached. No marked point.',
              'script_sha256':sha(Path(__file__))}
    OUT.mkdir(parents=True,exist_ok=True);retain(OUT/'protocol.json',protocol)
    source=json.loads(paths[i].read_text());R=PolynomialRing(QQ,'T');T=R.gen();q=R(source['q']).monic()
    r=int(pari(q).polsturm());dim=max(0,r//2-1);assert dim<=1,'HIGHER_DIMENSION_NOT_IMPLEMENTED'
    x0=QQ(source['base_x']);assert (QQ(source['scale'])*R(source['q'])(x0)).is_square()
    prefix={'classification':'exact real-place descent; independent replay required','case_index':i,
            'source':str(paths[i].relative_to(ROOT)),'source_sha256':sha(paths[i]),
            'real_branch_points':r,'real_Kummer_dimension':dim,'limits':protocol['limits']}
    retain(OUT/('case-%02d-input.json'%i),prefix)
    tests=[];rank=0
    if dim:
        A=matrix(QQ,6,6)
        for j in range(5):A[j+1,j]=1
        for j in range(6):A[j,5]=-q[j]
        powers=[identity_matrix(QQ,6)]
        for j in range(15):powers.append(powers[-1]*A)
        moments=[m.trace() for m in powers]
        for j,row in enumerate(source['generic_divisors']):
            g=R(row);degree=g.degree()
            # Division by the real anchor has the same sign as multiplication.
            h=R(((-1)**degree*g*R([x0,-1])**degree)%q)
            assert h.gcd(q)==1
            H=matrix(QQ,6,6,[sum(h[k]*moments[k+a+b] for k in range(6)) for a in range(6) for b in range(6)])
            V,D=congruence(H);signature=sum(d.sign() for d in D.diagonal())
            mixed=abs(signature)<r
            item={'generic_divisor_index':j,'h_mod_q':list(map(str,h.list())),
                  'congruence':[list(map(str,row)) for row in V.rows()],
                  'diagonal':list(map(str,D.diagonal())),'signature':int(signature),'mixed_real_signs':mixed}
            tests.append(item);retain(OUT/('case-%02d-divisor-%02d.json'%(i,j)),item)
            if mixed:rank=1;break
    assert rank==dim
    result={**prefix,'status':'PASS_COMPLETE_REAL_KUMMER_IMAGE_SPANNED_BY_INHERITED_DIVISORS',
            'generic_image_rank':rank,'divisor_tests':tests,
            'scope':'Complete real local image only. No new rational point or full global Selmer computation.',
            'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [paths[i],OUT/'protocol.json',Path(__file__)]}}
    retain(OUT/('case-%02d.json'%i),result)
    print('case',i,'real dimension',dim,'divisors tested',len(tests),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--case',type=int,required=True,choices=range(10));args=ap.parse_args()
    signal.alarm(25);construct(args.case)

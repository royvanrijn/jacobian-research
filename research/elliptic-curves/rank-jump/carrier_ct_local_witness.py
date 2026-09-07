#!/usr/bin/env python3
"""Observe and independently replay one Fisher/Cassels local-symbol certificate."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import subprocess
import retrospective as r
import labelled_carrier_ct as base
from cover_experiment import evaluate

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'CARRIER_CT_LOCAL_WITNESS_PROTOCOL.json'
INPUT=r.OUT/'rank_jump_carrier_ct_local_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_carrier_ct_local_witness_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-carrier-ct-local-witness-v1'


def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    raw=(base.WORK/'upstream_ellrank.c').read_bytes()
    assert r.digest(raw)==r.read(base.PROTOCOL)['upstream']['ellrank_c_sha256']
    text=raw.decode();replacements=[]
    def insert(old,new):
        nonlocal text
        assert text.count(old)==1,(old,text.count(old));text=text.replace(old,new);replacements.append([old,new])
    insert(base.MARKER,base.MARKER+'\n'+base.OBSERVATIONS)
    old='      GEN Mjj = gcoeff(M,j,j);'
    insert(old,old+'\n      rank_jump_observe_pair = (i == 3 && j == 1);')
    old='  long i, e = 0, lF = lg(F);'
    observe='\n  if (rank_jump_observe_pair) {\n'+''.join('    err_printf("RJ_PAIR_'+name+' %Ps\\n", '+name+');\n' for name in ['q1','q2','q3','T','m','gam','a','F'])+'  }'
    insert(old,old+observe)
    old='  if (!gequal0(c)) return c;'
    insert(old,'  if (!gequal0(c)) { if (rank_jump_observe_pair) err_printf("RJ_LOCAL %Ps %Ps %Ps\\n",p,a,c); return c; }')
    old='      if (!gequal0(c)) return gerepileupto(av,c);'
    insert(old,'      if (!gequal0(c)) { if (rank_jump_observe_pair) err_printf("RJ_LOCAL %Ps %Ps %Ps\\n",p,b,c); return gerepileupto(av,c); }')
    old='  return gc_long(av,e);'
    insert(old,'  if (rank_jump_observe_pair) err_printf("RJ_PAIR_PARITY %ld\\n",e);\n'+old)
    prefix=''.join('#define '+n+' rank_jump_'+n+'\n' for n in base.EXPORTS)+'static int rank_jump_observe_pair = 0;\n'
    modified=WORK/'observed_ellrank.c';binary=WORK/'observed_ellrank'
    code=(prefix+text+base.MAIN).encode()
    if not modified.exists():base.new_bytes(modified,code)
    assert modified.read_bytes()==code
    cmd=['cc','-O2','-I'+str(base.RUNTIME/'include/pari'),str(modified),'-L'+str(base.RUNTIME/'lib'),
         '-Wl,-rpath,'+str(base.RUNTIME/'lib'),'-lpari','-lm','-o',str(binary)]
    if not (WORK/'compile.json').exists():
        p=subprocess.run(cmd,capture_output=True,text=True,timeout=30)
        r.write_new(WORK/'compile.json',{'command':cmd,'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr})
    assert r.read(WORK/'compile.json')['returncode']==0
    prior=r.read(base.INPUT)
    if not (WORK/'descent.json').exists():
        p=subprocess.run([str(binary),'['+','.join(prior['original_target_model'])+']'],capture_output=True,text=True,timeout=30)
        r.write_new(WORK/'descent.json',{'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr})
    d=r.read(WORK/'descent.json');assert d['returncode']==0;lines=d['stderr'].splitlines()
    def one(prefix):
        a=[x[len(prefix):] for x in lines if x.startswith(prefix)];assert len(a)==1;return a[0]
    assert one('Cassels Pairing: ')==prior['Cassels_matrix_GP']
    assert [x[len('RJ_BASIS '):].split(' ',1)[1] for x in lines if x.startswith('RJ_BASIS ')]==prior['basis_GP']
    pair={k:one('RJ_PAIR_'+k+' ') for k in ['q1','q2','q3','T','m','gam','a','F']}
    local=[x[len('RJ_LOCAL '):].split(' ') for x in lines if x.startswith('RJ_LOCAL ')]
    assert all(len(x)==3 for x in local)
    r.write_new(INPUT,{'schema':'rank-jump.carrier-ct-local-inputs.v1',
                       'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (PROTOCOL,base.INPUT,base.OUTPUT,Path(base.__file__),Path(__file__))},
                       'source_sha256':r.digest(raw),'instrumented_source_sha256':r.digest(code),'observational_replacements':replacements,
                       'pair_GP':pair,'local_GP':local,'pairing_bit':int(one('RJ_PAIR_PARITY ')),
                       'compile':r.read(WORK/'compile.json'),'transcript':d})
    print('Captured',len(local),'local evaluations; pairing bit',one('RJ_PAIR_PARITY '))


def valuation_unit(q,p):
    q=F(q);n,d=q.numerator,q.denominator;v=0
    while n%p==0:n//=p;v+=1
    while d%p==0:d//=p;v-=1
    return v,F(n,d)


def hilbert(a,b,p):
    alpha,u=valuation_unit(a,p);beta,v=valuation_unit(b,p)
    if p==2:
        u,v=r.mod(u,8),r.mod(v,8)
        exponent=((u-1)//2)*((v-1)//2)+alpha*((v*v-1)//8)+beta*((u*u-1)//8)
        return -1 if exponent%2 else 1
    sign=-1 if (alpha*beta*((p-1)//2))%2 else 1
    if beta%2:sign*=1 if pow(r.mod(u,p),(p-1)//2,p)==1 else -1
    if alpha%2:sign*=1 if pow(r.mod(v,p),(p-1)//2,p)==1 else -1
    return sign


def compute():
    from sage.all import QQ,ZZ,PolynomialRing,pari
    inp=r.read(INPUT)
    for path,sha in inp['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    R=PolynomialRing(QQ,'v');v=R.gen()
    def polynomial(text):
        f=pari(text);return R([QQ(f.polcoef(i)) for i in range(int(f.poldegree())+1)])
    pair=inp['pair_GP'];qs=[polynomial(pair[k]) for k in ('q1','q2','q3')]
    invariants=[]
    for q in qs:
        e,d,c,b,a=[q[i] for i in range(5)]
        invariants.append((12*a*e-3*b*d+c*c,72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c**3))
    assert invariants[0]==invariants[1]==invariants[2];I,J=invariants[0]
    T=polynomial(pair['T']);assert T==v**3-3*I*v+J
    A=R.quotient(T,'phi');phi=A.gen();m=A(polynomial(pair['m']))
    zs=[4*q[4]*phi/3+q[3]**2-8*q[4]*q[2]/3 for q in qs]
    assert m*m==zs[0]*zs[1]*zs[2] and m.is_unit()
    V=PolynomialRing(A,'z');z=V.gen();q=qs[0];e,d,c,b,a=[q[i] for i in range(5)]
    hdd=12*(3*b*b-8*a*c)*z*z+24*(b*c-6*a*d)*z+8*c*c-12*(b*d+8*a*e)
    ddq=12*a*z*z+6*b*z+2*c
    H=-8*phi*phi+4*ddq*phi+hdd+8*I
    HM=(m/zs[0])*H
    gamma_raw=R([x.lift()[2] for x in HM.list()])
    gamma=polynomial(pair['gam']);ratio=gamma_raw/gamma;assert ratio in QQ and ratio!=0
    aa=F(pair['a']);assert aa==F(str(qs[1][4])) and aa>0
    support=[int(p) for p in pari(pair['F'])];assert all(ZZ(p).is_prime(proof=True) for p in support)
    assert [int(p) for p,x,c in inp['local_GP']]==support
    evaluations=[];product=1
    for pp,xx,cc in inp['local_GP']:
        p=int(pp);x=F(xx);c=F(cc)
        assert evaluate(list(map(F,map(str,gamma.list()))),x)==c and c
        y2=evaluate(list(map(F,map(str,qs[0].list()))),x);assert y2
        val,unit=valuation_unit(y2,p)
        assert val%2==0 and (r.mod(unit,8)==1 if p==2 else pow(r.mod(unit,p),(p-1)//2,p)==1)
        symbol=hilbert(aa,c,p);product*=symbol
        assert int(pari.hilbert(QQ(str(aa)),QQ(str(c)),p))==symbol
        evaluations.append({'prime':p,'local_coordinate':str(x),'quartic_value':str(y2),'gamma_value':str(c),'Hilbert_symbol':symbol})
    assert product==-1 and inp['pairing_bit']==1
    return {'schema':'rank-jump.carrier-ct-local-witness.v1','status':'PASS',
            'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,PROTOCOL,Path(__file__),HERE/'retrospective.py',HERE/'cover_experiment.py')},
            'quartics':[list(map(str,q.list())) for q in qs],'invariants_I_J':list(map(str,invariants[0])),
            'cubic_square_root':list(map(str,m.lift().list())),'auxiliary_quadratic':list(map(str,gamma.list())),
            'auxiliary_quadratic_scalar':str(ratio),'a':str(aa),'local_evaluations':evaluations,
            'negative_symbol_primes':[x['prime'] for x in evaluations if x['Hilbert_symbol']==-1],
            'real_Hilbert_symbol':1,'Hilbert_product':product,
            'boundary':r.read(PROTOCOL)['boundary']}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','build','check']);a=p.parse_args()
    if a.mode=='capture':capture()
    else:
        out=compute()
        if a.mode=='build':r.write_new(OUTPUT,out)
        else:assert r.read(OUTPUT)==out
        print('PASS; negative Hilbert symbols at',out['negative_symbol_primes'])

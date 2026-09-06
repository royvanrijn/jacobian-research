#!/usr/bin/env python3
"""Exact trace-form geometry of fifty retained cubic 2-cover classes."""
import argparse
from pathlib import Path
import subprocess
import retrospective as r
import bad_prime_support as bad
import strict_class_blocks as strict
import small_quotient_covers as small
import scalar_cup

PROTOCOL=Path(__file__).with_name('TWO_COVER_PENCIL_GEOMETRY_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_two_cover_pencil_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_two_cover_pencil_geometry_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-two-cover-pencil-geometry-v1'


def mul(x,y,f):
    out=[r.F(0)]*5
    for i in range(3):
        for j in range(3):out[i+j]+=x[i]*y[j]
    for k in [4,3]:
        for j in range(3):out[k-3+j]-=out[k]*f[j]
    return out[:3]


def export():
    rows=[]
    for index in r.read(PROTOCOL)['production_cases']:
        old=r.read(bad.INPUT)['cases'][index];block=r.read(strict.OUTPUT)['rows'][index]
        source=bad.cases()[index]
        _,points=r.short(source['model'],source['generic_points']+source['points'])
        d=r.F(old['elliptic_scaling_d']);f=list(map(r.F,old['integral_cubic_ascending']))
        gammas=[];ys=[]
        for i in old['selected_input_indices']:
            x,y=map(r.F,points[i]);x*=d*d;y*=d**3
            assert y*y==x**3+f[2]*x*x+f[1]*x+f[0]
            gammas.append([x,r.F(-1),r.F(0)]);ys.append(y)
        classes=[]
        for kind,masks in [('generic',block['generic_strict_kernel_masks']),('relative',block['relative_strict_lift_masks'])]:
            for mask in masks:
                beta=[r.F(1),r.F(0),r.F(0)];root=r.F(1)
                for i,g in enumerate(gammas):
                    if mask>>i&1:beta=mul(beta,g,f);root*=ys[i]
                classes.append({'kind':kind,'witness_mask':mask,'beta':list(map(str,beta)),'norm_root':str(root)})
        rows.append({'case_index':index,'id':old['id'],'cubic_ascending':list(map(str,f)),
            'classes':classes,'all_bad_places_complete':block['all_bad_places_complete'],
            'generic_rank':old['generic_dimension'],'known_independent_rank':old['witness_dimension']})
    control=r.read(scalar_cup.CONTROL)
    rows.append({'case_index':6,'id':'small-rational-versus-Sha','cubic_ascending':control['polynomial_ascending'],
        'classes':[{'kind':'small_strict','beta_index':x['beta_index'],'beta':x['beta'],'norm_root':'25'}
                   for x in r.read(small.OUTPUT)['records']]})
    r.write_new(INPUT,{'schema':'rank-jump.two-cover-pencil-inputs.v1',
        'source_hashes':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (r.INPUT,bad.INPUT,strict.OUTPUT,small.OUTPUT,scalar_cup.CONTROL)},'rows':rows})


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),PROTOCOL,INPUT)}


def gram_forms(f,beta):
    from sage.all import matrix
    ring=f.base_ring();K=f.parent().quotient(f,'theta');th=K.gen()
    beta=K(beta);basis=[K(1),th,th**2]
    forms=[matrix(ring,3) for _ in range(3)]
    for i in range(3):
        for j in range(3):
            value=(beta*basis[i]*basis[j]).lift()
            for k in range(3):forms[k][i,j]=value[k]
    # Columns are coordinates of multiplication by theta+a.
    T=matrix(ring,3,3,lambda i,j:((th+f[2])*basis[j]).lift()[i])
    M=matrix(ring,3,3,lambda i,j:(beta*basis[j]).lift()[i])
    return forms,T,M


def compute(index):
    from sage.all import QQ,PolynomialRing,block_diagonal_matrix,matrix
    from sage.version import version
    row=next(x for x in r.read(INPUT)['rows'] if x['case_index']==index)
    R=PolynomialRing(QQ,'z');f=R(list(map(QQ,row['cubic_ascending'])))
    assert f.degree()==3 and f[3]==1 and f.discriminant()!=0
    P=PolynomialRing(QQ,['lam','mu']);lam,mu=P.gens()
    a,b,c=f[2],f[1],f[0]
    normalized=mu*(lam**3+2*a*lam**2*mu+(a*a+b)*lam*mu**2+(a*b-c)*mu**3)
    one,T,_=gram_forms(f,[1,0,0]);records=[]
    for item in row['classes']:
        forms,T1,M=gram_forms(f,list(map(QQ,item['beta'])));Q0,Q1,Q2=forms
        norm=M.det();assert norm==QQ(item['norm_root'])**2 and norm
        assert T1==T and Q2==one[2]*M and Q1==Q2*T and Q2.det()==-norm
        assert Q2.inverse()*Q1==T
        A4=block_diagonal_matrix(Q2,matrix(QQ,1,[0]))
        signs=[]
        for sign in r.read(PROTOCOL)['limits']['signs']:
            B4=block_diagonal_matrix(Q1,matrix(QQ,1,[sign]))
            determinant=(lam*A4.change_ring(P)+mu*B4.change_ring(P)).det()
            assert -determinant/(sign*norm)==normalized
            signs.append({'sigma':sign,'pencil_determinant_coefficients':
                [str(determinant.monomial_coefficient(lam**(4-i)*mu**i)) for i in range(5)],
                'four_geometric_simple_pencil_roots':True,'smooth_genus_one':True})
        records.append({**item,'norm':str(norm),'Q0_Gram':[list(map(str,x)) for x in Q0.rows()],
            'Q1_Gram':[list(map(str,x)) for x in Q1.rows()],
            'Q2_Gram':[list(map(str,x)) for x in Q2.rows()],
            'first_conic_determinant':str(Q2.det()),'signs':signs})
    return {'bindings':bindings(),'case_index':index,'id':row['id'],'status':'PASS','software':{'sage':version},
        'cubic_discriminant':str(f.discriminant()),'common_operator':[list(map(str,x)) for x in T.rows()],
        'normalized_pencil_coefficients':[str(normalized.monomial_coefficient(lam**(4-i)*mu**i)) for i in range(5)],
        'class_count':len(records),'classes':records,
        'boundary':'Smoothness and common determinant are automatic. No new solubility, full-Selmer or rank assertion.'}


def capture(check=False):
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for i in range(7):
        if check:
            assert compute(i)==next(x for x in r.read(OUTPUT)['rows'] if x['case_index']==i)
            print('PASS pencil geometry',i,flush=True);continue
        path=WORK/f'case-{i}.json'
        if not path.exists():
            with (WORK/f'case-{i}.log').open('x') as log:
                try:
                    p=subprocess.run(['sage','-python',str(Path(__file__).resolve()),'worker','--index',str(i),
                        '--destination',str(path)],cwd=r.ROOT,stdout=log,stderr=log,timeout=30)
                    reason=None if p.returncode==0 else 'worker failure'
                except subprocess.TimeoutExpired:reason='30-second timeout'
                if reason and not path.exists():r.write_new(path,{'bindings':bindings(),'case_index':i,'status':'UNKNOWN','reason':reason})
        row=r.read(path);assert row['bindings']==bindings();rows.append(row)
        print('checkpoint',i,row['status'],row.get('class_count'),flush=True)
    if not check:r.write_new(OUTPUT,{'schema':'rank-jump.two-cover-pencil-geometry.v1','bindings':bindings(),'rows':rows})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','capture','worker','check'])
    p.add_argument('--index',type=int);p.add_argument('--destination',type=Path);args=p.parse_args()
    if args.mode=='export':export()
    elif args.mode=='worker':r.write_new(args.destination,compute(args.index))
    else:capture(args.mode=='check')

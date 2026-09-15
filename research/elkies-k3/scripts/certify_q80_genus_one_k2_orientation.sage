#!/usr/bin/env sage
"""Finite hypotheses and exact controls for the written Q80 k2 orientation proof."""
from sage.all import GF, PolynomialRing, QQ, Qp, ZZ
from pathlib import Path
from hashlib import sha256
import argparse
import json
import resource
import time

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'artifacts/generated-results/elkies-k3-q80-genus-one-k2-orientation-v1'
SOURCE='artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
BOUNDARY='artifacts/generated-results/elkies-k3-q80-genus-one-two-disagreement-v1'
DEGREE='artifacts/generated-results/elkies-k3-q80-degree-two-reciprocity-v1'
CHECKER='elkies-k3/scripts/verify_q80_genus_one_k2_orientation.py'
TEST='tests/test_q80_genus_one_k2_orientation.py'
P=131
def read(path):return json.loads(path.read_text())
def digest(path):return sha256(path.read_bytes()).hexdigest()
def write(name,data):
    with (OUT/name).open('x') as f:json.dump(data,f,indent=2,sort_keys=True);f.write('\n')
def files():
    return [SOURCE,CHECKER,TEST,str(Path(__file__).relative_to(ROOT)),
            'elkies-k3/scripts/verify_q80_genus_one_k2_nodal_boundary.py',
            *[BOUNDARY+'/'+n for n in ['input.json','result.json','character-gate.json','independent-replay.json']],
            *[DEGREE+'/'+n for n in ['input.json','result.json','quadratic-fibres.json','independent-replay.json']]]
def freeze():
    packet={'schema':1,'prime':P,'cpu_seconds':20,'memory_bytes':2*1024**3,
            'control_z':[38,13793,30954],'hensel_unit_digits':16,
            'bindings':{s:digest(ROOT/s) for s in files()},
            'preserved_preflight':{str(f.relative_to(ROOT)):digest(f) for f in (OUT/'preflight').rglob('*') if f.is_file()},
            'scope':'Finite nodal hypotheses and local controls only; the all-denominator Q80 genus1 k2 exclusion is the written analytic proof.'}
    write('input.json',packet);print(json.dumps({'status':'FROZEN','input_sha256':digest(OUT/'input.json')}),flush=True)
def run():
    resource.setrlimit(resource.RLIMIT_CPU,(20,25))
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    started=time.process_time();wall=time.monotonic();packet=read(OUT/'input.json')
    assert set(packet['bindings'])==set(files())
    for s,h in {**packet['bindings'],**packet['preserved_preflight']}.items():assert digest(ROOT/s)==h
    source=read(ROOT/SOURCE);F=GF(P);R=PolynomialRing(F,'t');t=R.gen()
    coeff=lambda key:R([F(QQ(v)) for v in source['weierstrass_model'][key+'_coefficients_low_to_high']])
    A,B=coeff('A'),coeff('B');q=t*t+62*t+88;assert q.is_irreducible()
    K=GF(P*P,name='theta',modulus=q);theta=K.gen()
    encode=lambda e:[int(e[i]) for i in range(2)]
    rho=73+6*theta;node=29+128*theta
    aa,bb=A(theta),B(theta);ap,bp=A.derivative()(theta),B.derivative()(theta)
    derivative=3*rho*rho+aa
    rho_prime=-(ap*rho+bp)/derivative
    vprime=-ap-3*rho*rho_prime/2
    delta=4*A**3+27*B**2
    assert not delta(theta) and delta.derivative()(theta)
    assert rho==-2*node and rho and derivative and vprime
    assert not rho**3+aa*rho+bb and not 3*node*node+aa
    assert not (-aa-3*rho*rho/4)
    assert -delta.derivative()(theta)==4*vprime*derivative**2
    nodal={'q':[88,62,1],'A_mod_p':[int(v) for v in A.list()],'B_mod_p':[int(v) for v in B.list()],
           'node':encode(node),'rho':encode(rho),'rho_prime':encode(rho_prime),
           'simple_root_derivative':encode(derivative),'V_prime':encode(vprime),
           'delta_prime':encode(delta.derivative()(theta))}
    codes=[]
    for root in [node,rho]:
        code=0
        for j,section in enumerate(source['sections']['records']):
            x=section['X']
            num=R([F(QQ(v)) for v in x['numerator_coefficients_low_to_high']])
            den=R([F(QQ(v)) for v in x['denominator_coefficients_low_to_high']])
            assert den(theta) # All literal basis abscissas are regular here.
            value=num(theta)/den(theta)-root
            assert value # No regularization is required at either selected root.
            norm=value**(P+1);assert norm in F
            code |= int(not F(norm).is_square())<<j
        codes.append(code)
    assert codes==[13412,0];nodal['norm_codes_node_simple']=codes
    old=read(ROOT/BOUNDARY/'character-gate.json')
    assert old['node_root']==[29,128] and old['simple_root']==[73,6]
    assert old['quadratic_zero_rows'][0]['norm_codes']==[13412,13412,0]
    write('nodal-hypotheses.json',nodal)
    # Identities are polynomial identities, with no evaluation at sampled parameters.
    M=PolynomialRing(QQ,names=('r','a','x'));r,a,x=M.gens();b=-r**3-a*r;V=-a-3*r*r/4
    assert x**3+a*x+b==(x-r)*((x+r/2)**2-V)
    assert -4*a**3-27*b*b==4*V*(3*r*r+a)**2
    Z=PolynomialRing(QQ,'z');z=Z.gen();T=PolynomialRing(Z,'t');t=T.gen()
    c=z*z-3;s=4*z*c;beta=c*c;gamma=-3*c*(z*z+1)
    Q=(t-beta)**2+4*s*z*(t-beta)+6*s*s*(z*z-1)
    Hnum=((t-gamma)**2-3*s*s)**2-s**4*t
    assert Hnum==(t-beta)**2*Q
    assert Q[1]**2-4*Q[0]==-8*s*s*c
    write('symbolic-identities.json',{'cubic_factorization':True,'nodal_discriminant_identity':True,
          'toy_factorization':True,'toy_discriminant_identity':True,'checked_as_polynomial_identities':True,
          'variable_order':['t','z'],'toy_Q_coefficients':[[str(v) for v in Q[k].list()] for k in range(3)]})
    controls=[]
    Rq=PolynomialRing(QQ,'t');t=Rq.gen()
    for z in packet['control_z']:
        z=ZZ(z);c=z*z-3;s=4*z*c;beta=c*c;gamma=-3*c*(z*z+1)
        Q=(t-beta)**2+4*s*z*(t-beta)+6*s*s*(z*z-1)
        X=-2+(t-gamma)**2/s**2;Y=(t-gamma)*(t-beta)/s**3
        assert (X+2)*((X-1)**2-t)==Q*Y**2
        disc=Q[1]**2-4*Q[0];vd=int(disc.valuation(P));unit=disc/P**vd
        square=bool(GF(P)(unit).is_square())
        row={'z':int(z),'contact_valuation':int(c.valuation(P)),'Q':[str(v) for v in Q.list()],
             'disc_valuation':vd,'disc_unit_mod_p':int(GF(P)(unit)),'disc_unit_square_mod_p':square,
             'branch_field_over_Qp':'ramified_quadratic' if vd%2 else ('split' if square else 'unramified_quadratic')}
        if not vd%2 and square:
            digits=packet['hensel_unit_digits'];pp=Qp(P,digits+5);root=pp(unit).sqrt()
            seeds=GF(P)(unit).sqrt(all=True);seed=min(int(v) for v in seeds)
            if int(root.residue())!=seed:root=-root
            root=ZZ(root.lift())%(P**digits)
            assert (root*root-unit)%(P**digits)==0
            centers=[(-Q[1]+sign*P**(vd//2)*root)/2 for sign in [1,-1]]
            radius=min(int(Q(b).valuation(P)-Q.derivative()(b).valuation(P)) for b in centers)
            xv=[X(b) for b in centers]
            values={'beta':int(beta.valuation(P)),'gamma':int(gamma.valuation(P)),
                    'branch1':int(centers[0].valuation(P)),'branch2':int(centers[1].valuation(P)),
                    'branch_separation':int((centers[0]-centers[1]).valuation(P)),
                    'x_difference':int((xv[0]-xv[1]).valuation(P)),
                    'node_offsets_sum':int((xv[0]+xv[1]-2).valuation(P)),
                    'node_offset1':int((xv[0]-1).valuation(P))}
            assert values=={'beta':4,'gamma':2,'branch1':2,'branch2':2,'branch_separation':3,
                            'x_difference':1,'node_offsets_sum':2,'node_offset1':1}
            assert radius>=19
            row.update({'sqrt_unit_mod_p_digits':str(root),'sqrt_unit_seed':seed,
                        'branch_centers':[str(v) for v in centers],'branch_error_valuation_lower_bound':radius,
                        'x_error_valuation_lower_bound':radius-2,'valuations':values})
        controls.append(row)
    assert [r['branch_field_over_Qp'] for r in controls]==['ramified_quadratic','unramified_quadratic','split']
    write('local-controls.json',{'rows':controls})
    result={'status':'PASS','input_sha256':digest(OUT/'input.json'),
            'records':{s:digest(OUT/s) for s in ['nodal-hypotheses.json','symbolic-identities.json','local-controls.json']},
            'finite_replay_scope':'Nodal units/derivatives, two inherited norm codes, four symbolic identities, three exact local controls.',
            'written_theorem':'All three polynomial Q80 genus1 k2 allocations excluded over Q_131; analytic proof required.',
            'genus1_allocations_remaining':9,'positive_mw17_target_complete':False,
            'formal_verification':False,'old_norm8_replay_upgraded':False}
    write('result.json',result)
    write('execution.json',{'status':'PASS','input_sha256':digest(OUT/'input.json'),'result_sha256':digest(OUT/'result.json'),
                           'cpu_seconds':time.process_time()-started,'wall_seconds':time.monotonic()-wall,
                           'cpu_limit':20,'memory_limit_bytes':2*1024**3})
    print(json.dumps(result,sort_keys=True),flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['freeze','run']);args=parser.parse_args()
    freeze() if args.action=='freeze' else run()

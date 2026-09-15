#!/usr/bin/env sage
"""Exact finite contacts and p-adic unit controls for the Q80 mixed-valuation proof."""
from sage.all import GF, PolynomialRing, QQ, Qp, EllipticCurve, matrix, vector
from hashlib import sha256
from pathlib import Path
import argparse,json,resource,time
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'artifacts/generated-results/elkies-k3-q80-genus-one-k0-orthogonal-contact-v1'
SOURCE='artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
NORM4='artifacts/generated-results/elkies-k3-q80-genus-zero-integral-reduction-v1'
NODE='artifacts/generated-results/elkies-k3-q80-genus-one-k2-orientation-v1'
POLE='artifacts/generated-results/elkies-k3-q80-genus-one-k0-pole-boundary-v1'
CHECKER='elkies-k3/scripts/verify_q80_k0_mixed_exclusion.py'
HELPER='elkies-k3/scripts/verify_q80_genus_one_k2_nodal_boundary.py'
P=131
def read(p):return json.loads(p.read_text())
def digest(p):return sha256(p.read_bytes()).hexdigest()
def write(n,o):
    with (OUT/n).open('x') as f:json.dump(o,f,indent=2,sort_keys=True);f.write('\n')
def paths():
    return [SOURCE,CHECKER,HELPER,str(Path(__file__).relative_to(ROOT)),
      *[NORM4+'/'+n for n in ['input.json','result.json','norm4-sections.json','complete-replay-input.json','complete-independent-replay.json']],
      *[NODE+'/'+n for n in ['input.json','result.json','nodal-hypotheses.json','independent-replay.json']],
      *[POLE+'/'+n for n in ['input.json','receipt.json']],
      'elkies-k3/Q80_GENUS_ONE_K0_POLE_BOUNDARY_2026-09-14.md']
def freeze():
    write('input.json',{'schema':1,'prime':P,'cpu_seconds':30,'memory_bytes':2*1024**3,'producer_padic_precision':20,
      'scope':'All131 constant-twist contact interpolants; a rank13 formal-deformation Jacobian; the exact norm4 word at the lifted nodal fibre modulo131 squared; a nonrational necessary leading scalar. The exclusion needs the written analytic deformation proof.',
      'bindings':{p:digest(ROOT/p) for p in paths()},
      'preserved_preflight':{str(p.relative_to(ROOT)):digest(p) for p in (OUT/'preflight').glob('*') if p.is_file()}})
def run():
    started=time.process_time();packet=read(OUT/'input.json');assert packet['bindings']=={p:digest(ROOT/p) for p in paths()}
    for p,h in packet['preserved_preflight'].items():assert digest(ROOT/p)==h
    source=read(ROOT/SOURCE);model=source['weierstrass_model'];point=read(ROOT/NORM4/'norm4-sections.json')['records'][732]
    F=GF(P);R=PolynomialRing(F,'t');t=R.gen();poly=lambda a:R([F(QQ(c)) for c in a]);coeff=lambda f:[int(c) for c in f.list()]
    A=poly(model['A_coefficients_low_to_high']);B=poly(model['B_coefficients_low_to_high']);q=t*t+62*t+88;D=q*q;S=73+6*t
    x0=(S-(S**3+A*S+B)*(3*S*S+A).inverse_mod(D))%D;rows=[]
    for c in F:
        x=x0+c*D;f,rem=(x**3+A*x+B).quo_rem(D);assert not rem
        lead=f.leading_coefficient();monic=f/lead
        if monic.is_square():rows.append({'c':int(c),'x':coeff(x),'r_monic':coeff(monic.sqrt()),'scalar':int(lead),'scalar_square':bool(lead.is_square())})
    assert len(rows)==1 and rows[0]['c']==82 and rows[0]['x']==point['x'] and rows[0]['scalar']==9
    write('constant-contacts.json',{'x0':coeff(x0),'q0':coeff(q),'tried':131,'rows':rows})
    XT,YT=R(point['x']),R(point['y']);r,rem=YT.quo_rem(q);assert not rem
    E=EllipticCurve(R.fraction_field(),[A,B]);T=E(0)
    for i,w in enumerate(point['word']):
        if not w:continue
        row=source['sections']['records'][i]
        xx=poly(row['X']['numerator_coefficients_low_to_high'])/poly(row['X']['denominator_coefficients_low_to_high'])
        yy=poly(row['Y']['numerator_coefficients_low_to_high'])/poly(row['Y']['denominator_coefficients_low_to_high'])
        T+=int(w)*E(xx,yy)
    assert T==E(XT,YT)
    w=point['word'];G=source['sections']['height_gram'];assert sum(w[i]*QQ(G[i][j])*w[j] for i in range(17) for j in range(17))==4
    fp=3*XT*XT+A;J=matrix(F,13,14)
    columns=[fp*t**k for k in range(5)]+[-2*D*r*t**k for k in range(5)]+[-r*r*t**k for k in range(4)]
    for j,f in enumerate(columns):
        for k in range(13):J[k,j]=f[k]
    minor=J.matrix_from_columns([i for i in range(14) if i!=4]);dd=(fp*r.inverse_mod(D))%D;dr=(fp-r*dd)//(2*D)
    kernel=vector(F,[r[k] for k in range(5)]+[dr[k] for k in range(5)]+[dd[k] for k in range(4)])
    assert J.rank()==13 and minor.det()!=0 and J*kernel==0 and r[4]==3
    write('formal-deformation.json',{'point_index':732,'word':w,'word_height':4,'x':coeff(XT),'y':coeff(YT),'r':coeff(r),'q0':coeff(q),
        'jacobian_rank':13,'fixed_x4_minor_determinant':int(minor.det()),'normalized_tangent':{'dx':coeff(r),'dr':coeff(dr),'dD':coeff(dd)},
        'r_mod_q0':coeff(r%q),'fprime_mod_q0':coeff(fp%q),'dD_mod_q0':coeff(dd%q),'word_identity_over_F131_t':True})
    base=Qp(P,prec=20);Z=PolynomialRing(base,'z');z=Z.gen();K=base.extension(z*z+62*z+88,names='theta');theta=K.gen();Q=PolynomialRing(K,'t')
    pk=lambda a:Q([K(QQ(c)) for c in a]);AK=pk(model['A_coefficients_low_to_high']);BK=pk(model['B_coefficients_low_to_high']);delta=4*AK**3+27*BK**2;alpha=theta
    for _ in range(6):alpha-=delta(alpha)/delta.derivative()(alpha)
    e=-3*BK(alpha)/(2*AK(alpha));s=(3*e).sqrt();torus=K(1);basis=[]
    def mod2(a):
        assert a.valuation()>=0 and a.precision_absolute()>=3
        f=a.polynomial();return [int(f[i].lift())%(P*P) for i in range(2)]
    for i,wi in enumerate(w):
        if not wi:continue
        row=source['sections']['records'][i]
        x=pk(row['X']['numerator_coefficients_low_to_high'])(alpha)/pk(row['X']['denominator_coefficients_low_to_high'])(alpha)
        y=pk(row['Y']['numerator_coefficients_low_to_high'])(alpha)/pk(row['Y']['denominator_coefficients_low_to_high'])(alpha)
        torus*=((y+s*(x-e))/(y-s*(x-e)))**int(wi);basis.append({'index':i,'word':int(wi),'x_mod_p2':mod2(x),'y_mod_p2':mod2(y)})
    u=s*(torus+1)/(torus-1);yt=u*(u*u-3*e);assert u.valuation()==yt.valuation()==1
    up=mod2(u);yp=mod2(yt);assert all(c%P==0 for c in up+yp)
    residue=GF(P**2,modulus=q,names='tau');tau=residue.gen();enc=lambda a:[int(a.polynomial()[i]) for i in range(2)]
    Yover=residue(yp[0]//P)+residue(yp[1]//P)*tau;Yprime=residue(YT.derivative()(tau));g=-Yover/Yprime;qp=2*tau+62
    required=12*(29+128*tau)/(qp*qp*g*g)
    assert enc(g)==[91,113] and enc(required)==[69,40] and required**P!=required
    write('nodal-contact.json',{'precision':20,'alpha_mod_p2':mod2(alpha),'node_mod_p2':mod2(e),'basis_rows':basis,
        'T_u_mod_p2':up,'T_y_mod_p2':yp,'T_y_valuation':1,'T_x_minus_simple_valuation':2,
        'T_y_over_p_mod_q0':enc(Yover),'T_yprime_mod_q0':enc(Yprime),'contact_shift_over_p':enc(g),
        'necessary_negative_abscissa_valuation':-2,'required_primitive_scalar':enc(required),'required_scalar_in_F131':False})
    write('result.json',{'status':'PASS','input_sha256':digest(OUT/'input.json'),
        'records':{n:digest(OUT/n) for n in ['constant-contacts.json','formal-deformation.json','nodal-contact.json']},
        'constant_twist_contacts':1,'nonsquare_contacts':0,'mixed_valuation_excluded_by_written_proof':True,
        'remaining_genus1_k0_j':[2],'remaining_valuation_pattern':'both negative','remaining_boundary':'UNKNOWN',
        'written_deformation_and_valuation_argument_required':True,'new_analytic_argument_independently_verified':False,
        'formal_verification':False,'external_review':False,'positive_mw17_target_complete':False})
    write('execution.json',{'status':'PASS','cpu_seconds':time.process_time()-started,'cpu_limit':30,'memory_bytes':2*1024**3})
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['freeze','run']);args=p.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(30,35));resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    freeze() if args.mode=='freeze' else run()

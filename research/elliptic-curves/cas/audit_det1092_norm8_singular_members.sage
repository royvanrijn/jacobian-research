#!/usr/bin/env sage-python
"""One existing generic pencil: exact singular-member obstruction, 25s cap.

No member/point sweep, exceptional point, new CVP or polynomial factorization.
Extract repeated support by gcd; reuse the ten saved vertical section words.
"""
import hashlib,json,signal,time,traceback
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector,prime_range
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_norm8_singular_members_v1'
GEN=ART/'det1092_norm8_seed_cover_v2/generic.json'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
FRAME=ART/'det1092_genus1_picard_image_v1/frame.json'
PRIMES=list(map(int,prime_range(3,212)))
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    p=OUT/name
    if p.exists():assert read(p)==data
    else:
        with p.open('x') as stream:json.dump(data,stream,indent=2,sort_keys=True);stream.write('\n')
def coeff(f):return list(map(str,f.list()))
def rec(f):return {'numerator':coeff(f.numerator()),'denominator':coeff(f.denominator())}
def main():
    begun=time.monotonic();generic=read(GEN);parent=read(PARENT);frame=read(FRAME)
    R=PolynomialRing(QQ,'z');z=R.gen();K=R.fraction_field()
    A=R(generic['Jacobian_A']);B=R(generic['Jacobian_B'])
    delta=-16*(4*A**3+27*B**2)
    repeated=delta.gcd(delta.derivative()).monic()
    nodal,remainder=delta.quo_rem(repeated**2);assert not remainder
    scalar=nodal.leading_coefficient();nodal=nodal.monic()
    assert A.degree()==8 and B.degree()==12
    assert repeated.is_squarefree() and nodal.is_squarefree()
    assert repeated.gcd(nodal)==1 and A.gcd(delta)==1
    save('geometry.json',{'status':'PASS_EXACT_DISCRIMINANT_DECOMPOSITION',
        'A':coeff(A),'B':coeff(B),'delta_scalar':str(scalar),
        'repeated_support':coeff(repeated),'simple_support':coeff(nodal),
        'finite_discriminant_degree':int(delta.degree()),
        'infinity_discriminant_order':24-int(delta.degree()),
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [GEN,PARENT,FRAME]}})
    print('CHECKPOINT_SINGULAR_SUPPORT',repeated.degree(),nodal.degree(),
          'infinity',24-delta.degree(),flush=True)
    attempts=[];witness=None
    for p in PRIMES:
        if any(c.denominator()%p==0 for c in nodal):
            attempts.append({'p':p,'status':'denominator_bad'});continue
        Rp=PolynomialRing(GF(p),'z');f=Rp(nodal)
        residues=[int(f(a)) for a in GF(p)]
        attempts.append({'p':p,'status':'checked','rational_residue_roots':residues.count(0)})
        if 0 not in residues:
            witness={'p':p,'polynomial':list(map(int,f.list())),
                     'value_residues':residues};break
    save('no-rational-nodal-member.json',{'status':'PASS_NO_RATIONAL_NODAL_PARAMETER' if witness else 'UNKNOWN',
        'prime_pool':PRIMES,'attempts':attempts,'witness':witness,
        'geometry_sha256':sha(OUT/'geometry.json')})
    assert witness is not None
    T=PolynomialRing(QQ,'t');t=T.gen();F=T.fraction_field()
    def dec(d):return F(T(d['numerator']))/T(d['denominator'])
    E=EllipticCurve(F,[dec(d) for d in parent['a_invariants']])
    basis=[E([dec(d) for d in P]) for P in parent['basis_weierstrass_coordinates']]
    w=vector(ZZ,[0]*14+[1,-1,0]);G=matrix(ZZ,parent['generic_height_gram'])
    center=basis[14]-basis[15]
    cx=center[0]+E.b2()/12;cy=center[1]+(E.a1()*center[0]+E.a3())/2
    h=T(generic['pole_h']);shift=T(generic['shift'])
    sections=[]
    for row in frame['vertical_old_sections']:
        word=vector(ZZ,row['word']);P=sum((a*Q for a,Q in zip(word,basis)),E(0))
        assert word*G*word==word*G*w
        if P.is_zero() or P==center:
            sections.append({'word':list(map(int,word)),'parameter':'infinity'});continue
        xx=P[0]+E.b2()/12;yy=P[1]+(E.a1()*P[0]+E.a3())/2
        m=(yy+cy)/(xx-cx);zz=(h*m+shift)/(h*h)
        assert zz.numerator().degree()<=0 and zz.denominator().degree()==0
        address=QQ(zz);WW=(2*xx+cx-m*m)/h
        branch=T([R(row)(address) for row in generic['quartic_t_coefficients_in_z']])
        assert WW*WW==branch and WW.denominator().degree()==0
        sections.append({'word':list(map(int,word)),'parameter':str(address),
                         'W':coeff(T(WW))})
    pairs=[]
    for address in sorted(set(row['parameter'] for row in sections)):
        rows=[row for row in sections if row['parameter']==address]
        assert len(rows)==2
        v1,v2=[vector(ZZ,row['word']) for row in rows]
        assert v1+v2==w and v1*G*v1+v2*G*v2==8
        if address!='infinity':
            assert repeated(QQ(address))==0 and T(rows[0]['W'])==-T(rows[1]['W'])
        pairs.append({'parameter':address,'components':rows})
    finite=[QQ(row['parameter']) for row in pairs if row['parameter']!='infinity']
    product=R.one()
    for address in finite:product*=z-address
    assert product==repeated and len(pairs)==5 and repeated.degree()==4 and nodal.degree()==14
    save('rational-singular-members.json',{'status':'PASS_ALL_RATIONAL_SINGULAR_MEMBERS_INHERITED',
        'classification':'generic-only exact obstruction for the fixed norm8 pencil',
        'members':pairs,'Kodaira_configuration':{'I2':5,'I1':14},
        'generic_geometric_MW_rank':12,
        'conclusion':'All rational singular parameters are the five I2 fibres, each the sum of two original generic sections. The14 irreducible nodal fibres have no rational parameter. No rational-normalization seed cover arises from singular members of this pencil over Q.',
        'boundary':'Does not exclude smooth positive-genus members, other pencils, higher arithmetic-genus rational curves, or seed incidence on302.',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [GEN,PARENT,FRAME,OUT/'geometry.json',OUT/'no-rational-nodal-member.json']}})
    print('PASS_ALL_RATIONAL_SINGULAR_MEMBERS_INHERITED','5 I2,14 I1',
          'no-root prime',witness['p'],'seconds',round(time.monotonic()-begun,3),flush=True)
if __name__=='__main__':
    OUT.mkdir(exist_ok=True)
    save('protocol.json',{'classification':'one generic pencil singular-member theorem test',
        'rule':'Use the already constructed norm8 orbit20124 pencil and ten saved vertical words only. Compute one discriminant and repeated support by polynomial gcd. Test the first prime in the fixed pool certifying no rational roots of the simple support. Reconstruct all saved vertical components, including infinity. No substitution of another pencil on failure.',
        'limits':{'seconds':25,'pencils':1,'new_member_searches':0,'point_searches':0,'CVPs':0,
                  'polynomial_factorizations':0,'integer_factorizations':0,'exceptional_point_inputs':0},
        'prime_pool':PRIMES,'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [GEN,PARENT,FRAME,Path(__file__)]}})
    def expired(signum,frame):raise TimeoutError('FROZEN_25_SECOND_LIMIT')
    signal.signal(signal.SIGALRM,expired);signal.alarm(25)
    try:main()
    except BaseException as exc:
        save('failure.json',{'error':str(exc),'traceback':traceback.format_exc()});raise
    finally:signal.alarm(0)

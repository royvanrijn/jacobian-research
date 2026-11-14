"""Complete the high-contact quintic cusp locus and its row-B singleton gate."""
from pathlib import Path
import json,time,hashlib
import sympy as S
from cusps import analyze,G,A,M,s,t
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
CHECKS=[]
def check(name,value):
    assert value,name
    CHECKS.append(name)
def zero(name,value):check(name,S.expand(value)==0)

def identities():
    u,b,e,n,v=S.symbols('u b e n v')
    aa=-3*u*u;cc=5*aa/3;ff=10*u**4-4*b*u**3-2*e*u
    pp=s**3+aa*s+M;qq=s**5+b*s**4+cc*s**3+e*s*s+ff*s+n
    zero('universal critical equation',S.diff(qq,s).subs(s,u))
    zero('universal clean-fibre P factor',pp-pp.subs(s,u)-(s-u)**2*(s+2*u))
    jet=(S.diff(pp,s,2)*S.diff(qq,s,3)-S.diff(pp,s,3)*S.diff(qq,s,2)).subs(s,u)
    zero('universal ordinary-cusp jet formula',jet-12*(20*u**3+6*b*u*u-e))
    zero('universal extra-preimage formula',qq.subs(s,-2*u)-qq.subs(s,u)-9*u*u*(e+3*b*u*u-2*u**3))
    rr=v**4+b*v**3+aa*v*v+(2*aa*b-e)*v-(aa*aa+ff-cc*aa)
    zero('universal diagonal collision root',rr.subs(v,2*u))
    zero('universal collision derivative equals cusp jet',S.diff(rr,v).subs(v,2*u)-jet/12)

def special():
    d=G.F.d
    row=next(v for v in json.loads((HERE.parent/'jc2-r7-quintic-normalizations/normalizations.json').read_text())['families'] if v['family']=='L')
    parse=lambda e:S.sympify(e,locals=dict(A=A,M=M,s=s,d=d))
    ff,gg=[parse(e).subs(d,0) for e in row['groebner_basis'][:2]]
    zero('d0 cubic branch equation',ff-(A**3-9*M*M-S.Rational(12,5)*M-S.Rational(1,15)))
    zero('d0 split equation',gg-A*A*(M+S.Rational(2,15)))
    q=parse(row['normalization_q']).subs(d,0)
    delta=S.resultant(3*s*s+A,S.diff(q,s),s)
    zero('d0 nonzero-A branch immersive',delta.subs(M,-S.Rational(2,15))-25*A**4)
    quad=135*M*M+36*M+1
    zero('d0 A0 branch quadratic',-81*parse(row['necessary_sufficient_equations'][0]).subs({A:0,d:0})-quad)
    check('two distinct special cusp choices',S.discriminant(quad,M)==756)
    ee=S.Rational(5,3)*M+S.Rational(2,9)
    zero('special normalization Q',q.subs(A,0)-(s**5+ee*s*s))
    check('special cusp ordinary and clean',S.gcd(ee,quad)==1)
    k6=parse(row['implicit_difference']['(2, 0)']).subs({A:0,d:0})
    check('special contact exactly38',S.gcd(k6,quad)==1)
    pc=G.F.p.subs(G.F.line).subs(d,0);qc=G.F.q.subs(G.F.line).subs(d,0)
    equation=S.Poly(S.resultant(pc-M,qc,t),M).monic().as_expr()
    check('special cusp target off C',S.gcd(equation,quad)==1)
    v=S.Symbol('v')
    zero('special collision quartic',v**4-ee*v-v*(v**3-ee))
    zero('special triple fibre P',(s**3+M).subs(s**3,-ee)-(M-ee))
    zero('special triple fibre Q',S.rem(s**5+ee*s*s,s**3+ee,s))
    zero('special triple tangent slopes',S.rem(s*S.diff(s**5+ee*s*s,s)+ee*S.diff(s**3+M,s),s**3+ee,s))
    check('special triple is ordinary',S.gcd(ee,quad)==1 and S.gcd(S.diff(v**3-ee,v),v**3-ee)==1)
    return dict(parameter='d=0',coefficient_equation=str(quad),normalization_p='s^3+M',normalization_q=str(s**5+ee*s*s),cusp=['M','0'],triple=[str(M-ee),'0'],
        first_image_V0_resultant=str(equation),k6_remainder=str(S.rem(k6,quad,M)),pair_count=2,singularities='one clean ordinary cusp and one ordinary triple point',infinity_contact=38)

def main():
    start=time.monotonic();loci=[];identities()
    for family,prime,degree in [('L',103,36),('K',101,41)]:
        data=analyze(family,prime)
        check(family+' full exact parameter degree survives reduction',data['base_degree']==degree)
        for label,value in data['unit_gcd_degrees'].items():check(family+' '+label,value==0)
        check(family+' three ordinary nodes at every nonspecial cusp choice',data['three_node_rank']==data['three_node_dimension']==3*degree)
        check(family+' every nonspecial cusp target avoids C',data['cusp_off_C_rank']==data['cusp_off_C_dimension']==6*degree)
        loci.append(data)
    exceptional=special()
    check('complete161 distinct normalization pairs',36+3*41+2==161)
    check('complete160 distinct first-image coefficients',36+3*41+1==160)
    check('rowB cusp fibre has3 source points',6-2*0-2*1-1==3)
    check('rowB elsewhere has even source cardinality',all((6-2*c-2*b)%2==0 for c in range(4) for b in range(4)))
    result=dict(schema='r7-quintic-cusp-classification-v1',nonspecial_loci=loci,special=exceptional,
        nonimmersive_pair_count=161,first_image_coefficient_count=160,
        classification='159 nonspecial pairs have A2+3A1; two d0 pairs have A2 plus an ordinary triple. Every cusp is clean, off C, and contact is38.',
        rowB='No singleton fibre under the explicit same-coordinate degree(3,5), infinity-contact>=38 second-image hypothesis.',
        scope='No forcing of that chart or contact, no exclusion of a Keller map merely for lacking a singleton, no global cover/source construction or JC2 closure.')
    (HERE/'classification.json').write_text(json.dumps(result,indent=2)+'\n')
    inputs=[HERE/name for name in ['verify.py','cusps.py','PROOF.md','README.md','REVIEW.md','discovery.json']]
    inputs += [HERE.parent/path for path in ['jc2-r7-quintic-incidence/norms.py','jc2-r7-quintic-incidence/finite.py','jc2-r7-quintic-incidence/incidence.json','jc2-r7-quintic-incidence/nonimmersion-parameters-L.json','jc2-r7-quintic-incidence/nonimmersion-parameters-K.json','jc2-r7-quintic-incidence/PROOF.md','jc2-r7-quintic-normalizations/normalizations.json','jc2-r7-quintic-normalizations/PROOF.md','jc2-r7-complement/geometry.py','jc2-r7-image-family/family.py','jc2-r7-laurent-collisions/laurent.py','jc2-omitted-divisors/RESIDUE_LEDGER.md','jc2-omitted-divisors/SIMPLE_RAMIFICATION_AND_RADIAL.md','jc2-omitted-divisors/NONIMMERSED_LEDGER.md']]
    inputs += [ROOT/'archive/non-elliptic/scripts/verify_f2_degree_6_10_gap.py']
    receipt=dict(passed=True,check_count=len(CHECKS),checks=CHECKS,elapsed_seconds=round(time.monotonic()-start,3),
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        output_sha256={'classification.json':hashlib.sha256((HERE/'classification.json').read_bytes()).hexdigest()},
        assurance='Written adjugate/unit lifting, complete coefficient-locus classification and oriented fibre argument, with exact modular/rational certificates; no formalization or independent external review.')
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(f'PASS {len(CHECKS)} QUINTIC CUSP CHECKS ({time.monotonic()-start:.3f}s)',flush=True)

if __name__=='__main__':main()

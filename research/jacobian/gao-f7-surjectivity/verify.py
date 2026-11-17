"""Exact F7 surjectivity certificate replay; no Groebner computation.

SymPy rederives the necessary coefficient equations. Independent FLINT
arithmetic replays the two exported Singular ideal identities.
"""
import hashlib
import json
from pathlib import Path
import platform
import sympy as S
from flint import fmpq_mpoly_ctx, fmpq

if not __debug__:
    raise RuntimeError('Run without -O')
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
checks=[]


def zero(name,expr):
    assert S.cancel(expr)==0,name
    checks.append(name)


def yes(name,condition):
    assert condition,name
    checks.append(name)


q=S.Rational
t,Y1,Y2,Y3,Y4=S.symbols('t Y1 Y2 Y3 Y4')
bt=Y4-3*Y2*t*t+2*Y1*t**3-t**4
at=Y3-2*Y2*t+Y1*t*t-t**3-t*t*bt
R=S.expand(at*at+at+t*bt+t**3+Y1*t-Y2)
yes('monic degree-twelve eliminant',S.degree(R,t)==12 and R.coeff(t,12)==1)
L=t*bt+t**3+Y1*t-Y2-q(1,4)
zero('square completion',R-(at+q(1,2))**2-L)
yes('square obstruction',S.degree(at,t)==6 and S.degree(L,t)==5 and S.expand(L).coeff(t,5)==-1)


def partitions(n,maximum=None):
    if n==0:
        yield []
        return
    for first in range(min(n,maximum or n),0,-1):
        for rest in partitions(n-first,first):yield [first]+rest


emptytypes=[p for p in partitions(12) if 1 not in p and any(m%2 for m in p)]
yes('ten possible nonsquare empty types',len(emptytypes)==10)
for part in emptytypes:
    odd=[m%2 for m in part]
    if sum(odd)==2:
        exps=[(m-3*o)//2 for m,o in zip(part,odd)]
        yes('quadratic-cube cubic-square type '+str(part),
            sum(exps)==3 and min(exps)>=0 and all(3*o+2*k==m for m,o,k in zip(part,odd,exps)))
    else:
        yes('quartic-cube type '+str(part),part==[3,3,3,3])

# Quartic-cube branch. The explicit solutions below follow from the
# coefficients11,10,9,8,6 respectively; no solve/GB call is used here.
u,v,z=S.symbols('u v z')
w4=u**3/24-u*v/2-q(2,3)
quartic=t**4+u*t**3+v*t*t+w4*t+z
target=S.Poly(S.expand(quartic**3),t)
y1=-target.nth(11)/4
y2=(target.nth(10)-4*y1*y1)/6
y4=(6*y1+9*y2*y2-target.nth(8))/2
y3=(target.nth(6)-14*y1*y2+6*y2*y4-2)/2
res=S.Poly(S.expand(R.subs({Y1:y1,Y2:y2,Y3:y3,Y4:y4})-target.as_expr()),t)
for j in [12,11,10,9,8,6]:zero('quartic eliminated coefficient '+str(j),res.nth(j))
quartic_data=json.loads((HERE/'quartic_cube-input.json').read_text())
for row in quartic_data['equations']:
    saved=S.sympify(row['polynomial'])
    ratio=S.cancel(res.nth(row['coefficient'])/saved)
    yes('quartic input coefficient '+str(row['coefficient']),ratio.is_Rational and ratio!=0)
yes('quartic remaining coefficients exhausted',
    set(row['coefficient'] for row in quartic_data['equations'])=={0,1,2,3,4,5,7})

# Quadratic-cube times cubic-square branch, centered at the quadratic.
d,alpha,beta,kappa=S.symbols('d alpha beta kappa')
s=t+u/2
Q=s*s-d/4
T=s**3+alpha*s*s+beta*s+kappa
F=s**6+alpha*s**5+(beta-3*d/8)*s**4+(kappa-3*alpha*d/8)*s**3 \
    +(-3*beta*d/8+3*d*d/128)*s*s \
    +(-3*kappa*d/8+3*alpha*d*d/128)*s+3*beta*d*d/128+d**3/1024
fp=S.Poly(S.expand(F),t)
lp=S.Poly(S.expand(Q**3*T*T-F*F),t)
yes('polynomial part of square root',lp.degree()<=5)
tail1=3*kappa*d*d/128+alpha*d**3/1024
tail2=beta*d**3/1024+3*d**4/32768
tail3=kappa*d**3/1024+3*alpha*d**4/32768
c1=tail1+q(1,2)
c2=tail2+u/4
c3=tail3-q(1,2)+u*u/8
f4s=beta-3*d/8
zero('leading remainder constraint',lp.nth(5)+1-2*c1)
zero('next remainder constraint',lp.nth(4)+fp.nth(5)-2*c2-(2*alpha+5*u)*c1)
zero('third remainder constraint',lp.nth(3)+fp.nth(4)-1-2*c3
     -(2*alpha+4*u)*c2-(2*f4s+4*alpha*u+5*u*u)*c1)
det=S.Matrix([c1,c2,c3]).jacobian([alpha,beta,kappa]).det()
yes('linear reduction invertible for d nonzero',S.cancel(det/d**9).is_Rational and S.cancel(det/d**9)!=0)
fit={alpha:q(2048,5)*(24-6*u*u+d)/d**4,
     beta:-256*u/d**3-3*d/32}
fit[kappa]=-q(64,3)/d**2-fit[alpha]*d/24
for i,c in enumerate([c1,c2,c3],1):zero('unique fitted parameter '+str(i),c.subs(fit))
fittedF=S.Poly(S.expand(F.subs(fit)),t)
fittedL=S.Poly(S.expand(lp.as_expr().subs(fit)),t)
necessary=[fittedF.nth(3)+1,3*fittedF.nth(1)+2*fittedF.nth(4),
           fittedL.nth(2),fittedL.nth(1)+fittedF.nth(5)+fittedF.nth(2),
           fittedL.nth(0)+fittedF.nth(4)/3+q(1,4)]
centered_data=json.loads((HERE/'centered-quadratic-input.json').read_text())
yes('five centered residual equations',len(centered_data['equations'])==len(necessary)==5)
for expr,row in zip(necessary,centered_data['equations']):
    num,den=S.cancel(expr).as_numer_denom()
    denpoly=S.Poly(den,u,d,domain=S.QQ)
    yes('only invert d in '+row['name'],len(denpoly.terms())==1 and denpoly.monoms()[0][0]==0)
    ratio=S.cancel(num/S.sympify(row['polynomial']))
    yes('centered input '+row['name'],ratio.is_Rational and ratio!=0)

# Literal multiplications/additions over Q: the replay never asks whether
# an element belongs to an ideal and never invokes Singular.
certs=json.loads((HERE/'certificates.json').read_text())
stats={}
for name,data in certs.items():
    vs=data['variables'];symbols=[S.Symbol(x) for x in vs]
    source=json.loads((HERE/data['source']).read_text())
    context=fmpq_mpoly_ctx.get(vs)
    total=context.constant(0)
    yes(name+' generators exhausted',sorted(data['generator_indices'])==list(range(len(source['equations']))))
    yes(name+' multiplier count',len(data['multipliers'])==len(data['generator_indices']))
    for row_index,terms in zip(data['generator_indices'],data['multipliers']):
        dictionary={}
        for term in terms:
            exps=tuple(term['exponents'])
            yes(name+' valid sparse term '+str(row_index)+' '+str(exps),
                len(exps)==len(vs) and min(exps)>=0 and exps not in dictionary)
            dictionary[exps]=fmpq(term['coefficient'])
        multiplier=context.from_dict(dictionary)
        poly=S.Poly(S.sympify(source['equations'][row_index]['polynomial']),*symbols,domain=S.QQ)
        generator=context.from_dict({m:fmpq(int(c.p),int(c.q)) for m,c in poly.terms()})
        total+=multiplier*generator
    expected=context.constant(1) if name=='quartic' else context.gen(vs.index('d'))
    yes(name+' exact ideal identity',total==expected)
    stats[name]={'multipliers':len(data['multipliers']),
                 'multiplier_terms':sum(len(t) for t in data['multipliers']),
                 'target':data['target'],'identity':'PASS'}

manifest=json.loads((HERE/'sources/manifest.json').read_text())
for row in manifest['sources']:
    yes('source hash '+row['repository_path'],hashlib.sha256((ROOT/row['repository_path']).read_bytes()).hexdigest()==row['sha256'])
receipt={'result':'PASS','claim_id':'GAOF7SURJ1','checks':len(checks),'check_names':checks,
         'python':platform.python_version(),'sympy':S.__version__,
         'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'certificates_sha256':hashlib.sha256((HERE/'certificates.json').read_bytes()).hexdigest(),
         'certificate_replays':stats,
         'input_hashes':{n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in ['quartic_cube-input.json','centered-quadratic-input.json']},
         'conclusion':'Every F7 eliminant has a simple root; together with the complete boundary count, F7 is surjective over C.',
         'remaining_global_cardinality_question':3,
         'assurance':'Independent FLINT arithmetic replay of Singular identities plus written reduction/exhaustion proof; no independent end-to-end proof review, formal verification, external review or novelty priority.'}
(HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['result','checks','certificate_replays','conclusion']}))

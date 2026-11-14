"""Assemble a finite-algebra candidate. Only verify.py can accept the theorem."""
from pathlib import Path
import json,time,math
from fractions import Fraction
import sympy as S
from flint import fmpz_poly,fmpz_mpoly_ctx
HERE=Path(__file__).resolve().parent;START=time.monotonic()
rur=json.loads((HERE/'rur-candidate.json').read_text());cofactor=json.loads((HERE/'cofactor-candidate.json').read_text())['records'][0]
old=json.loads((HERE.parent/'gao-f7-linear-cube/high-root-certificate.json').read_text());record=old['resultants'][1]
sf=fmpz_poly(rur['squarefree_resultant']);pr=fmpz_poly(record['primitive_reduced_coefficients']);cover,remainder=divmod(sf*sf,pr);assert not remainder
out,remainder=divmod(fmpz_poly(rur['generic_initial']),fmpz_poly(rur['generic_retained']));assert not remainder
exceptions=[]
for name in ['quintuple-szero','quintuple-cubic-two','quintuple-cubic-plus','quintuple-cubic-minus']:
 source=json.loads((HERE/(name+'-input.json')).read_text());text=(HERE/(name+'.log')).read_text();assert '?' not in text,name
 lines=text.split('BEGIN_MULTIPLIERS\n')[1].split('END_MULTIPLIERS')[0].strip().splitlines();V=S.symbols(' '.join(source['variables']));cols=2 if name=='quintuple-szero' else 1;assert len(lines)==cols*len(source['equations'])
 multipliers=[[] for _ in range(cols)]
 for i in range(len(source['equations'])):
  for j in range(cols):
   q=S.Poly(S.sympify(lines[cols*i+j].replace('^','**')),*V,domain=S.QQ)
   multipliers[j].append([{'exponents':list(ex),'coefficient':str(c)} for ex,c in q.terms() if c])
 ctx=fmpz_mpoly_ctx.get(source['variables']);originals=[ctx.from_dict({ex:int(c) for ex,c in S.Poly(S.sympify(row['polynomial']),*V).terms()}) for row in source['equations']]
 for j,identity in enumerate(multipliers):
  denominator=math.lcm(*(Fraction(term['coefficient']).denominator for terms in identity for term in terms));total=ctx.constant(0)
  for terms,original in zip(identity,originals):
   q={tuple(term['exponents']):int(Fraction(term['coefficient'])*denominator) for term in terms};total+=ctx.from_dict(q)*original
  target=ctx.gens()[j] if cols==2 else ctx.constant(1);assert total==denominator*target,name
 exceptions.append({'name':name,'source':name+'-input.json','targets':['a','b'] if cols==2 else ['1'],'identities':multipliers})
 print(json.dumps({'exception':name,'literal_candidate_identity':'PASS','terms':[sum(map(len,row)) for row in multipliers]}),flush=True)
certificate={'status':'CANDIDATE_ONLY; complete independent replay pending','source_resultant':record,'cofactor':cofactor,'squarefree_projection':rur['squarefree_resultant'],'projection_cover_quotient':list(map(int,cover.coeffs())),'cofactor_exception':rur['cofactor_exception'],'generic_initial':rur['generic_initial'],'generic_retained':rur['generic_retained'],'generic_out':list(map(int,out.coeffs())),'prime':1000000007,'exceptions':exceptions}
(HERE/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
print(json.dumps({'status':certificate['status'],'certificate_bytes':(HERE/'certificate.json').stat().st_size,'seconds':time.monotonic()-START}),flush=True)

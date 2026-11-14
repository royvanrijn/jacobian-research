"""Assemble the existing resultant and exceptional-identity candidate; no proof decision."""
from pathlib import Path
import json,sympy as S,time
from flint import fmpz_poly,nmod_mat
HERE=Path(__file__).resolve().parent;START=time.monotonic()
rows=[r for r in json.loads((HERE/'resultants-candidate.json').read_text())['records'] if r['pair']!=[1,3]]
assert [r['pair'] for r in rows]==[[0,1],[0,3],[0,2]]
prime=1000000007
f=fmpz_poly(rows[0]['primitive_reduced_coefficients']);g=fmpz_poly(rows[1]['primitive_reduced_coefficients'])+fmpz_poly(rows[2]['primitive_reduced_coefficients'])
def sylvester_coefficients(f,g):
 m,n=f.degree(),g.degree();N=m+n;matrix=[[0]*N for _ in range(N)]
 for i in range(n):
  for j,c in enumerate(reversed(list(f.coeffs()))):matrix[i][i+j]=int(c)
 for i in range(m):
  for j,c in enumerate(reversed(list(g.coeffs()))):matrix[n+i][i+j]=int(c)
 return matrix
matrix=sylvester_coefficients(f,g);det=int(nmod_mat([[c%prime for c in row] for row in matrix],prime).det());assert det!=0
identities=[]
for name in ['high-root-szero-integral','high-root-cubic-exception-integral']:
 source=json.loads((HERE/(name+'-input.json')).read_text());log=(HERE/(name+'.log')).read_text();assert '?' not in log
 content=log.split('BEGIN_MULTIPLIERS\n')[1].split('END_MULTIPLIERS')[0].strip().splitlines();assert len(content)==len(source['equations'])
 variables=S.symbols(' '.join(source['variables']))
 multipliers=[]
 for line in content:
  poly=S.Poly(S.sympify(line.replace('^','**')),*variables,domain=S.QQ)
  multipliers.append([{'exponents':list(e),'coefficient':str(c)} for e,c in poly.terms()])
 identities.append({'source':name+'-input.json','multipliers':multipliers})
certificate={'status':'CANDIDATE_ONLY; separate replay required','resultants':rows,'combination':[1,1],'prime':prime,'coprimality_sylvester_shape':[len(matrix),len(matrix)],'det_mod':det,'exceptional_identities':identities}
(HERE/'high-root-certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
print(json.dumps({'status':certificate['status'],'resultant_degrees':[r['original_degree'] for r in rows],'coprimality_shape':certificate['coprimality_sylvester_shape'],'det_mod':det,'exceptional_multiplier_terms':[sum(map(len,row['multipliers'])) for row in identities],'certificate_bytes':(HERE/'high-root-certificate.json').stat().st_size,'seconds':time.monotonic()-START}))

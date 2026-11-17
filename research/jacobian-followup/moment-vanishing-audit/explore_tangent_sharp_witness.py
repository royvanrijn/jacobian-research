from pathlib import Path
import json,sympy as S
p=Path(__file__).resolve().parent;e,f,h=S.symbols('e f h')
data=json.loads((p/'line_conic_charts_input.json').read_text())[0]
polys=[]
for q in data['remaining_moments'][:3]:
 q=S.Poly(S.sympify(q,locals={'e':e,'f':f}),e,f)
 assert all(j%2==0 for (i,j),v in q.terms())
 polys.append(S.expand(sum(v*e**i*h**(j//2) for (i,j),v in q.terms())))
res=S.factor(S.resultant(polys[0],polys[1],h));print('RESULTANT',res)
gb=S.groebner(polys[:2],h,e);print('LEX_GB',*[S.factor(q.as_expr()) for q in gb.polys],sep='\n')
(p/'tangent_sharp_witness_exploration.json').write_text(json.dumps({'M4_M6_resultant_in_f_squared':str(res),'lex_basis':list(map(lambda q:str(q.as_expr()),gb.polys)),'M4_M6_M8_in_e_h':list(map(str,polys))},indent=2)+'\n')

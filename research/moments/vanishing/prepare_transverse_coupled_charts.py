#!/usr/bin/env python3
from pathlib import Path
import sympy as S,json
p=Path(__file__).resolve().parent
a,b,c,d,e,f=S.symbols('a b c d e f');syms={str(x):x for x in [a,b,c,d,e,f]}
mom=[S.sympify(x,locals=syms) for x in json.loads((p/'line_conic_moments_input.json').read_text())['families']['transverse']['primitive_moments']]
det=S.det(S.Matrix([[a,d,e],[d,b,f],[e,f,c]]))
def prim(expr,vs):
 n,de=S.fraction(S.cancel(expr));return S.Poly(n,*vs).clear_denoms()[1].primitive()[1].as_expr()
records=[]
for label,vs,sub in [('TRANSVERSE_COUPLED_OPEN',[b,e,f],{a:1,d:1}),('TRANSVERSE_COUPLED_BOUNDARY',[e,f,c],{a:1,d:1,b:-2})]:
 base=[S.expand(m.subs(sub)) for m in mom]
 if label.endswith('OPEN'):
  sol=S.solve(base[0],c)[0];ps=[prim(q.subs(c,sol),vs) for q in base[1:]];delta=prim(det.subs(sub).subs(c,sol),vs)*(b+2)
  info={'eliminated':'c','solution':str(sol),'localizer_extra':'b+2'}
 else:
  relation=S.Poly(base[0],e,domain=S.QQ.frac_field(f,c));ps=[base[0]]+[S.Poly(q,e,domain=S.QQ.frac_field(f,c)).rem(relation).as_expr() for q in base[1:]]
  ps=[prim(q,vs) for q in ps];delta=prim(det.subs(sub),vs);info={'reduction':'Every higher moment is reduced by M2, a quadratic in e with constant nonzero leading coefficient.'}
 source='LIB "elim.lib";\nring r=0,('+','.join(map(str,vs))+'),dp;\noption(redSB);\nideal I='+','.join(str(q).replace('**','^') for q in ps)+';\npoly delta='+str(S.expand(delta)).replace('**','^')+';\n'
 source+='print("'+label+'_STANDARD_START");ideal G=std(I);\n'
 source+='print("'+label+'_BASIS_SIZE");print(size(G));print("'+label+'_BASIS");print(G);\n'
 source+='for(int kk=1;kk<=4;kk++){poly remainderPoly=reduce(delta^kk,G);print("'+label+'_POWER");print(kk);print(size(remainderPoly));if(remainderPoly==0){break;}}\n'
 source+='print("END");quit;\n'
 stem=label.lower();runner=(p/'run_unmixed_depth_audit.py').read_text();i=runner.index("text=(p/");j=runner.index("fpath=p/",i)
 runner=runner[:i]+'source='+repr(source)+'\n'+runner[j:];runner=runner.replace('unmixed_depth_audit',stem).replace('(20,20)','(30,30)').replace('timeout=20','timeout=30').replace("'wall_limit_seconds':20","'wall_limit_seconds':30")
 runner=runner.replace('Minimal free resolution/depth of the retained exact isolated-primary intersection U, solely to test whether the quaternary generic-section torsion lower bound is exact.','Exact coupled transverse line-conic first-ten moment chart, with explicit treatment of the M2 elimination denominator; no theorem from a timeout or nonzero target remainder.')
 (p/('run_'+stem+'.py')).write_text(runner)
 records.append({'label':label,'variables':list(map(str,vs)),'normalization':{str(k):str(v) for k,v in sub.items()},'polynomials':list(map(str,ps)),'localizer':str(S.expand(delta)),**info});print('Prepared',stem,flush=True)
(p/'transverse_coupled_charts_input.json').write_text(json.dumps(records,indent=2)+'\n')

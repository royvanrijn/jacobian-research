#!/usr/bin/env python3
from pathlib import Path
import sympy as S,json
p=Path(__file__).resolve().parent
d1,d2,d3,h0,h1,h2,t=S.symbols('d1 d2 d3 h0 h1 h2 t');vs=[d1,d2,d3,h0,h1,h2]
coeffs=[]
for k in range(1,6):
 alpha=-k-S.Rational(1,2);ser=[S.Integer(1)]
 for n in range(1,2*k+1):
  ser.append(S.Poly(sum(((alpha+1)*i-n)*[d1,d2,d3][i-1]*ser[n-i] for i in range(1,min(3,n)+1))/n,*vs).as_expr())
 hp=S.Poly((h0+h1*t+h2*t*t)**k,t)
 result=S.Poly(sum(hp.nth(j)*ser[2*k-j] for j in range(2*k+1)),*vs)
 coeffs.append(result.as_expr());print('C',k,'terms',len(result.terms()),flush=True)
(p/'linear_quadratic_invariants_input.json').write_text(json.dumps({'D':'1+d1*t+d2*t^2+d3*t^3','H':'h0+h1*t+h2*t^2','indices':[2,4,6,8,10],'normalized_moment_coefficients':list(map(str,coeffs))},indent=2)+'\n')
def prim(q,variables):return S.Poly(S.expand(q),*variables).clear_denoms()[1].primitive()[1].as_expr()
records=[];blocks=[]
cases=[('LQ_NONISO_D1',[d2,d3,h1],{h0:1,d1:1}),('LQ_NONISO_D1ZERO_H1',[d2,d3],{h0:1,d1:0,h1:1}),('LQ_NONISO_D1H1ZERO_D2',[d3],{h0:1,d1:0,h1:0,d2:1})]
for label,v,sub in cases:
 sol=S.solve(coeffs[0].subs(sub),h2)[0];ps=[prim(q.subs(sub).subs(h2,sol),v) for q in coeffs[1:]]
 source='ring r'+label+'=0,('+','.join(map(str,v))+'),dp;\noption(redSB);\nideal All='+','.join(str(q).replace('**','^') for q in ps)+';\n'
 source+='for(int kk=2;kk<=4;kk++){ideal I=0;for(int jj=1;jj<=kk;jj++){I=I,All[jj];}\n'
 source+='print("'+label+'_GATE");print(2*(kk+1));ideal G=std(I);print("'+label+'_BASIS_SIZE");print(size(G));print("'+label+'_BASIS");print(G);\n'
 source+='for(int pp=1;pp<=4;pp++){poly remainderPoly=reduce(d3^pp,G);print("'+label+'_POWER");print(pp);print(size(remainderPoly));if(remainderPoly==0){break;}}\n'
 source+='ideal Sat=sat(G,ideal(d3));print("'+label+'_SATURATED");print(Sat);if(reduce(1,Sat)==0){break;}\n}\n'
 records.append({'label':label,'variables':list(map(str,v)),'normalization':{str(k):str(v) for k,v in sub.items()},'h2_solution':str(sol),'remaining_moments':list(map(str,ps))});blocks.append(source)
# The isotropic L branch with h1 nonzero is normalized by rescaling L.
sub={h0:0,h1:1};sol=S.solve(coeffs[0].subs(sub),h2)[0];ps=[prim(q.subs(sub).subs(h2,sol),[d1,d2,d3]) for q in coeffs[1:]]
source='ring rLQ_ISO_H1=0,(d1,d2,d3),dp;\noption(redSB);\nideal I='+','.join(str(q).replace('**','^') for q in ps[:3])+';\nideal G=std(I);print("LQ_ISO_H1_BASIS");print(G);print("LQ_ISO_H1_D3_POWER");print(reduce(d3^4,G));\n'
records.append({'label':'LQ_ISO_H1','variables':['d1','d2','d3'],'normalization':{'h0':'0','h1':'1'},'h2_solution':str(sol),'remaining_moments':list(map(str,ps))})
(p/'linear_quadratic_invariant_charts_input.json').write_text(json.dumps(records,indent=2)+'\n')
for stem,selected in [('linear_quadratic_noniso_generic',blocks[:1]),('linear_quadratic_other_invariant_charts',blocks[1:]+[source])]:
 src='LIB "elim.lib";\n'+''.join(selected)+'print("END");quit;\n'
 runner=(p/'run_unmixed_depth_audit.py').read_text();i=runner.index("text=(p/");j=runner.index("fpath=p/",i)
 runner=runner[:i]+'source='+repr(src)+'\n'+runner[j:];runner=runner.replace('unmixed_depth_audit',stem).replace('(20,20)','(30,30)').replace('timeout=20','timeout=30').replace("'wall_limit_seconds':20","'wall_limit_seconds':30")
 runner=runner.replace('Minimal free resolution/depth of the retained exact isolated-primary intersection U, solely to test whether the quaternary generic-section torsion lower bound is exact.','Exact invariant D/H coefficient charts for a linear form times a nondegenerate quadratic in three Gaussian variables; finite exclusions require the written Gaussian generating-function identity and exhaustive rescaling cover.')
 (p/('run_'+stem+'.py')).write_text(runner)
 print('Prepared',stem)

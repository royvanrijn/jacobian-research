#!/usr/bin/env python3
from pathlib import Path
import sympy as S,json
p=Path(__file__).resolve().parent
a,b,c,d,e,f,u=S.symbols('a b c d e f u');syms={str(x):x for x in [a,b,c,d,e,f,u]}
data=json.loads((p/'line_conic_moments_input.json').read_text());mom={k:[S.sympify(x,locals=syms) for x in v['primitive_moments']] for k,v in data['families'].items()}
det=S.det(S.Matrix([[a,d,e],[d,b,f],[e,f,c]]))
def prim(expr,vs):return S.Poly(S.cancel(expr),*vs).clear_denoms()[1].primitive()[1].as_expr()
def block(label,vs,polys,delta):
 source='ring r'+label+'=0,('+','.join(map(str,vs))+'),dp;\noption(redSB);\n'
 source+='ideal All='+','.join(str(x).replace('**','^') for x in polys)+';\npoly delta='+str(delta).replace('**','^')+';\n'
 source+='for(int kk=2;kk<=4;kk++){\n ideal I=0;for(int jj=1;jj<=kk;jj++){I=I,All[jj];}\n'
 source+=' ideal G=std(I);print("'+label+'_GATE");print(2*(kk+1));\n'
 source+=' print("'+label+'_BASIS");print(G);\n'
 source+=' poly remainderPoly=reduce(delta^12,G);print("'+label+'_DET_REMAINDER_SIZE");print(size(remainderPoly));\n'
 source+=' ideal Sat=sat(G,ideal(delta));print("'+label+'_SATURATED");print(Sat);\n'
 source+=' if(reduce(1,Sat)==0){break;}\n}\n'
 return source
records=[];blocks=[]
for label,vs,sub,elim in [('TANGENT_A1B1D0',[e,f],{a:1,b:1,d:0},c),('TANGENT_A1B0D0',[e,f],{a:1,b:0,d:0},c),('TANGENT_A0B0D1',[c,f],{a:0,b:0,d:1},e),('TRANSVERSE_A0D1E1',[b,f],{a:0,d:1,e:1},c)]:
 kind='tangent' if label.startswith('TANGENT') else 'transverse';eq=S.expand(mom[kind][0].subs(sub));sol=S.solve(eq,elim)[0]
 ps=[prim(m.subs(sub).subs(elim,sol),vs) for m in mom[kind][1:]];dd=prim(det.subs(sub).subs(elim,sol),vs)
 records.append({'label':label,'kind':kind,'variables':list(map(str,vs)),'normalization':{str(x):str(v) for x,v in sub.items()},'eliminated':str(elim),'solution':str(sol),'remaining_moments':list(map(str,ps)),'determinant':str(dd)})
 blocks.append(block(label,vs,ps,dd))
# In the decoupled transverse chart every monomial has equal b,c powers.
ps=[]
for m in mom['transverse']:
 q=S.Poly(m.subs({a:1,d:0,e:0}),b,c,f);assert all(i==j for (i,j,k),v in q.terms())
 ps.append(S.expand(sum(v*u**i*f**k for (i,j,k),v in q.terms())))
sol=S.solve(ps[0],u)[0];dd=prim((u-f*f).subs(u,sol),[f]);other=[prim(q.subs(u,sol),[f]) for q in ps[1:]]
records.append({'label':'TRANSVERSE_DECOUPLED','kind':'transverse','variables':['f'],'normalization':{'a':'1','d':'0','e':'0'},'eliminated':'u=bc','solution':str(sol),'remaining_moments':list(map(str,other)),'determinant':str(dd)})
blocks.append(block('TRANSVERSE_DECOUPLED',[f],other,dd))
(p/'line_conic_charts_input.json').write_text(json.dumps(records,indent=2)+'\n')
for stem,selected in [('tangent_line_conic_charts',blocks[:3]),('transverse_easy_charts',blocks[3:])]:
 source='LIB "elim.lib";\n'+''.join(selected)+'print("END");quit;\n'
 runner=(p/'run_unmixed_depth_audit.py').read_text();i=runner.index("text=(p/");j=runner.index("fpath=p/",i)
 runner=runner[:i]+'source='+repr(source)+'\n'+runner[j:]
 runner=runner.replace('unmixed_depth_audit',stem)
 runner=runner.replace('Minimal free resolution/depth of the retained exact isolated-primary intersection U, solely to test whether the quaternary generic-section torsion lower bound is exact.','Exact line-conic covariance charts, beginning with the first-six pure moments and extending only if the determinant-saturated ideal is nonunit; interpretation requires the complete written chart cover.')
 (p/('run_'+stem+'.py')).write_text(runner)
 print('Prepared',stem)

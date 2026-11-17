#!/usr/bin/env python3
from pathlib import Path
import sympy as S,json
p=Path(__file__).resolve().parent
a,b,c,d,e,f=S.symbols('a b c d e f');syms={str(v):v for v in [a,b,c,d,e,f]}
data=json.loads((p/'cuspidal_moments_ten_input.json').read_text());mom=[S.sympify(x,locals=syms) for x in data['primitive_moments']]
det=S.det(S.Matrix([[a,d,e],[d,b,f],[e,f,c]]))
cases=[('AB11',[d,e,f],{a:1,b:1}),('A0B1D1',[c,e],{a:0,b:1,d:1}),('A1B0D1',[e,f],{a:1,b:0,d:1}),('A1B0D0',[e,f],{a:1,b:0,d:0}),('A0B0D1',[c,e],{a:0,b:0,d:1})]
blocks=[];records=[]
for label,vs,sub in cases:
 eliminated=f if label in ['A0B1D1','A0B0D1'] else c
 eq=S.expand(mom[0].subs(sub));solution=S.solve(eq,eliminated)[0]
 polys=[]
 for m in mom[1:]:
  poly=S.Poly(S.expand(m.subs(sub).subs(eliminated,solution)),*vs,domain=S.QQ)
  poly=poly.clear_denoms()[1].primitive()[1];polys.append(poly.as_expr())
 dd=S.expand(det.subs(sub).subs(eliminated,solution));den,ddp=S.Poly(dd,*vs).clear_denoms();dd=ddp.primitive()[1].as_expr()
 block='ring r'+label+'=0,('+','.join(map(str,vs))+'),dp;\noption(redSB);\n'
 block+='ideal I='+','.join(str(x).replace('**','^') for x in polys)+';\n'
 block+='poly delta='+str(dd).replace('**','^')+';\n'
 block+='print("'+label+'_START");ideal G=std(I);\n'
 block+='print("'+label+'_BASIS_SIZE");print(size(G));print("'+label+'_BASIS");print(G);\n'
 block+='print("'+label+'_DET_POWER");print(reduce(delta^12,G));\n'
 block+='ideal Sat=sat(G,ideal(delta));print("'+label+'_SATURATED");print(Sat);\n'
 blocks.append(block)
 records.append({'label':label,'variables':list(map(str,vs)),'normalization':{str(k):str(v) for k,v in sub.items()},'eliminated':str(eliminated),'solution':str(solution),'primitive_remaining_moments':list(map(str,polys)),'determinant_primitive':str(dd)})
(p/'cuspidal_covariance_charts_input.json').write_text(json.dumps(records,indent=2)+'\n')
for stem,selected,cap in [('cuspidal_generic_chart',blocks[:1],30),('cuspidal_boundary_charts',blocks[1:],20)]:
 source='LIB "elim.lib";\n'+''.join(selected)+'print("END");quit;\n'
 runner=(p/'run_unmixed_depth_audit.py').read_text();i=runner.index("text=(p/");j=runner.index("fpath=p/",i)
 runner=runner[:i]+'source='+repr(source)+'\n'+runner[j:]
 runner=runner.replace('unmixed_depth_audit',stem).replace('(20,20)',f'({cap},{cap})').replace('timeout=20',f'timeout={cap}').replace("'wall_limit_seconds':20",f"'wall_limit_seconds':{cap}")
 runner=runner.replace('Minimal free resolution/depth of the retained exact isolated-primary intersection U, solely to test whether the quaternary generic-section torsion lower bound is exact.','Exact normalized covariance chart moment prefix through10 for the cuspidal cubic; conclusions require the exhaustive written torus chart cover and exact input reconstruction.')
 (p/('run_'+stem+'.py')).write_text(runner)
 print('Prepared',stem,'cap',cap)

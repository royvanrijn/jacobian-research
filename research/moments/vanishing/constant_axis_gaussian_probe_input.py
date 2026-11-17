from pathlib import Path
import sys,json,time
p=Path(__file__).resolve().parents[1]/'research/moments/vanishing';sys.path.insert(0,str(p))
import sympy as s
import verify_ordinary_rank_one_lattice as old
import verify_gaussian_regular_projection as dets
t,u,v=old.t,old.u,old.v
original=old.parameters

def params(kind,mode):
 den=1-2*t if mode=='linear' else s.Integer(1);m=1;A=11;B=13
 F=u*v-6+t*t-t;xx=u+3;yy=v+2
 G=(xx**3 if kind=='transverse' else 0)+7*t*t*xx+11*t*t*yy+2*t**3-(A*xx*xx+2*xx*yy+B*yy*yy)/2-5*t*xx-7*t*yy-2*t*t
 W=s.expand(den*G+F*F/2);metric=s.Matrix([[A,1,2],[1,B,3],[2,3,m]]).det();return W,den,metric
old.parameters=params;rows=[]
for kind in ['transverse','triangle']:
 for mode in ['linear','constant']:
  print('Starting',kind,mode,flush=True);m=old.model(kind,mode);n=len(m['seed_columns'])
  cols=[[s.Poly(int(i==0),t) for i in range(n)]]+[[s.Poly(s.sympify(x,locals={'t':t}),t) for x in z['N']] for z in m['seed_columns']]
  co=dets.determinant_cofactors(cols,n);den=s.Poly(s.sympify(m['D'],locals={'t':t}),t);axis=s.Poly(s.sympify(m['Lf'],locals={'t':t}),t)
  C0=co[0].exquo(den**(n*(n-1))*axis**(20 if kind=='transverse' else 16));print(kind,mode,'degree',C0.degree(),'Lf',axis.as_expr(),flush=True)
  assert C0 and C0.degree()<=(77 if kind=='transverse' else 82)
  m['new_C0_degree']=int(C0.degree());rows.append(m)
Path(__file__).with_suffix('.json').write_text(json.dumps(rows,indent=2)+'\n')

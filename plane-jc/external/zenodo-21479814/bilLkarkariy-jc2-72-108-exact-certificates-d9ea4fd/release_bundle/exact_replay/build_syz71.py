from pathlib import Path
import pickle,sympy as sp
p=71;S=pickle.load(open('hne0_polred.pkl','rb'));w=sp.symbols('w');assert sp.Poly(w**5-w**4+3*w**3+3*w**2+26,w,modulus=p).is_irreducible
def cm(c):
 t=[]
 for i,(n,d) in sorted(c.items()):
  a=n%p*pow(d%p,-1,p)%p
  if a:
   m='' if i==0 else ('w' if i==1 else f'w^{i}');t.append(str(a) if not m else (m if a==1 else f'{a}*{m}'))
 return '+'.join(t) or '0'
def es(e):
 t=[]
 for m,c in sorted(e.items(),key=lambda z:(sum(z[0]),z[0]),reverse=True):
  a=cm(c);mon='*'.join(n if q==1 else f'{n}^{q}' for n,q in zip(('h','u1','u2'),m) if q);t.append(a if not mon else (mon if a=='1' else f'({a})*{mon}'))
 return '+'.join(t) or '0'
L=[f'ring R=({p},w),(h,u1,u2),dp;','minpoly=w^5-w^4+3*w^3+3*w^2+26;','ideal I='+',\n'.join(es(e) for e in S)+';','ideal J0=I[1],I[2],I[3],I[4];','option(redSB);','matrix T;ideal J=liftstd(J0,T);','int jj=0;int k;for(k=1;k<=size(J);k++){if(J[k]==h||J[k]==-h){jj=k;}}','module SY=syz(J0);module GS=std(SY);vector vv;int i;for(i=1;i<=4;i++){vv=vv+T[i,jj]*gen(i);}vector vr=reduce(vv,GS);','intvec ee;poly q;number c;poly chk=0;for(i=1;i<=4;i++){chk=chk+J0[i]*vr[i];q=vr[i];while(q!=0){ee=leadexp(q);c=leadcoef(q);print("C|"+string(i)+"|"+string(ee[1])+"|"+string(ee[2])+"|"+string(ee[3])+"|"+string(c));q=q-lead(q);}}','if(chk==h){print("CHECK_PASS");}else{print("CHECK_FAIL");}','quit;']
Path('hcert71.sing').write_text('\n'.join(L)+'\n')

import pickle,re
p=71;S=pickle.load(open('hne0_polred.pkl','rb'));pat=re.compile(r'^C\|(\d+)\|(\d+)\|(\d+)\|(\d+)\|')
cols=[]
for l in open('hcert71.out'):
 m=pat.match(l)
 if m:cols.append(tuple(map(int,m.groups())))
rows={(1,0,0)}
for i,a,b,c in cols:
 for e in S[i-1]:rows.add((a+e[0],b+e[1],c+e[2]))
rows=sorted(rows,key=lambda m:(sum(m),m));ri={m:i+1 for i,m in enumerate(rows)}
def cm(c):
 t=[]
 for j,(n,d) in sorted(c.items()):
  a=n%p*pow(d%p,-1,p)%p
  if a:
   mon='' if j==0 else ('w' if j==1 else f'w^{j}');t.append(str(a) if not mon else (mon if a==1 else f'{a}*{mon}'))
 return '+'.join(t) or '0'
L=[f'ring R=({p},w),x,dp;','minpoly=w^5-w^4+3*w^3+3*w^2+26;',f'matrix M[{len(rows)}][{len(cols)}];']
for j,(i,a,b,c) in enumerate(cols,1):
 for e,co in S[i-1].items():L.append(f'M[{ri[(a+e[0],b+e[1],c+e[2])]}, {j}]={cm(co)};')
L += ['print("ROWS='+str(len(rows))+' COLS='+str(len(cols))+'");','print("RANK="+string(rank(M)));','quit;']
open('rank71.sing','w').write('\n'.join(L)+'\n')
print(len(rows),len(cols),len(L))

"""Exact local nilpotence orders of the fourteen cyclic global targets.

Internal research computation; no global ideal membership is inferred.
"""
import sys,json,itertools,time
from pathlib import Path
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0,str(Path(__file__).parent/'rodrigues_local_length17_certificate'))
import rod_local_taylor_septic as c

data=json.loads((Path(__file__).parent/'rodrigues_local_length17_certificate/certificate_data.json').read_text())
N=17
Z=(F(0),)*N
ONE=(F(1),)+(F(0),)*(N-1)
Ma=sp.Matrix([[sp.Rational(x) for x in row] for row in data['multiplication_a']])
Mb=sp.Matrix([[sp.Rational(x) for x in row] for row in data['multiplication_b']])
basis=[tuple(x) for x in data['quotient_basis']]
MM=[Ma**i*Mb**j for i,j in basis]
table={}
for u in range(N):
    for v in range(N):
        col=MM[u][:,v]
        table[u,v]=tuple((k,F(int(x.p),int(x.q))) for k,x in enumerate(col) if x)

def add(*xs):
    return tuple(sum(x[i] for x in xs) for i in range(N)) if xs else Z
def sc(x,n):
    n=F(n)
    return tuple(v*n for v in x)
def neg(x): return sc(x,-1)
def sub(x,y): return add(x,neg(y))
def mul(x,y):
    out=[F(0)]*N
    for i,xi in enumerate(x):
        if not xi: continue
        for j,yj in enumerate(y):
            if not yj: continue
            xy=xi*yj
            for k,v in table[i,j]:out[k]+=xy*v
    return tuple(out)
def prod(*xs):
    out=ONE
    for x in xs:out=mul(out,x)
    return out
def power(x,n):
    return prod(*([x]*n))
def unit_power(x,exponent):
    assert x[0] != 0
    assert x[0] == 1, 'constant-one units only'
    delta=sub(x,ONE)
    total=ONE
    term=ONE
    coefficient=F(1)
    for k in range(1,7):
        term=mul(term,delta)
        coefficient*=F(exponent-(k-1),k)
        total=add(total,sc(term,coefficient))
    assert power(delta,7)==Z
    return total
def from_expr(expr):
    out=Z
    for (i,j),coef in sp.Poly(sp.expand(expr),c.a,c.b).as_dict().items():
        col=(Ma**i*Mb**j)[:,0]
        v=tuple(F(int(x.p),int(x.q)) for x in col)
        out=add(out,sc(v,F(int(coef.p),int(coef.q))))
    return out
def dot(x,y):return add(*(mul(a,b) for a,b in zip(x,y)))
def mv(M,v):return [dot(row,v) for row in M]
def mm(A,B):return [[dot(A[i],[B[k][j] for k in range(3)]) for j in range(3)] for i in range(3)]
def transpose(A):return [[A[j][i] for j in range(3)] for i in range(3)]
def det(A):
    return add(prod(A[0][0],A[1][1],A[2][2]),prod(A[0][1],A[1][2],A[2][0]),prod(A[0][2],A[1][0],A[2][1]),neg(prod(A[0][2],A[1][1],A[2][0])),neg(prod(A[0][1],A[1][0],A[2][2])),neg(prod(A[0][0],A[1][2],A[2][1])))
def adj(A):
    out=[[None]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            rows=[k for k in range(3) if k!=j]
            cols=[k for k in range(3) if k!=i]
            out[i][j]=sc(sub(mul(A[rows[0]][cols[0]],A[rows[1]][cols[1]]),mul(A[rows[0]][cols[1]],A[rows[1]][cols[0]])),(-1)**(i+j))
    return out

phase=c.F0+c.T1+c.T2+c.T3+c.T4+c.T5+c.T6
coefs={e:from_expr(x) for e,x in sp.Poly(sp.expand(phase),c.z,c.t,c.w).as_dict().items()}
get=lambda e:coefs.get(tuple(e),Z)
v=[get((0,0,1)),get((0,1,0)),get((1,0,0))]
C=[[get((2,0,0)),sc(get((1,1,0)),F(1,2)),sc(get((1,0,1)),F(1,2))],
   [sc(get((1,1,0)),F(1,2)),get((0,2,0)),sc(get((0,1,1)),F(1,2))],
   [sc(get((1,0,1)),F(1,2)),sc(get((0,1,1)),F(1,2)),get((0,0,2))]]
B=[C[2],C[1],C[0]]
Bv=mv(B,v);BBv=mv(B,Bv)
E=transpose([v,Bv,BBv])
D=det(E)
assert D[0]==F(1,8)
assert add(*(B[i][i] for i in range(3)))==Z
B2=mm(B,B)
p=sc(add(*(B2[i][i] for i in range(3))),F(1,2))
r=det(B)
pair=lambda x,y:dot(x,[y[2],y[1],y[0]])
s0=pair(v,v);s1=pair(v,Bv);s2=pair(v,BBv)
assert s2[0]==F(1,4)
G=[[pair(Ecol,Eother) for Eother in (v,Bv,BBv)] for Ecol in (v,Bv,BBv)]
assert G==[[s0,s1,s2],[s1,s2,add(mul(p,s1),mul(r,s0))],[s2,add(mul(p,s1),mul(r,s0)),add(mul(p,s2),mul(r,s1))]]
assert add(det(G),power(D,2))==Z
B3=mm(B2,B)
assert all(B3[i][j]==add(mul(p,B[i][j]),r if i==j else Z) for i in range(3) for j in range(3))
EA=mm(E,adj(E))
assert all(EA[i][j]==(D if i==j else Z) for i in range(3) for j in range(3))

# x=J*adj(E)^T*u/D.  Transform the cubic without dividing by D.
M=transpose(adj(E));M=[M[2],M[1],M[0]]
def polmul(P,Q):
    out={}
    for e,x in P.items():
        for f,y in Q.items():
            ef=tuple(e[i]+f[i] for i in range(3))
            out[ef]=add(out.get(ef,Z),mul(x,y))
    return {e:v for e,v in out.items() if v!=Z}
linear=[{tuple(int(j==i) for j in range(3)):M[k][i] for i in range(3)} for k in range(3)]
upowers=[]
for L in linear:
    ps=[{(0,0,0):ONE}]
    for k in range(3):ps.append(polmul(ps[-1],L))
    upowers.append(ps)
Hbar={}
for e,coef in coefs.items():
    if sum(e)!=3:continue
    terms=polmul(polmul(upowers[0][e[0]],upowers[1][e[1]]),upowers[2][e[2]])
    for f,x in terms.items():Hbar[f]=add(Hbar.get(f,Z),mul(coef,x))

targets={'p':p,'r':r,'s0':s0,'s1':s1}
Ubar={}
for indices in itertools.combinations_with_replacement(range(3),3):
    e=tuple(indices.count(i) for i in range(3))
    multiplicity=sp.factorial(3)
    for i in e:multiplicity/=sp.factorial(i)
    val=sc(Hbar.get(e,Z),F(1,int(multiplicity)))
    Ubar[indices]=val
    if indices==(2,2,2):val=add(val,sc(mul(D,s2),2))
    targets['U'+''.join(map(str,indices))]=val
assert len(targets)==14
for k in range(3):
    assert add(*(mul(G[i][j],Ubar[tuple(sorted((i,j,k)))]) for i in range(3) for j in range(3)))==Z

# Substitution back checks that Hbar(u)=D^3*H3(x) under u=E^T*J*x.
back_matrix=[[E[2][i],E[1][i],E[0][i]] for i in range(3)]
back_linear=[{tuple(int(j==i) for j in range(3)):back_matrix[k][i] for i in range(3)} for k in range(3)]
back_powers=[]
for L in back_linear:
    ps=[{(0,0,0):ONE}]
    for k in range(3):ps.append(polmul(ps[-1],L))
    back_powers.append(ps)
back={}
for e,coef in Hbar.items():
    terms=polmul(polmul(back_powers[0][e[0]],back_powers[1][e[1]]),back_powers[2][e[2]])
    for f,x in terms.items():back[f]=add(back.get(f,Z),mul(coef,x))
for e in itertools.product(range(4),repeat=3):
    if sum(e)==3: assert back.get(e,Z)==mul(power(D,3),get(e))

# Normalize D by lambda=(D0/D)^(1/6), and compare the exact target vectors.
lambda_unit=unit_power(sc(D,8),-F(1,6))
assert mul(power(lambda_unit,6),D)==sc(ONE,F(1,8))
inverseD=sc(unit_power(sc(D,8),-F(1)),8)
normalized={'p':mul(power(lambda_unit,2),p),'r':mul(power(lambda_unit,3),r),
            's0':mul(power(lambda_unit,2),s0),'s1':mul(power(lambda_unit,3),s1)}
normalized_s2=mul(power(lambda_unit,4),s2)
for indices,val in Ubar.items():
    weight=-(sum(indices)+2)
    actual=prod(unit_power(lambda_unit,F(weight)),power(inverseD,3),val)
    name='U'+''.join(map(str,indices))
    if indices==(2,2,2):
        actual=add(actual,sc(normalized_s2,128))
    normalized[name]=actual
    assert actual==prod(unit_power(lambda_unit,F(weight)),power(inverseD,3),targets[name])
print('ADAPTER CHECKS: Gram, Cayley-Hamilton, adjugate, all harmonic traces, full cubic reconstruction, D normalization PASSED.',flush=True)

def ni(x):
    y=ONE
    for k in range(1,8):
        y=mul(y,x)
        if y==Z:return k
    raise AssertionError('nilpotence7 failure')

out={}
for name,x in targets.items():
    assert x[0]==0,(name,x[0])
    n=ni(x)
    assert ni(normalized[name])==n
    out[name]={'nilpotence_index':n,'vector':list(map(str,x)),'last_nonzero_power':list(map(str,power(x,n-1)))}
    print('TARGET',name,'nilindex',n,'nonzero coefficients',[(basis[i],str(v)) for i,v in enumerate(x) if v],flush=True)
payload={'target_profile':out,'normalized_target_vectors':{name:list(map(str,v)) for name,v in normalized.items()},
         'normalized_s2':list(map(str,normalized_s2)),
         'raw_unit_equivalent_coordinates':'p,r,s0,s1; D^3 U_ijk for ijk !=222; D^3 U222+2D s2 for final target',
         'normalization':'lambda=(1/(8D))^(1/6); final normalized target=U222+128s2',
         'adapter_checks_passed':True}
Path(__file__).with_suffix('.json').write_text(json.dumps(payload,indent=2)+'\n')
print('NILINDICES', {name:x['nilpotence_index'] for name,x in out.items()},flush=True)

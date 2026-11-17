"""Global polynomial coefficient reconstruction for normalized tautological sweeps."""
import sympy as S
def gradient(f,ps):return S.Matrix([[S.diff(f,p) for p in ps]])
def divergence(g,ps):return S.expand(sum(S.diff(g[i],ps[i]) for i in range(len(ps))))
def homogeneous_parts(f,ps):
    out={}
    for mon,co in S.Poly(S.expand(f),*ps).terms():
        term=co*S.prod(p**j for p,j in zip(ps,mon));degree=sum(mon)
        out[degree]=out.get(degree,0)+term
    return out
def from_g(g,ps,constant=0):
    g=S.Matrix(g);ps=list(ps);div=divergence(g,ps)
    parts=homogeneous_parts(div,ps)
    if S.expand(parts.get(0,0))!=0:raise ValueError('linear g part must have trace zero')
    f=S.expand(constant+sum(value/S.Integer(n) for n,value in parts.items() if n>0))
    P=(g.jacobian(ps)-S.Matrix(ps)*gradient(f,ps)).applyfunc(S.expand)
    return dict(X=S.Matrix([f,*g]),f=f,g=g,P=P,df=gradient(f,ps))
def divergence_columns(M,ps):
    return S.Matrix([[S.expand(sum(S.diff(M[j,i],ps[j]) for j in range(len(ps))))
                      for i in range(len(ps))]])
def compatibility(P,ps):
    r=len(ps)
    if r<2:raise ValueError('matrix reconstruction needs r>=2')
    b=divergence_columns(P,ps);a=-b/S.Integer(r-1);curl=[];closed=[]
    for i in range(r):
     for j in range(r):
      for k in range(j+1,r):
       curl.append(S.expand((r-1)*(S.diff(P[i,k],ps[j])-S.diff(P[i,j],ps[k]))
                  -(b[k] if i==j else 0)+(b[j] if i==k else 0)))
    for j in range(r):
     for k in range(j+1,r):
      closed.append(S.expand(S.diff(a[k],ps[j])-S.diff(a[j],ps[k])))
    return dict(a=a,b=b,curl=curl,closed=closed,trace=S.trace(P))
def primitive_form(row,ps):
    answer=0
    for i,entry in enumerate(row):
     for mon,co in S.Poly(S.expand(entry),*ps).terms():
      answer+=co*ps[i]*S.prod(p**j for p,j in zip(ps,mon))/S.Integer(sum(mon)+1)
    return S.expand(answer)
def from_matrix(P,ps):
    ps=list(ps);data=compatibility(P,ps);f=primitive_form(data['a'],ps)
    g=S.Matrix([primitive_form(list(P[i,:]+ps[i]*data['a']),ps) for i in range(len(ps))])
    rebuilt=(g.jacobian(ps)-S.Matrix(ps)*gradient(f,ps)).applyfunc(S.expand)
    data.update(f=f,g=g,X=S.Matrix([f,*g]),rebuilt=rebuilt)
    return data
def coframe(f,P,ps):
    rows=[gradient(f,ps)]
    for i in range(1,len(ps)):rows.append((rows[-1]*P).applyfunc(S.expand))
    return S.Matrix.vstack(*rows)
def nilpotent_residuals(P):
    return [S.expand(S.trace(P**j)) for j in range(1,P.rows+1)]
def directional_matrix(P,ps,u):
    return sum((S.diff(P,p)*u[i] for i,p in enumerate(ps)),S.zeros(len(ps))).applyfunc(S.expand)
def nij_component(P,ps,i,j,k):
    return S.expand(sum(P[l,j]*S.diff(P[i,k],ps[l])-P[l,k]*S.diff(P[i,j],ps[l])
        -P[i,l]*(S.diff(P[l,k],ps[j])-S.diff(P[l,j],ps[k])) for l in range(len(ps))))
def haantjes_component(P,ps,i,j,k):
    r=len(ps);P2=P*P
    return S.expand(sum(P2[i,a]*nij_component(P,ps,a,j,k) for a in range(r))+
      sum(P[a,j]*P[b,k]*nij_component(P,ps,i,a,b) for a in range(r) for b in range(r))-
      sum(P[i,a]*P[b,k]*nij_component(P,ps,a,j,b)+P[i,a]*P[b,j]*nij_component(P,ps,a,b,k)
           for a in range(r) for b in range(r)))

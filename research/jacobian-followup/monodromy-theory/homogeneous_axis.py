"""Homogeneous-axis controls and exact exceptional11 subgroup quotients."""
import sympy as S
from psl11_curve_spectrum import closure
from local_image import fixed_sheets
def norm_control(r):
    T=S.Symbol('T');poly=S.Poly(S.expand((T**(r+1)-(T-1)**(r+1))/S.Integer(r+1)),T)
    co=poly.all_coeffs();M=S.zeros(r)
    for i in range(1,r):M[i,i-1]=1
    for i in range(r):M[i,r-1]=-co[-1-i]
    ps=list(S.symbols('p1:'+str(r+1)))
    Q=sum((ps[j]*M**j for j in range(r)),S.zeros(r))
    F=S.expand(Q.det(method='domain-ge'));f=-F
    g=(F*(M-S.eye(r))*S.Matrix(ps)).applyfunc(S.expand)
    return dict(r=r,T=T,d=poly.as_expr(),M=M,parameters=ps,Q=Q,F=F,f=f,g=g,
                X=S.Matrix([f,*g]))
def subgroups_inside(H,one,table):
    found={frozenset([one]):()};todo=list(found)
    while todo:
        D=todo.pop()
        for g in sorted(H-D):
            E=closure((*found[D],g),one,table)
            if E not in found:found[E]=(*found[D],g);todo.append(E)
    return found
def quotient_rows(H,one,table,inverses,permutations):
    groups=subgroups_inside(H,one,table);pairs=[];accepted=[]
    for D in sorted(groups,key=lambda G:(len(G),tuple(sorted(G)))):
     for I in sorted(groups,key=lambda G:(len(G),tuple(sorted(G)))):
      if not I<D:continue
      if any(table[table[g][i]][inverses[g]] not in I for g in D for i in I):continue
      gen=next((g for g in sorted(D-I) if closure((*groups[I],g),one,table)==D),None)
      if gen is None:continue
      n=len(D)//len(I);fixed=fixed_sheets(I,permutations);left=set(fixed);orbits=[]
      while left:
          a=min(left);orb={permutations[g][a] for g in D};orbits.append(len(orb));left-=orb
      row={'D':sorted(D),'I':sorted(I),'D_order':len(D),'I_order':len(I),'n':n,
           'quotient_generator':gen,'I_fixed_sheets':list(fixed),'D_orbits_on_I_fixed':sorted(orbits)}
      pairs.append(row)
      if len(fixed)>=n+1 and n in orbits:
          row=dict(row);row['allowed_actual_k']=[1+n*m for m in range(1,orbits.count(n)+1)]
          accepted.append(row)
    return groups,pairs,accepted

"""Symbolic controls for generic ramification-surface completeness."""
import sympy as S

def line_plucker(X,w,z):
    point=S.Matrix([1,*X]);direction=S.Matrix([0,1,w,z])
    return {(i,j):S.expand(point[i]*direction[j]-point[j]*direction[i])
            for i in range(4) for j in range(i+1,4)}

def ruled_data(f,r,s,t,v,w):
    fp=S.diff(f,t);rp=S.diff(r,t)
    A=S.expand(fp*S.diff(r,t,2)-S.diff(f,t,2)*rp)
    B=S.expand(fp*S.diff(s,t,2)-S.diff(f,t,2)*S.diff(s,t))
    X=S.Matrix([f,v,r*v+s])
    normal=S.Matrix([-(rp*v+S.diff(s,t))/fp,-r,1])
    return dict(X=X,normal=normal,A=A,B=B,
       asymptotic_numerator=A*v+B+2*rp*fp**2*w,
       p=S.cancel(r-2*rp**2*fp/A) if A!=0 else None,
       q=S.cancel(S.diff(s,t)/fp-rp*B/(fp*A)) if A!=0 else None)

def nonfinite_control(w,z,g,k):
    a=-k-2
    X=S.Matrix([a*w*z,(a+1)*w*w*z,-w*z*z])
    sweep=X+g*S.Matrix([1,w,z])
    return dict(a=a,X=X,S=sweep)

def binary_control(w,z,g):
    sig=w**3-z
    X=S.Matrix([sig,w**4/2,w**6/3-sig*sig/2])
    return dict(X=X,S=X+g*S.Matrix([1,w,z]))

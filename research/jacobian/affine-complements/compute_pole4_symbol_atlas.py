"""Complete repeated-root cubic-symbol atlas for the pole-(3,4) normal form.

Produces exact necessary pole-four coefficient equations. It does not
infer emptiness from a bounded parameter specialization.
"""
from pathlib import Path
from time import monotonic
import sympy as S,json,hashlib

HERE=Path(__file__).resolve().parent
s,t,w=S.symbols("s t w")
n0,n1,d0,d1,d2=S.symbols("n0 n1 d0 d1 d2")
t0,x0,y0,k,ell=S.symbols("t0 x0 y0 k ell")
D=d0*s*s+d1*s*t+d2*t*t

def equations(kind,L,M=None):
    # Choose an independent linear form K; nonproportionality is exactly n1!=0.
    K=t if L==s else s
    if L==s-t:
        K=t
    N=n0*L+n1*K
    if kind=="triple":
        f=[L**3,L*L*N,L*D,3*N*D-N**3]
    else:
        f=[L*L*M,L*M*N,S.Rational(1,3)*M*N*N+L*D,3*N*D]
    f=[S.expand(a.subs({s:w,t:1})) for a in f]
    p=[w*a for a in f]
    q=[S.Poly(a,w).coeff_monomial(w*w)*w
       +2*S.Poly(a,w).coeff_monomial(w**3)*w*w+c
       for a,c in zip(f,[t0,x0,y0,0])]
    def W3(i,j):return S.diff(p[i],w)*p[j]-p[i]*S.diff(p[j],w)
    def W2(i,j):return S.diff(q[i],w)*p[j]-p[i]*S.diff(q[j],w)
    U={(0,1):12*p[1]**2-12*p[0]*p[2],
       (0,2):-12*p[0]*p[1],(0,3):4*p[0]**2}
    V={(0,1):24*p[1]*q[1]-12*(p[0]*q[2]+q[0]*p[2])+6*k*p[1]-3*ell*p[0],
       (0,2):-12*(p[0]*q[1]+q[0]*p[1])-12*p[2]-6*k*p[0],
       (0,3):8*p[0]*q[0]+4*p[1],
       (1,2):-12*p[1],(1,3):4*p[0]}
    pole5=S.expand(sum(value*W3(*pair) for pair,value in U.items()))
    assert pole5==0
    pole4=S.Poly(S.expand(sum(value*W2(*pair) for pair,value in U.items())
                         +sum(value*W3(*pair) for pair,value in V.items())),w)
    coeff=[(degree[0],S.factor(value)) for degree,value in pole4.terms() if value!=0]
    return dict(kind=kind,L=str(L),M=None if M is None else str(M),N=str(N),
                leading_sections=f,coefficient_equations=coeff)

def main():
    start=monotonic()
    cases=[]
    for name,kind,L,M in [
        ("triple_at_marked","triple",s,None),
        ("triple_away_marked","triple",t,None),
        ("double_at_marked","double_simple",s,t),
        ("simple_at_marked","double_simple",t,s),
        ("neither_at_marked","double_simple",s-t,t)]:
        row=equations(kind,L,M)
        row["name"]=name
        row["leading_sections"]=list(map(str,row["leading_sections"]))
        row["coefficient_equations"]=[[degree,str(eq)] for degree,eq in row["coefficient_equations"]]
        cases.append(row)
        print(name,"equations",len(row["coefficient_equations"]))
        for degree,eq in row["coefficient_equations"]:
            if not S.sympify(eq).free_symbols.intersection({t0,x0,y0,k,ell}):
                print("  pure leading condition",degree,eq)
    payload=dict(as_of="2026-11-14",outcome="EQUATIONS_GENERATED",
                 nonproportionality_open="n1 != 0",cases=cases,
                 source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 elapsed_seconds=round(monotonic()-start,6))
    (HERE/"pole4_symbol_atlas.json").write_text(json.dumps(payload,indent=2)+"\n")
    print("seconds",payload["elapsed_seconds"])

if __name__=="__main__":main()


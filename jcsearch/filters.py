"""Exact candidate metrics; conservative flags are not presented as proofs."""
import sympy as sp
from sympy import isprime

x,y,s,t=sp.symbols("x y s t")

def generic_resultant_degree(P,Q,samples=((2,3),(5,7))):
    degrees=[]
    for a,b in samples:
        resultant=sp.resultant(P-a,Q-b,y)
        degrees.append(sp.Poly(resultant,x).degree() if resultant else None)
    return degrees[0] if degrees and len(set(degrees))==1 else None

def generic_degree_flags(P,Q):
    degree=generic_resultant_degree(P,Q)
    return {"generic_degree":degree,"birational":degree==1,
            "prime_degree":bool(degree and isprime(degree)),
            "below_borisov_six":bool(degree and degree<6)}

def collision_groebner(P,Q,p1=(0,0),p2=(1,0)):
    return all(sp.expand(f.subs({x:p1[0],y:p1[1]}))==0 and
               sp.expand(f.subs({x:p2[0],y:p2[1]}))==0 for f in (P,Q))

def normalized(P,Q):
    checks=[P.subs({x:0,y:0}),Q.subs({x:0,y:0}),
            sp.diff(P,x).subs({x:0,y:0})-1,sp.diff(P,y).subs({x:0,y:0}),
            sp.diff(Q,x).subs({x:0,y:0}),sp.diff(Q,y).subs({x:0,y:0})-1]
    return all(sp.expand(c)==0 for c in checks)

def line_restriction_metrics(P,Q):
    result={}
    for name,sub,var in (("y=0",{y:0},x),("x=0",{x:0},y)):
        p,q=sp.Poly(P.subs(sub),var),sp.Poly(Q.subs(sub),var)
        result[name]={"degrees":(p.degree(),q.degree()),"linear_coordinate":p.degree()==1 or q.degree()==1}
    return result

def newton_edges(f):
    support=[m for m,_ in sp.Poly(f,x,y).terms()]
    return {"support_size":len(support),"total_degree":max(map(sum,support),default=-sp.oo),
            "top_edge":[m for m in support if sum(m)==max(map(sum,support),default=-1)]}

def puiseux_leading_scan(P,Q,mmax=6,nmax=6):
    """One-term bad-branch screen x=s^-m, y=c*s^-n.

    This reports valuation candidates; it is not a full Newton-Puiseux expansion.
    """
    s,c=sp.symbols("s c")
    records=[]
    for m in range(1,mmax+1):
      for n in range(1,nmax+1):
        vals=[]
        for f in (P,Q):
          terms=sp.Add.make_args(sp.expand(f.subs({x:s**-m,y:c*s**-n})))
          exps=[int(term.as_powers_dict().get(s,0)) for term in terms]
          low=min(exps);lead=sp.expand(sum(term/s**low for term,e in zip(terms,exps) if e==low))
          vals.append((low,lead))
        if vals[0][0]<0 and vals[1][0]<0:
          common=sp.gcd(sp.Poly(vals[0][1],c),sp.Poly(vals[1][1],c))
          # The Puiseux leading coefficient c is required nonzero. Remove the
          # spurious component c=0 before reporting a common leading root.
          while common.degree()>0 and common.eval(0)==0:
              common=sp.div(common,sp.Poly(c,c))[0]
          if common.degree()>0: records.append({"m":m,"n":n,"valuations":(vals[0][0],vals[1][0]),"common_lead":sp.sstr(common.as_expr())})
    return records

def zero_fiber_collisions(P,Q):
    """Exact zero-fiber solve for moderate explicit candidates."""
    basis=sp.groebner([P,Q],x,y,order="lex")
    roots=[]
    if basis.is_zero_dimensional:
      try: roots=sp.solve_poly_system([P,Q],x,y) or []
      except NotImplementedError: roots=[]
    return {"zero_dimensional":basis.is_zero_dimensional,"known_roots":[tuple(map(sp.sstr,r)) for r in roots],
            "collision":len(roots)>1}

def candidate_metrics(P,Q):
    determinant=sp.factor(sp.Matrix([P,Q]).jacobian((x,y)).det())
    return {"jacobian":sp.sstr(determinant),"keller":bool(determinant.is_number and determinant!=0),
      "degree_flags":generic_degree_flags(P,Q),"zero_fiber":zero_fiber_collisions(P,Q),
      "lines":line_restriction_metrics(P,Q),"newton":{"P":newton_edges(P),"Q":newton_edges(Q)},
      "puiseux_leads":puiseux_leading_scan(P,Q)}

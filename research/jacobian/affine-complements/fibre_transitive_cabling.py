"""PF104: full-shear transfer predictions and exact small pattern controls.

Universal validity is proved in the source note. No returned necessary
source budget constructs a surface or certifies a Keller map.
"""
import math
import sympy as S
import full_peripheral_linking as N


def order(p):
    seen=set();lengths=[]
    for s in range(len(p)):
        if s in seen:continue
        j=s;length=0
        while j not in seen:seen.add(j);length+=1;j=p[j]
        lengths.append(length)
    return math.lcm(*lengths)


def fibre_transitive(e,f,a):
    if e<1 or f<1 or not 0<=a<e:
        raise ValueError("normalized joint peripheral profile required")
    return math.gcd(e,f-a)==1


def transfer_matrix(matrix,profiles,winding):
    F=S.Matrix(matrix)
    if F.rows!=F.cols or F.rows!=len(profiles):
        raise ValueError("one full peripheral profile per matrix row required")
    diagonal=[]
    for e,f,a in profiles:
        if not fibre_transitive(e,f,a):
            raise ValueError("regular fibre must act transitively on each orbit")
        c=S.Rational(a*(1-winding*winding),e)
        if not c.is_Integer:
            raise ValueError("stable profiles require integral shear correction")
        diagonal.append(c)
    out=winding*winding*F+S.diag(*diagonal)
    return [[str(out[i,j]) for j in range(out.cols)] for i in range(out.rows)]


def spin_budget(degree,known_pairs=()):
    """Necessary ONLY: actual smooth whole-plane source, ALL e odd, H1end=0.

    known_pairs is any list of DISTINCT omitted primes, not necessarily
    the entire list. Every nonautomorphic all-odd packet has a ramified
    prime of index at least three. All unknown costs/excess are nonnegative.
    """
    if not isinstance(degree,int) or degree<2:
        raise ValueError("nontrivial integral field degree required")
    if any(not isinstance(e,int) or not isinstance(f,int) or e<1 or f<1 or e%2==0
           for e,f in known_pairs):
        raise ValueError("positive pairs with odd transverse indices required")
    upper=min(degree-3,degree-1-sum(e*f-1 for e,f in known_pairs))
    lower=max(1,len(known_pairs));first=8*((lower+7)//8);last=8*(upper//8)
    return dict(maximum_omitted_primes=upper,first_possible_rank=first,
                last_possible_rank=last,necessary_rank_interval_nonempty=first<=last,
                required_rank_mod8=0,assumptions=["actual whole marked A2 source",
                "smooth finite normalization","all omitted indices odd","integral homology sphere end"],
                necessary_only=True)


def pattern_packet(e,f,a,r,k):
    """Exact integral H1/ports of one finite cyclic-or-abelian pattern cover.

    Does NOT assume fibre transitivity: retained failure controls need the
    extra torsion. Parameters must preserve the supplied profile.
    """
    if e<1 or f<1 or not 0<=a<e or r<=1 or k<=1 or math.gcd(r,k)!=1:
        raise ValueError("valid profile and genuine coprime ordinary cable required")
    count=e*f
    def idx(i,j):return i+e*j
    mu=tuple(idx((i+1)%e,j) for j in range(f) for i in range(e))
    lam=tuple(idx(i,j+1) if j+1<f else idx((i-a)%e,0) for j in range(f) for i in range(e))
    zperm=tuple(mu[lam[i]] for i in range(count))
    F=N.Flows([mu,lam,zperm]);m,l,z=F.generators
    u,v,g=S.gcdex(r,k);assert g==1;u=int(u);v=-int(v)
    h=F.product(F.power(m,k),F.power(l,r))
    mn=F.product(F.power(m,u),F.power(l,v),F.power(z,-v))
    ln=F.product(F.power(mn,-r*k),h)
    rel=[F.product(m,l,F.inv(m),F.inv(l)),F.product(F.power(z,r),F.inv(h))]
    if any(x[0]!=F.ident for x in rel) or mn[0]!=m[0] or ln[0]!=l[0]:
        raise ValueError("root relation or full peripheral action not preserved")
    par=list(range(count));tree=set()
    def root(i):
        while par[i]!=i:i=par[i]
        return i
    for g,A in enumerate(F.generators):
        for s in range(count):
            t=A[0][s]
            if root(s)!=root(t):par[root(s)]=root(t);tree.add(g*count+s)
    assert len(tree)==count-1
    keep=[j for j in range(3*count) if j not in tree]
    project=lambda row:[row[j] for j in keep]
    rows=[project(R[1][s]) for R in rel for s in range(count)]
    oldL=F.product(F.power(m,a),F.power(l,f));newL=F.product(F.power(mn,a),F.power(ln,f))
    ports=[project(F.power(m,e)[1][0]),project(oldL[1][0]),
           project(F.power(mn,e)[1][0]),project(newL[1][0])]
    c=S.Rational(a*(1-r*r),e);lattice=N.RowLattice(rows)
    differences=[[ports[0][i]-r*ports[2][i] for i in range(len(keep))],
                 [ports[3][i]-r*ports[1][i]-c*ports[2][i] for i in range(len(keep))]]
    return dict(profile=[e,f,a],parameters=[r,k],fibre_order=order(zperm),
                fibre_transitive=fibre_transitive(e,f,a),H1=N.reduced_h1(rows),
                quotient_by_proposed_basis=N.reduced_h1(rows+[ports[2],ports[1]]),
                shear_correction=str(c),port_difference_orders=[lattice.order(row) for row in differences],
                rows=rows,ports=ports,word_expansion_used=False)

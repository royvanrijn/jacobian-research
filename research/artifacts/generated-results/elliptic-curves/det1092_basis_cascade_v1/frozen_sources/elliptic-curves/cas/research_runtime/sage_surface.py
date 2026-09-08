"""Generic exact fibre arithmetic for separated quadratic surface twists.

Supports a nonisotrivial short elliptic surface with squarefree discriminant,
smooth infinity, and any even-degree squarefree twist divisor disjoint from
the discriminant. Other configurations use an explicit fibre-data replayer
with research_runtime.regulator.replay_reduction and remain UNKNOWN here.
"""

from .regulator import Surface, replay_reduction
from .store import FiniteFieldFacts, default_store


class UnsupportedFibres(ValueError):
    pass


def separated_twist_fibres(surface:Surface,prime:int):
    from sage.all import GF,PolynomialRing,QQ,ZZ
    if not ZZ(prime).is_prime(proof=True) or prime<5:
        raise ValueError("good reduction requires a prime >=5")
    rational=PolynomialRing(QQ,"t")
    A0,B0,d0=(rational(list(row)) for row in (surface.A,surface.B,surface.d))
    delta0=-16*(4*A0**3+27*B0**2)
    if delta0==0 or delta0.degree()<12 or delta0.degree()%12:
        raise UnsupportedFibres("need integral elliptic surface with smooth infinity")
    chi0=int(delta0.degree())//12
    degree=int(d0.degree())
    if A0.degree()>4*chi0 or B0.degree()>6*chi0 or degree%2:
        raise UnsupportedFibres("unverified infinity fibre")
    if not delta0.is_squarefree() or not d0.is_squarefree() or delta0.gcd(d0).degree()>0:
        raise UnsupportedFibres("twist and discriminant must be squarefree and disjoint")
    if (A0**3/delta0).derivative()==0:
        raise UnsupportedFibres("height gate requires nonisotriviality")
    field=GF(prime);ring=PolynomialRing(field,"t")
    A,B,d=(ring([field(c) for c in pol]) for pol in (A0,B0,d0))
    delta=-16*(4*A**3+27*B**2)
    if (A.degree(),B.degree(),d.degree(),delta.degree())!=(A0.degree(),B0.degree(),d0.degree(),delta0.degree()):
        raise ValueError("surface has degree loss at this prime")
    if not delta.is_squarefree() or not d.is_squarefree() or delta.gcd(d).degree()>0:
        raise ValueError("surface does not have the asserted good reduction")
    branches=[]
    for factor,multiplicity in d.factor():
        if multiplicity!=1:raise ArithmeticError("repeated twist branch")
        residue=field.extension(factor,"v") if factor.degree()>1 else field
        v=residue.gen() if factor.degree()>1 else -factor[0]/factor[1]
        cubic_ring=PolynomialRing(residue,"x");x=cubic_ring.gen()
        cubic=x**3+A(v)*x+B(v)
        factors=list(cubic.factor())
        if any(m!=1 for _,m in factors):raise ArithmeticError("singular branch cubic")
        degrees=[int(f.degree()) for f,_ in factors]
        tamagawa=1+degrees.count(1)
        if tamagawa not in (1,2,4):raise ArithmeticError("invalid I0* splitting")
        branches.append({"base_factor_coefficients":[str(c) for c in factor],
                        "place_degree":int(factor.degree()),"residual_factor_degrees":degrees,
                        "kodaira":"I0*","tamagawa":tamagawa})
    chi=chi0+degree//2
    tamagawa_product=1
    for row in branches:tamagawa_product*=row["tamagawa"]
    # Closed-place c_v appears ONCE in the BSD product, not once for each
    # geometric point above it. Degrees do count in Euler and root ranks.
    euler=int(delta.degree())+6*degree
    if euler!=12*chi:raise ArithmeticError("incomplete geometric fibre data")
    return {"surface_key":surface.key,"prime":prime,"chi":chi,
            "tamagawa_product":tamagawa_product,"branches":branches,
            "geometric_I1_count":int(delta.degree()),"geometric_I0star_count":degree,
            "infinity":"smooth","trivial_lattice_rank":2+4*degree,
            "expected_L_degree":int(delta.degree())+2*degree-4,
            "nonisotrivial":True,"height_specialization_isometry":True,
            "good_reduction":True}


class SurfaceProofEngine:
    def __init__(self,store=None):
        self.facts=FiniteFieldFacts(store or default_store())

    def fibres(self,surface,prime,*,discover=False):
        return self.facts.query(surface,prime,"surface-fibres",
            build=(lambda:separated_twist_fibres(surface,prime)) if discover else None,
            version="separated-even-twist-1")

    def reduction(self,surface,frobenius,*,verify_frobenius,discover=False):
        """Produce verified rank/regulator facts from an exact polynomial proof.

        The injected verifier checks the complete Frobenius witness, including
        its binding to surface and prime. No point/section/BNF search is invoked.
        Fibre arithmetic is cached in discovery and cheaply replayed on demand.
        """
        if frobenius.get("surface_key")!=surface.key:
            raise ValueError("Frobenius witness model mismatch")
        prime=frobenius["prime"]
        fibre=self.fibres(surface,prime,discover=discover)
        if fibre!=separated_twist_fibres(surface,prime):
            raise ArithmeticError("retained fibre proof does not replay")
        witness={**fibre,"frobenius_coefficients":frobenius["coefficients"],
                 "moments":frobenius.get("moments"),"frobenius_proof":frobenius}
        return replay_reduction(surface,witness,
            replay=lambda model,row:verify_frobenius(model,row["frobenius_proof"]))

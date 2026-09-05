"""Cached adapters for retained PARI callers using a polynomial variable.

Variable names are presentation details; the labelled cubic generator remains
theta. New curve-aware code should pass ArithmeticContext/TwoTorsionContext
directly. These adapters let older exact replayers share the same field facts.
"""
from .arithmetic import TwoTorsionContext

_arithmetic = None


def arithmetic():
    global _arithmetic
    if _arithmetic is None:
        from .sage_arithmetic import SageArithmetic
        _arithmetic = SageArithmetic()
    return _arithmetic


def identify(polynomial):
    from sage.all import pari
    polynomial = pari(polynomial)
    if int(pari.poldegree(polynomial)) != 3:
        raise ValueError("prepared cubic adapter received a different degree")
    context = TwoTorsionContext(tuple(str(pari.polcoef(polynomial, i)) for i in range(4)))
    return context, str(pari.variable(polynomial))


def rename(value, variable):
    from sage.all import pari
    if variable == "y":
        return value
    if value.type() == "t_POLMOD":
        # subst deliberately leaves the bound variable of a POLMOD alone.
        return pari.Mod(pari.subst(value.lift(), "y", pari(variable)),
                        pari.subst(value.mod(), "y", pari(variable)))
    if value.type() in ("t_VEC", "t_COL", "t_MAT"):
        entries = [rename(entry, variable) for entry in value]
        if value.type() == "t_COL":
            return pari.Col(entries)
        if value.type() == "t_MAT":
            return pari.Mat(entries)
        return pari(entries)
    return pari.subst(value, "y", pari(variable))


def prepared_nf(polynomial, factor_primes=(), *, discover=True):
    """Maximal order and reduced maps are built once, then restored exactly."""
    context, variable = identify(polynomial)
    adapter = arithmetic()
    row = adapter.field(context, factor_primes=factor_primes, discover=discover)
    if len(row["components"]) != 1 or row["components"][0]["degree"] != 3:
        raise ValueError("this legacy adapter requires an irreducible cubic field")
    return rename(adapter.nf(context), variable)


def prepared_prime_ideals(nf, prime, *, discover=True):
    """Prime decomposition in the caller's exact integral-basis frame.

    A retained legacy BNF may use a different maximal basis. Reusing ideal
    coordinates from a merely isomorphic order would silently corrupt local
    computations, so those frames receive their own cached decomposition.
    """
    from sage.all import pari
    context, variable = identify(nf.nf_get_pol())
    adapter = arithmetic()
    supplied_basis = tuple(str(pari.subst(z, variable, pari("y"))) for z in nf.nf_get_zk())
    try:
        canonical = adapter.nf(context)
    except FileNotFoundError:
        canonical = None
    if canonical is not None and supplied_basis == tuple(map(str, canonical.nf_get_zk())):
        return [rename(ideal, variable) for ideal in
                adapter.prime_ideals(context, int(prime), discover=discover)]
    identity = {"algebra": context.key, "integral_basis": supplied_basis,
                "variable": variable, "prime": int(prime)}
    value = adapter._fact("two-torsion/prime-ideals-in-retained-basis", identity,
        lambda: {"binary": adapter._blob(pari.idealprimedec(nf, int(prime)))}, discover)
    return list(adapter._restore(value["binary"]))


def prepared_polredabs(polynomial, factor_primes=(), *, discover=True):
    """Retained generator maps, plus a transported maximal order for the target."""
    from sage.all import pari
    context, variable = identify(polynomial)
    adapter = arithmetic()
    row = adapter.field(context, factor_primes=factor_primes, discover=discover)
    adapter.reduced_field(context, discover=discover)
    item = row["components"][0]
    return (rename(pari(item["reduced_polynomial"]), variable),
            rename(pari(item["original_generator_in_reduced"]), variable))


def prepared_bnf(nf, flag=1, tech=(), *, certify_flag=0, discover=True):
    """Explicit completeness or 2-primary-upper-bound BNF, certified once."""
    context, variable = identify(nf.nf_get_pol())
    value = arithmetic().bnf(context, flag=flag, tech=tech, certify_flag=certify_flag,
        requirement="two-primary-upper-bound" if certify_flag else "complete-selmer", discover=discover)
    return rename(value, variable)


def prepared_factor(value, *, discover=True):
    """Certified rational integer factors, with PARI's matrix interface."""
    from sage.all import pari
    factors = arithmetic().factor_integer(int(value), discover=discover)
    return pari.matrix(len(factors), 2, [entry for pair in factors for entry in pair])


def certified_bnf_checkpoint(path, expected_hash, *, discover=True):
    """Certify a retained binary once, in its owning PARI build and basis."""
    from pathlib import Path
    from hashlib import sha256
    from sage.all import pari
    path = Path(path)
    if sha256(path.read_bytes()).hexdigest() != expected_hash:
        raise ArithmeticError("BNF checkpoint hash mismatch")
    adapter = arithmetic()
    def build():
        bnf = pari.read(str(path))
        if not bool(pari.bnfcertify(bnf)):
            raise ArithmeticError("retained BNF failed unconditional certification")
        return {"binary": adapter._blob(bnf), "unconditionally_certified": True}
    row = adapter._fact("two-torsion/certified-bnf-checkpoint", expected_hash, build, discover)
    if row["unconditionally_certified"] is not True:
        raise ArithmeticError("BNF checkpoint lacks full certification")
    return adapter._restore(row["binary"])

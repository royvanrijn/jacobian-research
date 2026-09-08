from sage.all import *
from newfamily_rank11_common import build_family


def valuation_poly(F,q):
    if F == 0:
        return Infinity

    n = 0

    while True:
        h,r = F.quo_rem(q)

        if r:
            return n

        F = h
        n += 1


def build_finite_minimal_family(roots):

    D = build_family(roots)

    Araw = D["A"]
    Braw = D["B"]

    A = Araw
    B = Braw

    RT = A.parent()
    scale = RT(1)

    while True:

        G = gcd(A,B)

        if G.degree() <= 0:
            break

        changed = False

        for q,_ in factor(G):

            va = valuation_poly(A,q)
            vb = valuation_poly(B,q)

            k = min(
                va//4,
                vb//6,
            )

            if k <= 0:
                continue

            qk = q**k

            A = A.quo_rem(qk**4)[0]
            B = B.quo_rem(qk**6)[0]

            scale *= qk

            changed = True

        if not changed:
            break

    Delta = -16*(
        4*A**3
        + 27*B**2
    )

    return {
        **D,
        "Araw":Araw,
        "Braw":Braw,
        "Amin":A,
        "Bmin":B,
        "Deltamin":Delta,
        "finite_scale":scale,
    }

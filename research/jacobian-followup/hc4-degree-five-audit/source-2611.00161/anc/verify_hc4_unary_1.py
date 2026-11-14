# Verificateur cone unaire, partie 1 : A4 (partie de f4 sans x1) binaire non degeneree.
# Lancer : python3 verify_hc4_unary_1.py   (aucun argument)
exec(open(__file__.replace('verify_hc4_unary_1.py', 'unary_core.py')).read())
say("=" * 72)
say("Cone unaire, partie 1 : A4 quartique binaire en (x2, x3) a hessien non nul")
say("=" * 72)
Z4 = dict(ZERO_A4)
say("\n[U1] couche 9 : det Hess A4 = 0 (donc A4 cone, Gordan-Noether)")
# generique : A4 = x2**4 + l x2**2 x3**2 + x3**4, l**2 != 4 ; autres types par multiplicites de racines
TYPES = [("generique x2**4 + l x2**2 x3**2 + x3**4", dict(Z4, a14=1, a4=1, a11=None), ['a11**2-4']),
         ("x2**2 x3 (x2 + x3)", dict(Z4, a11=1, a13=1), []),
         ("x2**2 x3**2", dict(Z4, a11=1), []),
         ("x2**3 x3", dict(Z4, a13=1), [])]
for name, sub, nz in TYPES:
    sub = {k: v for k, v in sub.items() if v is not None}
    nzz = [eval(z.replace('**', '**'), {}, {'a11': RV['a11']}) for z in nz]
    say("\n  type", name)
    inc, fin = peel(sub, nzz)
    check("    peel initial coherent", not inc)
    leaf_incoherent("    a31 != 0 (normalise a31 = 1)", dict(fin, a31=1), nzz)
    leaf_dead_x4("    a31 = 0", dict(fin, a31=0), nzz, tot=(9, 8, 7), kmax=6)
bilan()

"""Symbolic identities and fail-closed mutation checks for projection certificates."""
import copy
import unittest
try:
    from sage.all import QQ, PolynomialRing, matrix, EllipticCurve
    HAVE_SAGE = True
except ImportError:
    HAVE_SAGE = False


@unittest.skipUnless(HAVE_SAGE, "requires Sage for exact polynomial arithmetic")
class ProjectionIdentities(unittest.TestCase):
    def test_generic_mumford_and_conic_determinant(self):
        C = PolynomialRing(QQ, ["A","p","q","u"])
        A,p,q,u = C.gens()
        K = C.fraction_field()
        A,p,q,u = map(K,(A,p,q,u))
        B = q*q-p**3-A*p
        R = PolynomialRing(K,["s","t","v"],order="lex")
        s,t,v = R.gens()
        kappa = 1+u*p+u*u*(A+p*p)
        Q1 = v*v+u*s*s-2*u*t-u*(u*p+2)
        Q2 = u*t*t+(1-u*p)*v*v+2*u*u*q*v-u*kappa
        ideal = R.ideal([Q1,Q2])
        Z = PolynomialRing(R,"z")
        z = Z.gen()
        D = 1+A*u*u+B*u**3
        g = u*D-(3*u+A*u**3)*z*z+3*u*z**4-u*z**6
        n = v*(u*p-1-t)-u*u*q
        remainder = (g-(s*v*z+n)**2) % (z*z-s*z+t)
        self.assertTrue(all(ideal.reduce(c) == 0 for c in remainder))
        G = matrix(K,[[u,0,0],[0,1-u*p,u*u*q],[0,u*u*q,-u*kappa]])
        self.assertEqual(G.det(),-u*u*D)
        anchor_x = v*v/(u*u)-(2-s*s+2*t)/u
        anchor_y = -n/(u*u)-v*(1+t)/(u*u)+v*anchor_x/u
        self.assertEqual(ideal.reduce(R(anchor_x-p)),0)
        self.assertEqual(ideal.reduce(R(anchor_y-q)),0)

    def test_norm_maps_at_other_parameters(self):
        import projection_fibres as p
        for u in (-3,-2,1,2,3):
            result = p.verify_geometry(QQ(-1),QQ(1),QQ(0),QQ(1),QQ(u))
            self.assertTrue(all(result.values()))

    def test_corrupt_cover_and_label_are_rejected(self):
        import projection_fibres as p
        import verify_projection_fibres as v
        data = p.r.read(p.OUTPUT)
        source = p.r.read(p.lc.INPUT)
        A,B = map(QQ,source["anchor"]["short_model_ainvariants"][3:])
        E = EllipticCurve([A,B])
        basis = [E(list(map(QQ,P))) for P in source["anchor"]["known_points_on_short_model"]]
        backend = v.RetainedBackend(v.ArithmeticContext.from_record(data["covers"][0]["context"]))
        protocol = p.r.read(p.PROTOCOL)
        mutations = [
            ("quartic",lambda c:c["quartic"].__setitem__(0,str(QQ(c["quartic"][0])+1))),
            ("square root",lambda c:c["cubic_invariant_over_beta_square_root"].__setitem__(0,"0")),
            ("anchor label",lambda c:c.__setitem__("anchor_mask",6)),
            ("conic map",lambda c:c["conic_parameter_matrix"][0].__setitem__(0,"0")),
        ]
        for name,mutate in mutations:
            row = copy.deepcopy(data["covers"][0])
            mutate(row)
            with self.subTest(name=name),self.assertRaises(AssertionError):
                v.verify_cover(row,backend,protocol,E,basis)


if __name__ == "__main__":
    unittest.main()


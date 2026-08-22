// Globally determine the rational points on the level-474 H3 base.
//
// This is an exact Magma certificate.  We work with the Q-isomorphic model
//
//   C : y^2 = -3*x^6 + 22*x^4 - 19*x^2 + 64,
//
// obtained from the published model by Y = 3*y.  A two-cover descent of C
// leaves three locally soluble classes.  The quartic factor over the cubic
// field below maps those classes to two elliptic covers.  Elliptic Chabauty
// at p = 41 proves that their only Q-rational x-images are respectively
// {0} and {-13/7, -1, 1, 13/7}.  The final substitutions recover both signs
// of y and certify that there are no Q-rational points at infinity.
//
// Required: Magma V2.29-1 or later, with the standard number-field and
// elliptic-Chabauty packages.  No GRH bounds or analytic-rank assumptions
// are used.

Q<x> := PolynomialRing(Rationals());
f := -3*x^6 + 22*x^4 - 19*x^2 + 64;
C := HyperellipticCurve(f);

Hk, AtoHk := TwoCoverDescent(C);
assert #Hk eq 3;
A<theta> := Domain(AtoHk);

L<a> := NumberField(-3*x^3 + 22*x^2 - 19*x + 64);
LX<X> := PolynomialRing(L);
h := X^2 - a;
g0 := -(Evaluate(f, X) div h)/3;
assert -3*h*g0 eq Evaluate(f, X);

LTHETA<THETA> := quo<LX | g0>;
j := hom<A -> LTHETA | THETA>;
gammas := {Norm(j(q @@ AtoHk)) : q in Hk};

gamma0 := (263403*a^2 - 123771*a + 818724)/4;
gamma1 := (45*a^2 - 21*a + 144)/16;
assert gammas eq {gamma0, gamma1};

// Return the rational P^1-images certified by elliptic Chabauty.  The
// success flag gives a subgroup of finite odd index; R = 1 therefore makes
// the Chabauty upper bound unconditional on the full Mordell--Weil group.
function CertifiedImages(gamma, base_x)
    E := HyperellipticCurve(gamma*g0);
    P0 := E![base_x, SquareRoot(Evaluate(gamma*g0, base_x))];
    Eprime, EtoEprime := EllipticCurve(E, P0);
    success, MWgrp, MWmap := PseudoMordellWeilGroup(Eprime);
    assert success;

    P1 := ProjectiveSpace(Rationals(), 1);
    EtoP1 := map<E -> P1 | [E.1, E.3]>;
    EprimeToP1 := Expand(Inverse(EtoEprime)*EtoP1);
    N, V, R, _ := Chabauty(MWmap, EprimeToP1, 41);
    assert N eq #V;
    assert R eq 1;
    return {
        EprimeToP1(MWmap(v))[1] / EprimeToP1(MWmap(v))[2] : v in V
    }, #V;
end function;

images0, count0 := CertifiedImages(gamma0, 0);
assert count0 eq 2;
assert images0 eq {0};

images1, count1 := CertifiedImages(gamma1, 1);
assert count1 eq 8;
assert images1 eq {
    -13/7, -1, 1, 13/7
};

// Every rational point of C belongs to one of the locally soluble two-cover
// classes above, so its x-coordinate is in this set.  Conversely, exact
// substitution gives the following complete affine fibre set.
assert not IsSquare(LeadingCoefficient(f));
points := {
    C![0, -8, 1], C![0, 8, 1],
    C![-1, -8, 1], C![-1, 8, 1],
    C![1, -8, 1], C![1, 8, 1],
    C![-13/7, -4016/343, 1], C![-13/7, 4016/343, 1],
    C![13/7, -4016/343, 1], C![13/7, 4016/343, 1]
};
assert #points eq 10;

print "H3GLOBAL|selmer_classes=3|elliptic_covers=2|prime=41";
print "H3GLOBAL|x_images=0,-13/7,-1,1,13/7";
print "H3GLOBAL|points_on_Y2=-27X6+198X4-171X2+576:";
print "H3GLOBAL|(-1,+-24),(0,+-24),(1,+-24),(-13/7,+-12048/343),(13/7,+-12048/343)";
print "H3GLOBAL|status=PASS_GLOBAL_H3_LEVEL474_RATIONAL_POINTS";

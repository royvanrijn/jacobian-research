"""Replay the icosahedral face-height identities and good-reduction counts."""

import ast
from math import gcd
import subprocess

import sympy as sp


def check_gram():
    coords = sp.symbols("a b c")
    field = sp.QQ.algebraic_field(sp.sqrt(5))
    a, b, c = (sp.Poly(x, *coords, domain=field) for x in coords)
    root5 = field.convert(sp.sqrt(5))
    k = field.convert((3 + sp.sqrt(5)) / 2)
    kinv = field.one / k
    q = a + b + c
    p = ((b - c.mul_ground(k)) * (a - b.mul_ground(k))
         * (a.mul_ground(k) - c)).mul_ground(field.one / (field.convert(2) + root5))

    def inner(left, right):
        return 4 * sum(coord * left.diff(symbol) * right.diff(symbol)
                       for coord, symbol in zip((a, b, c), coords))

    h = inner(p, p)
    r = h + 7*q**2*p
    face = ((q**2 - 4*(a*b + b*c + c*a))
            * (b.mul_ground(k) - c.mul_ground(kinv))
            * (a.mul_ground(k) - b.mul_ground(kinv))
            * (a.mul_ground(kinv) - c.mul_ground(k))).mul_ground(root5)
    assert r == face
    assert inner(p, r) == 3*q**2*r - 5*q**4*p - 40*q*p**2
    assert inner(r, r) == (
        (25*q**4 + 240*q*p)*r - 15*p*(q**3 + 8*p)**2
    )


def check_face_curves():
    p, r, w, e = sp.symbols("p r w e")
    h = r - 7*p
    f = sp.expand(
        -64*p**2 + 2240*p**3 - 21760*p**4 + 8640*p**5
        + 16*p*h + 80*p**2*h - 3600*p**3*h + 4*h**2 - 200*p*h**2 - 25*h**3
    )
    face = 100*r**2 - w*(25 + 240*p)*r + 15*w*p*(1+8*p)**2
    discriminant = sp.discriminant(face, r)
    assert sp.factor(sp.discriminant(discriminant, p)) == (
        -276480000000000*w**5*(48*w**2 + 33*w + 40)
    )
    half = sp.Rational(1, 2)
    specialized = face.subs(w, half)
    sextic = (233280000*p**6 + 119880000*p**5 + 21261800*p**4
              + 1456955*p**3 + 26342*p**2 - 277*p - 7)
    assert sp.expand(sp.resultant(f, specialized, r)
                     - 5000*p**2*(8*p+1)**2*sextic) == 0

    at_zero = {p: e, r: sp.Rational(3, 5)*e + sp.Rational(168, 25)*e**2}
    at_other = {p: e - sp.Rational(1, 8), r: 24*e**2}
    for location, expected in ((at_zero, -sp.Rational(64, 25)),
                                (at_other, -sp.Rational(25, 4))):
        model = sp.expand(specialized.subs(location, simultaneous=True))
        value = sp.expand(f.subs(location, simultaneous=True))
        assert all(model.coeff(e, i) == 0 for i in range(3))
        assert value.coeff(e, 0) == value.coeff(e, 1) == 0
        assert value.coeff(e, 2) == expected

    for prime in (11, 13, 17, 19):
        branch = sp.Poly(sextic, p, modulus=prime)
        elliptic = sp.Poly(4*discriminant.subs(w, half), p, modulus=prime)
        assert branch.degree() == 6 and elliptic.degree() == 3
        assert sp.gcd(branch, branch.diff()).degree() == 0
        assert sp.gcd(elliptic, elliptic.diff()).degree() == 0
        assert sp.gcd(branch, elliptic).degree() == 0
        assert branch.eval(0) and branch.eval(-pow(8, -1, prime))
        assert elliptic.eval(0) and elliptic.eval(-pow(8, -1, prime))


GP_COUNTS = r"""
value(x,r)={
    my(h=r-7*x);
    return(-64*x^2+2240*x^3-21760*x^4+8640*x^5+16*x*h
           +80*x^2*h-3600*x^3*h+4*h^2-200*x*h^2-25*h^3);
}
facecount(p,n)={
    my(z=ffgen(ffinit(p,n,'u),'u), Q=p^n, z1=z^0);
    my(c=0,x,disc,d,roots,r,f);
    for(j=0,Q-1,
        x=sum(i=0,n-1,(j\p^i%p)*z^i);
        disc=(25+240*x)^2/4-3000*x*(1+8*x)^2;
        if(issquare(disc,&d),
            roots=if(d==0,[d],[d,-d]);
            for(k=1,#roots,
                r=((25+240*x)/2+roots[k])/200;
                f=if(r==0 && x==0,-64*z1/25,
                     if(r==0 && x==-z1/8,-25*z1/4,value(x,r)));
                c+=if(f==0,1,2*issquare(f));
            );
        );
    );
    return(c+2*issquare(-2*z1));
}
forprime(p=11,19,print([p,vector(4,n,facecount(p,n))]));
"""


def jacobian_counts(prime, counts):
    sums = [prime**degree + 1 - count
            for degree, count in enumerate(counts, start=1)]
    coeff = [1] + [0]*8
    for degree in range(1, 5):
        scaled = -sum(sums[j-1] * coeff[degree-j]
                      for j in range(1, degree+1))
        coeff[degree], rem = divmod(scaled, degree)
        assert rem == 0
    for j in range(4):
        coeff[8-j] = prime**(4-j) * coeff[j]
    return sum(coeff), sum((-1)**j * term for j, term in enumerate(coeff))


def brute_prime_count(prime):
    squares = [0]*prime
    for root in range(prime):
        squares[root*root % prime] += 1
    other = -pow(8, -1, prime) % prime
    count = squares[-2 % prime]
    for p in range(prime):
        for r in range(prime):
            if (200*r*r - (25+240*p)*r + 15*p*(1+8*p)**2) % prime:
                continue
            h = r - 7*p
            value = (-64*p**2 + 2240*p**3 - 21760*p**4 + 8640*p**5
                     + 16*p*h + 80*p**2*h - 3600*p**3*h + 4*h**2
                     - 200*p*h**2 - 25*h**3)
            if r == 0 and p == 0:
                value = -64 * pow(25, -1, prime)
            elif r == 0 and p == other:
                value = -25 * pow(4, -1, prime)
            count += squares[value % prime]
    return count


def check_point_counts():
    gp = subprocess.run(
        ["gp", "-fq", "-s", "64000000"], input=GP_COUNTS, text=True,
        capture_output=True, check=True,
    )
    assert not gp.stderr.strip(), gp.stderr
    expected = {
        11: ([20, 156, 1304, 14544], 34580, 8180),
        13: ([12, 208, 2466, 28704], 28512, 35904),
        17: ([20, 326, 4814, 84506], 99600, 79800),
        19: ([20, 414, 6518, 130254], 137775, 142335),
    }
    divisors = {}
    for line in gp.stdout.splitlines():
        if not line.startswith("["):
            continue
        prime, counts = ast.literal_eval(line)
        assert counts[0] == brute_prime_count(prime)
        plus, minus = jacobian_counts(prime, counts)
        assert (counts, plus, minus) == expected[prime], (prime, counts, plus, minus)
        splits = pow(-2, (prime-1)//2, prime) == 1
        divisors[prime] = plus if splits else plus*minus
        print(prime, counts, plus, minus)
    assert set(divisors) == set(expected)
    assert gcd(*divisors.values()) == 1
    for excluded in divisors:
        assert any(divisor % excluded
                   for prime, divisor in divisors.items() if prime != excluded)


if __name__ == "__main__":
    check_gram()
    check_face_curves()
    check_point_counts()
    print("PASS: Gram identities, smooth normalization data and torsion bound.")

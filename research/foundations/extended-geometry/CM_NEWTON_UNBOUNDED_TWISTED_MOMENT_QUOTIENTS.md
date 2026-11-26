# Unbounded non-CM sphere moments from twisted cyclic quotients

A second deck action of order seven gives a family of quotients of the
same cyclic cover. One gives an exact congruence criterion for a sphere
moment and another proves that the primitive factor is non-CM.
Together the quotient signatures force nonordinary reduction at every
good nonsplit prime.

## Statement

Let \(M\ge5\) be odd with \(M\equiv5\pmod7\). Set
\[
 n=7M,\qquad e=M+1,\qquad f_0=6M,\qquad
 P(U)=(U-1)^e(U-c)^{f_0},
\]
where \(c\in\mathbb Q\) and \(v_2(c)=1\). Define \(C\) as the
smooth normalization of \(Y^n=P(z^7)\). Write \(J^{[n]}(C)\)
for the primitive Jacobian factor for the order-\(n\) action on \(Y\).

**Theorem.** The cover is connected. Its primitive factor has
dimension \(13\varphi(n)/2\) and is not potentially CM. At every
prime \(p\nmid n\) with \(c,1-c\) units, the factor is nonordinary
if \(p\not\equiv1\pmod n\).

Let \(F\) have highest literal weight \(j+2i=n\) with angular
polynomial \(z^{-7n}P(z^7)\) on that face. At a prime \(p\) as above,
let \(d\) be the order of \(p\) modulo \(n\). For each \(k\ge1\) put
\(q=p^{kd}\) and \(m=(q-1)/n\). If the coefficients are integral at
a place over \(p\) of their number field, \(q\mathcal A(F^m)\) is
integral and its residue is nonzero exactly when
\[
           i+M(i\bmod7)<n\quad\text{for every }i\in\langle p\rangle
\]
using least positive representatives modulo \(n\) and modulo seven.
For every fixed phase over characteristic zero, all sufficiently large
primes in these classes certify all the indicated moments. At least five
nonsplit classes always qualify, and exactly seven qualify at \(n=35\),
giving natural density \(7/24\) in that case. All compatible strict lower
terms are arbitrary. The phases are absolute nonpowers in the sphere
function field, and there is a compatible phase of least sphere degree
\(8n\).

## Cyclic quotients and the non-CM factor

Let \(U=z^7\). A root of \(z^7=1\) has exponent \(e\)
in the Kummer equation over \(z\), with \(\gcd(n,e)=1\).
Thus \(C\) is connected even geometrically with degree \(n\)
over the \(z\)-line. Its geometric deck group over the \(U\)-line
is \(\mu_n\times\mu_7\), acting independently on \(Y,z\).
For \(0\le k\le6\), the degree-seven quotient with \(Z_k=Yz^k\)
has model
\[
               D_k:\quad Z_k^n=U^{kM}P(U).
\]
The maps are defined over \(\mathbb Q\) and equivariant for
the \(n\)-th root action. Pullback on Jacobians is injective
up to isogeny because its composition with pushforward is
multiplication by seven. Thus each primitive factor of \(D_k\)
is an isogeny factor of \(J^{[n]}(C)\).

For \(D_1\) the branches at \(0,c\) have exponents \(M,6M\).
They collide at 2 and their sum is \(n\). Blow up the marked
base at the collision. At the node the base parameters
satisfy \(Us=c\), and
\[
                   (Z_1/U)^n=(U-1)^e(1-s)^{6M}.
\]
The right side is a unit and \(n\) is odd. The normalization
is therefore étale at the node with \(n\) geometric nodes.
The main component is the \(n\)-th root of \((U-1)^e\),
which is connected and rational with two fully ramified
points. On the exceptional component use \(u=U/c\) and
\(Z_1/c\). Its reduced right side is a nonzero constant times
\(u^M(u-1)^{6M}\). It has \(M\) geometric components, each a
degree-seven cover branched at two points and hence rational.
Away from the nodes the normalization is smooth by tame
ramification along the separated sections.

The graph has one main vertex, \(M\) bubble vertices and
\(n\) edges. The deck group acts regularly on the edges and
through its quotient of order \(M\) on the bubble vertices.
The vertex chain group contains no primitive order-\(n\)
characters, so graph homology contains each such character
once from the edge chain group. This homology is the character
lattice of the Jacobian torus. Every primitive character of
\(D_1\) sees its four branches, so its first-cohomology rank
is two and the primitive factor has dimension \(\varphi(n)\).
The factor is totally toric at 2. Its toric rank persists
after finite extension, precluding potentially good reduction.
CM abelian varieties and their isogeny factors have potentially
good reduction everywhere. This proves the non-CM assertion
for \(J^{[n]}(C)\).
[SGA 7 I, Exposé IX, (12.3.7) and §12.4](https://grothendiecksga.com/read/sga7/en/9-12.html)
· [Serre–Tate, Theorem 6](https://annals.math.princeton.edu/1968/88-3/p05).

The cover \(C\) has fourteen distinct finite branch points
and one at infinity. Its exponents are \(e,f_0,-7\bmod n\),
all nonzero. Every primitive character thus sees all fifteen
points and has first-cohomology rank thirteen, giving the
stated dimension.
[Bouw, Lemma 4.3](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/6D1BEC9200F6849339665AC4371AB7AF/S0010437X01000343a.pdf).

## Nonordinary reduction at nonsplit primes

Coherent Frobenius on \(D_k\)
sends character \(i\) to \(pi\bmod n\). For a primitive label
put \(x=i/n\), \(t=i\bmod7\) with \(1\le t\le6\), and
\(\epsilon_j(x)=\mathbf 1_{x+j/7>1}\) for \(0\le j\le6\).
The branch exponents \(kM,e,f_0,-kM-e-f_0\) give dimension
\[
        h_k(i)=1-\epsilon_t(x)+\epsilon_{kt\bmod7}(x).
\]
Indeed for \(v=(kt\bmod7)/7\), their fractional parts are
\(v\), \(x+t/7-\epsilon_t(x)\), \(1-t/7\), and
\(1-x-v+\epsilon_{kt\bmod7}(x)\). Their sum minus one gives
the dimension. An endpoint equality would make \(i\) a
multiple of \(M\), impossible for a primitive label.
[Bouw, Lemma 4.3 and §5](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/6D1BEC9200F6849339665AC4371AB7AF/S0010437X01000343a.pdf).

For a primitive \(i\) the sum of \(h_k(i)-h_0(i)\) over \(0\le k\le6\)
is \(\lfloor i/M\rfloor\), since \(t\) is a unit modulo seven and the
number of residues \(j\) with \(x+j/7>1\) equals that floor. If the
primitive factor of \(C\) were ordinary, each primitive quotient factor
would be ordinary. Its coherent dimensions would be constant on
Frobenius orbits, so multiplication by \(a=p\bmod n\) would preserve
\(\lfloor i/M\rfloor\). At \(i=1\) this gives \(a<M\). For \(a>1\)
choose the first power of two \(j\) for which \(aj\ge M\). Because
\(aj<2M\) and \(j<M\), the two floors differ. Both labels are primitive
since \(n\) is odd. Thus the only possible ordinary class is \(a=1\).

## The nonzero moment

In the sphere chart \(2xy+t^2=1\), use \(z=x/t^2\) for the
weight \(j+2i\) on \(x^it^j\). At a good prime let \(d\) be the order
of \(p\) modulo \(n\). First put \(q=p^d\) and \(m=(q-1)/n\).
Angular-constant heights in \(F^m\) are at most \(nm=q-1\).
The denominator \(q\) is the only interval denominator of
valuation \(d\), and its height is
even. At that height only the prescribed top face contributes.
Writing \(\mathcal A\) for the normalized sphere functional gives
\[
    q\mathcal A(F^m)\equiv
          [z^{7nm}]P(z^7)^m=[U^{q-1}]P(U)^m\pmod p .
\]
For characteristic-zero coefficients in an arbitrary field,
take one number-field specialization of the finitely
generated coefficient algebra, keeping the rational face
fixed. Its coefficients are integral away from finitely
many primes and its nonzero moments certify the original ones.
[Sphere pole calculation](CM_NEWTON_SELECTED_EIGENLINES.md#a-moment-detecting-a-transition-rather-than-a-return).

The selected coefficient is the \(q\)-Frobenius return on
the coherent character-one line of \(D_0\). Its class \(Y/U\)
pulls back to \(W/z\), where \(W=Y/z^6\) is the equation
adapted to the angular face. For direct nonvanishing
translate \(U=V+1\). For \(q\le d\le 2q-2\), the identity
\((V+1)^d=(V^q+1)(V+1)^{d-q}\) in characteristic \(p\)
gives coefficient zero at \(q-1\). Since
\(\deg P^m=(n+1)m=q-1+m\), this proves
\[
      [U^{q-1}]P(U)^m
         \equiv (1-c)^m\binom{6Mm}{m}\pmod p .
\]
The base-\(p\) digits of \(m\) and \(6Mm\) come from long division of
\(1/n\) and \(6/7\), respectively. At label \(i=p^j\bmod n\), the
two normalized remainders are \(i/n\) and \(1-(i\bmod7)/7\).
Under the orbit criterion the second is larger, so its digit is at
least the corresponding digit of \(m\). The base-\(p\) factorization
of \((1+X)^{6Mm}\) makes the binomial coefficient a unit. If the
criterion fails, the corresponding \(D_0\) dimension \(h_0(i)\) is
zero and Frobenius has zero return. This is the classical three-point
orbit condition for an endpoint slope, now expressed as the residue
of the specified sphere moment.
[Li–Mantovan–Pries–Tang, Theorem 3.2](https://arxiv.org/abs/1805.04598).
In either case the residue of the \(kd\)-fold return is the \(k\)-th
power of the \(d\)-fold coefficient, which lies in \(\mathbb F_p\).
Denominator isolation proves the corresponding congruence at order
\((p^{kd}-1)/n\).

The order-six subgroup generated by \(M-1\) consists of
\(1,M-1,1+3M,6M-1,1+2M,2M-1\). All six labels meet the inequality,
so its five nonidentity classes certify at nonordinary primes.
Checking the orbit criterion over the units modulo 35 gives
\(1,4,8,9,11,16,22,29\), of which seven are nonsplit. The prime
number theorem in arithmetic progressions gives the densities.

For an explicit regular phase, lift each angular term of
angle \(-v<0\) by the same coefficient times
\((-2)^v y^v t^n\). Lift angle zero by \(t^n\) and
angle seven by \(x^7t^{n-14}\), with their prescribed
coefficients. These give the required face since
\(2xy=1-t^2\). The term \(t^ny^{7n}\) has nonzero coefficient
and degree \(8n\). Eliminating \(xy\) never raises degree,
and the resulting monomials \(x^it^j\), \(y^vt^j\)
with \(i\ge0,v>0\) form a basis. Thus this lift has least
sphere degree \(8n\).

In the rational function field with coordinates \(z,t\),
the order of \(F\) at \(t=\infty\) is \(-n\). A proper
power would have exponent dividing \(n\) and would make
the leading angular function \(z^{-7n}P(z^7)\) a power of
the same order. At a root of \(z^7=1\) that function has
multiplicity \(e\), coprime to \(n\). This rules out proper
powers even over an algebraic closure of the coefficient field.

The exact replay computes the coefficient before translation by
expanding the digit blocks and carrying their exponents. At \(n=35\)
it checks all 24 unit classes and obtains a quotient signature mismatch
in each nonidentity class. At \(c=2\), its direct and translated
calculations agree at \((p,d)=(29,2),(43,4),(79,6),(3,12)\), with
residues \(1,14,67,0\), respectively. The order \((29^4-1)/35=20208\)
also gives residue one. The original order-three examples and their
inverse classes are included in the replay.
[Exact modular replay](verify_cm_newton_unbounded_twisted_moment_quotients.py).

# Rank-one binary quadratics with unequal linear coefficients: the six-window

**Theorem 1.** Let K be any field of characteristic zero. If a,d,e,k are
nonzero elements of K and d!=e, then

    f(x,y)=a(x+y)^2+d x+e y+k

has a nonzero factorial moment among every six consecutive positive
indices. In detail, with Lambda(x^r y^s)=r!s! and M_j=Lambda(f^j), for
every integer n>=1 at least one of M_n,...,M_(n+5) is nonzero. The phase
has exactly six monomials, so this proves its Strong Factorial window.
This is a restricted binary quadratic theorem, not a proof for every
binary quadratic, full SFC(2), or SFC(3).

The argument below is self-contained. It uses both actual boundary
moment sequences, rather than arbitrary solutions of their recurrence.
No bounded computation, generic-rank assumption, irreducibility theorem,
or unproved prime-existence assertion is a premise.

## 1. Native identities

Write Lambda_1(X^r)=r! and

    H_j(a,t,k)=Lambda_1((aX^2+tX+k)^j),  H_0=1.

For every j>=1, omitting the H_(-1) term at j=1, one has the polynomial
identity

    H_j=2a j(2j-1) H_(j-1)
        +j(j-1)(t^2-4ak) H_(j-2)+k^j+t j k^(j-1).       (1)

In particular the forcing terms have plus signs. Here and below a term
with a zero integer prefactor is omitted rather than interpreted as an
undefined negative-index moment.

For completeness, let E(T)=sum_j H_j T^j/j!, and let J(T) and L(T)
be the corresponding series with an extra factor X and X^2 inside
Lambda_1. All series are formal. The monomial Stein identity
Lambda_1(P')=Lambda_1(P)-P(0) gives

    2a T J+t T E=E-exp(kT),
    E+T(2a L+t J)=J,  E'=a L+t J+kE.

Eliminating J,L without division by a yields

    4a T^2 E'=[1-2aT+(4ak-t^2)T^2]E-(1+tT)exp(kT).

Coefficient comparison proves (1) for all a,t,k.

The two-variable monomial Stein identity also gives, for every r>=0,

    (r+1)(d-e)M_r=H_(r+1)(a,d,k)-H_(r+1)(a,e,k).        (2)

Indeed f_x-f_y=d-e and the difference of the x- and y-boundary terms
is the displayed right-hand side. Thus (2) is a polynomial identity in
the coefficients over Q; it does not require an analytic integral.

## 2. What an actual six-zero block would force

Suppose, for a contradiction, that M_n,...,M_(n+5) all vanish. Put
N=n+4, so N>=5, and S=d+e. Equation (2) gives common boundary values

    R_l=H_l(a,d,k)=H_l(a,e,k),  l=n+1,...,n+6.

Subtract (1) for the two linear coefficients at index l+2. For
l=n+1,...,n+4 this yields

    (l+1) S R_l=-k^(l+1).                                (3)

Since k!=0, this first proves S!=0. Substitute (3) into (1) at
j=n+3 and j=n+4. After multiplying the nonzero rational denominators,

    P(j)=j(j+1)de+(j+1)k(2a+S)+k^2=0
                  at j=N-1,N.

The difference of these equations, followed by substitution, gives

    k(2a+S)=-2Nde,    k^2=N(N+1)de.                       (4)

Scaling f by k^(-1) preserves its zero moments. From now on all
coefficients mean the normalized coefficients, so k=1. Put

    P=1/[N(N+1)].

The necessary equations become

    de=P,  S=d+e!=0,  a=-S/2-1/(N+1),                    (5)
    H_j(a,d,1)=H_j(a,e,1)=-1/[S(j+1)]
                              for j=N-3,...,N.          (6)

Only the first two indices in (6) will be needed. In particular we
exclude the actual simultaneous boundary equations. We do not assert
that one arbitrarily prescribed boundary sequence alone is impossible.

## 3. Reduction of the normalized coefficients to a number field

There is a small algebraicity point before using p-adic valuations.
For the fixed integer N parameterize (5) over Q by a nonzero variable D:

    e(D)=P/D,  S(D)=D+P/D,
    a(D)=-S(D)/2-1/(N+1).

Choose the even integer j among N-3,N-2 and consider the rational
function

    B_j(D)=H_j(a(D),D,1)+1/[S(D)(j+1)].                  (7)

This is not identically zero. To see this, evaluate at the rational
point D=1. Then e=P>0, S=1+P>0, and a=-(1+P)/2-1/(N+1)<0. The factorial
integral identity for a real polynomial gives

    H_j(a,1,1)=integral_0^infinity (aX^2+X+1)^j exp(-X)dX>0.

Here j is even and the integrand is strictly positive near X=0.
On the other hand the required value -1/[S(j+1)] is negative. This
positivity is used only to exhibit a rational nonidentity witness;
no realness is assumed for a hypothetical counterexample.

Clearing the D and S(D) denominators in (7) therefore gives a nonzero
polynomial in Q[D]. At a hypothetical solution of (5)--(6), neither
denominator vanishes and that polynomial vanishes at d. Hence d is
algebraic over Q; so are e and a by (5). Their field is a number field.
For any rational prime p it embeds into an algebraic closure of Q_p.
An embedding preserves all displayed algebraic equalities and nonzero
conditions. The p-adic inequalities are deduced after choosing the embedding.
We may consequently work in a finite valued extension of Q_p, with the
normalized valuation v(p)=1.

This also handles arbitrary original fields K: only the normalized
coefficients are forced algebraic. The original scale k may remain
transcendental, which does not affect any zero moment.

## 4. A literal normalized native expansion

Let d denote either one of the two linear coefficients, and put

    A=a/d,  zeta=1/d.

Direct multinomial expansion gives, for every j>=0,

    H_j(a,d,1)/(j! d^j)
      =sum_(i,c>=0,i+c<=j) [(2i)!/i!]
         binom(j+i-c,2i) A^i zeta^c/c!.                  (8)

No division by a or by a moment occurs. We will choose the boundary
whose linear coefficient has the smaller p-adic valuation. After that
choice every term needed to reduce (8) modulo the maximal ideal is
integral; the terms with c>=1 will be strictly positive in valuation.

## 5. Any odd prime divisor of N is impossible

Suppose an odd prime p divides N, and write t=v(N)>=1. Since de=P and
N+1 is a p-unit, choose d from the two boundary coefficients so that

    alpha=v(d)=min(v(d),v(e))<=-t/2<0.

Write s=v(S); the ultrametric inequality gives s>=alpha. By (5), since
2 and N+1 are p-units,

    v(a)>=min(s,0)>=alpha,  hence v(A)>=0.               (9)

For every c>=1, Legendre's formula and -alpha>=t/2 show

    v(zeta^c/c!)>=ct/2-(c-s_p(c))/(p-1)>0,              (10)

where s_p is the base-p digit sum. Strictness includes p=3,t=1,
when the right-hand side equals s_3(c)/2>0. Thus (8) reduces to its
c=0 terms. Their coefficients are the integers

    L_(j,i)=(j+i)!/[i!(j-i)!]
           =[(2i)!/i!] binom(j+i,2i).

For every i>=ceil(p/2), the product (i+1)...(2i) contains a multiple
of p. For i<p this is the multiple p itself; for i>=p the interval
has length at least p. All those coefficients vanish modulo p.
For the remaining smaller i, denominators i! are p-units. At j=N-3,
every i>=3 in that smaller range has the factor N in the factorial
interval; at j=N-2 the same holds for every i>=2.

It follows that, in the residue field,

    H_(N-3)/[(N-3)! d^(N-3)] = 1+6A+12A^2,
    H_(N-2)/[(N-2)! d^(N-2)] = 1+2A,                    (11)

when p>=5. The first polynomial takes the value 1 at A=-1/2, so these
two residue polynomials have no common root. For p=3 the first residue
polynomial is simply 1, and the second is again 1+2A. Thus in every
case at least one of the two normalized moments is a unit.

For j=N-3,N-2 define its baseline valuation

    B_j=v(j!)+j alpha.

Both are strictly negative, since

    v(j!)=(j-s_p(j))/(p-1)<j/(p-1)<=jt/2<=-j alpha.

At least one actual H_j has exactly its baseline valuation. But both
required targets in (6) have valuation -s, because N-2 and N-1 are
p-units. Therefore s>0. Now (5) says that a is a p-unit: the term
-S/2 has positive valuation and -1/(N+1) is a unit. In particular
v(A)=-alpha>0. Both residues in (11) are then 1, so both moments have
exact baseline valuations. Their difference is

    B_(N-2)-B_(N-3)=v(N-2)+alpha=alpha<0,

whereas their required targets have the same valuation -s. This is a
contradiction. Hence N has no odd prime divisor.

## 6. The remaining powers of two are impossible

The remaining possibility is N=2^t with t>=3. Again choose the smaller
linear-coefficient valuation and write

    alpha=v_2(d)<=-t/2<=-3/2,  s=v_2(S)>=alpha,
    T=2a/d=-S/d-2/[d(N+1)],  zeta=1/d.

These formulas imply v_2(T)>=0 and v_2(zeta)>=3/2. Rewrite (8) as

    H_j(a,d,1)/(j!d^j)
      =sum_(i,c>=0,i+c<=j) (2i-1)!!
        binom(j+i-c,2i) T^i zeta^c/c!,                  (12)

with (-1)!!=1. The odd double factorial is a unit. For c>=1,

    v_2(zeta^c/c!)>=3c/2-(c-s_2(c))>0.

Consequently the residue of (12) is

    F_j(T)=sum_(i=0)^j binom(j+i,2i) T^i  in F_2[T].

These polynomials satisfy

    F_0=1, F_1=1+T, F_j=T F_(j-1)+F_(j-2) (j>=2).       (13)

For example (13) follows from their generating function
(1+u)/(1+Tu+u^2) over F_2. The recurrence shows that consecutive F_j
are coprime: a common zero would propagate back to F_0=1.
Thus at least one normalized moment at j=N-3,N-2 is a unit.

Both baselines B_j=v_2(j!)+j alpha are negative. More specifically
B_(N-3)<-1: j=N-3>=5 and

    B_j<=j-1-3j/2=-j/2-1.

The two target valuations in (6), in this order, are -s-1 and -s,
because v_2(N-2)=1 and v_2(N-1)=0. If the first normalized moment is
a unit, its exact baseline below -1 forces s>0; if the second is a
unit, its negative baseline also forces s>0. With s>0 the two terms
in the formula for T both have strictly positive valuation. Therefore
T has positive valuation and all F_j(T) reduce to F_j(0)=1. Both
moments consequently have exact baseline valuations.

Their baseline difference is

    B_(N-2)-B_(N-3)=v_2(N-2)+alpha=1+alpha<0,

whereas the target-valuation difference is +1. This is impossible.
Thus no power of two N>=8 can occur either.

## 7. Conclusion and scope

Every integer N>=5 either has an odd prime divisor or is a power of
two at least eight. Sections 5 and 6 therefore contradict all possible
critical systems (5)--(6). The six-zero block assumed in Section 2
cannot exist, proving Theorem 1 over every characteristic-zero field.

The result concerns rank-one homogeneous quadratic part a(x+y)^2,
nonzero a,d,e,k, and unequal d,e. The argument uses only two of the
four actual common boundary values, but it uses both boundary
sequences to select the smaller valuation. It does not establish an
unqualified one-boundary exclusion. The aligned case d=e, sparse
strata, general rank-two quadratic parts, full binary quadratic SFC,
and SFC in three variables are not conclusions of this note.


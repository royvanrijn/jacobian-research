# Central governing fields cannot create the missing jump classes

The degree-192 fields in the frozen sixteen-fibre governing panel contain
**no additional cubic Kummer class beyond their two generic inputs**.
Moreover, none of their three nonzero input combinations is strict on
**any** of the sixteen fibres. All 48 exclusions have explicit retained
local witnesses, including the rows whose full factorization is incomplete.

There is a general closure theorem behind this result. Adjoining all
pair-governing fields of the entire generic subgroup, or further central
obstruction extensions, still cannot add a new S4 Kummer class. The
additional central bits carry obstruction information; a new incidence
direction requires a new standard S3 module in the extension kernel.

This sharpens the earlier observation that governing degree did not
correlate with rank. The retained construction has **zero capacity for
new classes**, regardless of how its central CT data vary.

Three direct finite-prime witnesses verify the separation on the matched
103b2 high/low pair and a historic control. These use the third **generic**
section, not an exceptional direction:

| Frozen case | Witness prime | First-pair governing octic | Third generic class quartic |
|---|---:|---|---|
| 103b2 high, 3726/881 | 2551 | eight linear factors | two irreducible quadratics |
| 103b2 low, -1049/2296 | 971 | eight linear factors | two irreducible quadratics |
| ICARM356 | 853 | eight linear factors | two irreducible quadratics |

Thus even the third inherited class field is absent from the first-pair
governing field. The high/low difference requires additional arithmetic,
not a larger interpretation of those same two inputs.

## The relevant class-field capacity

Let V=E[2] be the standard two-dimensional F2 representation of S3,
and M the cubic splitting field. For a subspace W of H1(Q,V), let L_W
be the compositum over M of the normal closures of its S4 class fields.
If dim W=k, then

\[
 \operatorname{Gal}(L_W/M)\simeq V^k,
 \qquad [L_W:\mathbb Q]=6\cdot4^k.
\tag{1}
\]

This applies to arbitrary norm-square Kummer classes, before imposing
Selmer or strict conditions. Restriction to G_M is injective because
H1(S3,V)=0. Independence of the classes makes the joint image all of V^k:
a proper S3 submodule would be annihilated by a nonzero map to V,
contradicting independence. Here V is simple and End_S3(V)=F2.
These finite representation facts are independently checked in the
[standard S3 calculation](STRICT_CLASS_CREATION_IS_A_STANDARD_S3_BLOCK.md).

The retained pair field N contains L_W of degree96 and has degree192.
A third independent class would require a joint class field of degree384
inside N, which is impossible. Hence the only nonzero cubic Kummer
classes with normal closure contained in N are alpha, beta and alpha+beta.
All three are inherited rational classes. There can be other subfields;
the assertion concerns S4 classes with this specified cubic resolvent.

## Why central obstruction data cannot enlarge incidence

Here is a more general statement which does not rely on the number192.
Suppose Gamma acts on V through the full S3 quotient, and C is a central
subgroup acting trivially on V. Then

\[
 \boxed{H^1(\Gamma/C,V)\ \xrightarrow{\mathrm{inflation}}\
        H^1(\Gamma,V)\text{ is an isomorphism}.}
\tag{2}
\]

For a cocycle a and c∈C, its restriction to C is a homomorphism. The
cocycle identity gives

    a(g*c*g^(-1)) = g*a(c),

because C acts trivially on V. Centrality makes the left side a(c), so
a(c) belongs to V^S3=0. Therefore every cocycle vanishes on C and descends
to Gamma/C. Coboundaries descend as well. This proves (2), including
surjectivity; no arithmetic solubility assumption is used.

For the pair-governing group, write an element as (a,b,g,z), with
 a,b∈V, g∈S3 and z∈F2. Its multiplication is

    (a,b,g,z)*(a',b',g',z')
      =(a+g*a', b+g*b', gg', z+z'+e2(a,g*b')).

The bit z is central. The kernel over M has order32, commutator and
Frattini subgroup of order2, and elementary abelian quotient V⊕V.
An exact cochain calculation gives

    dim Z1=4,       dim B1=2,       dim H1=2.

Representatives of H1 are exactly the two coordinate cocycles a and b.
All cocycles vanish on the central bit. The worker solves the twelve
binary generator-value variables through all Cayley-graph constraints.
The independent verifier instead computes the kernel's commutators,
enumerates its 256 linear maps to V, and checks S3 equivariance. Exactly
four maps survive. An independent enumeration of all4096 S3 cochains
checks H1(S3,V)=0.

Thus the concrete finite-group calculation and the degree argument
independently give the same class capacity.

## The full generic governing compositum is also closed

Let G be the m-dimensional marked generic subgroup. Let L_G be its joint
S4 class field, and N_G the compositum of all pair-governing fields for
a generic basis. This is a mathematical compositum; the enormous full
field is not constructed numerically here.

Each pair field is a central C2 extension of its corresponding joint
class field. Consequently Gal(N_G/L_G) embeds in the product of those
central C2 groups and is central in Gal(N_G/Q): conjugation is trivial
in every factor, and restriction to the factors is injective. Equation
(2) therefore gives

\[
 \boxed{\{\text{cubic Kummer classes split by }N_G\}=G.}
\tag{3}
\]

Here a class is split by N_G when its cocycle restricts to zero there,
equivalently its S4 normal closure is contained in N_G. Arbitrary further
central extensions over Q preserve (3), as do their finite composita.
This includes additional scalar quadratic choices in governing lifts.
The theorem does not cover an extension kernel on which S3 acts through
a new copy of V; that is exactly the missing incidence structure.

In particular the strict part contained in N_G is exactly G∩U. Collecting
more CT bits for the same generic subgroup cannot produce classes in U
outside G.

## An additional strict block needs a disjoint standard layer

There is also a relative-degree consequence. Suppose r independent
strict classes modulo G have been independently constructed. Their joint
class field with G has

    Gal(L_(G+W)/L_G)=V^r.

Its intersection with N_G is exactly L_G. Indeed, the intersection would
give a quotient of V^r which also comes from the central kernel of N_G.
That quotient has trivial S3 action, whereas V^r has no nonzero trivial
quotient. Thus

\[
 \boxed{[N_G L_{G+W}:N_G]=4^r.}
\tag{4}
\]

The new extension is unramified above M and split above the strict places,
by the strict-class correspondence; those properties persist after base
change to N_G. Over Q its new kernel has standard S3 action, rather than
the trivial action of central CT bits.

For context, joining the **already retained** rank lower bounds to the
[point-independent strict-boundary bounds](FRESH_STRICT_BLOCK_NECESSITIES.md)
forces the following necessities. These are deductions from labels, not
new independent class constructions:

| Frozen high | Additional strict dimension beyond G required at least r | New relative degree over N_G required at least 4^r |
|---|---:|---:|
| 103b2, 3726/881 | 9 | 262144 |
| 11952, -2448/11 | 5 | 1024 |
| 11952, 110314/102227, retained local bound27 | 6 | 4096 |
| 11952, 2828/2015 | 8 | 65536 |
| ICARM356 | 11 | 4194304 |
| ICARM385 | 4 | 256 |
| ICARM398 | 5 | 1024 |

These are relative to **G and its governing compositum**, not to the
larger global root-curve pool used in the separate capacity theorem.
The later public bound28 for inventory188 is not substituted into the
frozen panel or used as an arithmetic input here.

Equation (4) identifies the kind of extension that must appear. It does
not explain which t creates it or construct the required classes.

## Exact strict exclusions on all sixteen fibres

For each selected pair we evaluate all three words alpha, beta,
alpha+beta using the retained local squareclass signatures. The table
lists a witness place for each word, in that order. Every signature is
nonzero there, and each place belongs to the verified strict set S.

| Token | Frozen fibre | Witness places for alpha, beta, alpha+beta | Strict dimension in all-generic governing compositum |
|---|---|---|---:|
| 00 | 074d9 high, 2818/1535 | 2,2,2 | 0–5 |
| 01 | 074d9 low, 2824/885 | 3,3,13 | 9 |
| 02 | 103b2 high, 3726/881 | 2,2,3 | 0 |
| 03 | 103b2 low, -1049/2296 | 2,2,3 | 0 |
| 04 | 11952 high, -2448/11 | 2,5,2 | 0 |
| 05 | 11952 low, -1171/1683 | 2,13,2 | 10 |
| 06 | 11952 high, 110314/102227 | 2,2,2 | 0 |
| 07 | 11952 low, 130349/28916 | 3,2,2 | 5 |
| 08 | 11952 high, 2012/211 | 2,2,2 | 0–5 |
| 09 | 11952 high, 2828/2015 | 2,2,2 | 0 |
| 10 | 11952 high, 4286/1881 | 2,2,2 | 0–5 |
| 11 | MW16 high, -1867/270 | 2,41,2 | 0–6 |
| 12 | MW16 low, -3187/3697 | 2,2,2 | 0–7 |
| 13 | ICARM356 | 2,2,2 | 1 |
| 14 | ICARM385 | 5,2,2 | 0 |
| 15 | ICARM398 | 2,2,2 | 0 |

The strict dimension in **each selected pair field** is exactly zero
in every row. The last column concerns the much larger all-generic
compositum and equals the already computed inherited strict dimension
or its retained interval. Partial factorization is sufficient for the
pair exclusions: one verified nonsquare place excludes a strict class.
An interval in the last column is not replaced by zero.

Three same-family contrasts now have a precise interpretation. The
103b2 high and low both have zero strict classes in their entire generic
governing composita, while the high requires at least nine additional
strict directions. The compact11952 low has ten inherited strict classes
against zero on its matched high. The larger11952 low has five against
zero on its matched high. Those inherited CT structures measure neither
the missing class-group excess nor its rational solubility.

## Three arithmetic separation witnesses

The [small separation protocol](GOVERNING_THIRD_CLASS_SEPARATION_PROTOCOL.json)
uses only generic section index2 on tokens02,03,13. If that section is
(x,y), its quartic is

    Z^4 - 3*x/2*Z^2 - y*Z - 3*x^2/16 - A/4.

Its cubic field is K and its class is x-theta, by the
[fixed-resolvent calculation](FIXING_THE_RESOLVENT_POLYNOMIAL_TESTS_SOLUBILITY.md).
At the first witnessing prime at most10009, the retained governing octic
splits completely and this quartic has factor degrees(2,2). The cubic
has three roots and the third Kummer class has exactly two nonsquare
components. These patterns prove noncontainment directly.

The independent replay checks eight distinct octic roots. It proves
the quartic has no linear factor and divides Z^(p^2)-Z, hence is a product
of two distinct irreducible quadratics, without calling a factorization
routine. All three cases pass at the primes in the opening table.

The third section remains generic. This is a falsifiable calibration of
the field-capacity statement, not an exceptional-class discovery.

## Consequence for the next incidence experiment

The useful new restriction is structural: new strict directions must
enter through a **noncentral standard-module extension** outside the
whole generic governing compositum. A central cochain or a new value of
an existing CT bit cannot supply that extension.

Governing bits can supply second-descent obstruction information under
the appropriate local twist hypotheses. Neither their values nor the
governing degrees supply newly constructed incidence directions. The missing computation remains an independent
unramified cubic class or corresponding S4 field outside G, with the
local and independence certificates. This turn supplies no such class,
no additional-quotient CT matrix and no parameter-selection rule.

## Artifacts and replay

The [incidence protocol](GOVERNING_FIELD_INCIDENCE_CLOSURE_PROTOCOL.json),
[main certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_governing_field_incidence_closure_v1.json)
and [independent verification](../../artifacts/generated-results/elliptic-curves/rank_jump_governing_field_incidence_closure_verification_v1.json)
cover the finite group and all48 local exclusions. Their inputs project
only generic arithmetic and bind the completed local verification; no
new number-field setup or local-power computation is performed.

The [separation certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_governing_third_class_separation_v1.json)
and [independent polynomial replay](../../artifacts/generated-results/elliptic-curves/rank_jump_governing_third_class_separation_verification_v1.json)
record the three finite-prime witnesses.

```sh
timeout 30 python3 elliptic-curves/rank-jump/governing_field_incidence_closure.py check
timeout 30 python3 elliptic-curves/rank-jump/verify_governing_field_incidence_closure.py check
timeout 30 sage -python elliptic-curves/rank-jump/governing_third_class_separation.py check
timeout 30 python3 elliptic-curves/rank-jump/verify_governing_third_class_separation.py check
```

Only new rank-jump-specific files are added. Active search protocols,
scoring, candidates, worker limits and mathematical-status entries are
untouched.

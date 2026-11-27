# Literature screen: equal-mass all-k Gaussian first moments

Screen date: 8 October 2026. This records source scope only; it is not a proof
review, exhaustive citation search, or priority assessment.

## Exact target and direct prior conjecture

For k >= 2, d >= k-1, and any measurable (or fractional) partition with
gamma_d(C_i)=1/k, the candidate bound is

    sum_i || integral_{C_i} x d gamma_d ||^2
       <= (E max_{1<=i<=k} Z_i)^2/(k-1),

where the Z_i are independent standard normals, as in the paper.
Equivalently one can subtract their common average: the centered score
vector has covariance P = I - 11^T/k (diagonal (k-1)/k), and its
expected maximum is unchanged. A unit-variance simplex score vector
requires an additional normalization and is not used in this formula. The proposed equality case is the
regular central simplex fan,
cylindrically extended, up to null sets and symmetries.

This is directly covered as a conjecture by Steven Heilman, “Stable Gaussian
Minimal Bubbles,” arXiv:1901.03934v1 (13 January 2019),
https://arxiv.org/abs/1901.03934. Problem 1.15 and Conjecture 1.16, printed
p. 8, formulate the prescribed-positive-mass first-moment optimization for
m>3 and m-1 <= ambient dimension, and conjecture regular-simplex simplicial
cones are maximizers. At equal masses, the prescribed centering shift is
zero. Thus the present all-k statement is an equal-mass specialization of
that conjecture; the closed-form right side is the candidate fan's value.
The introduction, printed p. 4, calls the m>3 first-moment problem open at
the time. This does not establish its status as of this screen date.

An explicit low-dimensional predecessor is Steven Heilman, “Euclidean
Partitions Optimizing Noise Stability,” *Electronic Journal of Probability*
19 (2014), no. 71, 1–37, arXiv:1211.7138v2,
https://arxiv.org/abs/1211.7138. Conjecture 3, printed pp. 36–37, treats
k=4,d=3 and says a maximizer of the first-moment functional psi_0 among
equal-mass fractional partitions is simplicial conical. It distinguishes
the unrestricted fractional class (cited as known) from the equal-measure
constraint. This directly overlaps the k=4 case, not the all-k theorem.

These Heilman sources were read in full from `/tmp/1901.03934v1.pdf/.txt`
and `/tmp/1211.7138v2.pdf/.txt`. The candidate all-k statement settles the
equal-mass subcase of Conjecture 1.16 if its proof is validated; it is not a
previously unstated conjecture.

## Nearby results checked and distinguished

- Emanuel Milman and Joe Neeman, “The Gaussian Double-Bubble and Multi-Bubble
  Conjectures,” *Annals of Mathematics* 195 (2022), no. 1, 89–206,
  DOI 10.4007/annals.2022.195.1.2,
  https://doi.org/10.4007/annals.2022.195.1.2, arXiv:1805.10961v3,
  https://arxiv.org/abs/1805.10961. Full text checked at
  `/tmp/1805.10961.pdf/.txt`. Theorems 1.1–1.2 prove unique minimization of
  Gaussian perimeter by the simplicial cluster for 2 <= q <= d+1 and every
  positive prescribed mass vector. This includes the regular-simplex fan,
  but the optimized functional is perimeter, not squared first moments.
- Emanuel Milman, “Multi-Bubble Isoperimetric Problems,” arXiv:2510.07078v2
  (31 July 2026), https://arxiv.org/abs/2510.07078. Full text checked at
  `/tmp/2510.07078v2.pdf/.txt`. Section 8.5, Proposition 8.4, printed
  pp. 17–18, gives the Gaussian model-profile PDE
  `tr((-nabla^2 I_model)^(-1)) = 2 I_model`. Sections 8.5–8.6 derive profile
  inequalities by variations of perimeter-minimizing clusters. This is
  profile/perimeter theory, not an all-partition centroid bound.
- Abhijeet Mulgund, “Stochastic Domination of Gaussian Maxima by the Regular
  Simplex,” arXiv:2609.28452v2 (27 September 2026),
  https://arxiv.org/abs/2609.28452. Full text checked at
  `/tmp/2609.28452v2.pdf/.txt`. Theorem 1.1 compares CDFs of maxima of
  centered unit-variance Gaussian vectors for all thresholds; Corollary 1.2
  is a circumscribed-simplex Gaussian-volume bound; Theorem 6.2 treats
  signal identification with an inactive state. It does not state a theorem
  for arbitrary equal-mass measurable partitions. In particular, it does
  not directly cover a general centroid score matrix whose row variances
  need not be equal.
- OpenAI Mathematics, “The Gaussian Propeller Bound in Every Dimension,”
  24 September 2026, local source at
  `../openai-math/preprints/The-Gaussian-Propeller-Bound-in-Every-Dimension-September-24-2026/main.pdf`;
  checked repository snapshot `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`.
  Theorem 1.1 permits arbitrary cell probabilities and empty cells, and
  gives the sharp unrestricted value 9/(8 pi), attained by three planar
  sectors. It does not fix k positive equal masses. For k=3 the regular
  simplex fan is the propeller; for k>=4 this unrestricted result is not the
  claimed equal-mass sharp bound.
- Steven Heilman, Elchanan Mossel, and Joe Neeman, “Standard Simplices and
  Pluralities are Not the Most Noise Stable,” *Israel Journal of Mathematics*
  213 (2016), no. 1, 33–53, arXiv:1403.0885v3,
  https://arxiv.org/abs/1403.0885. Theorem 2.6 gives positive-noise,
  unequal-mass counterexamples. This refutes a neighboring shifted-simplex
  noise-stability claim for unequal masses; it is not a counterexample to
  the equal-mass first-Hermite endpoint.

## Multi-bubble/profile implications checked

In the 2022 Milman–Neeman paper, equations (1.2)–(1.4) define the
interface-area Laplacian `L_A` and identify the model Gaussian multi-bubble
profile Hessian as `nabla^2 I_model(v) = -L_A_model(v)^(-1)`. Lemmas 11.1
and 11.2 (printed pp. 56–57) apply to an **isoperimetric minimizing Gaussian
cluster** at a prescribed interior mass vector. They identify its profile
Hessian and interface areas with those of the model simplex cluster and
derive the second-variation inequality

    -(delta_X V)^T L_A^(-1) (delta_X V) <= Q(X).

For translations, the paper writes `delta_w V = M w`, where
`M = sum_{i<j} A_ij (e_i-e_j) n_ij^T`, and computes the Gaussian translation
second variation. Its Cauchy–Schwarz argument classifies perimeter minimizers
as model simplicial clusters. These statements assume perimeter minimality;
they are not bounds for the centroid matrix of an arbitrary partition.

Milman–Neeman §12.4 mentions a possible functional version of the Gaussian
multi-bubble inequality for simplex-valued functions as “in preparation”; no
paper or explicit first-moment theorem is supplied there. The 2026 survey
does not identify such a result. The searched theorem statements do not give
the proposed all-partition anisotropic bound
`tr(B^T L_p^+ B) <= c_p` for arbitrary cell-moment matrix B at quota p.

## Limited status and adjacent resolved/refuted claims

The exact all-k first-moment assertion is an established prior conjecture.
The sources above do not state a proof or refutation of its equal-mass
specialization. The directly solved low-cardinality cases are k=2, by the
Gaussian isoperimetric/halfspace first-moment bound, and the all-partition
upper bound for k=3, by OpenAI result 096 (the equal-mass propeller fan
attains it). Other nearby resolved claims are: (i) Gaussian multi-bubble
perimeter minimization in the Milman–Neeman range; (ii) the unrestricted
propeller first-moment maximum in OpenAI result 096; and (iii) unit-variance
Gaussian maximum domination in Mulgund’s Theorem 1.1. The nearby refuted
claim is the unequal-mass positive-noise shifted-simplex claim in
Heilman–Mossel–Neeman Theorem 2.6. None is the target theorem.

No equivalent all-partition theorem or counterexample was located in this
focused check of the cited primary sources, local Gaussian partition records,
quantization/principal-points search results, or the named OpenAI snapshot.
This is a limited non-hit, not proof that the conjecture remains open and not
a priority claim. If validated, the candidate result would prove the
equal-mass all-k subcase of Heilman Conjecture 1.16.

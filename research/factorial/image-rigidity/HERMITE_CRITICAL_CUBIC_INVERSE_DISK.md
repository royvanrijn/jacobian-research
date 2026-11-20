# An explicit full cubic inverse disk at the Hermite critical point

This lemma concerns the exact limiting phase at delta0,z=2i. It does not certify a full inverse domain on a nonzero parameter window and does not give a numerical Strong Factorial cutoff.

**Lemma CriticalHermiteFullInverseDisk.**

Put v=-i+u and

S(u)=F_0(-i+u,2i)-F_0(-i,2i)
    =log(1-i*u)+i*u-u^2/2
    =i*u^3*I(u),
I(u)=integral_0^1 s^2/(1-i*u*s) ds.

For |u|<=r<1, the reciprocal-disk inequality gives Re I(u)>=1/[3*(1+r)]>0. Hence the principal cube root in

t(u)=-i*u*[3*I(u)]^(1/3)

is analytic on |u|<1, has a single simple zero u0, and satisfies t^3=3*S and t'(0)=-i. Moreover

|t(u)|^3>=|u|^3/(1+r) on |u|=r.

At r=3/5 this gives |t|>51/100, since

(51/100)^3*(8/5)<(3/5)^3.

For every |tau|<51/100, Rouche applied to t(u)-tau on |u|=3/5 therefore gives exactly one zero, counting multiplicity. It is simple, so the inverse u(tau) is holomorphic on the entire disk |tau|<51/100 and |u|<3/5 there. In particular the closed radius1/2 disk lies strictly inside its analytic domain. This is a proved fixed analytic inverse disk, not an inference from the stationary interpolation function H.

Let h(t)=dv/dt=du/dt. The Cauchy derivative bound, using the inverse disk radius1/2, gives

|h(t)|<=6 for |t|<=2/5.

The critical interpolation values are h(0)=i and h'(0)=i/2. Write

h(t)=i+(i/2)*t+t^2*G(t).

The exact Cauchy divided-amplitude derivative formula on the circle |tau|=2/5 gives, for |t|<=3/10,

|G'(t)|<= (2/5)*6 / [((2/5)-(3/10))^2*(2/5)^2]=1500.

This controls the TRUE inverse Jacobian and its first divided amplitude at the exact critical base point. Uniform persistence is not established by this basepoint argument; the full thin window is treated in [the uniform cubic inverse theorem](GAP_TWO_FULL_ENDPOINT_UNIFORM_CUBIC_INVERSE.md), without asserting the whole broader complex disk. The corresponding complete-period vector and nonlocal constants are separate quantitative obligations. The auxiliary stationary extension in [the matched-data note](QUANTITATIVE_GAP_TWO_ENDPOINT_DATA.md) cannot supply that missing persistence, as its exact quadratic mismatch shows.


The three rational margins are replayed by the shared verifier:

    python3 research/factorial/image-rigidity/verify_quantitative_gap_two_endpoint_data.py \
      --repo-root . --output /tmp/hermite_critical_inverse_review.json

Updated 20 November 2026.

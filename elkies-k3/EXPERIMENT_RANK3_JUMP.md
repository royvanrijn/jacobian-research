# Explicit E8+A2^3 rank-jump system

From the exact NS computation the inherited generic fibration has:

    fibers: II* + IV + IV + IV + II
    ADE: E8 + A2^3
    MW rank on the target X(6,79) locus: 3
    MW height Gram, after reduction:
        (1/3) * [[8,-1,0],[-1,10,0],[0,0,12]]

Normalize the singular fibers to

    infinity : II*
    0,1,lambda : IV
    mu : II

giving

    y^2 = x^3 + [t(t-1)(t-lambda)]^2 (t-mu).

At the discriminant-3 CM point lambda=mu=0 this becomes

    y^2 = x^3 + t^5(t-1)^2,

Utsumi No.1.

A polynomial section with deg x<=4, deg y<=6 and deg(y^2-x^3)<=7 can be
parameterized generically by

    x=q^2+r
    y=q^3+s

with deg q=2 and deg r,s<=1.

Therefore the existence of one extra section is 8 coefficient equations in only
9 variables (lambda,mu plus 7 section parameters). Its solution should contain
the desired one-dimensional rank-jump curve.

First inspect:

    sage elkies-k3/scripts/show_rank3_jump_system.sage

Then finite-field sliced samples:

    python3 elkies-k3/scripts/run_rank3_jump_msolve.py \
      --primes 101,103,107 --seeds 1,2,3 --threads 8 --timeout 600

Unlike the previous 48/167-variable systems, this is only 9 variables / 9 equations
after one generic slice. If msolve still struggles, the next step is direct
elimination using the explicit coefficient triangularity, not more compute.

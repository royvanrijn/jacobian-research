# MW3 local-Tate construction

This refines MW3_CONSTRUCTION.md by replacing the huge Delta(lambda) equations
with local multiplicative-fiber equations.

At a multiplicative fiber t=a with singular cubic double root s:

    A(a) = -3 s^2
    B(a) =  2 s^3

and ord Delta >=2 is equivalent (on the multiplicative chart) to

    B'(a) + s A'(a) = 0.

For I3, the next condition is

    A'(a)^2 = 6 s (s A''(a) + B''(a)).

This dramatically lowers the degree/term count of the lambda-fiber equations.

The known component profiles give first-blowup section constraints. A section on
a nonidentity component specializes to the singular point:

    x(a)=s, y(a)=0.

For sections with denominator z=t-r:

    X(a)=s z(a)^2, Y(a)=0.

At infinity our normalization A8=-3, B12=2 gives singular point (X,Y)=(1,0),
so every nonidentity section has the simple linear conditions on leading terms:

    x_lead=1, y_lead=0.

Run:

    sage elkies-k3/scripts/build_mw3_local_tate_system.sage --stage p1
    sage elkies-k3/scripts/build_mw3_local_tate_system.sage --stage p12
    sage elkies-k3/scripts/build_mw3_local_tate_system.sage --stage all

Only after inspecting counts/complexity should you try:

    python3 elkies-k3/scripts/run_mw3_local_probe.py --stage p1 --p 101

The deeper I11 component indices 2,6,10 are not yet encoded; this first layer
only distinguishes identity vs nonidentity. Those deeper blow-up equations are
the next refinement if needed.

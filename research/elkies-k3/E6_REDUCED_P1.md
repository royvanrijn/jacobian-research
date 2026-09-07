# Reduced E6 P1 construction

The first E6 section scaffold had 32 variables. Before solving, eliminate
variables forced linearly by local fiber and section conditions:

At I4@0:

    a0 = -3 s0^2
    b0 =  2 s0^3
    b1 = -s0 a1

At I4@1:

    A(1)=-3s1^2  -> solve a5
    B(1)= 2s1^3  -> solve b7 in terms of b6
    B'(1)+s1 A'(1)=0 -> solve b6

For P1=(1,0,1,0,0):

    x4=0, y6=0             (nonidentity at IV*)
    X(1)=s1, Y(1)=0        (nonidentity at I4@1)

so solve x3 and y5.

This removes ten variables before any Gröbner work.

Run:

    sage elkies-k3/scripts/build_e6_reduced_p1_system.sage

Then modular probe:

    python3 elkies-k3/scripts/run_e6_reduced_p1_probe.py \
      --p 101 --threads 8 --timeout 300

If still hard, the next reduction is to derive the deeper I4 component-1
blow-up condition for P1 and use the low-order I4 equations to eliminate b2,b3.

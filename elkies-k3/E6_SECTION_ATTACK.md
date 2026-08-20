# E6 MW3 section construction

Canonical component profile chosen from the 32 Shioda-compatible triples:

    P1=(1,0,1,0,0), P1.O=0
    P2=(1,1,2,1,1), P2.O=1
    P3=(2,3,0,0,0), P3.O=0

Tuple order:

    (IV*, I4@0, I4@1, I2@lambda, I2@mu)

All mutual section intersections are 2.

The preferred fiber model is:

    IV* + I4 + I4 + I2 + I2 + 4 I1.

Two generators, P1 and P3, are polynomial sections. This motivates solving:

    fiber family -> P1 -> P1+P3 -> all three

rather than adding the rational P2 second.

Run metadata + one P1 probe:

    bash elkies-k3/scripts/start_e6_section_attack.sh

Or manually:

    sage elkies-k3/scripts/build_e6_mw3_section_system.sage --stage p1
    sage elkies-k3/scripts/build_e6_mw3_section_system.sage --stage p13
    sage elkies-k3/scripts/build_e6_mw3_section_system.sage --stage all

Then:

    python3 elkies-k3/scripts/run_e6_mw3_probe.py \
      --stage p1 --p 101 --threads 8 --timeout 300

Only first singular-hit component conditions are encoded initially. Exact IV*
classes 1/2 and I4 components 1/2/3 require deeper blow-up equations if the
first system remains too large.

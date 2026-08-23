# H92 q24 Hensel/SRR dead end

Status: **HISTORICAL_DIAGNOSTIC**.

The q24 direct Hensel route reached p-adic precision 512 at prime `100003`, but
the checkpoint remained `complete=false` and `exact_verified=false`.
`exactify_h92_q24_from_padic_srr.sage` then tried simultaneous rational
reconstruction on the 156 residues. Every block/chunk strategy reconstructed a
candidate vector, but none passed the exact Weierstrass identity and modular
regression gate.

Current instruction for the fixed corridor:

- do not increase q24 Hensel precision;
- do not use SRR as the next exactification step;
- use `replay_h92_q24_d12_discovery_and_theta.sage` and the explicit
  `D24 = Theta + sum m_i C_i` construction to build the exact resolved-RR
  q24/orbit85 pencil;
- keep orbit42 on the corrected profile `height=7`, `correction=3`, `P.O=3`.

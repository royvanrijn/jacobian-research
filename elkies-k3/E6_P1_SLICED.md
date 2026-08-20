# Frozen four-elimination sliced E6 P1 sampler

This freezes the four safe triangular eliminations already observed:

    P1_7 -> b4
    P1_6 -> b5
    P1_5 -> y2
    P1_4 -> a3

This leaves 13 variables on the expected 3-dimensional P1 locus.

Instead of further symbolic elimination, add 3 generic affine linear slices over
GF(p), producing a zero-dimensional sample suitable for msolve.

Run several seeds:

    python3 elkies-k3/scripts/run_e6_p1_sliced_batch.py \
      --p 101 \
      --seeds 1,2,3,4 \
      --slices 3 \
      --workers 2 \
      --threads 4 \
      --timeout 300

Each seed generates a different generic 3-plane slice through the P1 locus.

If solutions are returned, use the accompanying .meta.txt file to reconstruct
b4,b5,y2,a3 and then the complete A(t),B(t),P1. Those sampled surfaces become
starting points for adding P3.

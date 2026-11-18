# Hyperbolicity certificates for Prime Knots with Augmentation-Rank Gaps at Unbounded Bridge Number

Extract the archive into an ordinary working directory and run the commands below
from the extracted `prime-knot-certificates` directory. No research checkout
is needed. Each hyperbolicity input contains the literal braid, planar diagram,
component marking, decorated ideal triangulation, face pairings, logarithmic
gluing equations and rational interval bounds. The diagram-move input contains
the two finite disk moves and the marked band of the fixed six-strand braid.

## Exact diagram checks

Run with Python 3.12 or later. The diagram checker verifies the two oriented
disk moves with their boundary endpoints fixed and compares the literal source
and target diagrams at bridge numbers 5 through 64. The uniform conclusion
follows from the induction in Appendix A of the article.

```sh
python3 checks/verify_threaded_curl_uniform_transfer.py
```

The three scripts below also run without `--sage` for their exact braid,
diagram and marking checks. Those checks use permutations, free-group reduction
or oriented port pairings. They reject modified detours and incorrect component
markings. The general coloring and rank-one arguments are written in the article.

## Interval verification

Use the Linux x86-64 SageMath 10.9 environment with Python 3.14.7, SnapPy
3.3.2 and Spherogram 2.4.1. Use the Python executable in that environment in
the commands below. If SnapPy and Spherogram are installed in a
separate package directory, put that directory on `PYTHONPATH`.
The input JSONs pin the analytic and topology modules and
the SnapPy kernel binary by SHA-256. The replay checks those hashes in addition
to the package versions.

```sh
python checks/verify_protected_five_bridge_uniform_family.py --sage
python checks/verify_protected_six_bridge_uniform_family.py --sage
python checks/verify_protected_universal_parent_candidate.py --sage
```

The scripts certify `M_5`, `M_6` and `N`, in that order. Each run should
end in `PASS` and a successful rigorous geometry check. All three use 100-bit
interval arithmetic and zero triangulation randomizations. They bind the
triangulation to the exact diagram, verify gluing and orientation conditions
on the retained boxes, and rerun the interval existence test. The article
deduces the prime-family conclusion by two successive Dehn fillings.

# Bridge index of thirteen-crossing knots

Both results of this repository are written up in one paper, *The bridge index and the nonorientable four-genus of knots with at most thirteen crossings*: [`paper/paper1/knot-tables.pdf`](paper/paper1/knot-tables.pdf), with the TeX source in [`paper/paper2/knot-tables-tex.zip`](paper/paper2/knot-tables-tex.zip). Part I is the bridge index, Part II the nonorientable four-genus. Release 1.3.0 is archived on Zenodo, doi:[10.5281/zenodo.23561157](https://doi.org/10.5281/zenodo.23561157); all versions doi:10.5281/zenodo.23521966.

KnotInfo (version 2026.10.5) lists the bridge index of 536 prime knots with 13 crossings as
`[3,4]`: it is known to be 3 or 4, and the exact value is not known.

**Result.** 243 of these knots have bridge index 4 and meridional rank 4, so the meridional
rank conjecture of Cappell and Shaneson holds for them (Kirby's problem list, Problem 1.11).
The list is in [`paper/bridge13.pdf`](paper/bridge13.pdf), Appendix A, and in
[`data/certificates.json`](data/certificates.json). The other 293 knots, listed in
[`data/still_open.json`](data/still_open.json), stay open.

| group | image of a meridian | knots |
|---|---|---|
| Sp(4, F₃) | transvection | 67 |
| Sp(4, F₅) | transvection | 161 |
| O(4, F₅) | reflection | 15 |

## Second result: the nonorientable four-genus

The folder [`crosscap4d/`](crosscap4d/) settles the smooth nonorientable four-genus γ₄ of
2,537 of the 3,525 knots with 11 to 13 crossings that KnotInfo leaves open (520 at γ₄ = 1,
1,995 at γ₄ = 2, 22 at γ₄ = 3) and narrows the interval of 376 more. The upper bounds are
pinch moves checked in Regina; the new lower bounds come from Donaldson's theorem through
the Goeritz-lattice obstructions of Jabuka and Kelly. The paper is
[`crosscap4d/paper/crosscap4d.pdf`](crosscap4d/paper/crosscap4d.pdf).

## The proof

Write μ(K) for the meridional rank and b(K) for the bridge index; always μ(K) ≤ b(K).

- **Upper bound.** The KnotInfo diagram of each knot has Wirtinger number 4. By Blair,
  Kjuchukova, Velazquez and Villanueva, the bridge index is at most the Wirtinger number of
  any diagram, so b(K) ≤ 4.
- **Lower bound.** Each certificate is a homomorphism ρ from the knot group to GL(4, F_q)
  that sends the meridian to a pseudo-reflection: a matrix fixing a hyperplane.
  - Each meridian is sent to a pseudo-reflection, and the hyperplanes they fix meet only
    in 0, so no nonzero vector is fixed by the whole image.
  - If k meridians generated the group, their images would fix the intersection of k
    hyperplanes. That intersection is nonzero when k ≤ 3, so μ(K) ≥ 4.

The pseudo-reflections used are transvections x ↦ x + ω(x,v)v of a symplectic form, or
reflections of a quadratic form. Earlier lower bounds came from Coxeter quotients, where
meridians go to order-2 reflections. The order-3 and order-5 transvections settle most of
these knots.

## Checking it

```
python3 verify/verify.py              # standard library only
pip install database_knotinfo
python3 verify/verify.py --knotinfo   # also compares PD codes and [3,4] entries with KnotInfo
```

For every knot the verifier checks each of the following:

- the matrices are constant along arcs;
- the Wirtinger relation holds at every crossing;
- each arc matrix minus the identity has rank 1;
- the stacked matrices minus the identity have rank 4;
- the four seed arcs color the whole diagram by the Wirtinger coloring moves.

All 243 pass in a few seconds. The verifier rejects four kinds of corrupted certificate:

- a changed matrix entry;
- seeds on three arcs;
- one crossing with its overstrand reversed;
- identity matrices.

## Layout

- `paper/`: the joint paper in `paper1/` (PDF) and `paper2/` (TeX zip), its source in
  `knot-tables/` (`make_generated.py` copies the two `generated.tex` files after
  `paper/make_tables.py` and `crosscap4d/paper/make_tables.py`), and the earlier paper
  `bridge13.tex` with its PDF. `make_tables.py` writes `generated.tex` from the data.
- `data/certificates.json`: for each knot, the KnotInfo PD code, the prime q, one 4×4 matrix
  per edge, and four seed edges.
- `data/still_open.json`: the 293 knots with entry `[3,4]` that remain open.
- `verify/verify.py`: the verifier.
- `search/`: the programs that found the certificates (`symp.py`, `run_sp2.py`, `generic.py`
  for searches over colorings; `gram.py`, `ogram.py`, `glsearch.py`, `gram_gf.py` for searches
  over Gram matrices of the seeds; `w3.py` for 3-bridge diagrams; `alexq.py`, `alexq_exact.py` for Alexander-quandle colourings and
  `braid3.py` for closed 3-braids, both of which settle none of the open knots). `make_certificates.py`
  writes `data/certificates.json` from their output. These programs need numpy, SnapPy and
  database_knotinfo, and are run from `search/`; the proof does not depend on them.
- `CITATION.cff` and `.zenodo.json`: how to cite this repository, and the metadata of its Zenodo
  archive. `.github/workflows/release.yml`, run by hand, checks every certificate against KnotInfo
  and publishes a release with the paper and a zip of the repository; Zenodo archives it.

## Citing

Release 1.1.0, with both results:
[github.com/creelie/conjecture-to-choose/releases/tag/v1.1.0](https://github.com/creelie/conjecture-to-choose/releases/tag/v1.1.0).
Release 1.0.0 (bridge index only) is archived on Zenodo:
[10.5281/zenodo.23521967](https://doi.org/10.5281/zenodo.23521967). The DOI
[10.5281/zenodo.23521966](https://doi.org/10.5281/zenodo.23521966) always resolves to the latest
version. [`CITATION.cff`](CITATION.cff) has the full citation.

## References

- R. Blair, A. Kjuchukova, R. Velazquez, P. Villanueva, Wirtinger systems of generators of
  knot groups, Comm. Anal. Geom. 28 (2020), 243–262, arXiv:1705.03108.
- R. Blair, A. Kjuchukova, N. Morrison, Coxeter quotients of knot groups through 16
  crossings, Experimental Mathematics 34 (2025), 27–37, arXiv:2208.09032.
- R. Blair, A. Kjuchukova, E. Pfaff, The plain sphere number of a link, arXiv:2504.10517.
- C. Livingston, A. H. Moore, KnotInfo, https://knotinfo.org.

Author: Deep Bhattacharjee.

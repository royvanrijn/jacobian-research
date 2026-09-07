# Current ICARM data

[`icarm_current.json`](icarm_current.json) selects the complete compressed
public database and its current publication overlay on the frozen V22 inventory.
The full raw snapshot preserves ICARM attribution, metadata and points.

From Python with `elliptic-curves/cas` on the import path:

```python
from refresh_icarm_local_database import load_catalogue, load_inventory
public_database = load_catalogue()
research_curves = load_inventory()
```

For the same research rows with certified conductors and unresolved bounds:

```python
from local_conductor_database import load_conductor_inventory
curves_with_conductors = load_conductor_inventory()
```

The added `conductor_information` field distinguishes `EXACT`, `UNKNOWN`, and
publicly `REPORTED` values. Exact rows carry complete bad-prime lists and proof
paths; unresolved rows carry a proved divisor and an upper bound. The
[`conductor_screen_current.json`](conductor_screen_current.json) manifest pins
the arithmetic certificates. Record comparisons retain their dated snapshot.

`load_inventory()` updates publication status only. Discovery point sets, rank
certificates and the distinction between local search and public reproduction
remain unchanged. Public rank and conductor fields are reported metadata unless
separately certified. See the
[audit and six recovery controls](../notes/CURRENT_ICARM_DATABASE_AND_CONTROLS_2026-09-07.md).

Use this view for post-discovery deduplication and declared retrospective work.
Do not substitute it into frozen candidate selection or treat its public points
as prospective discoveries.

The current [research database export](research_curves/database.json) and
[expanded inventory](../INVENTORY.md) contain all201 minimal equations, metrics
and point sets, with129 exact conductors and72 unresolved values. Rebuild these
views with `render_main_readme_curves.py` after selecting a new certified summary.

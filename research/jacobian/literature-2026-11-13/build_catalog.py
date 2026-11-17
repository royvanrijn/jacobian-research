"""Build the dated literature inventory from retained metadata, without network or CAS."""
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
rows = {}
for filename in ('all-years-records.json', 'api-additional-records.json', 'related-records.json'):
    for row in json.loads((HERE / filename).read_text()):
        key = row['id'].rsplit('/', 1)[-1]
        if row['updated'][:10] < '2026-01-01':
            continue
        if filename == 'related-records.json' and not key.startswith('2608.07338'):
            continue
        if key.startswith('2610.25831'):
            continue  # Mean-curvature self-shrinkers: an irrelevant nonproperness hit.
        rows[key] = {k: v for k, v in row.items() if k != 'query'}
        rows[key]['source'] = row['id'].replace('http:', 'https:')
        rows[key]['review_depth'] = 'metadata and abstract screening'
        rows[key]['independent_proof_audit_in_this_pass'] = False
        rows[key]['priority'] = 'adjacent'

def annotate(key, priority, depth, relevance):
    rows[key].update(priority=priority, review_depth=depth, repository_relevance=relevance)

annotate('2609.33491v1', 'high', 'main results and sections 5.1 through 5.4', 'PF2R7CUSP1 / OP-PLANE: complete missing-divisor ledger includes unramified components in dimension at least three; no plane descent follows.')
annotate('2608.27341v2', 'high', 'main theorem and introductory proof strategy', 'Global affine-space recognition requires a projective bundle over P1 and hyperplane boundary on every fibre. Version 2 corrects the main theorem. Resolved plane boundaries need not have this form.')
annotate('2608.12543v3', 'high', 'main theorem, sections 5 and 6, source attribution', 'Directly credits Roy van Rijn. The displayed collision fills the 19-dimensional polarization algebra, excluding only invariant linear compression preserving that collision, before the unchanged symmetric lift.')
annotate('2610.24889v1', 'high', 'full concrete-map audit in the sibling de-bondt directory; broader operator statements only partly inspected', 'BCRDB1 / BCRDBGRAPH1: 12/18/90 construction, economic stabilization, and VC witnesses; not JC2 descent.')
annotate('2610.36495v1', 'high', 'main theorem and proof outline', 'Claims ordinary total degree at most four for three-variable Keller maps. Not a bound on geometric degree; degree-five and degree-six cases remain proposed questions.')
annotate('2610.36736v1', 'high', 'abstract and stated criteria', 'Single Legendre-potential membership / finite derivative criterion. Doubling supplies no independent Nambu conditions; a criterion is not a termination proof.')
annotate('2610.33919v1', 'high', 'abstract and introductory main implications and hypotheses', 'JC2 implies generalized DC1. Converse needs characteristic-zero graph-module existence, equivalently the specified uniform BW boundedness condition. Ordinary DC1 alone does not close this route.')
annotate('2607.20597v2', 'medium', 'main closedness argument and component conclusion', 'Parameter-family arguments require irreducible-component control. Does not transport the single r7 braid witness to the entire polynomial-gap family.')
annotate('2610.02961v2', 'medium', 'abstract and arithmetic setup, sections 2 and 3 opening', 'Integral coefficients and Jacobian one are material assumptions. Integral-point conditions do not exclude complex r7 curves without an arithmetic bridge.')
annotate('2609.10180v1', 'medium', 'abstract and characteristic-p setup', 'Differential-operator extension preserving holonomicity in positive characteristic is distinct from the missing characteristic-zero lifting theorem.')
annotate('2609.05746v1', 'medium', 'abstract and introduction; not the full 175-page proof', 'Positive-characteristic counterexamples, moduli, and special projection-factorization results; do not assert unrestricted JC2 in characteristic zero.')
annotate('2608.14217v1', 'high', 'Theorem 1.2, structural proposition 6.1, scope discussion', 'Independent quartic HC4 claim overlaps HC4CQ1. Arbitrary-degree P is allowed in a special constant-direction form with degree Q at most two.')
annotate('2608.19112v2', 'medium', 'abstract and signature hypotheses', 'Degree-five HC4 and pivot criteria under real Lorentzian index-one hypotheses; not all complex potentials.')
annotate('2608.00222v2', 'medium', 'abstract and organization; follow-up inspection of sections5.5,5.6 and11.1', 'Revised4 November, still v2. Finite-cover one-point support and local fibre length differ from fibres of the original etale source. Complete F4-F7 fibre stratifications remain open. Relevant to the singleton-zero/properness distinction in the local HC4MYOR2 application; no full external proof replay is claimed.')
annotate('2611.03572v1', 'adjacent', 'metadata and abstract', 'Families of higher-dimensional counterexamples; no plane invertibility theorem.')
annotate('2607.22198v2', 'medium', 'metadata and abstract', 'Five-variable HC counterexample; distinguish constant-Hessian potentials from homogeneous Hessian-nilpotent VC quartics.')
annotate('2608.07338v1', 'medium', 'metadata and explicit identity in abstract', 'Generalized VC counterexample uses a third-order operator and an auxiliary multiplier. It is not a five-variable quartic Laplacian VC witness.')
annotate('2607.18186v1', 'medium', 'metadata and abstract', 'Explicit Gaussian moment witnesses in dimension three are not derived through the much larger standard Keller-reduction chain.')
annotate('2609.16178v1', 'medium', 'metadata and abstract', 'Mathieu property classified by torus condition; published 31 July despite September identifier. Moment constructions may be more economical than Keller reductions.')
annotate('1203.1691v40', 'claim requiring audit', 'updated abstract and version history', 'Claims JC2 through an affine open-embedding theorem. Not accepted as input: the generalized simply-connected assertion requires particular care in light of Kulikov-type surfaces.')
annotate('2304.14675v5', 'incompatible with verified higher-dimensional examples', 'updated abstract', 'Claims all-dimensional JC; cannot be imported alongside the explicit noninjective Keller maps verified in this repository.')
annotate('2410.06959v5', 'claim requiring audit', 'abstract and comparison with Kanel-Belov--Fedorov introduction', 'Ordinary first-Weyl Dixmier proof claim. The October paper says its authors could not fully understand the proof; no generalized DC1 or JC2 consequence is imported.')

manual = [
 ('hexagon:2611.00156v1', 'Equivalence of the Four-Variable Hessian and Planar Jacobian Conjectures', ['Qiyue Tang'], '2026-11-06', 'https://hexagonmath.org/pdf/2611.00156v1', 'reproduced structural theorem', 'Full written proof review of Theorems2.1/2.3, De Poi import, Lemmas3.2-3.8 and Proposition4.1; exact supporting algebra replayed', 'HCTANG2 reproduces the constant-linear cotangent alternative and HC4-JC2 equivalence. HCTANG3 gives potential degree101, marked-collision Meng-Yang graph degree48, degree7 polynomial pullbacks with collision jet match, and a finite constant-covector certificate. HCTANG1 separately adds sharp rational pole bounds. The later local HC4MYGALL1 corollary excludes the entire standard polynomial-graph collision chart in every degree via the elementary HC4MYGAP1 obstruction, without a plane degree import. HC4MYGT2 additionally excludes both arbitrary-degree triangular dual orders with base-dependent coefficients and every base-independent plane Keller dual substitution, using the elementary HC4MYGT1 theorem. Neither conjecture is settled; no formal verification or external review.'),
 ('hexagon:2611.00161v1', 'The four-dimensional Hessian conjecture for potentials of degree at most five', ['David Niedbala Giraudin'], '2026-11-06', 'https://hexagonmath.org/pdf/2611.00161v1', 'reproduced follow-up', 'Full manuscript and source inspected; 247 main checks replayed, normalization and stabilizer coverage reviewed after initial literature triage', 'HC4D5EXT1 closes OP-HC4-D5 by reproduction, together with 29 dependency checks. HC4D5Z1 adds an independently constructed complete 35-coefficient unit ideal. Explicitly unrefereed; no full independent reimplementation, formal verification or external review. Not arbitrary-degree HC4 or JC2.'),
 ('hexagon:2611.00134v1', 'The four-dimensional Hessian conjecture for quintic potentials with a ternary cone top', ['David Niedbala Giraudin'], '2026-11-05', 'https://hexagonmath.org/pdf/2611.00134v1', 'reproduced dependency', 'Full manuscript and source inspected; all 29 checks replayed and written ternary normalizations reviewed after initial triage', 'Ternary quintic-top dependency used by HC4D5EXT1. The next-day paper supplies binary and unary cases. These are dependent parts of one reproduction, not independent complete proofs.'),
 ('orevkov-germs-2026', 'On germs of mappings C2 to C2', ['S. Yu. Orevkov'], '2026', 'https://www.math.univ-toulouse.fr/~orevkov/k-en.pdf', 'high', 'section 1 hypotheses and author publication record', 'Published Sbornik Math. 217(2026), 572-582. Requires a smooth analytic source germ and bijectivity on every component over the branch curve. Does not identify which unramified components are deleted from an affine source.'),
 ('2503.07099v1', 'On irreducible germs of generic morphisms', ['Vik. S. Kulikov'], '2025-03-10', 'https://arxiv.org/abs/2503.07099v1', 'supporting prior work', 'metadata and abstract', 'Local branch-germ classification underlying the 2026 Orevkov addendum; not a complete global boundary theorem.'),
 ('2407.13795v1', 'There are no Keller maps having prime degree field extensions', ['Vered Moskowicz'], '2024-07-16', 'https://arxiv.org/abs/2407.13795v1', 'claim requiring audit', 'metadata and abstract; no newer revision listed', 'Prime extension-degree exclusion is a preprint claim, not an independently audited theorem in this pass; not a new 2026 result.'),
 ('zenodo:21479814', 'Exact Computer-Assisted Exclusion of the (72,108) Frontier in the Two-Dimensional Jacobian Problem', ['Billel Helali'], '2026-07-21', 'https://doi.org/10.5281/zenodo.21479814', 'retained baseline', 'existing pinned repository provenance and fresh discovery search; live Zenodo unavailable', 'Local exact replay and transcription audit already retained in plane-jc/PROVENANCE.md. External record scopes its degree-125 consequence to the imported GGHV reduction. No higher frontier established by this scan.'),
 ('2204.14178', 'Increasing the degree of a possible counterexample to the Jacobian Conjecture from 100 to 108', ['Jorge Alberto Guccione', 'Juan Jose Guccione', 'Rodrigo Horruitiner', 'Christian Valqui'], '2022', 'https://arxiv.org/abs/2204.14178', 'retained baseline', 'repository pinned provenance; not freshly re-audited', 'External reduction behind the retained (72,108) exclusions. Keep its exact hypotheses when comparing bounds.'),
 ('doi:10.2989/16073606.2026.2701437', 'The lower side of the Newton polygon of hypothetical counterexamples to the plane Jacobian conjecture', ['Jorge A. Guccione', 'Juan J. Guccione', 'Christian Valqui'], '2026', 'https://cris.pucp.edu.pe/entities/publication/4ac4fae4-49ed-421c-a246-4d9534888bd6', 'medium', 'institutional bibliographic record and abstract; manuscript not audited', '2026 journal appearance of Newton-polygon restrictions related to the older 1605.09430 work. Publication year alone does not establish a new degree frontier or an unchanged text.')
]
for key, title, authors, published, source, priority, depth, relevance in manual:
    rows[key] = dict(id=key, title=title, authors=authors, published=published, source=source,
                     priority=priority, review_depth=depth, repository_relevance=relevance,
                     independent_proof_audit_in_this_pass=False)

for key in ('hexagon:2611.00161v1', 'hexagon:2611.00134v1'):
    rows[key]['followup_audit'] = '../../hessian/hc4-degree-five-audit/README.md'
    rows[key]['followup_status'] = 'Reproduced with written coverage review; HC4D5EXT1 closes OP-HC4-D5'

rows['hexagon:2611.00156v1']['followup_audit'] = '../../hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md'
rows['hexagon:2611.00156v1']['followup_status'] = 'Written structural reproduction HCTANG2 and quantitative corollary HCTANG3; geometric coverage is not certified by the finite checker'

for row in rows.values():
    row.setdefault('repository_relevance', 'Screened as adjacent background; no result imported into the mathematical-status ledger.')
catalog = {'as_of': '2026-11-13', 'purpose': 'Literature triage before further Jacobian calculations',
           'mathematical_status_authority': '../../MATH_STATUS.json',
           'independent_audit_note': 'The initial literature triage promoted no new external theorem. Follow-up reproduction HC4D5EXT1 closes OP-HC4-D5, and the complete written Tang audit registers HCTANG2/HCTANG3. Geometric coverage is written, bulk Giraudin algorithms are not independently reimplemented, and no formal or external review is claimed. Earlier de Bondt and Tang scoped additions are linked separately.',
           'records': sorted(rows.values(), key=lambda r: r.get('updated', r['published']), reverse=True)}
(HERE / 'catalog.json').write_text(json.dumps(catalog, indent=2, ensure_ascii=False) + '\n')
ns = {'a': 'http://www.w3.org/2005/Atom', 'o': 'http://a9.com/-/spec/opensearch/1.1/'}
snapshots = []
for path in sorted(HERE.glob('*.xml')):
    root = ET.fromstring(path.read_bytes())
    entries = root.findall('a:entry', ns)
    snapshots.append(dict(file=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                          returned=len(entries), total_results=root.findtext('o:totalResults', None, ns),
                          oldest_returned_update=min(e.findtext('a:updated', '', ns) for e in entries),
                          newest_returned_update=max(e.findtext('a:updated', '', ns) for e in entries),
                          self_links=[l.attrib.get('href') for l in root.findall('a:link', ns) if l.attrib.get('rel') == 'self']))
(HERE / 'coverage.json').write_text(json.dumps(dict(as_of='2026-11-13', catalog_records=len(rows),
    search_terms=['Jacobian Conjecture', 'Keller map', 'Hessian Conjecture', 'Dixmier Conjecture',
                  'Druzkowski / Druzkowski with diacritics', 'nonproperness and polynomial',
                  'etale and endomorphisms', 'vanishing conjecture and Zhao',
                  'Generalized Vanishing Conjecture'],
    date_strategy='2026 submissions; separately all-years records sorted by last update to catch older papers revised in 2026; targeted November primary-source searches outside arXiv.',
    snapshots=snapshots, exclusions=['2610.25831: mean-curvature self-shrinkers', 'Related-query generic-vanishing and nilpotent-orbit false positives'],
    limitations=['Search coverage is not a proof of bibliographic exhaustiveness.',
                 'The all-years query returns the 200 most recently updated records; 2026 updates are filtered locally.',
                 'One compound adjacent-topic API query returned HTTP 500 twice; narrower generalized-VC query and primary web searches supplied the fallback.',
                 'Live Zenodo record/API access failed; existing pinned provenance remains the baseline.',
                 'Initial Hexagon abstract routes failed; follow-up browser-UA downloads supplied HTML and complete source archives, pinned in hc4-degree-five-audit.',
                 'Initial triage executed no ancillary programs; the follow-up audit actually replayed all 276 degree-five and ternary-dependency checks.']), indent=2) + '\n')
print(f'Catalog: {len(rows)} records; {len(snapshots)} hashed API snapshots')

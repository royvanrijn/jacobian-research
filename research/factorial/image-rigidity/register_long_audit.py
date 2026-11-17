"""Record the completed internal Long audit without rewriting unrelated rows."""
from pathlib import Path
import hashlib
import json
import os

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = 'factorial/image-rigidity/'
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()


def update_array(path, key, changes, additions=()):
    before = path.read_bytes()
    text = before.decode()
    data = json.loads(text)
    known = {r['id']:r for r in data[key]}
    def serialized(row):
        return '\n'.join('    '+line for line in json.dumps(row,indent=2).splitlines())
    for id, change in changes.items():
        old = known[id]
        row = change(dict(old))
        if row != old:
            needle = serialized(old)
            assert text.count(needle) == 1, ('noncanonical formatting',id)
            text = text.replace(needle,serialized(row),1)
            known[id] = row
    new = []
    for row in additions:
        if row['id'] in known:
            assert known[row['id']] == row, ('changed existing audit row',row['id'])
        else:
            new.append(row)
    if new:
        head,end = text.rsplit('\n  ]\n}',1)
        assert not end.strip()
        text = head+',\n'+',\n'.join(serialized(row) for row in new)+'\n  ]\n}'+end
    if text.encode() != before:
        assert path.read_bytes() == before, 'Concurrent registry edit'
        temporary = path.with_suffix(path.suffix+'.long-audit-research/.scratch')
        temporary.write_bytes(text.encode())
        os.replace(temporary,path)


if any(row['id']=='FIRPENCILGAMMA20261116' for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']):
    raise SystemExit('Historical binary-audit registration: later pencil scopes are maintained by register_simplex_pencils.py; refusing to overwrite them.')

ledger_path = HERE/'source_ledger.json'
ledger = json.loads(ledger_path.read_text())
for row in ledger['cached_primary_sources']:
    if row['file']=='sources/long-fc2.tex':
        row['review'] = 'Complete main-text internal proof audit, including corrected leaf extraction/transport, rational path module, all-block traces, conjugate arithmetic pairing, endpoints, flat limit and algebraic specialization. No external or formal verification.'
    if row['file']=='sources/long-eit-pinned.tex':
        row['review'] = 'Complete internal audit of the algebraic-value rigidity and single-integral nonvanishing proof: moment E-functions, Fourier semisimplicity, constant rational projectors, Beukers lifting, exponential-polynomial reduction and tree Cauchy jumps. Published inputs matched. No amplitude generalization.'
additional = [
    ('long-eit-vi-current.tex','https://raw.githubusercontent.com/long-mathematics/exponential-integral-theorem/92db79eedeb064694ff20e9201d73fc2bc0857f4/papers/06-exponential-factorial-moments/exponential-factorial-moments.tex',
     'Scope and theorem inventory read; higher-dimensional family nonvanishing inputs from Paper V remain unaudited. No theorem promotion.'),
    ('beukers-math0405549.html','https://arxiv.org/html/math/0405549',
     'Theorem 3.2 and Lemma 3.1 read and matched to EIT ordinary-point linear lifting with prescribed specialization; no independence hypothesis added.'),
    ('pakovich-muzychuk-07104085v2.html','https://arxiv.org/html/0710.4085v2',
     'Theorem 3.1 and Section 4.1/Proposition 4.1 read and matched to rational top-component irreducibility and actual endpoint path-module use. No absolute irreducibility assumed.'),
]
existing = {r['file'] for r in ledger['cached_primary_sources']}
for name,url,review in additional:
    if 'sources/'+name not in existing:
        ledger['cached_primary_sources'].append(dict(file='sources/'+name,url=url,sha256=digest(HERE/'sources'/name),review=review))
for url in ['https://arxiv.org/abs/0712.0349','https://arxiv.org/html/math/0405549','https://arxiv.org/html/0710.4085v2']:
    if url not in ledger['additional_primary_references']:
        ledger['additional_primary_references'].append(url)
ledger['published_input_pinpoints'] = {
    'Beukers2006':'Theorem3.2 + Lemma3.1: all linear relations at a nonzero ordinary algebraic point lift with exact prescribed specialization; dependent coordinate lists allowed.',
    'deCataldoMigliorini2009':'Remark4.2.4: small proper smooth-source map gives the intermediate extension of its covering local system. Finite curve map has zero-dimensional fibers and satisfies the strict nondense-stratum inequality.',
    'PakovichMuzychuk2009v2':'Theorem3.1, Section4.1, Proposition4.1: rational irreducible top periodic component and actual endpoint path module. The result over Q is not substituted by a false absolute-complex irreducibility claim.',
    'standard_module_theory':'Riemann-Hilbert for regular holonomic modules, exact Fourier Weyl automorphism and generic localization, checked for the displayed finite polynomial direct image.',
}
ledger_path.write_text(json.dumps(ledger,indent=2)+'\n')

audit = {
    'date':'2026-11-16',
    'previous_goal_turn':'progress',
    'previous_turn_evidence':'Revalidated existing rigidity proof, adapters, registered entries and 152-control receipt against current files before continuing.',
    'current_turn':'progress',
    'canonical_proof':BASE+'LONG_EIT_AND_FC2_AUDIT.tex',
    'source_sha256':{r['file']:r['sha256'] for r in ledger['cached_primary_sources']},
    'eit_main_labels_checked':['lem:D-spectrum','prop:moment-system','thm:Beukers','prop:semisimplicity','prop:endomorphism-rigidity','cor:constant-splitting','prop:exp-poly','lem:tree-representation','thm:moment-rigidity','thm:rigidity','cor:single'],
    'fc2_main_labels_checked':['lem:specialization','lem:kernel','lem:normality','lem:local-extraction','thm:transport','lem:zero-leaf','thm:rational-input','lem:dichotomy','lem:averaging','lem:separation','prop:path-module','cor:annihilation','lem:trace','prop:arithmetic','lem:pure-power','lem:endpoints','thm:nonflat','lem:radial-limit','cor:fc2'],
    'conclusion':'Complete internal verification of the pinned main proof chains, with published inputs treated as theorems and their hypotheses/coordinate conventions checked. No appendix or finite grid is used to replace an analytic argument.',
    'new_transfer':'Two-exponential/block-radial pure eventual radical is zero over fields; over every Q-algebra it equals coefficient nilpotence. Scaled actual-period normalization is exactly ordinary coefficient torsion Gamma_J.',
    'scope_boundaries':['General FC in dimension >=3 unresolved','Unrestricted SFC unresolved','Paper VI/Paper V higher-dimensional family audit incomplete','Complete bibliography audit incomplete','Arbitrary weighted Gamma FC not inferred','Nonzero limit at every transcendental coefficient choice not inferred'],
    'independent_replay':False,'formal_verification':False,'external_review':False,'novelty_claim':False,
}
(HERE/'long_source_audit.json').write_text(json.dumps(audit,indent=2)+'\n')

checker = BASE+'verify_long_audit.py'
proof = BASE+'LONG_EIT_AND_FC2_AUDIT.tex'
locks = [proof,checker,BASE+'long_audit_verification.json',BASE+'long_source_audit.json',BASE+'source_ledger.json',BASE+'sources/long-fc2.tex',BASE+'sources/long-eit-pinned.tex']
common = dict(state='proved',canonical_source=proof,checker=checker,proof_type='reproduction',
              independent_replay=False,formal_verification=False,external_review=False,
              artifact_hash='sha256:'+digest(HERE/'verify_long_audit.py'),software_lock=locks,
              supersedes=[],closes_problems=[],narrows_problems=['OP-FACTORIAL-GENERAL20261116'],
              consumers=[],invalidates_assumptions=[],replaced_by=[],priority='reference',programme_status='archived')
new = [
    dict(common,id='FIREIT20261116',kind='reproduction',title='The pinned Long--Le Blanc EIT I algebraic-value rigidity and single-integral nonvanishing proof passes internal verification',
         scope='All algebraic polynomial phases q on the affine line and finite algebraic degree-zero boundary cycles, fixed amplitude dx. Internal full main-proof audit shows I(1) algebraic implies I identically zero, hence every nonzero member is transcendental/nonzero at every nonzero algebraic scale. Verified moment E-function arithmetic, finite small direct-image/Fourier semisimplicity and signs, constant rational-projector spectra even for mixed multiple critical points sharing a value, Beukers lifting and constant relation spaces, rational-to-polynomial coefficient reduction and tree Cauchy jumps. Published inputs explicitly credited and matched, not reproved. Arbitrary amplitudes are excluded by an exact counterexample. Finite controls check transcription only; no external refereeing or formal assurance.',
         dependencies=['external: Long--Le Blanc EIT I pinned commit22d98c4f27e606b5b25fa673c33f9c56c53cc86b','external: Beukers2006 Theorem3.2 and Lemma3.1','external: deCataldo--Migliorini2009 Remark4.2.4 and standard Riemann-Hilbert/Fourier equivalences']),
    dict(common,id='FIRFC2NONFLAT20261116',kind='reproduction',title='Long binary factorial nonflat nonvanishing is internally verified with full leaf transport and rational block traces',
         scope='All nonzero complex binary polynomials whose top homogeneous part is not c*(x+y)^D, every total degree. The complete pinned main proof is checked: Mittag-Leffler kernel and normal convergence, correct positive Volterra correction, full root/path histories, absence of negative Puiseux history terms at nonzero critical values, zero-leaf redundancy, rational top path module, simultaneous unweighted block traces, maximal-ramification weighted coefficient and positive conjugate Lindemann pairing, and final endpoint pure-power exclusion. Algebraic specialization preserves actual total/projective degrees. Some factorial moment is nonzero and a power argument yields infinitely many. No genericity, absolute irreducibility, monodromy-preservation-under-specialization, finite-prefix argument or appendix multiplicity bound is substituted for the proof. Credited reproduction, no external/formal assurance.',
         dependencies=['external: Long FC2 commit8b1eb46ff530122ed25f588b78bf9fdc98cb9577','external: Pakovich--Muzychuk0710.4085v2 Theorem3.1 and Proposition4.1','external: Lindemann--Weierstrass theorem']),
    dict(common,id='FIRFC2FULL20261116',kind='reproduction',title='Long full binary Factorial proof is internally verified, including its pinned EIT dependency',
         scope='For every f in C[y1,y2], L(f^m)=0 for every positive m forces f=0. In fact eventual pure zero forces f=0, so every nonzero f has infinitely many nonzero factorial moments. Nonflat verification plus the all-layer dominated radial limit and audited EIT handle all binary top shapes; a constant phase is immediate. Algebraic specialization of a putative counterexample is essential: a transcendental flat coefficient choice can have zero leading radial limit without being a counterexample. The zero-integral kernel for product-exponential measure in two variables is consequently MZ with zero eventual radical. This is a complete internal written reproduction of Long, with standard published inputs, not an originality/peer-review/formalization claim. Unrestricted SFC and higher-dimensional FC are not implied.',
         dependencies=['FIRFC2NONFLAT20261116','FIREIT20261116','FIRFLATAFFINE20261116']),
    dict(common,id='FIRRADIAL2RINGS20261116',kind='corollary',proof_type='hybrid',title='Two independent exponential radial blocks have only nilpotent eventually pure-null phases over every Q-algebra',
         scope='Entire independent two-block Gaussian invariant polynomial subalgebra with blocks (X1^2+X2^2)/2 and (X3^2+X4^2)/2, arbitrary complex degree/support. Over every characteristic-zero field, eventual pure vanishing, all-marker safety and the zero phase coincide. Over ANY commutative Q-algebra B, eventual factorial pure-nullity, coefficientwise nilpotence, polynomial nilpotence and all-marker eventual safety are equivalent. Reduce at every prime, apply binary FC to the finite coefficient field, then use finite generation of the coefficient ideal to obtain nilpotence without any Noetherian assumption. All original Gaussian markers are retained; nilpotence supplies uniform-in-marker tails. All-positive pure-nullity over a nonreduced ring can still occur for epsilon*(y1-y2) with epsilon^2=0. This is not a theorem for arbitrary Gamma shapes, three radial blocks or general Gaussian MZ.',
         dependencies=['FIRFC2FULL20261116','FIRGAMMAADAPTER20261116']),
    dict(common,id='FIRRADIAL2NORM20261116',kind='corollary',proof_type='hybrid',title='Binary factorial and two-block radial actual-period normalization is exactly ordinary coefficient torsion',
         scope='For f in B[y1,y2], B ANY commutative Q-algebra, let J be its finite coefficient ideal. The scalars u making u*sum L(q*f^m)t^m/m! polynomial in t for EVERY fixed original marker q are exactly Gamma_J(B)=union_k Ann_B(J^k). Necessity already follows from the unit period after passing to B/Ann(u) and the coefficient-nilpotence theorem; sufficiency is uf^m=0 for m>=k. The same result retains every Gaussian marker after radial substitution. Ordinary, not reduced, torsion matters: in Q[s,t]/(s^2*t^2), f=s*y1 has normalization ideal(t^2), while t is not admissible. Finite nilpotence exponents yield a concrete uniform tail bound. This is a structural special-sector consequence of audited binary FC, and does not claim a new general Gaussian family theorem.',
         dependencies=['FIRRADIAL2RINGS20261116']),
]
new_ids = [r['id'] for r in new]
def change_open(row):
    row['title'] = 'Resolve higher-dimensional Factorial and unrestricted Strong Factorial and complete the literature audit'
    row['scope'] = 'The full user programme remains open: prove or refute general FC (dimensions>=3 remain) and unrestricted SFC; complete the current primary-source literature and higher-dimensional-family audits; exploit Gaussian/Gamma/Dirichlet certificates without a false purity-to-safety inference. Full rigidity and its special SFC family are internally verified. The complete pinned Long binary FC/EIT main proof chains are now internally verified by FIREIT20261116/FIRFC2NONFLAT20261116/FIRFC2FULL20261116, with external/formal assurance false. FIRRADIAL2RINGS20261116/FIRRADIAL2NORM20261116 classify the entire two-exponential sector over all Q-algebras, retaining ordinary nilpotents. General triangle integration is homogeneous FC3 and remains unresolved; newly located EIT VI/Paper V family inputs and the remaining bibliography require review. Unrestricted Image is already false by SIC2C4. No finite-prefix or diagonal-shadow shortcut, unrestricted Gamma-shape theorem, nonzero radial limit at all transcendental coefficients, bibliography completeness, peer review or full goal completion is claimed.'
    row['dependencies'] = list(dict.fromkeys(row['dependencies']+new_ids))
    return row
update_array(ROOT/'MATH_STATUS.json','entries',{'OP-FACTORIAL-GENERAL20261116':change_open},new)
known = {r['id']:r for r in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']}
fingerprint = lambda id:hashlib.sha256('\0'.join(known[id][k] for k in ('state','title','scope')).encode()).hexdigest()[:16]
def refresh(row):
    row['claims'] = list(dict.fromkeys(row['claims']+new_ids))
    row['reviewed_claims'] = {id:fingerprint(id) for id in row['claims']}
    row['sources'] = list(dict.fromkeys(row['sources']+[proof,BASE+'long_source_audit.json']))
    return row
def change_work(row):
    row = refresh(row)
    row['title'] = 'Resolve general higher-dimensional FC/SFC after the completed binary FC/EIT and rigidity audits'
    row['next_step'] = 'Audit the newly located EIT VI families and exact Paper V nonvanishing inputs, finish primary bibliography including Rocks, build original-unit period certificates for purity on unsafe phases in three or more radial blocks, and seek a universal radial-sector structure theorem. Full binary FC/EIT and full rigidity have complete internal written audits; do not repeat them or finite R3 windows without a specific mathematical defect.'
    row['prerequisites'] = 'Read the completed Long audit and source ledger, binary all-ring/ordinary-torsion corollaries, original Gaussian adapter, SIC2C4 and FM1. Retain analytic-proof versus finite-check assurance, dx-only EIT scope, ordinary nilpotents, absence of a universally nonzero flat limit for transcendental coefficients, and the nonmultiplicativity of diagonal extraction. No external messaging or publishing authorized.'
    return row
def change_lesson(row):
    row = refresh(row)
    row['use'] = 'Use the exact Gamma/Gaussian marker transfer and credited full rigidity proof. The completed internal Long FC2/EIT audit now classifies the entire two-exponential radial sector: zero over fields, coefficient nilpotence over all Q-algebras, and actual-period scalar normalization Gamma_J(B) with ordinary torsion retained.'
    caution = ' Do not assert that every transcendental flat coefficient choice has a nonzero leading radial limit, extend dx-only EIT to arbitrary amplitudes, or replace the ordinary normalization ideal by its radical.'
    if caution not in row['avoid']:
        row['avoid'] += caution
    row['boundary'] = 'General FC in dimension>=3, unrestricted SFC, full higher-dimensional-family/Paper V review and the complete bibliography remain open. Full binary FC/EIT main chains and rigidity now have complete internal written verification; external refereeing and formal verification remain false. Uniform triangles are homogeneous FC3, not the solved homogeneous binary class. Arbitrary Gamma shapes are not covered by the two-exponential classification.'
    row['implementation'] = list(dict.fromkeys(row['implementation']+[checker]))
    return row
update_array(ROOT/'knowledge/work_items.json','actions',{'WORK-FACTORIAL-GENERAL':change_work})
update_array(ROOT/'knowledge/lessons.json','lessons',{'METHOD-FACTORIAL-RADIAL-PURITY':change_lesson})
print(json.dumps({'registered':new_ids,'updated_open_scope':'OP-FACTORIAL-GENERAL20261116','goal_complete':False}))

#!/usr/bin/env python3
import hashlib, json, re, sys
from pathlib import Path

ROOT=Path('.')
STATUS=ROOT/'STATUS.md'; INDEX=ROOT/'MATH_STATUS.json'
idx=json.loads(INDEX.read_text()); entries=idx['entries']; by={e['id']:e for e in entries}

# Reconstruct the Fermigier entries from the already-merged generated status rows.
rows={}
for line in STATUS.read_text().splitlines():
    if line.startswith('| EC-F'):
        cols=[c.strip() for c in line.strip().strip('|').split(' | ')]
        if cols[0] in {'EC-FG12','EC-FXPT1','EC-FXPT2','EC-FXPT3','EC-FXPT4','EC-FXPT5','EC-FXPT6'}:
            rows[cols[0]]=cols
for id in ['EC-FG12','EC-FXPT1','EC-FXPT2','EC-FXPT3','EC-FXPT4','EC-FXPT5','EC-FXPT6']:
    c=rows[id]; src=re.search(r'\]\(([^)]+)\)',c[3]).group(1); checker=re.search(r'\]\(([^)]+)\)',c[5]).group(1)
    deps=[x.strip('` ') for x in c[4].split(', ')] if c[4]!='—' else []
    proof=c[6].split(';',1)[0]
    e={'id':id,'kind':'theorem' if id=='EC-FG12' else 'example','state':'proved','title':c[1],'scope':c[2],
       'canonical_source':src,'dependencies':deps,'checker':checker,'proof_type':proof,'independent_replay':False,
       'formal_verification':False,'external_review':False,'artifact_hash':'sha256:'+hashlib.sha256((ROOT/checker).read_bytes()).hexdigest(),
       'software_lock':['.python-version','requirements.txt'],'supersedes_notes':[],'supersedes':[],'closes_problems':[],
       'narrows_problems':[],'consumers':[],'invalidates_assumptions':[],'replaced_by':[],
       'priority':'primary' if id=='EC-FG12' else 'reference'}
    if id=='EC-FXPT3': e['replaced_by']=['EC-FXPT4']
    if id=='EC-FXPT4': e['supersedes']=['EC-FXPT3']
    by[id]=e

mr13=json.loads(r'''{"id":"EC-MR13","kind":"theorem","state":"proved","title":"A low-height six-root Mestre family of generic rank at least 13","scope":"For the six centers (0,23,93,128,133,175), exact square approximation gives a primitive Mestre quartic with fixed square content 2400^2, twelve paired-root sections, and six verified nonvisible linear-companion identities. The leading coefficient is 9*(T^2+14406), so the base change T=(14406-u^2)/(2u) supplies a rational split-infinity section. At u=1, the twelve paired-root images, one chosen companion, and split infinity have stacked good-reduction images of rank thirteen over F_3; reduction at p=7 excludes rational 3-torsion. Infinite descent proves the selected specialized points independent, and specialization therefore proves generic rank at least thirteen over Q(u). The primitive discriminant frontier is irreducible and squarefree of degree 20 in T and degree 40 after base change, compared with degree 398 for the reproduced Kihara rank-at-least-14 family. The degree comparison is conductor geometry, not an exact conductor formula, and neither generic rank equality nor independence of all six companions is claimed.","canonical_source":"elliptic-curves/notes/MESTRE_RANK13_02393128133175.md","dependencies":["external: good-reduction specialization homomorphism for elliptic curves","external: finite-reduction and infinite-descent independence criterion"],"checker":"elliptic-curves/cas/verify_mestre_rank13_02393128133175.py","proof_type":"exact_symbolic","independent_replay":true,"formal_verification":false,"external_review":false,"artifact_hash":"sha256:918756aa04d2f6716548ef9228125fbe07a30177cb8b8f161f53796d54701bde","software_lock":[".python-version","requirements.txt"],"supersedes_notes":[],"supersedes":[],"closes_problems":[],"narrows_problems":["OP-EC-RANK-CONDUCTOR"],"consumers":[],"invalidates_assumptions":[],"replaced_by":[],"priority":"primary"}''')
md17=json.loads(r'''{"id":"EC-MD17","kind":"example","state":"proved","title":"Split-infinity Mestre family specialization of rank at least 17 below the conductor cutoff","scope":"For the six-root Mestre family with centers (0,25,95,143,168,205), the primitive quartic has leading coefficient T^2+39146, so T=(39146-u^2)/(2u) supplies split infinity. An exact companion-section certificate proves this base-changed family has generic rank at least 13. Its primitive discriminant core is irreducible and squarefree of degree 20 in T and degree 40 after base change, compared with degree 398 for the reproduced Kihara family; this is conductor geometry, not an exact conductor formula. At u=197, hence T=337/394, seventeen pinned rational points on the specialized short Jacobian have full F_3 column rank in the stacked good-reduction quotients at p=37,41,61,67,79,83,101,103,137,139,149,163,167,173,181,193. The good-reduction order 20 at p=13 excludes rational 3-torsion, so infinite descent proves rank E(Q)>=17. PARI/GP independently reconstructs the minimal model [1,1,1,-1163348683373499147707371416562962,15227131493689013260364706485730874765958430844575], root number -1, and conductor 2462086522751621334987931952469307556796057284118717977320345864383117775914. The exact inequalities N<10^76 and log(10)<231/100 give log(N)<4389/25<182.72. This is a rank lower bound only, four points short of the operational rank-21 target; no saturation or upper bound is claimed. The bounded four-family discovery screen and its conductor timeouts are not upper-bound evidence.","canonical_source":"elliptic-curves/notes/MESTRE_DSQUARE_FOUR_SCREEN.md","dependencies":["external: finite-reduction and infinite-descent independence criterion","external: good-reduction specialization homomorphism for elliptic curves","external: PARI/GP 2.15.4 global minimization, conductor, and root-number routines"],"checker":"elliptic-curves/cas/verify_mestre_dsquare_four_u197.py","proof_type":"hybrid","independent_replay":false,"formal_verification":false,"external_review":false,"artifact_hash":"sha256:68bb4137f8fdd642df33833f790947560ed23c035d6e024c99ffa0266166e188","software_lock":[".python-version","requirements.txt"],"supersedes_notes":[],"supersedes":[],"closes_problems":[],"narrows_problems":["OP-EC-RANK-CONDUCTOR"],"consumers":[],"invalidates_assumptions":[],"replaced_by":[],"priority":"reference"}''')
by['EC-MR13']=mr13; by['EC-MD17']=md17

for id in ['EC-CRT1','EC-R20','EC-FSEED1']:
    p=ROOT/by[id]['checker']; by[id]['artifact_hash']='sha256:'+hashlib.sha256(p.read_bytes()).hexdigest()
by['EC-R20']['scope']='At canonical identity fermigier-mestre-v1:u=28917/20 (whose legacy literal-shift alias is T=28917/10), a bounded ratpoints search supplies a rational point cloud from which exact finite-reduction arithmetic selects 20 independent points. PARI/GP gives global minimal model [1,1,1,-4437412060110743641525245114305,3586842216822165612930264910099076801587288127] and exact conductor 2876153493562761211278364526603564191699143885403233935132057708367930, whose natural logarithm is 159.9348252255254533984...<182.72. The stored points and conductor are exact, but the search is bounded: no twenty-first point, global saturation, or unconditional rank upper bound is proved. This is one point short and is not a target solution.'
op=by['OP-EC-RANK-CONDUCTOR']
op['scope']='Find either an elliptic curve over Q with at least 21 rigorously independent rational points and natural-log conductor below the literal cutoff 182.72, or a curve with at least 30 rigorously independent rational points. EC-CRT1 establishes the exact one-parameter Hensel--CRT--Gauss mechanism; EC-FERM1 and EC-FRANK1 supply the canonical Fermigier adapter, all thirteen quartic points, and E22 lower-bound replay while retaining the unresolved source normalization. EC-FG12 proves that this adapter has arithmetic generic rank exactly twelve; EC-FXPT1 closes the declared affine, quadratic, Mobius, and direct pair-product transports between its rank-22 and rank-20 anchors; EC-FXPT2 closes every affine transport between signed support-at-most-two quotient representatives; EC-FXPT4 closes every rational discriminant component across all 80 bidegree-(2,1) pencils; EC-FXPT5 finds no third simultaneous-square parameter on any genuine pair cover through projective height 200000; and EC-FXPT6 finds no point through projective height 1024 on the pilot nonlinear discriminant component. EC-R20 remains the strongest new conductor-qualified survivor: rank at least 20 with log conductor 159.9348252255..., explicitly one point short. EC-MR13 supplies a new generic-rank-at-least-13 Mestre family with degree-40 discriminant geometry, and EC-MD17 supplies a distinct generic-rank-at-least-13 family together with a conductor-qualified rank-at-least-17 specialization at u=197. EC-E29 locally replays the public lower bound 29 but supplies no thirtieth point, while EC-KIH14 remains the strongest fully public higher-generic-rank fallback reproduced here. No target curve has been found. Immediate gates are rational points beyond the completed boxes on the nonlinear discriminant or genus-nine Fermigier cover curves, higher-support exceptional transport, four further independent directions on EC-MD17 or comparable low-frontier Mestre fibers, two further dimensions on the odd-root-number rank-19 survivors, and access to or reconstruction of the unpublished rank-17 fibration behind E29. Exact rank 21 would additionally require a rigorous upper bound; the operational lower-bound targets do not.'
for id in ['EC-FG12','EC-FXPT1','EC-FXPT2','EC-FXPT3','EC-FXPT4','EC-FXPT5','EC-FXPT6','EC-MR13','EC-MD17']:
    if id not in op['dependencies']: op['dependencies'].append(id)

order=[]
for e in entries:
    if e['id'] in {'EC-FG12','EC-FXPT1','EC-FXPT2','EC-FXPT3','EC-FXPT4','EC-FXPT5','EC-FXPT6','EC-MR13','EC-MD17'}: continue
    if e['id']=='EC-R20':
        order += [by[x] for x in ['EC-FG12','EC-FXPT1','EC-FXPT2','EC-FXPT3','EC-FXPT4','EC-FXPT5','EC-FXPT6']]
    if e['id']=='EC-E29': order += [by['EC-MR13'],by['EC-MD17']]
    order.append(by[e['id']])
idx['entries']=order
INDEX.write_text(json.dumps(idx,indent=2,ensure_ascii=False)+'\n')

note=ROOT/'elliptic-curves/notes/FERMIGIER_REPRODUCTION.md'; text=note.read_text()
marker='## Search and evidence boundary\n'
bridge=r'''## Constructive rank-gain bridge

A cover that merely makes an existing Mordell--Weil section divisible does
**not** increase rank. If a new section \(S\) satisfies \([n]S=R\) for an old
section \(R\), then \(S=R/n\) after tensoring the Mordell--Weil group with
\(\mathbf Q\); such a cover can improve saturation, but supplies no new
rational direction.

For a quadratic base change \(K=\mathbf Q(u)(\sqrt{d(u)})\), the useful
rank-gain test is instead the anti-invariant summand. Over \(\mathbf Q\),

\[
\operatorname{rank} E(K)=\operatorname{rank} E(\mathbf Q(u))
+\operatorname{rank} E^{(d)}(\mathbf Q(u)),
\]

where \(E^{(d)}\) is the quadratic twist. Thus a productive low-degree cover
must create a genuinely new non-torsion point on the twist, not a half of a
known section. The practical search order is: force low-height squareclass
conditions to genus zero or one; parameterize/solve them; reject candidates
whose good specialization lies in the old twelve-dimensional span; then
certify survivors by exact finite reductions and score conductor geometry
before expensive point searches.

This matches explicit high-rank constructions based on quadratic sections and
conic conditions, and the elliptic-surface/K3 strategy of treating the
Mordell--Weil/Neron--Severi lattice plus good specializations as primary search
objects. For the present target, conductor should remain part of the search
objective rather than a final filter.

References: N. D. Elkies, *Three lectures on elliptic surfaces and curves of
high rank*, arXiv:0709.2908; N. D. Elkies and M. Watkins, *Elliptic curves of
large rank and small conductor*, arXiv:math/0403374; A. Dujella,
M. Kazalicki, and J. C. Peral, *Elliptic curves with torsion groups Z/8Z and
Z/2Z x Z/6Z* (quadratic-section construction in Sec. 1).

'''
if '## Constructive rank-gain bridge\n' not in text:
    if marker not in text: raise SystemExit('note insertion marker missing')
    note.write_text(text.replace(marker,bridge+marker))

sys.path.insert(0,str(ROOT/'scripts')); import render_status
STATUS.write_text(render_status.render(idx))
print('ELLIPTIC_CONSOLIDATION_REPAIR_WRITTEN')

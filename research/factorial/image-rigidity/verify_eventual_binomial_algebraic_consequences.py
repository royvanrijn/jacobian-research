#!/usr/bin/env python3
"""Exact/source support for ordinary CI, nilpotence and adjunction consequences."""
from __future__ import annotations
if not __debug__:raise SystemExit("Refusing -O/-OO: exact assertions must execute.")
from fractions import Fraction as Q
from pathlib import Path
from math import factorial,comb,isqrt
import argparse,copy,hashlib,json,re,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
PN="EVENTUAL_BINOMIAL_ALGEBRAIC_CONSEQUENCES.md"
CN="verify_eventual_binomial_algebraic_consequences.py"
IN="eventual_binomial_algebraic_inputs.json"
LN="eventual_binomial_algebraic_source_ledger.json"
RN="eventual_binomial_algebraic_verification.json"
PREFIX="research/factorial/image-rigidity/"
FIXTURES=[[0, 1, 2], [0, 2, 3], [1, 2, 4], [2, 2, 5], [1, 3, 4], [3, 2, 3], [2, 4, 4], [4, 3, 3]]
ROLES={'checker': ('local', 'verify_eventual_binomial_algebraic_consequences.py'), 'eventual_adjacent_all_order_proof': ('repo', 'research/factorial/image-rigidity/UNIFORM_PADE_PARABOLA_SFC.tex'), 'eventual_constant_plus_power_analytic_proof': ('repo', 'research/factorial/image-rigidity/ALL_UNIVARIATE_TWO_TERM_SFC.tex'), 'eventual_gap405_coalescing_bridge_proof': ('repo', 'research/factorial/image-rigidity/COALESCING_BINOMIAL_SADDLE_BRIDGE.md'), 'eventual_gap405_complement_checker': ('repo', 'research/factorial/image-rigidity/verify_cusp_global_complement.py'), 'eventual_gap405_complement_inputs': ('repo', 'research/factorial/image-rigidity/cusp_global_complement_inputs.json'), 'eventual_gap405_complement_ledger': ('repo', 'research/factorial/image-rigidity/cusp_global_complement_source_ledger.json'), 'eventual_gap405_complement_proof': ('repo', 'research/factorial/image-rigidity/CUSP_GLOBAL_REGULAR_COMPLEMENT_AND_STOKES_DETOUR.md'), 'eventual_gap405_complement_receipt': ('repo', 'research/factorial/image-rigidity/cusp_global_complement_verification.json'), 'eventual_gap405_cusp_main_checker': ('repo', 'research/factorial/image-rigidity/verify_cusp_eventual_strong_factorial_completion.py'), 'eventual_gap405_cusp_main_inputs': ('repo', 'research/factorial/image-rigidity/cusp_eventual_strong_factorial_inputs.json'), 'eventual_gap405_cusp_main_ledger': ('repo', 'research/factorial/image-rigidity/cusp_eventual_strong_factorial_source_ledger.json'), 'eventual_gap405_cusp_main_proof': ('repo', 'research/factorial/image-rigidity/CUSP_EVENTUAL_STRONG_FACTORIAL_THEOREM.md'), 'eventual_gap405_cusp_main_receipt': ('repo', 'research/factorial/image-rigidity/cusp_eventual_strong_factorial_verification.json'), 'eventual_gap405_fixed_bounded_checker': ('repo', 'research/factorial/image-rigidity/verify_bounded_quartic_cusp.py'), 'eventual_gap405_fixed_bounded_proof': ('repo', 'research/factorial/image-rigidity/BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md'), 'eventual_gap405_fixed_bounded_receipt': ('repo', 'research/factorial/image-rigidity/bounded_quartic_cusp_verification.json'), 'eventual_gap405_fixed_completion_checker': ('repo', 'research/factorial/image-rigidity/verify_quartic_cusp_local_completion.py'), 'eventual_gap405_fixed_completion_inputs': ('repo', 'research/factorial/image-rigidity/quartic_cusp_local_completion_inputs.json'), 'eventual_gap405_fixed_completion_ledger': ('repo', 'research/factorial/image-rigidity/quartic_cusp_local_completion_source_ledger.json'), 'eventual_gap405_fixed_completion_receipt': ('repo', 'research/factorial/image-rigidity/quartic_cusp_local_completion_verification.json'), 'eventual_gap405_fixed_cusp_checker': ('repo', 'research/factorial/image-rigidity/verify_quartic_cusp_unfolding.py'), 'eventual_gap405_fixed_cusp_proof': ('repo', 'research/factorial/image-rigidity/QUARTIC_CUSP_UNFOLDING_AND_ACTUAL_ORDER_DRIFT.md'), 'eventual_gap405_fixed_cusp_receipt': ('repo', 'research/factorial/image-rigidity/quartic_cusp_unfolding_verification.json'), 'eventual_gap405_fixed_model_checker': ('repo', 'research/factorial/image-rigidity/verify_quartic_lee_yang_zero_and_robin_bridge.py'), 'eventual_gap405_fixed_model_inputs': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_zero_and_robin_inputs.json'), 'eventual_gap405_fixed_model_ledger': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_zero_and_robin_source_ledger.json'), 'eventual_gap405_fixed_model_proof': ('repo', 'research/factorial/image-rigidity/QUARTIC_LEE_YANG_ZERO_AND_ROBIN_BRIDGE.md'), 'eventual_gap405_fixed_model_receipt': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_zero_and_robin_verification.json'), 'eventual_gap405_fixed_primary_hjn': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_sources/hou_jiang_newman_2208_00917v1.pdf'), 'eventual_gap405_fixed_primary_newman': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_sources/newman_fourier_transforms_1976.pdf'), 'eventual_gap405_ledger': ('repo', 'research/factorial/image-rigidity/square_root_gap_eventual_source_ledger.json'), 'eventual_gap405_local_checker': ('repo', 'research/factorial/image-rigidity/verify_cusp_shrinking_window_log_ratio.py'), 'eventual_gap405_local_inputs': ('repo', 'research/factorial/image-rigidity/cusp_shrinking_window_inputs.json'), 'eventual_gap405_local_ledger': ('repo', 'research/factorial/image-rigidity/cusp_shrinking_window_source_ledger.json'), 'eventual_gap405_local_proof': ('repo', 'research/factorial/image-rigidity/CUSP_SHRINKING_WINDOW_LOG_RATIO_NONCOLLISION.md'), 'eventual_gap405_local_receipt': ('repo', 'research/factorial/image-rigidity/cusp_shrinking_window_verification.json'), 'eventual_gap405_main_checker': ('repo', 'research/factorial/image-rigidity/verify_square_root_gap_eventual_completion.py'), 'eventual_gap405_main_inputs': ('repo', 'research/factorial/image-rigidity/square_root_gap_eventual_inputs.json'), 'eventual_gap405_main_proof': ('repo', 'research/factorial/image-rigidity/SQUARE_ROOT_GAP_EVENTUAL_STRONG_FACTORIAL_THEOREM.md'), 'eventual_gap405_morse_checker': ('repo', 'research/factorial/image-rigidity/verify_cusp_two_thimble_positive_cycle.py'), 'eventual_gap405_morse_inputs': ('repo', 'research/factorial/image-rigidity/cusp_two_thimble_inputs.json'), 'eventual_gap405_morse_ledger': ('repo', 'research/factorial/image-rigidity/cusp_two_thimble_source_ledger.json'), 'eventual_gap405_morse_proof': ('repo', 'research/factorial/image-rigidity/CUSP_TWO_THIMBLE_POSITIVE_CYCLE.md'), 'eventual_gap405_morse_receipt': ('repo', 'research/factorial/image-rigidity/cusp_two_thimble_verification.json'), 'eventual_gap405_primary_vidunas': ('repo', 'research/factorial/image-rigidity/sources/vidunas-belyi-2025.pdf'), 'eventual_gap405_receipt': ('repo', 'research/factorial/image-rigidity/square_root_gap_eventual_verification.json'), 'eventual_gap405_retained_gamma_upper': ('repo', 'research/factorial/image-rigidity/GAMMA_PRODUCT_ROOT_MODULUS.tex'), 'eventual_gap405_stationary_checker': ('repo', 'research/factorial/image-rigidity/verify_cusp_stationary_belyi_monodromy.py'), 'eventual_gap405_stationary_inputs': ('repo', 'research/factorial/image-rigidity/cusp_stationary_belyi_inputs.json'), 'eventual_gap405_stationary_ledger': ('repo', 'research/factorial/image-rigidity/cusp_stationary_belyi_source_ledger.json'), 'eventual_gap405_stationary_primary_audit': ('repo', 'research/factorial/image-rigidity/cusp_stationary_belyi_primary_audit.md'), 'eventual_gap405_stationary_proof': ('repo', 'research/factorial/image-rigidity/CUSP_STATIONARY_BELYI_MONODROMY_AND_RADIAL_BARRIER.md'), 'eventual_gap405_stationary_receipt': ('repo', 'research/factorial/image-rigidity/cusp_stationary_belyi_verification.json'), 'eventual_gap405_strict_checker': ('repo', 'research/factorial/image-rigidity/verify_strict_subcritical_eventual_sfc.py'), 'eventual_gap405_strict_inputs': ('repo', 'research/factorial/image-rigidity/strict_subcritical_inputs.json'), 'eventual_gap405_strict_ledger': ('repo', 'research/factorial/image-rigidity/strict_subcritical_source_ledger.json'), 'eventual_gap405_strict_proof': ('repo', 'research/factorial/image-rigidity/SUBCRITICAL_ORDERED_PAIR_EVENTUAL_SFC.md'), 'eventual_gap405_strict_receipt': ('repo', 'research/factorial/image-rigidity/strict_subcritical_verification.json'), 'eventual_ledger': ('repo', 'research/factorial/image-rigidity/ordinary_binomial_eventual_source_ledger.json'), 'eventual_main_checker': ('repo', 'research/factorial/image-rigidity/verify_ordinary_binomial_eventual_completion.py'), 'eventual_main_inputs': ('repo', 'research/factorial/image-rigidity/ordinary_binomial_eventual_inputs.json'), 'eventual_main_proof': ('repo', 'research/factorial/image-rigidity/ALL_ORDINARY_BINOMIAL_EVENTUAL_STRONG_FACTORIAL.md'), 'eventual_novertical_checker': ('repo', 'research/factorial/image-rigidity/verify_no_vertical_physical_equimodular_pair.py'), 'eventual_novertical_inputs': ('repo', 'research/factorial/image-rigidity/no_vertical_pair_inputs.json'), 'eventual_novertical_ledger': ('repo', 'research/factorial/image-rigidity/no_vertical_pair_source_ledger.json'), 'eventual_novertical_proof': ('repo', 'research/factorial/image-rigidity/NO_VERTICAL_PHYSICAL_EQUIMODULAR_PAIR.md'), 'eventual_novertical_receipt': ('repo', 'research/factorial/image-rigidity/no_vertical_pair_verification.json'), 'eventual_receipt': ('repo', 'research/factorial/image-rigidity/ordinary_binomial_eventual_verification.json'), 'eventual_regular_checker': ('repo', 'research/factorial/image-rigidity/verify_above_critical_regular_complement.py'), 'eventual_regular_inputs': ('repo', 'research/factorial/image-rigidity/above_critical_regular_inputs.json'), 'eventual_regular_ledger': ('repo', 'research/factorial/image-rigidity/above_critical_regular_source_ledger.json'), 'eventual_regular_proof': ('repo', 'research/factorial/image-rigidity/ABOVE_CRITICAL_REGULAR_COMPACT_COMPLEMENT.md'), 'eventual_regular_receipt': ('repo', 'research/factorial/image-rigidity/above_critical_regular_verification.json'), 'eventual_small_prime_power_proof': ('repo', 'research/factorial/image-rigidity/PRIME_POWER_TWO_TERM_SFC.tex'), 'eventual_turn_checker': ('repo', 'research/factorial/image-rigidity/verify_simple_turn_shrinking_window.py'), 'eventual_turn_inputs': ('repo', 'research/factorial/image-rigidity/simple_turn_shrinking_window_inputs.json'), 'eventual_turn_ledger': ('repo', 'research/factorial/image-rigidity/simple_turn_shrinking_window_source_ledger.json'), 'eventual_turn_proof': ('repo', 'research/factorial/image-rigidity/SIMPLE_TURN_SHRINKING_WINDOW_CONSECUTIVE_EXCLUSION.md'), 'eventual_turn_receipt': ('repo', 'research/factorial/image-rigidity/simple_turn_shrinking_window_verification.json'), 'inputs': ('local', 'eventual_binomial_algebraic_inputs.json'), 'primary_rocks2015': ('repo', 'research/factorial/image-rigidity/sources/rocks-2015-author-preprint.pdf'), 'proof': ('local', 'EVENTUAL_BINOMIAL_ALGEBRAIC_CONSEQUENCES.md')}
SCOPE={'schema': 1, 'support': {'a_min': 0, 'd_min': 1, 'N': 'fixed N>=N0(a,d), N0>=2'}, 'coefficient_scheme': 'Q[A,B]/(mu_(N-1),mu_N)', 'CI_degrees': ['N-1', 'N'], 'CI_length': 'N*(N-1)', 'top_degree': '2*N-3', 'sharp_ideal_nilpotence_index': '2*N-2', 'sharp_phase_nilpotence_index_in_universal_quotient': True, 'ordinary_ideal_containments': True, 'radical_only': False, 'ring_scope': 'commutative Q-algebras', 'reduced_ring_required': False, 'coefficients_forced_zero': False, 'exists_large_zero_pair_iff_nilpotent_coefficient_ideal': True, 'prescribed_same_N_iff_nilpotent_ideal': False, 'marked_moment_tail_from_actual_phase_nilpotence': True, 'adjunction_scope': 'independent exponential variables over every characteristic-zero field', 'adjunction_window': 't+2', 'adjunction_terminal': 'N>=max(N0,t+2)', 'monomial_base_window': 't+1, terminal N>=t+1', 'independent_monomial_factors': True, 'same_cutoff_for_independent_monomial_factors': True, 'dependent_additions_covered': False, 'dependent_monomial_factor_identity': False, 'effective_N0': False, 'uniform_N0_over_supports': False, 'all_terminal_orders': False, 'fixture_N_asserted_above_unknown_N0': False, 'new_adjunction_or_CI_mechanism_claim': False, 'formal_verification': False, 'external_review': False, 'CI_fixtures': [[0, 1, 2], [0, 2, 3], [1, 2, 4], [2, 2, 5], [1, 3, 4], [3, 2, 3], [2, 4, 4], [4, 3, 3]]}
ROUTE=['retained407_fixed_support_eventual_coprimality', 'homogeneous_regular_sequence', 'ordinary_Hilbert_complete_intersection', 'universal_Q_algebra_specialization', 'existential_large_pair_reverse', 'classical_independent_adjunction', 'independent_monomial_factors']
def require(ok,msg):
    if not ok:raise AssertionError(msg)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path):return json.loads(path.read_text(encoding="utf-8"))
def same_typed(a,b):
    if type(a) is not type(b):return False
    if type(b) is dict:return set(a)==set(b) and all(same_typed(a[k],b[k]) for k in b)
    if type(b) is list:return len(a)==len(b) and all(same_typed(x,y) for x,y in zip(a,b))
    return a==b
def check_inputs(data):require(same_typed(data,SCOPE),"Exact algebraic/ring/independence scope and types required")
def bases(repo,support):return {"local":HERE,"support":support,"repo":repo}
def check_sources(doc,repo,support):
    require(type(doc) is dict and set(doc)=={"schema","formal_verification","external_review","logical_route","sources"},"Source schema")
    require(type(doc["schema"]) is int and doc["schema"]==1,"Schema1")
    require(doc["formal_verification"] is False and doc["external_review"] is False,"Assurance flags")
    require(doc["logical_route"]==ROUTE,"Exact logical route")
    entries=doc["sources"];require(type(entries) is list and len(entries)==len(ROLES),"All mandatory source roles")
    out={};seen=set();roots=bases(repo,support)
    for row in entries:
        require(type(row) is dict and set(row)=={"role","kind","path","sha256"},"Source row schema")
        role=row["role"];require(type(role) is str and role in ROLES and role not in out,"Unique role")
        require((row["kind"],row["path"])==ROLES[role],"Exact role path")
        key=(row["kind"],row["path"]);require(key not in seen,"Unique path")
        require(type(row["sha256"]) is str and re.fullmatch("[0-9a-f]{64}",row["sha256"]),"SHA256")
        root=roots[row["kind"]].resolve();p=(root/row["path"]).resolve()
        require(p.is_relative_to(root) and p.is_file() and sha(p)==row["sha256"],"Missing/stale source:"+role)
        out[role]=row;seen.add(key)
    forbidden={"GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md",
               "UNIFORM_BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md",
               "all_exponent_prefix_certificates.json",
               "verify_all_univariate_two_term_sfc.py",
               "verify_ordinary_lower_one_finite_bridge.py",
               "ordinary_lower_one_finite_bridge_inputs.json"}
    require(not any(Path(path).name in forbidden for _,path in seen),"Unneeded stronger or finite-prefix route excluded")
    return out
def links(repo,support,rows):
    expected={(r["kind"],r["path"]):r["sha256"] for r in rows.values()}
    specs=[
      (PREFIX+"ordinary_binomial_eventual_source_ledger.json","sources","mapped",77),
      (PREFIX+"simple_turn_shrinking_window_source_ledger.json","sources","mapped",5),
      (PREFIX+"above_critical_regular_source_ledger.json","sources","mapped",8),
      (PREFIX+"no_vertical_pair_source_ledger.json","sources","mapped",4),
      (PREFIX+"square_root_gap_eventual_source_ledger.json","sources","mapped",54),
      (PREFIX+"cusp_eventual_strong_factorial_source_ledger.json","sources","mapped",44),
      (PREFIX+"strict_subcritical_source_ledger.json","sources","mapped",6),
      (PREFIX+"cusp_two_thimble_source_ledger.json","sources","mapped",8),
      (PREFIX+"cusp_shrinking_window_source_ledger.json","sources","mapped",9),
      (PREFIX+"cusp_global_complement_source_ledger.json","sources","mapped",6),
      (PREFIX+"quartic_cusp_local_completion_source_ledger.json","sources","canonical",15),
      (PREFIX+"quartic_lee_yang_zero_and_robin_source_ledger.json","sources","canonical",5),
      (PREFIX+"cusp_stationary_belyi_source_ledger.json","files","stationary",6)]
    out=[]
    for path,key,mode,count in specs:
        d=load(repo/path)
        require(type(d.get("schema")) is int and d["schema"]==1,"Imported schema1")
        entries=d.get(key);require(type(entries) is list and len(entries)==count,"Imported exact cardinality")
        seen=set()
        for row in entries:
            require(type(row) is dict and type(row.get("path")) is str and type(row.get("sha256")) is str,"Imported row")
            if mode=="canonical":target=("repo",PREFIX+row["path"])
            elif mode=="stationary":target=("repo",PREFIX+row["path"]) if row["root"]=="local" else ("repo",row["path"])
            else:
                kind=row.get("kind",row.get("root"));require(kind in {"local","support","repo"},"Imported root")
                target=("repo",PREFIX+row["path"]) if kind in {"local","support"} else ("repo",row["path"])
            require(target not in seen,"Unique imported source");seen.add(target)
            require(target in expected and expected[target]==row["sha256"],"Imported source must match direct mandatory pin")
        out.append({"kind":"repo","path":path,"sources":count})
    return out
def equal_receipt(a,b):require(a==b,"Fresh child receipt differs")
def replay(repo,support,rows):
    role="eventual_main_checker";receipt="eventual_receipt";proof="eventual_main_proof"
    with tempfile.TemporaryDirectory(prefix="eventual-binomial-algebraic-") as folder:
        out=Path(folder)/"actual.json"
        cmd=[sys.executable,"-B",str(repo/ROLES[role][1]),"--repo-root",str(repo),"--output",str(out)]
        r=subprocess.run(cmd,cwd=folder,capture_output=True,text=True,timeout=90)
        require(r.returncode==0,"Retained407 fresh replay failed:"+r.stdout+";"+r.stderr)
        require(out.is_file(),"Fresh407 receipt missing")
        raw=out.read_bytes();equal_receipt(raw,(repo/ROLES[receipt][1]).read_bytes())
        d=json.loads(raw);require(d.get("status")=="PASS","Retained407 PASS")
        require(d.get("proof_sha256")==rows[proof]["sha256"],"Retained407 proof identity")
        require(d.get("checker_sha256")==rows[role]["sha256"],"Retained407 checker identity")
        return [{"checker_role":role,"receipt_sha256":hashlib.sha256(raw).hexdigest(),"byte_identical":True}]
def padd(a,b):
    out=dict(a)
    for m,c in b.items():out[m]=out.get(m,Q(0))+c
    return {m:c for m,c in out.items() if c}
def pscale(a,c):return {m:Q(c)*v for m,v in a.items() if c*v}
def pmul(a,b):
    out={}
    for m,c in a.items():
        for n,d in b.items():
            z=tuple(x+y for x,y in zip(m,n))
            out[z]=out.get(z,Q(0))+c*d
    return {m:c for m,c in out.items() if c}
def ppow(a,n,variables):
    out={(0,)*variables:Q(1)}
    while n:
        if n&1:out=pmul(out,a)
        a=pmul(a,a);n//=2
    return out
def direct_moment(phase,m,variables):
    return sum((c*prod_factorials(ex) for ex,c in ppow(phase,m,variables).items()),Q(0))
def prod_factorials(ex):
    out=1
    for k in ex:out*=factorial(k)
    return out
def mu(a,d,m):
    return [Q(comb(m,j)*factorial(a*m+d*j)) for j in range(m+1)]
def macaulay(a,d,N,r):
    rows=[]
    for degree in [N-1,N]:
        if r>=degree:
            f=mu(a,d,degree)
            for shift in range(r-degree+1):
                row=[Q(0)]*(r+1)
                for j,c in enumerate(f):row[shift+j]=c
                rows.append(row)
    return rows
def rr(rows,ncols,track=False):
    a=[list(r) for r in rows];nr=len(a)
    tr=[[Q(int(i==j)) for j in range(nr)] for i in range(nr)] if track else None
    pivots=[];row=0
    for col in range(ncols):
        pivot=next((j for j in range(row,nr) if a[j][col]),None)
        if pivot is None:continue
        a[row],a[pivot]=a[pivot],a[row]
        if track:tr[row],tr[pivot]=tr[pivot],tr[row]
        value=a[row][col];a[row]=[x/value for x in a[row]]
        if track:tr[row]=[x/value for x in tr[row]]
        for j in range(nr):
            if j!=row and a[j][col]:
                c=a[j][col];a[j]=[x-c*y for x,y in zip(a[j],a[row])]
                if track:tr[j]=[x-c*y for x,y in zip(tr[j],tr[row])]
        pivots.append(col);row+=1
        if row==nr:break
    return a,pivots,tr
def ci_controls(fixtures):
    hilbert_rows=0;memberships=0;survivors=0;ordinary=[]
    for a,d,N in fixtures:
        quotient=[]
        for r in range(2*N):
            matrix=macaulay(a,d,N,r);reduced,pivots,_=rr(matrix,r+1)
            actual=r+1-len(pivots)
            conv=sum(1 for i in range(N-1) for j in range(N) if i+j==r)
            exact=(r+1)-max(0,r-N+2)-max(0,r-N+1)+max(0,r-2*N+2)
            require(actual==conv==exact,"Actual graded quotient versus independent Hilbert formulas")
            quotient.append(actual);hilbert_rows+=1
        require(sum(quotient)==N*(N-1),"Actual CI length")
        require(quotient[2*N-3]==1 and not any(quotient[2*N-2:]),"Actual top degree and sharp ideal index")
        L=2*N-2;matrix=macaulay(a,d,N,L)
        reduced,pivots,tr=rr(matrix,L+1,True)
        require(pivots==list(range(L+1)),"Full ordinary Macaulay rank")
        for i in range(L+1):
            witness=[sum((tr[i][j]*matrix[j][k] for j in range(len(matrix))),Q(0)) for k in range(L+1)]
            require(witness==[Q(int(i==k)) for k in range(L+1)],"Literal ordinary degree-L monomial membership")
            memberships+=1
        top=2*N-3;matrix=macaulay(a,d,N,top)
        reduced,pivots,_=rr(matrix,top+1)
        free=next(j for j in range(top+1) if j not in pivots)
        functional=[Q(int(j==free)) for j in range(top+1)]
        for i,pivot in enumerate(pivots):functional[pivot]=-reduced[i][free]
        require(all(sum((x*y for x,y in zip(row,functional)),Q(0))==0 for row in matrix),"Literal top functional kills the whole ordinary ideal component")
        require(functional[free]==1 and comb(top,free)*functional[free]!=0,"A phase coefficient survives in the universal quotient")
        survivors+=1
        # The full rank also independently certifies fixture coprimality.
        ordinary.append({"a":a,"d":d,"N":N,"length":sum(quotient),
                         "top_degree":top,"nilpotence_index":L,"surviving_B_exponent":free})
    return {"actual_CI_fixtures":len(fixtures),"actual_Hilbert_degree_controls":hilbert_rows,
            "literal_ordinary_membership_witnesses":memberships,"literal_top_degree_functionals":survivors},ordinary
def epmul(a,b,cut):
    out={}
    for i,c in a.items():
        for j,d in b.items():
            if i+j<cut:out[i+j]=out.get(i+j,Q(0))+c*d
    return {j:c for j,c in out.items() if c}
def eppow(a,n,cut):
    out={0:Q(1)}
    while n:
        if n&1:out=epmul(out,a,cut)
        a=epmul(a,a,cut);n//=2
    return out
def epadd(a,b):
    out=dict(a)
    for j,c in b.items():out[j]=out.get(j,Q(0))+c
    return {j:c for j,c in out.items() if c}
def epmoment(a,d,m,A,B,cut):
    out={}
    for j in range(m+1):
        c=comb(m,j)*factorial(a*m+d*j)
        z=epmul(eppow(A,m-j,cut),eppow(B,j,cut),cut)
        out=epadd(out,{k:c*v for k,v in z.items()})
    return out
def nilpotent_controls():
    e={1:Q(1)};neg={1:-Q(1)}
    require(epmoment(0,1,1,neg,e,2)=={} and epmoment(0,1,2,neg,e,2)=={},"Actual N2 ring moment pair vanishes")
    require(e and neg and not epmul(e,e,2),"Coefficient classes survive while their square ideal is zero")
    count=0
    for N in range(2,13):
        cut=N+1
        require(eppow(e,N,cut) and not eppow(e,N+1,cut),"Explicit coefficient ideal nilpotence index")
        require(epmoment(2,3,N-1,e,{},cut),"Nilpotence does not force a prescribed large pair to vanish")
        require(all(not epmoment(2,3,m,e,{},cut) for m in range(cut,cut+3)),"All sufficiently high moments do vanish")
        count+=1
    return count
def independent_controls():
    recurrence=0;factorization=0;index=0
    cases=[(0,2,Q(2,3),-Q(3,5)),(1,1,Q(1),-Q(1)),(2,3,-Q(2),Q(3)),(0,1,Q(3),Q(0))]
    for a,d,A,B in cases:
        for t in range(4):
            nvars=1+t
            phase={}
            ex=(a,)+(0,)*t;phase[ex]=phase.get(ex,Q(0))+A
            ex=(a+d,)+(0,)*t;phase[ex]=phase.get(ex,Q(0))+B
            lam=[Q(1,2),-Q(2,3),Q(3,4)][:t]
            for j,value in enumerate(lam):
                ex=[0]*nvars;ex[1+j]=1;ex=tuple(ex);phase[ex]=phase.get(ex,Q(0))+value
            phase={ex:c for ex,c in phase.items() if c}
            f=[sum((Q(comb(m,j)*factorial(a*m+d*j),factorial(m))*A**(m-j)*B**j for j in range(m+1)),Q(0)) for m in range(8)]
            g=[direct_moment(phase,m,nvars)/factorial(m) for m in range(8)]
            den=[Q(1)]
            for value in lam:
                z=[Q(0)]*(len(den)+1)
                for j,c in enumerate(den):z[j]+=c;z[j+1]-=value*c
                den=z
            for m in range(8):
                require(f[m]==sum((den[r]*g[m-r] for r in range(min(t,m)+1)),Q(0)),"Independent multivariate moment versus adjunction recurrence")
                recurrence+=1
            for c in [(0,1),(1,2)]:
                multiplied={ex+c:v for ex,v in phase.items()}
                for m in range(1,6):
                    direct=direct_moment(multiplied,m,nvars+2)
                    predicted=factorial(c[0]*m)*factorial(c[1]*m)*direct_moment(phase,m,nvars)
                    require(direct==predicted,"Independent monomial factor uses actual disjoint coordinates")
                    factorization+=1
    for t in range(8):
        for N0 in [2,3,7,15]:
            N=max(N0,t+2);block=set(range(N-t-1,N+1))
            require(min(block)>=1 and len(block)==t+2,"Positive t+2 window")
            for m in [N-1,N]:
                require(all(m-r in block for r in range(t+1)),"Both base equations use only the prescribed block")
            M=t+1;short=set(range(M-t,M+1))
            require(min(short)==1 and all(M-r in short for r in range(t+1)),"Monomial-base t+1 positive window")
            index+=1
    return recurrence,factorization,index
def universal_controls():
    # Exact formal polynomial identities in parameters p,q and coefficient ring variable.
    one={(0,0):Q(1)};p={(1,0):Q(1)};q={(0,1):Q(1)}
    sub=lambda a,b:padd(a,pscale(b,-1))
    lhs=sub(sub(pmul(padd(p,q),sub(padd(p,q),one)),pmul(p,sub(p,one))),pmul(q,sub(q,one)))
    require(lhs==pscale(pmul(p,q),2),"Second numerator derivative at1 gives2pq")
    N={(1,):Q(1)};O={(0,):Q(1)}
    Pn=sub(N,O);Qn=N
    require(sub(padd(Pn,Qn),{(0,):Q(2)})==sub(pscale(N,2),{(0,):Q(3)}),"CI top degree2N-3")
    require(sub(padd(Pn,Qn),O)==sub(pscale(N,2),{(0,):Q(2)}),"First ideal nilpotence index2N-2")
    geometric=0
    for M in range(10):
        z={(1,1):Q(1)}
        series={(j,j):Q(1) for j in range(M+1)}
        result=pmul(sub(one,z),series)
        require(result=={(0,0):Q(1),(M+1,M+1):-Q(1)},"Universal finite geometric inverse identity")
        geometric+=1
    return 3,geometric
def mathematical_negatives():
    names=[]
    # Two homogeneous equations without coprimality do not have this CI length.
    # At N3 choose A²,A³: every vector has zero pure-B coefficient.
    matrix=[]
    for degree,form in [(2,[Q(1),Q(0),Q(0)]),(3,[Q(1),Q(0),Q(0),Q(0)])]:
        for shift in range(4-degree+1):
            row=[Q(0)]*5
            for j,c in enumerate(form):row[shift+j]=c
            matrix.append(row)
    reduced,pivots,_=rr(matrix,5)
    require(len(pivots)<5 and all(row[-1]==0 for row in matrix),"Common-factor CI counterexample")
    names.append("common-factor-equations-not-a-zero-dimensional-CI")
    # In characteristic5, unit coefficient phase X has both moments at5,6 zero.
    require(factorial(5)%5==factorial(6)%5==0 and pow(1,100,5)==1,"Q-algebra hypothesis is necessary")
    names.append("positive-characteristic-nilpotence-reverse-fails")
    # For dependent g=X+X, f2=1 but g2-g1=4-2=2 after m! normalization.
    require(Q(factorial(2),factorial(2))!=Q(4*factorial(2),factorial(2))-2,"Dependent adjunction recurrence genuinely fails")
    names.append("dependent-linear-adjunction-false")
    require(factorial(4)!=factorial(2)**2,"Dependent monomial factor identity genuinely fails")
    names.append("dependent-monomial-factor-false")
    # This is already checked by actual truncated-Q-algebra computations.
    names.append("nilpotent-does-not-imply-the-same-prescribed-zero-pair")
    names.append("two-ring-moments-do-not-force-zero-coefficients")
    return names
def algebra():
    fixtures=FIXTURES
    ci,details=ci_controls(fixtures)
    ring=nilpotent_controls()
    rec,fac,idx=independent_controls()
    uni,geo=universal_controls()
    bad=mathematical_negatives()
    ci.update({"nonreduced_fixed_N_countercontrols":ring,
       "independent_adjunction_moment_controls":rec,"independent_monomial_moment_controls":fac,
       "positive_window_index_controls":idx,"universal_Hilbert_parameter_identities":uni,
       "universal_geometric_inverse_controls":geo,"mathematical_negative_controls":len(bad)})
    return ci,details,bad


def negatives(ledger,repo,support):
    out=[]
    def reject(name,fn):
        try:fn()
        except (AssertionError,TypeError,KeyError,ValueError):out.append(name);return
        raise AssertionError("Negative accepted:"+name)
    for role in sorted(ROLES):
        d=copy.deepcopy(ledger);d["sources"]=[r for r in d["sources"] if r["role"]!=role]
        reject("omit:"+role,lambda d=d:check_sources(d,repo,support))
        d=copy.deepcopy(ledger);next(r for r in d["sources"] if r["role"]==role)["sha256"]="0"*64
        reject("badpin:"+role,lambda d=d:check_sources(d,repo,support))
    for key in ["radical_only","reduced_ring_required","coefficients_forced_zero",
                "prescribed_same_N_iff_nilpotent_ideal","dependent_additions_covered",
                "dependent_monomial_factor_identity","effective_N0","uniform_N0_over_supports",
                "all_terminal_orders","fixture_N_asserted_above_unknown_N0","new_adjunction_or_CI_mechanism_claim"]:
        d=copy.deepcopy(SCOPE);d[key]=True
        reject("false-scope:"+key,lambda d=d:check_inputs(d))
    for key in ["ordinary_ideal_containments","exists_large_zero_pair_iff_nilpotent_coefficient_ideal",
                "independent_monomial_factors","same_cutoff_for_independent_monomial_factors"]:
        d=copy.deepcopy(SCOPE);d[key]=False
        reject("missing-scope:"+key,lambda d=d:check_inputs(d))
    for key,val in [("CI_length","2*N"),("top_degree","2*N-2"),("sharp_ideal_nilpotence_index","2*N-3"),
                    ("ring_scope","all commutative rings"),("adjunction_window","t+1 for all two-term bases"),
                    ("adjunction_terminal","every N>=N0")]:
        d=copy.deepcopy(SCOPE);d[key]=val
        reject("wrong-quantity:"+key,lambda d=d:check_inputs(d))
    d=copy.deepcopy(SCOPE);d["support"]["a_min"]=False
    reject("typed-support",lambda:check_inputs(d))
    for key in ["schema","empty","duplicate"]:
        d=copy.deepcopy(ledger)
        if key=="schema":d["schema"]=True
        if key=="empty":d["sources"]=[]
        if key=="duplicate":d["sources"][-1]=copy.deepcopy(d["sources"][0])
        reject("ledger:"+key,lambda d=d:check_sources(d,repo,support))
    reject("forged-child",lambda:equal_receipt(b'{"status":"PASS"}',b'{"status":"PASS","actual":true}'))
    return out
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--repo-root",required=True,type=Path)
    ap.add_argument("--output",type=Path,default=HERE/RN)
    ns=ap.parse_args();repo=ns.repo_root.resolve()
    require(repo.is_dir(),"Explicit repository root")
    inp=load(HERE/IN);check_inputs(inp)
    ledger=load(HERE/LN);rows=check_sources(ledger,repo,HERE)
    linked=links(repo,HERE,rows);children=replay(repo,HERE,rows)
    counts,details,mathbad=algebra();bad=negatives(ledger,repo,HERE);check_sources(ledger,repo,HERE)
    counts.update({"fresh_retained407_replays":len(children),"cross_linked_source_ledgers":len(linked),
                   "source_scope_controls":len(bad)})
    payload={"status":"PASS","schema":1,"scope":inp,"formal_verification":False,"external_review":False,
      "proof_sha256":sha(HERE/PN),"checker_sha256":sha(HERE/CN),"input_sha256":sha(HERE/IN),
      "source_ledger_sha256":sha(HERE/LN),"mandatory_source_roles":len(rows),
      "sources":[rows[r] for r in sorted(rows)],"verified_child_ledgers":linked,"fresh_children":children,
      "counts":counts,"actual_fixture_invariants":details,"mathematical_negative_names":mathbad,
      "source_scope_negative_names":bad,
      "meaning":"Retained407 supplies the large-pair premise. The written regular-sequence/Hilbert, universal specialization and formal independence arguments prove the consequences. Exact ordinary witnesses and finite/source replay support the proof; they are not formal verification or a finite-prefix inference."}
    ns.output.parent.mkdir(parents=True,exist_ok=True)
    ns.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":"PASS","mandatory_sources":len(rows),"counts":counts,"output":str(ns.output)}))
if __name__=="__main__":main()

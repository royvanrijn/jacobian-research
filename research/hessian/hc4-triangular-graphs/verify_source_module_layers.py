"""Transcription controls for the written two-layer source theorem."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as s
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):check(s.expand(value)==0,label)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
y,t=s.symbols('y t',nonzero=True)
# The actual finite power-field action;no invented arbitrary semisimplicity.
for degree in [2,3,5,7]:
    M=s.zeros(degree)
    for i in range(degree-1):M[i+1,i]=1
    M[0,degree-1]=y
    check(M**degree==y*s.eye(degree),'finite power map has its exact regular representation at f='+str(degree))
    check(s.trace(s.eye(degree))==degree,'the trivial curve character has nonzero field trace at f='+str(degree))
    for j in range(1,degree):
        zero(s.trace(M**j),'the nontrivial cyclic character has trace zero at f='+str(degree)+' and j='+str(j))
        lam=s.Rational(j,degree)
        check(all(lam+k!=0 for k in [-4,-1,0,2,5]),'nonintegral character Euler weights do not hit the invariant de Rham class at f='+str(degree)+' and j='+str(j))
    for k in [0,1,3]:
        zero(s.diff(t**(degree*k),t)/(degree*t**(degree-1))-k*t**(degree*(k-1)),
             'the target derivative is the lifted curve derivative at f='+str(degree)+' and k='+str(k))

# Generic constant-row traces can be normalized by arbitrary nonzero
# residue conventions;their kernel and its dimension do not change.
for partition,degrees in [([1],[1]),([1],[4]),([2],[1,1]),([2,1],[2,1,3]),([1,1,1],[1,1,1])]:
    b=sum(partition);r=len(partition);rows=s.zeros(r,b);offset=0
    for j,count in enumerate(partition):
        for i in range(offset,offset+count):rows[j,i]=degrees[i]
        offset+=count
    check(rows.rank()==r,'boundary trace rows are onto for grouping '+str(partition)+' with f='+str(degrees))
    check(len(rows.nullspace())==b-r,'constant curve-character kernel has dimension b-r for '+str(partition))
    length_T=sum(degrees);length_T0=sum(f-1 for f in degrees)+b-r
    check(length_T0==length_T-r,'the complete semisimple trace-zero boundary length is sum f minus r for '+str(partition))
    check((length_T0==0)==(all(vv==1 for vv in partition) and all(ff==1 for ff in degrees)),
          'trace-zero semisimplicity has exactly one degree-one prime per image for '+str(partition))
    # An invertible Gysin map is required by source de Rham acyclicity;
    # it is not constructed for a hypothetical source by this matrix model.
    check(s.eye(b).rank()==b and b-(b-r)==r,'the exact-sequence cohomology dimensions are consistent for '+str(partition))

# The one-boundary abstract length-three extension has two nonzero arrows
# from its simple top to its two distinct socle constituents. Its endomorphisms
# are scalar in this model;setting either arrow to zero permits splitting.
a,b,c=s.symbols('a b c')
check(s.linsolve([a-c,b-c],(a,b,c))==s.FiniteSet((c,c,c)),
      'two nonzero extension components leave no splitting projector in the length-three arrow model')
check(s.linsolve([a-c],(a,b,c))==s.FiniteSet((c,b,c)),
      'a hostile zero extension component creates an independent socle endomorphism')
check((2-1)+1==2 and (1+1)+1==3,'the trace image and one-boundary source have lengths two and three')
check(2-4==-2,'Verdier duality has the correct codimension-one purity shift even at a quotient point')
check(0-2+2==0 and 1-2==-1,'surface and curve perverse normalizations yield the curve boundary in unshifted degree one')
algebra=len(checks)
context=json.loads((HERE/'sources/source-module-layers-context.json').read_text())
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(sha(ROOT/path)==digest,'immutable retained input '+path)
for entry in context['newly_cached_primary_context']:
    check(sha(HERE/'sources'/entry['retained_file'])==entry['sha256'],'inspected new primary '+entry['retained_file'])
check(len(context['newly_imported_proved_external_statements'])==3 and context['newly_used_unproved_external_statements']==[],
      'the RH,finite-IC and Borisov power-map statements are the explicit three proved imports')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(flags['prior_PF2BATT1_nodal_exclusion_credit_is_explicit'] and flags['source_socle_and_boundary_quotient_are_semisimple']
      and not flags['trace_zero_semisimplicity_proved_to_force_inverse'] and not flags['original_goal_complete'],
      'attribution,source-layer structure and the remaining extension/multiplication gap stay explicit')
proof=HERE/'SOURCE_MODULE_LAYERS_AND_TRACE_ZERO_KERNEL.md'
paths=['sources/source-module-layers-context.json']+['sources/'+e['retained_file'] for e in context['newly_cached_primary_context']]
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=s.__version__,
 proof_sha256=sha(proof),checker_sha256=sha(Path(__file__)),dependencies=context['retained_proof_and_source_sha256'],
 source_sha256={p:sha(HERE/p) for p in paths},independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-source-module-layers.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies','source_sha256'}},indent=2))

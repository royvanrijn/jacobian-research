"""Transcription controls for the all-degree canonical-energy criterion."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as s
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):check(s.expand(value)==0,label)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
z1,z2,z3,w1,w2,a,b,c,DD,shift=s.symbols('z1 z2 z3 w1 w2 a b c DD shift')
# Arbitrary source symbols,without imposing any small-degree Keller model.
for power in [1,2,3,5,7]:
    expanded=s.Poly((a*z1+b*z2+c*z3-DD)**power,a,b,c)
    for monomial,expected in [((power,0,0),z1**power),((0,power,0),z2**power),((0,0,power),z3**power)]:
        zero(expanded.coeff_monomial(monomial)-expected,'pure target coefficient recovers a full source power at k='+str(power)+' and '+str(monomial))
    zero((a*z1+b*z2-DD+shift)**power-sum(s.binomial(power,j)*shift**(power-j)*(a*z1+b*z2-DD)**j for j in range(power+1)),
         'adding a target polynomial shifts moments triangularly at k='+str(power))
K1=s.Function('K1')(z1,z2);K2=s.Function('K2')(z1,z2);D=s.Function('D')(z1,z2)
h=D+w1*K1+w2*K2
xi=[s.diff(h,z1),s.diff(h,z2)]
G=z1*xi[0]+z2*xi[1]+w1*s.diff(h,w1)+w2*s.diff(h,w2)-h
zero(G-(z1*xi[0]+z2*xi[1]-D),'the whole cotangent energy is linear in the dual gradient target')
zero(s.det(s.Matrix([s.diff(h,w1),s.diff(h,w2)]).jacobian([z1,z2]))-s.det(s.Matrix([K1,K2]).jacobian([z1,z2])),
     'the cotangent base outputs retain the actual Keller base Jacobian')
# Field trace base extension is represented by the same matrices with new
# independent target variables;trace is coefficientwise by linearity.
t=s.symbols('t');M=s.Matrix([[0,t],[1,0]]);I=s.eye(2)
for power in [1,2,4,6]:
    V=a*M+b*I
    zero(s.Poly(s.trace(V**power),a,b).coeff_monomial(a**power)-s.trace(M**power),
         'actual companion trace coefficient commutes with independent target extension at k='+str(power))
check(s.Poly(a*z1+b*z2-DD,a,b).coeff_monomial(a)==z1 and s.Poly(a*z1+b*z2-DD,a,b).coeff_monomial(b)==z2,
      'independent target residues prevent simultaneous cancellation of all negative-valuation source coefficients')
algebra=len(checks)
context=json.loads((HERE/'sources/canonical-energy-trace-context.json').read_text())
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(sha(ROOT/path)==digest,'retained immutable input '+path)
check(context['newly_imported_proved_external_statements']==context['newly_used_unproved_external_statements']==[],
      'the HC4 application explicitly uses retained Tang and the specified cotangent proof imports no new theorem')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(flags['universal_HC4_energy_criterion_uses_retained_Tang'] and not flags['canonical_energy_integrality_proved_automatic'] and not flags['original_goal_complete'],
      'the single-observable criterion retains its classification dependency and unsupplied integrality gate')
proof=HERE/'CANONICAL_ENERGY_INTEGRALITY_AND_POWER_TRACES.md'
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=s.__version__,
 proof_sha256=sha(proof),checker_sha256=sha(Path(__file__)),dependencies=context['retained_proof_and_source_sha256'],
 source_sha256={'sources/canonical-energy-trace-context.json':sha(HERE/'sources/canonical-energy-trace-context.json')},
 independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-canonical-energy-traces.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies','source_sha256'}},indent=2))

"""Generate and directly verify actual-power Bernstein data for finite QQ input.

Singular constructs graph identities. A separate sparse rational Weyl engine
verifies their witnesses against the explicit Fourier/graph generators.
No bounded moment prefix, actual-image kernel, or supplied guessed p is used.
"""
from fractions import Fraction
from itertools import product
from pathlib import Path
import argparse
import hashlib
import json
import math
import re
import shutil
import subprocess
import sympy as S
import bernstein_slope_producer as selector


class Weyl:
    def __init__(self,n):self.n=n;self.zero=(0,)*(2*n)
    def scalar(self,c):return {self.zero:Fraction(c)} if c else {}
    def var(self,i):
        e=list(self.zero);e[i]=1;return {tuple(e):Fraction(1)}
    def add(self,*ps):
        out={}
        for p in ps:
            for e,c in p.items():out[e]=out.get(e,Fraction(0))+c
        return {e:c for e,c in out.items() if c}
    def scale(self,p,c):return {e:a*c for e,a in p.items() if a*c}
    def sub(self,p,q):return self.add(p,self.scale(q,-1))
    def mul(self,p,q):
        n=self.n;out={}
        for e,c in p.items():
            for f,d in q.items():
                for k in product(*(range(min(e[n+i],f[i])+1) for i in range(n))):
                    co=c*d
                    for i,ki in enumerate(k):co*=math.comb(e[n+i],ki)*math.prod(range(f[i]-ki+1,f[i]+1))
                    mo=tuple(e[i]+f[i]-k[i] for i in range(n))+tuple(e[n+i]+f[n+i]-k[i] for i in range(n))
                    out[mo]=out.get(mo,Fraction(0))+co
        return {e:c for e,c in out.items() if c}
    def power(self,p,j):
        out=self.scalar(1)
        while j:
            if j&1:out=self.mul(out,p)
            j//=2
            if j:p=self.mul(p,p)
        return out
    def substitute(self,terms,images):
        out={}
        for exponents,c in terms:
            if len(exponents)!=len(images) or any(type(j)!=int or j<0 for j in exponents):raise ValueError('Invalid PBW exponent vector')
            v=self.scalar(c)
            for image,j in zip(images,exponents):
                if j:v=self.mul(v,self.power(image,j))
            out=self.add(out,v)
        return out
    def parse(self,terms):
        out={}
        for e,c in terms:
            if len(e)!=2*self.n or any(type(j)!=int or j<0 for j in e):raise ValueError('Wrong Weyl exponent vector')
            e=tuple(e);out[e]=out.get(e,Fraction(0))+Fraction(c)
        return {e:c for e,c in out.items() if c}
    def serialize(self,p):return [[list(e),str(c)] for e,c in sorted(p.items())]


def graph_data(phase_terms,symbol_terms,n):
    """Canonical variable order:(xi,eta,t_lambda,t_eta;their derivatives)."""
    N=n+3;A=Weyl(N);v=[A.var(i) for i in range(2*N)]
    xi=v[:n];eta=v[n];tl,te=v[n+1:n+3];dl,de=v[N+n+1:N+n+3]
    f=A.substitute(symbol_terms,xi)
    shifted=[]
    for i in range(n):
        diff=[]
        for e,c in symbol_terms:
            if e[i]:
                a=list(e);a[i]-=1;diff.append([a,str(Fraction(c)*e[i])])
        fi=A.substitute(diff,xi)
        shifted.append(A.add(v[N+i],A.mul(fi,dl)))
    minus=[A.scale(z,-1) for z in shifted]
    generators=[]
    for i in range(n):
        diff=[]
        for e,c in phase_terms:
            if e[i]:
                a=list(e);a[i]-=1;diff.append([a,str(Fraction(c)*e[i])])
        generators.append(A.sub(xi[i],A.mul(eta,A.substitute(diff,minus))))
    generators.append(A.sub(A.add(v[N+n],de),A.substitute(phase_terms,minus)))
    generators.extend([A.sub(tl,f),A.sub(te,eta)])
    ES=A.scale(A.mul(dl,tl),-1);EZ=A.scale(A.mul(de,te),-1)
    base_images=xi+[eta]+v[N:N+n]+[v[N+n]]
    lifted_images=xi+[eta]+v[N:N+n]+[A.add(v[N+n],de)]
    return A,generators,base_images,lifted_images,ES,EZ,f,eta


def rational_terms(terms,length):
    out={}
    for e,c in terms:
        if len(e)!=length or any(type(a)!=int or a<0 for a in e):raise ValueError('Invalid input exponent vector')
        e=tuple(e);out[e]=out.get(e,Fraction(0))+Fraction(c)
    return [[list(e),str(c)] for e,c in sorted(out.items()) if c]


def expression(terms,names):
    parts=[]
    for e,c in terms:
        co=Fraction(c);factors=[f'({co.numerator}/{co.denominator})']
        factors.extend(f'({name})^{a}' for name,a in zip(names,e) if a)
        parts.append('*'.join(factors))
    return '+'.join(parts) or '0'


EMITTER='''
proc emit(poly v,string label)
{
  print("BEGIN|"+label);
  while(v!=0)
  {
    print("TERM|"+label+"|"+string(leadexp(v))+"|"+string(leadcoef(v)));
    v=v-lead(v);
  }
  print("END|"+label);
}
'''


def parse_output(stdout):
    records={};checks={}
    u=S.Symbol('u')
    for line in stdout.splitlines():
        if line.startswith('BEGIN|'):records[line.split('|',1)[1]]=[]
        elif line.startswith('TERM|'):
            _,label,e,c=line.split('|',3)
            records[label].append([[int(a.strip()) for a in e.split(',')],str(S.sympify(c.replace('^','**'),locals={'u':u}))])
        elif line.startswith('CHECK|'):
            _,label,value=line.split('|',2);checks[label]=value.strip()
    if any(line.lstrip().startswith('?') for line in stdout.splitlines()):raise ValueError('Singular reported an error')
    return records,checks


def run_stage(script,workdir,name,timeout):
    path=workdir/(name+'.sing');outpath=workdir/(name+'.stdout')
    path.write_text(script)
    result=subprocess.run(['Singular','-q','--random=1',str(path)],capture_output=True,text=True,timeout=timeout)
    outpath.write_text(result.stdout+result.stderr)
    if result.returncode!=0:raise ValueError('Singular stage failed:'+name)
    if 'STAGE_COMPLETE|'+name not in result.stdout:raise ValueError('Singular stage incomplete:'+name)
    records,checks=parse_output(result.stdout)
    if any(value!='0' for value in checks.values()):raise ValueError('Nonzero certificate normal form:'+name)
    return records,checks,dict(stage=name,script_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),stdout_sha256=hashlib.sha256(outpath.read_bytes()).hexdigest())


def first_script(phase,symbol,n):
    m=n+1;names=[f'z{i+1}' for i in range(n)]+['h']+[f'D{i+1}' for i in range(n)]+['Dh']
    x=S.symbols('x0:'+str(n))
    P=sum(S.Rational(c)*S.prod(v**a for v,a in zip(x,e)) for e,c in phase)
    def from_sympy(v):return [[list(e),str(c)] for e,c in S.Poly(v,*x,domain=S.QQ).terms() if c]
    dp=[f'(-D{i+1})' for i in range(n)]
    I=[f'z{i+1}-h*({expression(from_sympy(S.diff(P,x[i])),dp)})' for i in range(n)]
    I.append('Dh-('+expression(phase,dp)+')')
    f=expression(symbol,names[:n]);allnames=','.join(names)
    code='LIB "dmodapp.lib";\n'+EMITTER
    code+=f'ring R=0,({allnames}),dp; def W=Weyl();setring W;\n'
    code+='ideal InputIdeal='+','.join(I)+';\n'
    code+=f'def A=SDLoc(InputIdeal,{f});setring A;\n'
    code+=f'ideal F=LD,{f};ideal C=std(F);ideal Bs=nselect(C,1..{2*m});poly B=Bs[1];\n'
    code+='matrix L=lift(F,ideal(B));poly Q=L[ncols(F),1];\n'
    code+=f'print("CHECK|ordinary|"+string(reduce(B-Q*({f}),std(LD))));\n'
    code+='emit(B,"ordinary_B");emit(Q,"ordinary_Q");\n'
    code+=f'ring V0=(0,u),({allnames}),dp;def V=Weyl();setring V;\n'
    code+='map Phi=A,'+','.join(names+['u'])+';ideal GenericIdeal=Phi(LD);\n'
    code+='def G=SDLoc(GenericIdeal,h);setring G;ideal F=LD,h;ideal C=std(F);\n'
    code+=f'ideal Bs=nselect(C,1..{2*m});poly B=Bs[1];matrix L=lift(F,ideal(B));poly Q=L[ncols(F),1];\n'
    code+='print("CHECK|eta_generic|"+string(reduce(B-Q*h,std(LD))));\n'
    code+='emit(B,"generic_eta_B");emit(Q,"generic_eta_Q");print("STAGE_COMPLETE|first");quit;\n'
    return code


def clear_generic(records):
    u=S.Symbol('u');den=S.Integer(1)
    for key in ['generic_eta_B','generic_eta_Q']:
        for e,c in records[key]:den=S.lcm(den,S.fraction(S.cancel(S.sympify(c,locals={'u':u})))[1])
    out={}
    for key in ['generic_eta_B','generic_eta_Q']:
        result=[]
        for e,c in records[key]:
            poly=S.Poly(S.cancel(den*S.sympify(c,locals={'u':u})),u,domain=S.QQ)
            for (a,),co in poly.terms():
                if co:result.append([e[:-1]+[a,e[-1]],str(co)])
        out[key]=result
    return out,str(den)


def second_script(phase,symbol,n,ordinary,cleared):
    m=n+1;N=n+3;A,gens,base,lifted,ES,EZ,f,eta=graph_data(phase,symbol,n)
    names=[f'z{i+1}' for i in range(n)]+['h','tl','te']+[f'D{i+1}' for i in range(n)]+['Dh','Dl','De']
    ob=A.substitute(ordinary['ordinary_B'],lifted+[ES])
    oq=A.substitute(ordinary['ordinary_Q'],lifted+[ES])
    eb=A.substitute(cleared['generic_eta_B'],base+[ES,EZ])
    eq=A.substitute(cleared['generic_eta_Q'],base+[ES,EZ])
    od=A.sub(ob,A.mul(oq,f));ed=A.sub(eb,A.mul(eq,eta))
    ex=lambda q:expression(A.serialize(q),names)
    code='LIB "bfun.lib";\n'+EMITTER
    code+='ring R=0,('+','.join(names)+'),dp;def W=Weyl();setring W;\n'
    code+='ideal J='+','.join('('+ex(g)+')' for g in gens)+';ideal G=std(J);\n'
    code+='poly ES='+ex(ES)+';poly OrdinaryDifference='+ex(od)+';poly EtaDifference='+ex(ed)+';\n'
    code+='print("CHECK|ordinary_graph|"+string(reduce(OrdinaryDifference,G)));\n'
    code+='ideal Rems;poly Remainder=reduce(EtaDifference,G);vector v;int k=0;poly H=0;\n'
    code+='while(1){Rems[k+1]=Remainder;v=linSyzSolve(Rems);if(v!=0){break;}k++;Remainder=reduce(ES*Remainder,G);}\n'
    code+='int j;for(j=1;j<=k+1;j++){H=H+v[j]*ES^(j-1);print("HCOEFF|"+string(j-1)+"|"+string(v[j]));}\n'
    code+='poly EtaFinal=H*EtaDifference;print("CHECK|eta_graph|"+string(reduce(EtaFinal,G)));\n'
    code+='matrix Witness=lift(J,ideal(OrdinaryDifference,EtaFinal));\n'
    code+=f'for(j=1;j<={len(gens)};j++){{emit(Witness[j,1],"ordinary_witness_"+string(j));emit(Witness[j,2],"eta_witness_"+string(j));}}\n'
    code+='print("STAGE_COMPLETE|second");quit;\n'
    return code


def verify_certificate(certificate):
    """Exact Weyl identities against input-derived generators;no Singular call."""
    phase=certificate['input']['phase_terms'];symbol=certificate['input']['symbol_terms'];n=certificate['input']['dimension']
    if type(n)!=int or n<1:raise ValueError('Invalid input dimension')
    if phase!=rational_terms(phase,n) or symbol!=rational_terms(symbol,n):raise ValueError('Noncanonical finite input')
    if not phase or max(sum(e) for e,c in phase)<1 or not symbol or any(sum(e)==0 for e,c in symbol):raise ValueError('Input outside positive-degree/zero-constant contract')
    A,gens,base,lifted,ES,EZ,f,eta=graph_data(phase,symbol,n)
    if certificate.get('graph_generators')!=[A.serialize(g) for g in gens]:raise ValueError('Recorded graph generators differ from the original input')
    one=certificate['ordinary'];two=certificate['eta']
    ob=A.substitute(one['B_terms'],lifted+[ES]);oq=A.substitute(one['Q_terms'],lifted+[ES])
    eb=A.substitute(two['cleared_B_terms'],base+[ES,EZ]);eq=A.substitute(two['cleared_Q_terms'],base+[ES,EZ])
    H=A.substitute(two['torsion_polynomial_terms'],[ES])
    od=A.sub(ob,A.mul(oq,f));ed=A.mul(H,A.sub(eb,A.mul(eq,eta)))
    for difference,witness in [(od,one['graph_witness']),(ed,two['graph_witness'])]:
        if len(witness)!=len(gens):raise ValueError('Witness length differs from the explicit graph presentation')
        total={}
        for w,g in zip(witness,gens):total=A.add(total,A.mul(A.parse(w),g))
        if total!=difference:raise ValueError('Exact Weyl certificate identity failed')
    m=n+1;s,t=S.symbols('s t')
    B1=sum(S.Rational(c)*s**e[-1] for e,c in one['B_terms'])
    if not B1 or any(any(e[:-1]) for e,c in one['B_terms']):raise ValueError('Invalid scalar ordinary Bernstein polynomial')
    J=max((e[2*m-1] for e,c in one['Q_terms']),default=0)
    roots=selector.integer_roots(S.Poly(B1,s,domain=S.QQ));k0=max([0]+roots)+1
    h=sum(S.Rational(c)*s**e[0] for e,c in two['torsion_polynomial_terms'])
    b2=sum(S.Rational(c)*s**e[-2]*t**e[-1] for e,c in two['cleared_B_terms'])
    if not h or not b2 or any(any(e[:-2]) for e,c in two['cleared_B_terms']):raise ValueError('Invalid two-parameter scalar polynomial')
    actual=S.Poly(S.expand((h*b2).subs(t,-selector.Y-1).subs(s,selector.X)),selector.X,selector.Y,domain=S.QQ)
    expected=dict(polynomial_terms=selector.terms(actual),validity_start=k0,step_bound=J)
    if expected!=certificate['certificate_summary']:raise ValueError('Incorrect specialization summary')
    starts=certificate.get('specialization_thresholds',{})
    if starts!=dict(ordinary_graph_relation_start=0,order_constraint_start=0,nonvanishing_and_steps_start=k0,ordinary_scalar_nonnegative_integer_roots=roots):raise ValueError('Incorrect uniform specialization thresholds')
    return True


def generate(input_data,workdir,timeout=30):
    n=input_data['dimension']
    if type(n)!=int or n<1:raise ValueError('Positive dimension required')
    phase=rational_terms(input_data['phase_terms'],n);symbol=rational_terms(input_data['symbol_terms'],n)
    if not phase or max(sum(e) for e,c in phase)<1:raise ValueError('This generator requires a nonzero positive-degree phase')
    if not symbol or any(sum(e)==0 for e,c in symbol):raise ValueError('Nonzero zero-constant polynomial symbol required')
    workdir=Path(workdir);workdir.mkdir(parents=True,exist_ok=True)
    one,c1,log1=run_stage(first_script(phase,symbol,n),workdir,'first',timeout)
    for name in ['ordinary_B','ordinary_Q','generic_eta_B','generic_eta_Q']:
        if name not in one:raise ValueError('Missing generated object:'+name)
    ordinary={key:[[e,str(S.Rational(c))] for e,c in one[key]] for key in ['ordinary_B','ordinary_Q']}
    cleared,den=clear_generic(one)
    code=second_script(phase,symbol,n,ordinary,cleared)
    two,c2,log2=run_stage(code,workdir,'second',timeout)
    text=(workdir/'second.stdout').read_text();hterms=[]
    for line in text.splitlines():
        if line.startswith('HCOEFF|'):
            _,i,c=line.split('|',2)
            if S.Rational(c):hterms.append([[int(i)],str(S.Rational(c))])
    s,t=S.symbols('s t');m=n+1
    B1=sum(S.Rational(c)*s**e[-1] for e,c in ordinary['ordinary_B'])
    roots=selector.integer_roots(S.Poly(B1,s,domain=S.QQ));k0=max([0]+roots)+1
    J=max((e[2*m-1] for e,c in ordinary['ordinary_Q']),default=0)
    h=sum(S.Rational(c)*s**e[0] for e,c in hterms)
    b2=sum(S.Rational(c)*s**e[-2]*t**e[-1] for e,c in cleared['generic_eta_B'])
    p=S.Poly(S.expand((h*b2).subs(t,-selector.Y-1).subs(s,selector.X)),selector.X,selector.Y,domain=S.QQ)
    algebra,generators,*_=graph_data(phase,symbol,n)
    certificate=dict(schema_version=1,completed=True,input=dict(dimension=n,phase_terms=phase,symbol_terms=symbol),
        ordinary=dict(B_terms=ordinary['ordinary_B'],Q_terms=ordinary['ordinary_Q'],
            graph_witness=[two[f'ordinary_witness_{i+1}'] for i in range(n+3)]),
        eta=dict(cleared_B_terms=cleared['generic_eta_B'],cleared_Q_terms=cleared['generic_eta_Q'],
            generic_coefficient_denominator_cleared=den,torsion_polynomial_terms=hterms,
            graph_witness=[two[f'eta_witness_{i+1}'] for i in range(n+3)]),
        certificate_summary=dict(polynomial_terms=selector.terms(p),validity_start=k0,step_bound=J),
        graph_generators=[algebra.serialize(g) for g in generators],
        specialization_thresholds=dict(ordinary_graph_relation_start=0,order_constraint_start=0,nonvanishing_and_steps_start=k0,ordinary_scalar_nonnegative_integer_roots=roots),
        specialization='Unlocalized graph evaluation at EVERY nonnegative first exponent;second exponent uses eta localization only.',
        stages=[log1,log2],normal_form_checks=dict(c1,**c2),
        exact_rational_Weyl_witness_verification=False,independent_replay=False,formal_verification=False,external_review=False)
    verify_certificate(certificate);certificate['exact_rational_Weyl_witness_verification']=True
    return certificate


def decide(input_data,workdir,timeout=30,*,all_positive_pure=False):
    n=input_data['dimension']
    if type(n)!=int or n<1:raise ValueError('Positive input dimension required')
    phase=rational_terms(input_data['phase_terms'],n);symbol=rational_terms(input_data['symbol_terms'],n)
    common=dict(completed=True,input=dict(dimension=n,phase_terms=phase,symbol_terms=symbol),
        certificate_source='Elementary exact finite-input case',independent_replay=False,formal_verification=False,external_review=False)
    if not phase or not symbol:
        return dict(common,decision=dict(alternative='zero phase' if not phase else 'zero operator',
            all_marker_eventual_safety=True,affine_cutoff=dict(C=0,E=0),pure_linear_cutoff=0))
    if any(sum(e)==0 for e,c in symbol):
        if all_positive_pure:raise ValueError('Nonzero phase and invertible operator contradict supplied purity')
        return dict(common,decision=dict(alternative='nonzero symbol constant makes the operator invertible',
            selected_rational_slope='0',all_marker_eventual_safety=False,bad_marker='1;every positive pure output is nonzero'))
    if max(sum(e) for e,c in phase)==0:
        return dict(common,decision=dict(alternative='nonzero constant phase and positive operator order',
            all_marker_eventual_safety=True,affine_cutoff=dict(C=1,E=0),pure_linear_cutoff=1))
    certificate=generate(input_data,workdir,timeout)
    plan=selector.make_plan(certificate['certificate_summary'])
    phase=certificate['input']['phase_terms'];symbol=certificate['input']['symbol_terms']
    D=max(sum(e) for e,c in phase)
    jet=max(D*plan['symbol_jet_order_required_per_phase_degree'],max(sum(e) for e,c in symbol))
    evaluations=selector.evaluate_rational(plan,phase,symbol,jet)
    decision=selector.resolve(plan,evaluations['evaluations'],all_positive_pure=all_positive_pure)
    return dict(completed=True,certificate_source='Generated from original finite rational input and verified by explicit exact Weyl identities',
        certificate=certificate,plan=plan,original_power_evaluations=evaluations,decision=decision,
        safety_decision_scope='qualitative eventuality for every fixed original polynomial marker',
        independent_replay=False,formal_verification=False,external_review=False)


def runtime_manifest():
    binary=Path(shutil.which('Singular'))
    version=subprocess.run([str(binary),'--version'],capture_output=True,text=True,check=True).stdout.splitlines()[0]
    libs=Path('/usr/share/singular/LIB')
    return dict(Singular_version=version,Singular_binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
        fixed_random_seed=1,primary_library_sha256={n:hashlib.sha256((libs/n).read_bytes()).hexdigest() for n in ['dmodapp.lib','bfun.lib','nctools.lib']},
        note='Primary graph/linear-syzygy library files pinned;this is not a complete transitive installation snapshot. Returned graph identities are separately checked without Singular.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--timeout',type=float,default=30)
    modes=parser.add_mutually_exclusive_group();modes.add_argument('--decide',action='store_true');modes.add_argument('--verify',action='store_true')
    parser.add_argument('--pure',action='store_true',help='Supply the all-positive purity premise;not inferred from finite tests')
    args=parser.parse_args();work=args.output.parent/(args.output.stem+'_work')
    data=json.loads(args.input.read_text())
    try:
        if args.verify:out=dict(completed=True,exact_rational_Weyl_witness_verification=verify_certificate(data),formal_verification=False)
        elif args.decide:out=decide(data,work,args.timeout or None,all_positive_pure=args.pure)
        else:out=generate(data,work,args.timeout or None)
    except subprocess.TimeoutExpired as e:
        out=dict(completed=False,classification='UNKNOWN',reason='Singular stage timeout',stage=str(e.cmd[-1]))
    args.output.write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()

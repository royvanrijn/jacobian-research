"""Cold, bounded MW16 transfer adapter; no BNF or historical relation inputs.

Reuses the successful small-representative Matrix and exact norm sieve. The
only arithmetic input is an equation and its generic subgroup. All field,
local, ideal and dependency coordinates are rebuilt for this input.
"""
import argparse
from math import gcd, prod
from pathlib import Path
import time

from transfer_common import *


def coefficients(a):
    from sage.all import pari
    return [str(pari.lift(a).polcoef(i)) for i in range(3)]


def setup(data):
    from sage.all import QQ, PolynomialRing, pari
    R = PolynomialRing(QQ, 'z')
    f = R(data['cubic_ascending'])
    nf = pari.nfinit([pari(f), data['S']])
    assert not len(pari.nfcertify(nf))
    assert str(nf.disc()) == data['field_discriminant']
    return R, f, nf


def local_generators(nf, P):
    """Complete squareclass generators, including residue degree >1 at two."""
    from sage.all import pari, GF, matrix
    p, e, f = int(P[0]), int(P[2]), int(P[3])
    pi = pari.nfbasistoalg(nf, pari.idealappr(nf, P))
    assert int(pari.idealval(nf, pi, P)) == 1
    if p == 2:
        bid = pari.idealstar(nf, pari.idealpow(nf, P, 2*e+1), 2)
        units = [pari.nfbasistoalg(nf, u) for u in bid.bid_get_gen()]
        expected = e*f+2
    else:
        mod = pari.nfmodprinit(nf, P)
        unit = None
        for b in nf.nf_get_zk():
            for k in range(64):
                u = b+k
                if u and pari.idealval(nf, u, P) == 0 and not pari.issquare(pari.nfmodpr(nf, u, mod)):
                    unit = u
                    break
            if unit is not None: break
        if unit is None: raise ArithmeticError('bounded local unit generation unresolved')
        units, expected = [unit], 2
    gens = [pi]+units
    assert all(pari.idealval(nf, u, P) == 0 for u in units)
    pairing = [[int(pari.nfhilbert(nf, a, b, P) == -1) for b in gens] for a in gens]
    assert matrix(GF(2), pairing).rank() == expected
    return gens, pairing


def prepare(folder):
    from sage.all import QQ, ZZ, EllipticCurve, PolynomialRing, pari
    from research_runtime.sage_arithmetic import SageArithmetic
    from research_runtime.memory_store import MemoryFactStore
    from prepare_small_conductor_norm_form import reduce_form, discriminant, hessian
    import sage.version
    packet = read(folder/'input.json')
    model, points = seed(packet)
    progress(folder, 'generic_rank_certificate')
    admitted, proof = admission(model, points)
    assert proof['rank_lower_bound'] == 16
    write(folder/'seed.json', dict(curve=list(map(str, model)), points=points, proof=proof))
    progress(folder, 'cold_field_preparation')
    arithmetic = SageArithmetic(MemoryFactStore())
    context = arithmetic.prepare(tuple(model), discover=True)
    E = EllipticCurve(QQ, model)
    minimal = EllipticCurve(QQ, context.minimal_model.coefficients)
    transport = E.isomorphism_to(minimal)
    aa = list(minimal.a_invariants())
    a1, a2, a3, a4, a6 = aa
    R = PolynomialRing(QQ, 'z')
    f = R([16*a3*a3+64*a6, 8*a1*a3+16*a4, a1*a1+4*a2, 1])
    S = sorted({2} | {int(p) for p, _ in ZZ(f.discriminant()).factor(proof=True)})
    assert f.is_irreducible()
    pari.addprimes(S)
    nf = pari.nfinit([pari(f), S])
    assert not len(pari.nfcertify(nf))
    gammas = []
    for point in points:
        x, y = transport(E(point)).xy()
        z, W = 4*x, 8*y+4*a1*x+4*a3
        assert W*W == f(z)
        denominator = ZZ(z.denominator()).sqrt()
        assert denominator**2 == z.denominator()
        beta = denominator**2*(z-R.gen())
        n = QQ(pari.nfeltnorm(nf, pari.Mod(pari(beta), pari(f))))
        assert n == (denominator**3*W)**2
        assert ZZ(beta[0]).gcd(denominator) == 1
        gammas.append(dict(beta_ascending=list(map(str, beta.list())), norm=str(n)))
    one, w0, t0 = nf.nf_get_zk()
    assert one == 1
    pair = pari.nfalgtobasis(nf, w0*t0)
    w, t = w0-pair[2], t0-pair[1]
    bs = pari.Mat([pari.nfalgtobasis(nf, v) for v in [one, w, t]])
    def coords(v): return list(bs**-1*pari.nfalgtobasis(nf, v))
    ww, wt, tt = coords(w*w), coords(w*t), coords(t*t)
    a, b, c, d = -ww[2], ww[1], -tt[2], tt[1]
    assert ww == [-a*c,b,-a] and wt == [-a*d,0,0] and tt == [-b*d,d,-c]
    form = list(map(int, [a,b,c,d]))
    assert discriminant(form) == nf.disc()
    reduced, M, steps = reduce_form(form)
    # Degree-three homogeneous identity is determined by four distinct slopes.
    for m, n in [(1,0),(0,1),(1,1),(-1,1)]:
        alpha = a*(M[0]*m+M[1]*n)+(M[2]*m+M[3]*n)*w
        assert pari.nfeltnorm(nf, alpha) == a*a*sum(v*m**(3-i)*n**i for i,v in enumerate(reduced))
    data = dict(cubic_ascending=list(map(str, f.list())), S=S,
        field_discriminant=str(nf.disc()), integral_basis=list(map(str,nf.nf_get_zk())),
        minimal_model=list(map(str,aa)), input_to_minimal=list(map(str,transport.tuple())),
        generic_classes=gammas, arithmetic_context=context.record(),
        form=dict(fixed_a=str(a), w_power_basis=coefficients(w),
                  binary_cubic_descending=reduced, matrix=M, hessian=hessian(reduced), reduction_steps=steps),
        sage=sage.version.version, pari=str(pari.version()), input_sha256=sha(folder/'input.json'),
        local_boundary='All primes dividing the integral two-division discriminant, 2 and every real embedding. Good odd primes outside S are unramified for point Kummer classes; their primitive linear form is recorded.')
    write(folder/'field.json', data)
    progress(folder, 'field_prepared', field_discriminant=str(nf.disc()), S=S)


class Collector:
    def __init__(self, folder):
        from sage.all import QQ, ZZ, AA, GF, pari, prime_range
        from research_runtime.local_kummer import LocalSquareclasses
        from reference_small_representative_class import Matrix
        self.folder, self.data = folder, read(folder/'field.json')
        self.R, self.f, self.nf = setup(self.data)
        self.pari, self.ZZ, self.GF = pari, ZZ, GF
        self.policy = read(folder/'policy.json')['constructor']
        self.S = set(self.data['S'])
        self.gammas = [pari.Mod(pari(self.R(g['beta_ascending'])),pari(self.f)) for g in self.data['generic_classes']]
        self.roots = self.f.roots(AA,multiplicities=False)
        self.locals = [LocalSquareclasses(self.nf,p) for p in sorted(self.S)]
        form = self.data['form']
        self.a = int(form['fixed_a'])
        self.w = pari.Mod(pari(self.R(form['w_power_basis'])),pari(self.f))
        self.M = list(map(int,form['matrix']))
        self.form = list(map(int,form['binary_cubic_descending']))
        self.normroots = self.R(list(reversed(self.form))).roots(AA,multiplicities=False)
        self.primes = list(map(int,prime_range(self.policy['smooth_bound']+1)))
        self.fixed_primes = [int(p) for p,e in ZZ(self.a).factor(proof=True)]
        self.support = prod(self.primes)
        pari.addprimes(self.primes+self.fixed_primes)
        self.cols, self.ideals, self.lookup, self.blocks = [], [], {}, {}
        for p in sorted(set(self.primes+self.fixed_primes)|self.S):
            block = []
            for P in pari.idealprimedec(self.nf,p):
                i = len(self.cols)
                h = str(pari.idealhnf(self.nf,P))
                self.cols.append(dict(p=p,hnf=h,e=int(P[2]),f=int(P[3])))
                self.ideals.append(P)
                self.lookup[(p,h)] = i
                block.append(i)
            assert sum(self.cols[i]['e']*self.cols[i]['f'] for i in block) == 3
            self.blocks[p] = block
        base = dict(columns=self.cols, canonical_rational_relations=[dict(parity_columns=[i for i in b if self.cols[i]['e']%2]) for b in self.blocks.values()])
        # Generic-only unramified, totally real quadratic characters protect
        # ideal columns. There is no inherited expected dimension.
        anchor_constraints = [[] for g in self.gammas]
        for p in sorted(self.S):
            for i in self.blocks[p]:
                P = self.ideals[i]
                for row,g in zip(anchor_constraints,self.gammas): row.append(int(pari.idealval(self.nf,g,P))%2)
                if p == 2:
                    gens, _ = local_generators(self.nf,P)
                    for row,g in zip(anchor_constraints,self.gammas):
                        row.extend(int(pari.nfhilbert(self.nf,g,u,P)==-1) for u in gens[1:])
        for row,g in zip(anchor_constraints,self.gammas):
            pol = self.R(coefficients(g))
            row.extend(int(pol(x)<0) for x in self.roots)
        self.anchor_masks = kernel(list(map(pack,anchor_constraints)))
        canonical = {max(r['parity_columns']) for r in base['canonical_rational_relations']}
        anchors, labels, records = [], [], []
        for i,c in enumerate(self.cols):
            p = c['p']
            if i in canonical or p<101 or p>self.policy['anchor_prime_bound'] or c['f']!=1 or c['e']!=1 or self.f.discriminant()%p==0: continue
            root = self.root(i)
            vals = [int(self.R(coefficients(g)).change_ring(GF(p))(root)) for g in self.gammas]
            if not all(vals): continue
            raw = pack(pow(v,(p-1)//2,p)==p-1 for v in vals)
            label = pack((raw&m).bit_count()%2 for m in self.anchor_masks)
            if rank(labels+[label]) == len(labels): continue
            anchors.append(i); labels.append(label); records.append(dict(column=i,p=p,root=root,label=str(label)))
            if len(anchors) == len(self.anchor_masks): break
        assert len(anchors) == len(self.anchor_masks), 'generic-only anchor certificate unresolved'
        self.state = Matrix(base,anchors,self.policy['smooth_bound'])
        for g in self.gammas:
            self.state.add([[i,int(pari.idealval(self.nf,g,P))] for i,P in enumerate(self.ideals) if self.cols[i]['p'] in self.S])
        self.atoms, self.keys, self.parities = [], set(), []
        self.local_cache, self.char_cache, self.history = {}, {}, {}
        write(folder/'factor-base.json', dict(**base, generic_anchor_masks=list(map(str,self.anchor_masks)), generic_anchors=records))
        progress(folder,'collector_ready',columns=len(self.cols),anchors=len(anchors),matrix=self.state.report())

    def root(self,i):
        p, h = self.cols[i]['p'], self.cols[i]['hnf']
        return next(int(x) for x in self.f.change_ring(self.GF(p)).roots(multiplicities=False)
                    if str(self.pari.idealhnf(self.nf,p,self.pari.Mod(self.pari(self.R.gen()-int(x)),self.pari(self.f)))) == h)

    def accept(self,m,n,value):
        pari, nf, M, a = self.pari,self.nf,self.M,self.a
        alpha = a*(M[0]*m+M[1]*n)+(M[2]*m+M[3]*n)*self.w
        assert pari.nfeltnorm(nf,alpha) == a*a*value
        v = pari.nfalgtobasis(nf,alpha)
        v *= pari.denominator(v)
        v /= pari.content(v)
        key = tuple(map(int,v))
        if next(x for x in key if x)<0: key = tuple(-x for x in key)
        if key in self.keys: return None
        alpha = pari.nfbasistoalg(nf,pari.Col(list(key)))
        N = int(pari.nfeltnorm(nf,alpha))
        fac = pari.idealfactor(nf,alpha)
        factors = [[self.lookup[(int(fac[j,0][0]),str(pari.idealhnf(nf,fac[j,0])))],int(fac[j,1])] for j in range(fac.nrows())]
        assert all(e>=0 for i,e in factors)
        assert prod(self.cols[i]['p']**(self.cols[i]['f']*e) for i,e in factors) == abs(N)
        self.state.add(factors)
        vals, parity = dict(factors), 0
        for p in {self.cols[i]['p'] for i in vals}-self.S:
            bit = sum(self.cols[i]['f']*vals.get(i,0) for i in self.blocks[p])%2
            for i in self.blocks[p]: parity ^= ((vals.get(i,0)+bit*self.cols[i]['e'])%2)<<i
        atom = dict(index=len(self.atoms), alpha_ascending=coefficients(alpha), norm=str(N), ideal_factorization=factors, source_pair=[m,n])
        self.atoms.append(atom); self.parities.append(parity); self.keys.add(key)
        return atom

    def batch(self,pairs,path):
        from reference_class_targeted_relations import residues, strip
        pairs = list(pairs)
        c0,c1,c2,c3 = self.form
        values = [((c0*m+c1*n)*m+c2*n*n)*m+c3*n*n*n for m,n in pairs]
        assert all(values)
        accepted = []
        rem = residues(list(map(abs,values)),self.support) if values else []
        for (m,n),value,r in zip(pairs,values,rem):
            if strip(abs(value),gcd(abs(value),r)) != 1: continue
            atom = self.accept(m,n,value)
            if atom is not None: accepted.append(atom)
        write(path,dict(candidate_count=len(pairs),relations=accepted,ending_matrix=self.state.report()))
        progress(self.folder,'collecting_relations',chunk=path.name,atoms=len(self.atoms),new_atoms=len(accepted),matrix=self.state.report())

    def monitor(self):
        from sage.all import prime_range
        kernels = kernel(self.parities)
        used = sorted({i for k in kernels for i in indices(k)})
        if not used: return False
        labels = [('generic',j) for j in range(16)]+[('atom',i) for i in used]
        amap = {i:16+j for j,i in enumerate(used)}
        def factor(label):
            kind,i = label
            if kind=='generic': return self.gammas[i]
            atom = self.atoms[i]
            return int(atom['norm'])*self.pari.Mod(self.pari(self.R(atom['alpha_ascending'])),self.pari(self.f))
        factors = [factor(label) for label in labels]
        locbits = []
        for label,g in zip(labels,factors):
            if label not in self.local_cache:
                pol = self.R(coefficients(g))
                self.local_cache[label] = pack([int(b) for L in self.locals for b in L.signature(g)]+[int(pol(x)<0) for x in self.roots])
            locbits.append(self.local_cache[label])
        candidates = [1<<j for j in range(16)]+[sum(1<<amap[i] for i in indices(k)) for k in kernels]
        strict = [xor(k,candidates) for k in kernel([xor(c,locbits) for c in candidates])]
        # Fixed support-disjoint character bank, independent of desired labels.
        finite = [0]*len(factors); coordinates=[]
        polys = [self.R(coefficients(g)) for g in factors]
        for p in prime_range(self.policy['proof_prime_min'],self.policy['proof_prime_max']+1):
            p=int(p)
            if p in self.S or self.f.discriminant()%p==0 or self.a%p==0: continue
            if any(v.denominator()%p==0 for pol in polys for v in pol): continue
            roots = sorted(map(int,self.f.change_ring(self.GF(p)).roots(multiplicities=False)))
            if len(roots)!=3: continue
            vals = [[int(pol.change_ring(self.GF(p))(x)) for x in roots] for pol in polys]
            if any(not v for row in vals for v in row): continue
            offset=len(coordinates)
            for j,row in enumerate(vals):
                finite[j] |= pack(pow(v,(p-1)//2,p)==p-1 for v in row)<<offset
            coordinates.extend([dict(p=p,root=x) for x in roots])
            if len(coordinates)>=self.policy['character_coordinates']: break
        assert rank(finite[:16])==16, 'generic character bank insufficient; independence UNKNOWN'
        selected, chars = [], finite[:16].copy()
        for ordinal,k in enumerate(strict):
            character=xor(k,finite)
            if rank(chars+[character]) == len(chars): continue
            selected.append(dict(ordinal=ordinal, factor_labels=[list(labels[j]) for j in indices(k)],proof_character=str(character)))
            chars.append(character)
            if len(selected)==2: break
        result=dict(status='CANDIDATES_PENDING_INDEPENDENT_REPLAY' if selected else 'NO_NEW_CLASS_YET',
            columns=self.cols,atoms=self.atoms,classes=selected,coordinates=coordinates,
            generic_strict_dimension=len(kernel(locbits[:16])),outside_kernel_count=len(kernels),
            generic_characters=list(map(str,finite[:16])),field_sha256=sha(self.folder/'field.json'))
        write(self.folder/'extraction-latest.json',result)
        progress(self.folder,'class_extraction',atoms=len(self.atoms),dependencies=len(kernels),new_classes=len(selected))
        if selected:
            write(self.folder/'classes.json',result)
            return True
        return False


def reduced_lattice(q,r,hessian):
    # Exact existing target_small_conductor_prime_ideals.reduce_lattice policy.
    A,B,C=map(int,hessian)
    def dot(u,v): return 2*A*u[0]*v[0]+B*(u[0]*v[1]+u[1]*v[0])+2*C*u[1]*v[1]
    u,v=[q,0],[r,1]
    for _ in range(1000):
        if dot(v,v)<dot(u,u):u,v=v,u
        length,inner=dot(u,u),dot(u,v)
        if 2*abs(inner)<=length:break
        k=(2*inner+length)//(2*length)
        v=[v[i]-k*u[i] for i in range(2)]
    else:raise ArithmeticError('lattice reduction cap')
    assert abs(u[0]*v[1]-u[1]*v[0])==q and all((t[0]-r*t[1])%q==0 for t in [u,v])
    return u,v


def construct(folder):
    C=Collector(folder)
    box=C.policy['initial_box']
    for lo in range(1,box+1,32):
        C.batch(((m,n) for n in range(lo,min(lo+32,box+1)) for m in range(-box,box+1) if gcd(m,n)==1),folder/f'box-{lo:04d}.json')
    if C.monitor():return
    done=0
    for level in C.policy['strip_levels']:
        eligible=[]
        for i,c in enumerate(C.cols):
            p=c['p']
            if i in C.state.canonical or C.state.anchors>>i&1 or p>C.policy['smooth_bound'] or c['f']!=1 or c['e']!=1 or C.history.get(i,0)>=level:continue
            if C.f.discriminant()%p==0 or C.a%p==0:continue
            projection=C.state.projection(i)
            if projection:eligible.append((i,projection))
        span,selected={},[]
        for i,v in eligible:
            r=reduce(v,span)
            if r:span[r.bit_length()-1]=r;selected.append(i)
            if len(selected)>=C.policy['targets_per_level']:break
        selected += [i for i,v in eligible if i not in selected][:C.policy['targets_per_level']-len(selected)]
        selected.sort()
        write(folder/f'level-{level}-selection.json',dict(columns=selected,eligible=len(eligible),projection_rank=len(span)))
        for i in selected:
            if not C.state.projection(i):continue
            p=C.cols[i]['p'];root=C.root(i)
            wbar=int(C.R(C.data['form']['w_power_basis']).change_ring(C.GF(p))(root))
            bm=(C.a*C.M[0]+wbar*C.M[2])%p;bn=(C.a*C.M[1]+wbar*C.M[3])%p
            if not bm:
                write(folder/f'level-{level}-column-{i}-skip.json',dict(reason='vertical lattice chart omitted by frozen policy',p=p))
                continue
            v1,v2=reduced_lattice(p,-bn*pow(bm,-1,p)%p,C.data['form']['hessian'])
            scale=2**96
            slopes=[int((((x*v2[1]-v2[0])/(v1[0]-x*v1[1]))*scale).floor()) for x in C.normroots]
            pairs=set()
            for v in range(C.history.get(i,0)+1,level+1):
                for slope in slopes:
                    center=slope*v//scale
                    for u in range(center-1,center+2):
                        if gcd(u,v)!=1:continue
                        m,n=u*v1[0]+v*v2[0],u*v1[1]+v*v2[1]
                        if gcd(m,n)==1:pairs.add((m,n))
            C.batch(sorted(pairs),folder/f'level-{level}-column-{i}.json')
            C.history[i]=level;done+=1
            if done%C.policy['monitor_interval']==0 and C.monitor():return
        if C.monitor():return
    write(folder/'construction-terminal.json',dict(status='NO_NEW_CLASS_WITHIN_FROZEN_BANK',atoms=len(C.atoms),matrix=C.state.report()))


def verify_classes(folder):
    """Independent local Hilbert and prime-valuation replay, not character-only."""
    from sage.all import QQ, ZZ, AA, GF, matrix, pari
    d=read(folder/'field.json'); c=read(folder/'classes.json')
    assert c['field_sha256']==sha(folder/'field.json')
    R,f,nf=setup(d);S=set(d['S']);cols=c['columns']
    roots=f.roots(AA,multiplicities=False)
    factors={('generic',i):pari.Mod(pari(R(g['beta_ascending'])),pari(f)) for i,g in enumerate(d['generic_classes'])}
    used={i for record in c['classes'] for kind,i in record['factor_labels'] if kind=='atom'}
    projected={};prime_cache={}
    for i in sorted(used):
        atom=c['atoms'][i];a=pari.Mod(pari(R(atom['alpha_ascending'])),pari(f));N=ZZ(pari.nfeltnorm(nf,a))
        assert str(N)==atom['norm']
        vals=dict(atom['ideal_factorization']);assert all(e>=0 for e in vals.values())
        assert prod(ZZ(cols[j]['p'])**(cols[j]['f']*e) for j,e in vals.items())==abs(N)
        nv={}
        for j,e in vals.items():
            p=cols[j]['p'];nv[p]=nv.get(p,0)+cols[j]['f']*e
            if p not in prime_cache:prime_cache[p]={str(pari.idealhnf(nf,P)):P for P in pari.idealprimedec(nf,p)}
            P=prime_cache[p][cols[j]['hnf']]
            assert int(P[2])==cols[j]['e'] and int(P[3])==cols[j]['f']
            assert int(pari.idealval(nf,a,P))==e
        projected[i]={j:vals.get(j,0)+col['e']*nv[col['p']] for j,col in enumerate(cols) if col['p'] in nv}
        factors[('atom',i)]=N*a
    local=[]
    for p in sorted(S):
        for P in pari.idealprimedec(nf,p):
            gens,pairing=local_generators(nf,P)
            signatures={k:[int(pari.nfhilbert(nf,a,b,P)==-1) for b in gens] for k,a in factors.items()}
            for record in c['classes']:
                assert all(sum(signatures[tuple(label)][j] for label in record['factor_labels'])%2==0 for j in range(len(gens)))
            local.append(dict(p=p,hnf=str(pari.idealhnf(nf,P)),generators=list(map(str,gens)),hilbert_gram=pairing))
        progress(folder,'independent_local_replay',p=p,used_atoms=len(used))
    for record in c['classes']:
        ids=[i for kind,i in record['factor_labels'] if kind=='atom']
        for j,col in enumerate(cols):
            if col['p'] not in S:assert sum(projected[i].get(j,0) for i in ids)%2==0
        for root in roots:assert sum(int(R(coefficients(factors[tuple(label)]))(root)<0) for label in record['factor_labels'])%2==0
    points=read(folder/'seed.json')['points'];model=list(map(F,read(folder/'seed.json')['curve']))
    from sage.all import EllipticCurve
    E=EllipticCurve(QQ,model);Emin=EllipticCurve(QQ,list(map(QQ,d['minimal_model'])));transport=E.isomorphism_to(Emin)
    for g,pt in zip(d['generic_classes'],points):
        x,y=transport(E(list(map(QQ,pt)))).xy();pol=R(g['beta_ascending']);den=ZZ(-pol[1]).sqrt()
        assert den**2==-pol[1] and pol==den**2*(4*x-R.gen()) and ZZ(pol[0]).gcd(den)==1
        assert QQ(pari.nfeltnorm(nf,pari.Mod(pari(pol),pari(f))))==QQ(g['norm']) and QQ(g['norm']).is_square()
    chars=[0]*(16+len(c['classes']))
    for k,coordinate in enumerate(c['coordinates']):
        p,root=coordinate['p'],coordinate['root']
        assert ZZ(p).is_prime(proof=True) and p not in S and f.discriminant()%p and f(root)%p==0
        bits={}
        for label,a in factors.items():
            v=int(R(coefficients(a)).change_ring(GF(p))(root));assert v
            bits[label]=int(pow(v,(p-1)//2,p)==p-1)
        row=[bits[('generic',i)] for i in range(16)]+[sum(bits[tuple(label)] for label in record['factor_labels'])%2 for record in c['classes']]
        for j,bit in enumerate(row):chars[j]|=bit<<k
    assert rank(chars)==len(chars)
    write(folder/'classes-verified.json',dict(status='PASS_EXACT_STRICT_CLASS_AND_INDEPENDENCE',additional_classes=len(c['classes']),combined_kummer_rank=len(chars),local_hilbert_certificates=local,classes_sha256=sha(folder/'classes.json'),ordinary_ideal_class_increment='UNKNOWN',rational_solubility='UNKNOWN'))


def compact(folder):
    from sage.all import QQ, ZZ, pari
    d=read(folder/'field.json');c=read(folder/'classes.json');v=read(folder/'classes-verified.json')
    assert v['classes_sha256']==sha(folder/'classes.json')
    R,f,nf=setup(d);cols=c['columns'];S=set(d['S'])
    gammas=[pari.Mod(pari(R(g['beta_ascending'])),pari(f)) for g in d['generic_classes']]
    gcds=[pari.idealadd(nf,ZZ(g['norm']).sqrt(),gamma) for gamma,g in zip(gammas,d['generic_classes'])]
    bad=[(p,str(pari.idealhnf(nf,P)),P) for p in sorted(S) for P in pari.idealprimedec(nf,p)]
    for number,record in enumerate(c['classes']):
        generic=[i for kind,i in record['factor_labels'] if kind=='generic']
        exponents={};bases=[];powers=[]
        for kind,i in record['factor_labels']:
            if kind=='generic':value=gammas[i]
            else:
                atom=c['atoms'][i];vals=dict(atom['ideal_factorization']);nv={}
                for j,e in vals.items():p=cols[j]['p'];nv[p]=nv.get(p,0)+cols[j]['f']*e
                for j,col in enumerate(cols):
                    if col['p'] in nv:exponents[j]=exponents.get(j,0)+vals.get(j,0)+col['e']*nv[col['p']]
                value=int(atom['norm'])*pari.Mod(pari(R(atom['alpha_ascending'])),pari(f))
            bases.append(pari.nfalgtobasis(nf,value));powers.append(1)
        ideal=pari.idealhnf(nf,1);steps=[]
        def multiply(I,source):
            nonlocal ideal
            product=pari.idealmul(nf,ideal,I);J,a=pari.idealred(nf,[product,1])
            assert pari.idealmul(nf,J,a)==product
            steps.append(dict(source=source,ideal_before=str(ideal),factor_ideal=str(I),ideal_after=str(J),multiplier=str(a)))
            ideal=J;bases.append(a);powers.append(-2)
        for i in generic:multiply(gcds[i],dict(generic=i))
        for j,e in sorted(exponents.items()):
            col=cols[j]
            if col['p'] in S:continue
            assert e%2==0
            if e:multiply(pari.idealpow(nf,pari(col['hnf']),e//2),dict(column=j,exponent=e//2))
        for p,h,P in bad:
            av=sum(e for j,e in exponents.items() if cols[j]['p']==p and cols[j]['hnf']==h)
            gv=sum(int(pari.idealval(nf,gammas[i],P)) for i in generic)
            gcdv=sum(int(pari.idealval(nf,gcds[i],P)) for i in generic)
            assert (av+gv)%2==0
            correction=(av+gv)//2-gcdv
            if correction:multiply(pari.idealpow(nf,P,correction),dict(p=p,hnf=h,exponent=correction))
        write(folder/f'compact-{number}-circuit.json',dict(steps=steps,bases=list(map(str,bases)),exponents=powers,class_index=number))
        progress(folder,'expanding_reduced_class',class_index=number,factors=len(bases))
        beta=pari.nfbasistoalg(nf,pari.nffactorback(nf,bases,powers))
        assert pari.idealhnf(nf,beta)==pari.idealpow(nf,ideal,2)
        N=QQ(pari.nfeltnorm(nf,beta));assert N>0 and N.is_square()
        bits=max(max(abs(QQ(x).numerator()).nbits(),QQ(x).denominator().nbits()) for x in coefficients(beta))
        write(folder/f'cover-{number}.json',dict(beta_ascending=coefficients(beta),norm=str(N),positive_norm_square_root=str(N.sqrt()),max_coefficient_bits=int(bits),half_ideal_hnf=str(ideal),class_index=number,circuit_sha256=sha(folder/f'compact-{number}-circuit.json'),status='COMPACT_REPRESENTATIVE_PENDING_REPLAY'))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['prepare','construct','verify_classes','compact']);ap.add_argument('--folder',type=Path,required=True)
    args=ap.parse_args();folder=args.folder.resolve();guard(folder)
    from sage.all import pari
    pari.allocatemem(64000000,268435456,silent=True)
    globals()[args.mode](folder)


if __name__=='__main__':main()

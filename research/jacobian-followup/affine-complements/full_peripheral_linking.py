"""Full peripheral homology and rational core linking for iterated knot covers.

Arbitrary commuting longitude action is allowed. The rational matrix requires
all filling cores rationally null; integral orders are reported separately.
This is not the complete torsion linking form of the filled manifold.
"""
from pathlib import Path
import sys,math
import sympy as S
from sympy.matrices.normalforms import hermite_normal_form

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P))
from peripheral_boundary_lattices import boundary_lattices
from composed_peripheral_cover import reduced_h1


def inverse(p):
    out=[0]*len(p)
    for i,j in enumerate(p):out[j]=i
    return tuple(out)


class Flows:
    """Destination/edge-flow semidirect product; binary signed powers."""
    def __init__(self,rho):
        self.degree=len(rho[0]);self.width=len(rho)*self.degree
        self.ident=tuple(range(self.degree))
        if any(tuple(sorted(p))!=self.ident for p in rho):
            raise ValueError('permutations on the same nonempty sheet set required')
        self.unit=(self.ident,[[0]*self.width for _ in self.ident])
        self.generators=[]
        for g,p in enumerate(rho):
            f=[[0]*self.width for _ in self.ident]
            for s in self.ident:f[s][g*self.degree+s]=1
            self.generators.append((inverse(p),f))

    def mul(self,a,b):
        da,fa=a;db,fb=b
        return tuple(db[da[s]] for s in self.ident),[
            [x+y for x,y in zip(fa[s],fb[da[s]])] for s in self.ident]

    def inv(self,a):
        d,f=a;back=inverse(d)
        return back,[[-x for x in f[back[s]]] for s in self.ident]

    def power(self,a,k):
        if k<0:return self.power(self.inv(a),-k)
        out=self.unit
        while k:
            if k&1:out=self.mul(out,a)
            a=self.mul(a,a);k//=2
        return out

    def product(self,*parts):
        out=self.unit
        for a in parts:out=self.mul(out,a)
        return out


def peripheral_cover(p,q,cables,rho,reverse=False):
    if not rho or not rho[0] or len(rho)!=2+len(cables):
        raise ValueError('one permutation per knot-presentation generator required')
    if p<1 or q<1 or math.gcd(p,q)!=1:
        raise ValueError('positive coprime initial torus parameters required')
    if any(r<1 or math.gcd(r,k)!=1 for r,k in cables):
        raise ValueError('positive winding and coprime signed cable parameters required')
    f=Flows(rho);a=f.generators
    u,v,g=S.gcdex(q,p);assert g==1
    mu=f.product(f.power(a[0],int(u)),f.power(a[1],int(v)))
    lam=f.product(f.power(a[0],p),f.power(mu,-p*q))
    relations=[f.product(f.power(a[0],p),f.power(a[1],-q))]
    for j,(r,k) in enumerate(cables,2):
        u,v,g=S.gcdex(r,k);assert g==1;u=int(u);v=-int(v)
        fibre=f.product(f.power(mu,k),f.power(lam,r))
        relations.append(f.product(f.power(a[j],r),f.inv(fibre)))
        mu=f.product(f.power(mu,u),f.power(lam,v),f.power(a[j],-v))
        lam=f.product(f.power(mu,-r*k),fibre)
    if any(rel[0]!=f.ident for rel in relations):raise ValueError('knot-group relation violated')
    # Use destinations consistently; inverting both commuting actions gives
    # the same stabilizer lattices, also returned as raw monodromy below.
    boundaries=boundary_lattices(mu[0],lam[0])
    edges=[(g,s) for g in range(len(a)) for s in f.ident]
    parent=list(f.ident);tree=set()
    def root(s):
        while parent[s]!=s:s=parent[s]
        return s
    for g,s in reversed(edges) if reverse else edges:
        t=a[g][0][s]
        if root(s)!=root(t):parent[root(s)]=root(t);tree.add((g,s))
    if len(tree)!=f.degree-1:raise ValueError('connected covering required')
    keep=[edge for edge in edges if edge not in tree]
    indices=[g*f.degree+s for g,s in keep]
    project=lambda row:[row[i] for i in indices]
    rows=[project(rel[1][s]) for rel in relations for s in f.ident]
    for b in boundaries:
        s=b['seed'];e=b['normal_index'];h=b['shear'];l=b['longitudinal_degree']
        m=f.power(mu,e);L=f.product(f.power(mu,h),f.power(lam,l))
        assert m[0][s]==L[0][s]==s
        b.update(meridian=project(m[1][s]),longitude=project(L[1][s]))
    return dict(degree=f.degree,relation_rows=rows,boundaries=boundaries,keep=keep,
                meridian_monodromy=inverse(mu[0]),longitude_monodromy=inverse(lam[0]),
                word_expansion_used=False)


class RowLattice:
    """Independent HNF row basis; exact membership and order in the quotient."""
    def __init__(self,rows):
        self.basis=hermite_normal_form(S.Matrix(rows).T).T

    def coordinates(self,row):
        sol=S.linsolve((self.basis.T,S.Matrix(row)))
        if sol is S.EmptySet or not sol:return None
        sol=next(iter(sol))
        assert not any(x.free_symbols for x in sol)
        return sol

    def order(self,row):
        c=self.coordinates(row)
        if c is None:return None  # infinite order
        return math.lcm(*(int(x.q) for x in c))

    def contains(self,row):return self.order(row)==1


def rational_core_linking(data):
    rows=data['relation_rows'];bs=data['boundaries'];count=len(bs)
    meridians=[b['meridian'] for b in bs];longitudes=[b['longitude'] for b in bs]
    closed=rows+meridians
    lattice=RowLattice(closed)
    orders=[lattice.order(row) for row in longitudes]
    if any(q is None for q in orders):
        raise ValueError('every filling core must be rationally null-homologous')
    null=S.Matrix(rows).nullspace()
    basis=S.Matrix.hstack(*null) if null else S.zeros(len(rows[0]),0)
    mu=S.Matrix(meridians)*basis;lam=S.Matrix(longitudes)*basis
    if mu.rank()!=count:raise ValueError('independent rational meridians required')
    coefficients=[]
    for i in range(count):
        sol=S.linsolve((mu.T,lam.row(i).T))
        if not sol:raise ValueError('longitude not a rational meridian combination')
        row=next(iter(sol));assert not any(x.free_symbols for x in row)
        coefficients.append(list(row))
    F=S.Matrix(coefficients)
    if F!=F.T:raise ValueError('rational linking symmetry failed')
    e=S.Matrix([b['normal_index'] for b in bs]);h=S.Matrix([b['shear'] for b in bs])
    if F*e!=h:raise ValueError('full peripheral surface balance failed')
    exterior_lattice=RowLattice(rows)
    row_denominators=[];defect_orders=[]
    for i,row in enumerate(coefficients):
        denom=math.lcm(*(int(x.q) for x in row))
        defect=[denom*longitudes[i][j]-sum(int(denom*row[k])*meridians[k][j] for k in range(count))
                for j in range(len(rows[0]))]
        order=exterior_lattice.order(defect)
        if orders[i]%denom or order!=orders[i]//denom:
            raise ValueError('integral torsion defect disagrees with filling core order')
        row_denominators.append(denom);defect_orders.append(order)
    balance=[sum(b['normal_index']*b['longitude'][j]-b['shear']*b['meridian'][j] for b in bs)
             for j in range(len(rows[0]))]
    if not exterior_lattice.contains(balance):raise ValueError('integral lifted-surface relation failed')
    return dict(matrix=[[str(x) for x in row] for row in coefficients],
                core_orders=orders,all_cores_integrally_zero=all(q==1 for q in orders),
                rational_row_denominators=row_denominators,cleared_relation_torsion_orders=defect_orders,
                closed_H1=reduced_h1(rows,meridians),
                normal_indices=list(map(int,e)),shears=list(map(int,h)),
                longitudinal_degrees=[b['longitudinal_degree'] for b in bs],
                weighted_balance='F*e=a',integral_lifted_surface_balance_verified=True,
                complete_filled_manifold_torsion_linking_form_computed=False)

"""Exact small matrix model and admissible affine-root equations over F3."""

from itertools import product
from verify_integer_lifts import nullspace


I3 = (1,0,0,0,1,0,0,0,1)


def multiply(a,b):
    n = 3 if len(a)==9 else 4
    assert len(a)==len(b)==n*n
    return tuple(sum(a[n*i+k]*b[n*k+j] for k in range(n))%3 for i in range(n) for j in range(n))


def identity(n):
    return tuple(int(i==j) for i in range(n) for j in range(n))


def power(a,n):
    value=identity(3 if len(a)==9 else 4)
    for _ in range(n):value=multiply(value,a)
    return value


def determinant(a):
    return (a[0]*(a[4]*a[8]-a[5]*a[7])-a[1]*(a[3]*a[8]-a[5]*a[6])+a[2]*(a[3]*a[7]-a[4]*a[6]))%3


def inverse(a):
    det=determinant(a);assert det
    return tuple(v*pow(det,-1,3)%3 for v in [
        a[4]*a[8]-a[5]*a[7],a[2]*a[7]-a[1]*a[8],a[1]*a[5]-a[2]*a[4],
        a[5]*a[6]-a[3]*a[8],a[0]*a[8]-a[2]*a[6],a[2]*a[3]-a[0]*a[5],
        a[3]*a[7]-a[4]*a[6],a[1]*a[6]-a[0]*a[7],a[0]*a[4]-a[1]*a[3]])


def image(a,v):
    return tuple(sum(a[3*i+j]*v[j] for j in range(3))%3 for i in range(3))


def rank(rows,n):
    return n-len(nullspace(rows,n,3))


def root(a):
    for j in range(3):
        column=tuple((a[3*i+j]-int(i==j))%3 for i in range(3))
        first=next((v for v in column if v),None)
        if first is not None:return tuple(v*pow(first,-1,3)%3 for v in column)
    raise AssertionError("identity is not a transvection")


def model():
    elements=set()
    for v in product(range(3),repeat=3):
        first=next((x for x in v if x),None)
        if first!=1:continue
        for f in product(range(3),repeat=3):
            if not any(f) or sum(a*b for a,b in zip(v,f))%3:continue
            elements.add(tuple((int(i==j)+v[i]*f[j])%3 for i in range(3) for j in range(3)))
    tau=(1,1,0,0,1,0,0,0,1)
    elements=[tau]+sorted(elements-{tau});assert len(elements)==104
    index={a:i for i,a in enumerate(elements)}
    inverses=[power(a,2) for a in elements]
    assert all(multiply(a,b)==I3 for a,b in zip(elements,inverses))
    forward=[[index[multiply(multiply(inverses[o],a),elements[o])] for a in elements] for o in range(104)]
    backward=[[index[multiply(multiply(elements[o],a),inverses[o])] for a in elements] for o in range(104)]
    general_linear=[a for a in product(range(3),repeat=9) if determinant(a)]
    assert len(general_linear)==11232
    centralizer=[a for a in general_linear if multiply(a,tau)==multiply(tau,a)]
    assert len(centralizer)==108
    actions=[[index[multiply(multiply(g,a),inverse(g))] for a in elements] for g in centralizer]
    return elements,forward,backward,centralizer,actions


def lift_matrix(crossings,labels,elements):
    n=len(labels);roots=[root(elements[i]) for i in labels];rows=[]
    for o,a,b,sign in crossings:
        word=[-(o+1),a+1,o+1,-(b+1)] if sign else [o+1,a+1,-(o+1),-(b+1)]
        prefix=I3;block=[[0]*n for _ in range(3)]
        for letter in word:
            i=abs(letter)-1;value=image(prefix,roots[i])
            for j in range(3):block[j][i]=(block[j][i]+(1 if letter>0 else -1)*value[j])%3
            base=elements[labels[i]]
            prefix=multiply(prefix,base if letter>0 else power(base,2))
        assert prefix==I3
        rows.extend(block)
    return rows,nullspace(rows,n,3),roots


def linear_group(generators,cap=5616):
    values=[I3];seen={I3}
    for a in values:
        for b in generators:
            value=multiply(a,b)
            if value not in seen:
                assert len(values)<cap
                seen.add(value);values.append(value)
    return values

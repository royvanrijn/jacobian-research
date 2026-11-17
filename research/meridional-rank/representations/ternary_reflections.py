"""The remaining rank-one class in GL3(F3) and its affine-root equations."""

from itertools import product
from ternary_transvections import determinant,identity,inverse,image,multiply,root
from verify_integer_lifts import nullspace


def model():
    elements=set()
    for v in product(range(3),repeat=3):
        first=next((x for x in v if x),None)
        if first!=1:continue
        for f in product(range(3),repeat=3):
            if sum(x*y for x,y in zip(v,f))%3!=1:continue
            elements.add(tuple((int(i==j)+v[i]*f[j])%3 for i in range(3) for j in range(3)))
    tau=(2,0,0,0,1,0,0,0,1)
    elements=[tau]+sorted(elements-{tau});assert len(elements)==117
    index={m:i for i,m in enumerate(elements)}
    table=[[index[multiply(multiply(a,b),a)] for b in elements] for a in elements]
    gl=[m for m in product(range(3),repeat=9) if determinant(m)]
    centralizer=[m for m in gl if multiply(m,tau)==multiply(tau,m)]
    assert len(centralizer)==96
    actions=[[index[multiply(multiply(g,a),inverse(g))] for a in elements] for g in centralizer]
    return elements,table,centralizer,actions


def lift_matrix(crossings,labels,elements):
    n=len(labels);roots=[root(elements[i]) for i in labels];rows=[]
    for o,a,b,positive in crossings:
        word=[-(o+1),a+1,o+1,-(b+1)] if positive else [o+1,a+1,-(o+1),-(b+1)]
        prefix=identity(3);block=[[0]*n for _ in range(3)]
        for letter in word:
            i=abs(letter)-1;base=elements[labels[i]]
            local=roots[i] if letter>0 else tuple(-v%3 for v in image(base,roots[i]))
            value=image(prefix,local)
            for j in range(3):block[j][i]=(block[j][i]+value[j])%3
            prefix=multiply(prefix,base)
        assert prefix==identity(3);rows.extend(block)
    return rows,nullspace(rows,n,3),roots

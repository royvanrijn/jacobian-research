"""Simon simultaneous Selmer linear algebra over explicitly prepared data.

This module never initializes or certifies a number field or BNF. Callers
supply a certified BNF, exact curve/labelled generator and complete known
rational discriminant support. Local structures and the result are shared
across requests; missing support fails before descent.
"""
from hashlib import sha256
from pathlib import Path


def prepared_simon(arithmetic, curve, bnf, theta, *, primes, class_data_id, discover):
    from sage.all import pari
    from sage.env import SAGE_EXTCODE
    from run_elkies_2026_relative_2selmer_checkpointed import SIMON_GP_FUNCTION
    from .pari_context import identify
    simon = Path(SAGE_EXTCODE)/'pari'/'simon'
    names = ('ellQ.gp', 'ell.gp', 'qfsolve.gp', 'resultant3.gp')
    implementation = {name: sha256((simon/name).read_bytes()).hexdigest() for name in names}
    implementation['function'] = sha256(SIMON_GP_FUNCTION.encode()).hexdigest()
    implementation['adapter'] = 'explicit-frame-and-support-1'
    context, variable = identify(bnf.nf_get_pol())
    frame = {'algebra': context.key, 'variable': variable, 'maximal_basis': list(map(str, bnf.nf_get_zk()))}
    primes = sorted({2, *map(int, primes)})
    identity = {'curve': [str(curve[i]) for i in range(5)], 'theta': str(theta),
                'frame': frame, 'class_data_id': class_data_id, 'implementation': implementation}

    def build():
        for name in names:
            pari.read(str(simon/name))
        pari('DEBUGLEVEL_ell=0;LIMBIGPRIME=0;LIM1=0;LIM3=0;LIMTRIV=0;')
        pari.addprimes(primes)
        # Verify that supplied support covers every local place before ppinit.
        derivative = 3*theta**2+2*curve[1]*theta+curve[3]
        ideal = abs(pari.idealadd(bnf, derivative, pari('b->b.index')(bnf)))
        support_integer = int(ideal[0, 0])*2*int(pari.numerator(pari('e->e.disc')(curve)))
        remaining = abs(support_integer)
        for prime in primes:
            while remaining and remaining % prime == 0:
                remaining //= prime
        if remaining != 1:
            raise ArithmeticError('unprepared descent discriminant support')
        places = [p for p in primes if support_integer % p == 0]
        pp = []
        for prime in places:
            row = arithmetic._fact('two-torsion/simon-local-frame',
                {'frame': frame, 'prime': prime, 'implementation': implementation},
                lambda prime=prime: {'binary': arithmetic._blob(pari('ppinit')(bnf, prime))}, True)
            pp.append(arithmetic._restore(row['binary']))
        pari('ECContextPrimes='+str(pari(places)))
        pari('ECContextPP='+str(pari(pp)))
        pari('eccontextpp(p)={for(i=1,#ECContextPrimes,if(ECContextPrimes[i]==p,return(ECContextPP[i])));error("unprepared descent place")}')
        pari('eccontextsupport(n)={my(m=abs(n));if(m==0,error("zero descent support"));for(i=1,#ECContextPrimes,my(p=ECContextPrimes[i]);m=m/p^valuation(m,p));if(m!=1,error("unprepared descent discriminant support"));select(p->valuation(n,p)!=0,ECContextPrimes)}')
        source = SIMON_GP_FUNCTION.replace('pp=ppinit(bnf.nf,p);', 'pp=eccontextpp(p);')
        source = source.replace('descentprimes=factorint(badideal[1,1]*2)[,1];',
                                'descentprimes=eccontextsupport(badideal[1,1]*2);')
        source = source.replace('badprimes=factorint(badideal[1,1]*2*numerator(ell.disc))[,1];',
                                'badprimes=eccontextsupport(badideal[1,1]*2*numerator(ell.disc));')
        if 'factorint(' in source or 'pp=ppinit(' in source:
            raise ArithmeticError('descent bypassed prepared arithmetic')
        for definition in source.split('/* ELKIES_R17_GP_DEFINITION_SPLIT */'):
            pari(definition)
        raw = pari('ell2selmer_basis_gen')(curve, bnf, 1, theta)
        return {'binary': arithmetic._blob(raw)}
    row = arithmetic._fact('arithmetic/simon-complete-matrix', identity, build, discover)
    return arithmetic._restore(row['binary'])

"""Compare a fresh replay with retained mathematics and explicit source history."""

from hashlib import sha256


SNAPSHOTS = {
    'elliptic-curves/cas/verify_mestre_fermigier_two_section_generic_rank13.py': (
        'archive/repository-cleanup-2026-09-12/research__elliptic-curves__cas__verify_mestre_fermigier_two_section_generic_rank13.py.txt',
        'c548a3c9febf37f3e8d8f92c857a50015ba87808c5a10815a63f3393c42498b8',
    ),
    'elliptic-curves/cas/verify_mestre_parent_and_label_audits.sage': (
        'archive/repository-cleanup-2026-09-12/research__elliptic-curves__cas__verify_mestre_parent_and_label_audits.sage.txt',
        '83720b2c122a0d715915a22ec4e1dd1811139b00f80ffada15df5dffb19aed63',
    ),
}


def check_retained_replay(actual, retained, root):
    """Require identical mathematical output; never rewrite historical hashes.

    The two permitted changes disable the rejected entry point and add this
    comparison. Every other source must retain its original bytes. Current
    checker source is separately pinned by MATH_STATUS.json.
    """
    if {k: v for k, v in actual.items() if k != 'sources'} != {
        k: v for k, v in retained.items() if k != 'sources'
    }:
        raise ValueError('retained mathematical payload differs')
    if actual['sources'].keys() != retained['sources'].keys():
        raise ValueError('retained source roster differs')
    for name, old_hash in retained['sources'].items():
        now = sha256((root / name).read_bytes()).hexdigest()
        if actual['sources'][name] != now:
            raise ValueError('current source binding differs: ' + name)
        if now == old_hash:
            continue
        if name not in SNAPSHOTS:
            raise ValueError('undeclared historical source change: ' + name)
        path, expected = SNAPSHOTS[name]
        if old_hash != expected or sha256((root / path).read_bytes()).hexdigest() != expected:
            raise ValueError('pinned historical source differs: ' + name)

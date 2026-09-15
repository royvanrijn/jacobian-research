"""The v1 exact replay with explicit already-compatible endpoint binding."""
from cancellation_basis_accessibility import PRIOR, read, need
from verify_cancellation_basis_accessibility import verify as replay


def verify(row, plan, folder):
    endpoints = next(r['packet'] for r in read(PRIOR/'endpoint-packets.json')['rows'] if r['case'] == row['id'])
    result = read(folder/'scan.json'); maps = read(folder/'maps.json')
    need(endpoints['curve'] == maps['seed']['curve'], 'endpoint equation differs from search model')
    for target in result['targets']:
        need(target['point'] == endpoints['points'][target['endpoint_index']], 'endpoint witness changed')
    need(len(result['endpoint_statuses']) == len(endpoints['points']), 'endpoint classification omitted rows')
    return replay(row, plan, folder)

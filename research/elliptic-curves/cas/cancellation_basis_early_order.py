"""A target-blind order that exposes the latest admitted generator block first."""


def visit_order(centres, initial_rank, fresh_start, refreshed):
    if not centres:
        raise ValueError('empty bank')
    dimension = len(centres[0]['representative'])
    if not 0 < initial_rank <= fresh_start <= dimension or any(
            len(row['representative']) != dimension for row in centres):
        raise ValueError('inconsistent basis dimensions')
    if not refreshed:
        if dimension != initial_rank:
            raise ValueError('initial bank changed dimension')
        return list(range(len(centres)))
    if fresh_start >= dimension:
        raise ValueError('refresh added no generator')
    # Stable within each group: the underlying verified bank is unchanged.
    def lane(i):
        word = centres[i]['representative']
        if any(word[fresh_start:]): return 0
        if any(word[initial_rank:]): return 1
        return 2
    return sorted(range(len(centres)), key=lane)

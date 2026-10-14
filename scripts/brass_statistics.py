"""Resolve BRASS metadata without silently inventing measured estimates."""
import math


def resolve_statistics(row, policy='native_defaults'):
    if policy not in ('require', 'native_defaults'):
        raise ValueError('Unknown missing_statistics_policy: ' + policy)
    if row.get('sex') not in ('female', 'male'):
        raise ValueError('BRASS requires explicit sex metadata for chromosome settings')
    values = [str(row.get(k, '') if row.get(k) is not None else '').strip() for k in ('purity', 'ploidy')]
    missing = [not value for value in values]
    if any(missing) and policy == 'require':
        raise ValueError('Missing purity/ploidy under require policy')
    purity, ploidy = [default if absent else float(value) for value, absent, default in zip(values, missing, (0.75, 2.0))]
    if not (math.isfinite(purity) and 0 < purity <= 1 and math.isfinite(ploidy) and ploidy > 0):
        raise ValueError('Invalid purity/ploidy; malformed values cannot trigger fallback')
    state = 'DEFAULT' if all(missing) else 'MIXED' if any(missing) else 'REAL'
    if state == 'DEFAULT':
        # The pinned parser recognizes '?' and sets native defaults + DEFAULT.
        # Supplying literal numeric defaults would incorrectly mark them REAL.
        text = 'rho ?\nPloidy ?\n'
    else:
        text = 'rho {}\nPloidy {}\n'.format(purity, ploidy)
        if state == 'MIXED':
            text += '# BRASS_PIPELINE_STATE MIXED\n'
    found = 'N' if row['sex'] == 'female' else 'Y'
    text += 'GenderChr chrY\nGenderChrFound {}\n'.format(found)
    return {'purity': purity, 'ploidy': ploidy, 'statistics_state': state,
            'gender_chr': 'chrY', 'gender_chr_found': found,
            'statistics_values_source': {k: 'BRASS 6.3.4 default (not measured)' if m else 'sample metadata' for k,m in zip(('purity','ploidy'),missing)}}, text

"""Frozen v3 frequency policy, independent of any previous output directory."""
import json
import math
from methanol_lamda.isotopologues.ch3oh13.cdms import catalog_file_to_overrides


def frequency_maps(root, xu_path):
    main = {}
    for row in json.loads(xu_path.read_text()):
        pair = frozenset((row['upper_key'], row['lower_key']))
        if pair in main:
            raise ValueError('duplicate accepted Xu observation endpoints')
        main[pair] = dict(row, frequency_policy='D24_raw_Xu_CALCULATION_identity_fitted_minimum_UNC')
    home = root/'isotopologues/13ch3oh'
    config13 = json.loads((home/'config/observed-frequencies.json').read_text())
    raw = catalog_file_to_overrides(home/config13['catalog'])
    isotope = {
        frozenset((x.identity.upper_key.encode(), x.identity.lower_key.encode())):
        dict(upper_key=x.identity.upper_key.encode(), lower_key=x.identity.lower_key.encode(),
             frequency_mhz=x.frequency_mhz, source=x.source_id,
             frequency_policy='original_negative_tag_CDMS_exact_identity')
        for x in raw.values()
    }
    maps = {'12': main, '13': isotope}
    for iso, config in [('12', root/'config/observed-frequencies.json'),
                        ('13', home/'config/observed-frequencies.json')]:
        for row in json.loads(config.read_text()).get('overrides', []):
            maps[iso][frozenset((row['upper_key'], row['lower_key']))] = dict(
                row, frequency_policy='explicit_owner_reviewed_override')
    for table in maps.values():
        for row in table.values():
            if not math.isfinite(float(row['frequency_mhz'])) or float(row['frequency_mhz']) <= 0:
                raise ValueError('invalid observed frequency')
    return maps

from TprParser.EdrReader import EdrReader
import numpy as np
import glob
import os

def read_xvg(fname:str) -> dict:
    data = {}
    keys = []
    with open(fname, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('#') or len(line.strip()) < 2:
                continue
            elif line.startswith('@ s') and 'legend' in line:
                keys.append(line.split('"')[1])
            else:
                if line.startswith('@'):
                    continue
                # includes Time column 
                values = list(map(float, line.strip().split()))
                assert len(keys)+1 == len(values)
                for k, v in zip(['Time']+keys, values):
                    if k not in data:
                        data[k] = []
                    data[k].append(v)
    return data


for edr in glob.glob('test_edr/*.edr'):
    print(f'do test {edr}', flush=True)
    xvg = edr.split('.edr')[0] + '.xvg'
    if not os.path.exists(xvg):
        raise FileNotFoundError(xvg)
    
    data = read_xvg(xvg)
    energies = EdrReader(edr).get_ene()
    for k in data:
        a = np.array(data[k])
        b = np.array(energies[k])
        # Constr. rmsd is not match to edr
        if k == 'Constr. rmsd':
            continue
        try:
            assert np.allclose(a, b, rtol=1E-3), f'{a} != {b} for {k} in {xvg}'
        except:
            assert k == 'Box-Vel-XX' or k == 'Box-Vel-YY' or k == 'Box-Vel-ZZ'
            assert np.allclose(a, b, rtol=1E-1), f'{a} != {b} for {k} in {xvg}'

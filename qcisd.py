'''
This is an example using cc-ade to perform a QCISD computation without DISS using
the HF data in a 6-31g basis (r=1.6 bohr). We are detecting diagrams which fall below 
the energy tolerence (this should be set a bit lower maybe), once detected the diagram 
is removed from the active ade list. Note the pruning, unlike CCSD, does not occur with
QCISD as these diagrams are ignored in the ansatz.
'''
import numpy as np
from cc_ade import expansion
from validation.validate import validate

cycle_limit = 50
data = np.load('../cc-ade/validation/hf.npz')
nvir, nocc = data['ts'].shape                               #get the occupied and virtual orbital numbers
o, v, n = slice(None, nocc), slice(nocc, None), np.newaxis  #make slices for occupied and virtual orbitals

#user_data are variable and values dictionary passed to eval
user_data = {}
user_data['gs'], user_data['fs'] = data['g'], data['f']     #spin eri and fock
user_data['ts'] = np.zeros_like(data['ts'])                 #singles amplitudes
user_data['td'] = np.zeros_like(data['td'])                 #doubles amplitudes, initialise to zero
user_data['o'], user_data['v'] = o, v                       #orbital slices

#fock denominator
eps = np.diag(user_data['fs'])
ds = -eps[v,n] + eps[n,o]
dd = -eps[v,n,n,n] - eps[n,v,n,n] + eps[n,n,o,n] + eps[n,n,n,o]

#mp2 initialisation
user_data['td'] = user_data['gs'][v,v,o,o] * np.reciprocal(dd)

#run parameters
max_cycle, tol, norm, converged, verbose = 50, 1e-8, 0, False, True

s = expansion(order='S', series='T_1 + T_2 + T_1 T_2')
s.order_by_reference()
vs = validate(s, 'hf')

d = expansion(order='D', series='1 + T_1 + T_2 + \\frac{1}{2} T_2 T_2')
d.order_by_reference()
vd = validate(d, 'hf')

e = expansion(order='E', series='1 + T_2')
ve = validate(e, 'hf')

cycle_energy = [0.0]
for cycle in range(cycle_limit):   

    #store pre-update amplitudes
    prune = []
    t1 = np.zeros_like(user_data['ts'])
    for diagram in s.ade:
        t =  vs.evaluate(diagram['reference'], user_data)

        if cycle >5 and np.linalg.norm(t) < tol: 
            prune.append(diagram['reference'])
            print('Removing ', diagram['reference'])

        t1 += t

    s.ade = [i for i in s.ade if not  i['reference']  in prune]

    prune = []
    t2 = np.zeros_like(user_data['td'])
    for diagram in d.ade:
        t = vd.evaluate(diagram['reference'], user_data)

        if cycle > 5 and np.linalg.norm(t) < tol: 

            prune.append(diagram['reference'])
            print('Removing ', diagram['reference'])
        
        t2 += t
 
    d.ade = [i for i in d.ade if not  i['reference']  in prune]

    user_data['ts'] = t1 * np.reciprocal(ds) + user_data['ts']
    user_data['td'] = t2 * np.reciprocal(dd) + user_data['td']

    energy = 0.0
    for diagram in e.ade:
        energy += ve.evaluate(diagram['reference'], user_data)

    #calculate current cycle energy
    cycle_energy.append(energy) 

    #test convergence
    delta_energy = np.abs(cycle_energy[-2] - cycle_energy[-1])
    if delta_energy < tol:
        converged = True
        break
    else:
        if verbose: print('cycle = {:>3d}  energy = {:>15.10f}   \u0394E = {:>12.10f} '.format(cycle, cycle_energy[-1], delta_energy))
        del cycle_energy[0]


if converged:
    print('final QCISD energy correction is ', cycle_energy[-1], 'au')


'''
cycle =   0  energy =   -0.1607374331   ΔE = 0.1607374331 
cycle =   1  energy =   -0.1704617873   ΔE = 0.0097243543 
cycle =   2  energy =   -0.1724974676   ΔE = 0.0020356802 
cycle =   3  energy =   -0.1745412891   ΔE = 0.0020438215 
cycle =   4  energy =   -0.1753135096   ΔE = 0.0007722206 
cycle =   5  energy =   -0.1758829631   ΔE = 0.0005694535 
cycle =   6  energy =   -0.1761493589   ΔE = 0.0002663957 
cycle =   7  energy =   -0.1763215299   ΔE = 0.0001721710 
cycle =   8  energy =   -0.1764106764   ΔE = 0.0000891465 
cycle =   9  energy =   -0.1764647146   ΔE = 0.0000540382 
cycle =  10  energy =   -0.1764941251   ΔE = 0.0000294105 
cycle =  11  energy =   -0.1765114081   ΔE = 0.0000172831 
cycle =  12  energy =   -0.1765210482   ΔE = 0.0000096401 
cycle =  13  energy =   -0.1765266294   ΔE = 0.0000055812 
cycle =  14  energy =   -0.1765297803   ΔE = 0.0000031509 
cycle =  15  energy =   -0.1765315916   ΔE = 0.0000018113 
cycle =  16  energy =   -0.1765326203   ΔE = 0.0000010287 
cycle =  17  energy =   -0.1765332097   ΔE = 0.0000005894 
cycle =  18  energy =   -0.1765335454   ΔE = 0.0000003357 
cycle =  19  energy =   -0.1765337375   ΔE = 0.0000001920 
cycle =  20  energy =   -0.1765338470   ΔE = 0.0000001096 
cycle =  21  energy =   -0.1765339096   ΔE = 0.0000000626 
cycle =  22  energy =   -0.1765339454   ΔE = 0.0000000357 
cycle =  23  energy =   -0.1765339658   ΔE = 0.0000000204 
cycle =  24  energy =   -0.1765339775   ΔE = 0.0000000117 
final QCISD energy correction is  -0.176533984148932 au
'''
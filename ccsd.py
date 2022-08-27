'''
This is an example using cc-ade to perform a CCSD computation without DISS using
the HF data in a 6-31g basis (r=1.6 bohr). We are detecting diagrams which fall below 
the energy tolerence (this should be set a bit lower maybe), once detected the diagram 
is removed from the active ade list. Note the pruning of diagrams doesn't start for a
few cycles to give things a chance to settle down. At cycle 5 S_{5a}, D_{5a} and {D_5b} 
are rejected at 1e-8 tolerance. At 1e-6 S_{1} and S_{2a} fall below tolerance too.
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

s = expansion('SD', 'S')
s.order_by_reference()
vs = validate(s, 'hf')

d = expansion('SD', 'D')
d.order_by_reference()
vd = validate(d, 'hf')

e = expansion('SD', 'E')
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
    print('final CCSD energy correction is ', cycle_energy[-1], 'au')


'''
cycle =   0  energy =   -0.1601976203   ΔE = 0.1601976203 
cycle =   1  energy =   -0.1697368248   ΔE = 0.0095392045 
cycle =   2  energy =   -0.1716299460   ΔE = 0.0018931212 
cycle =   3  energy =   -0.1737029874   ΔE = 0.0020730414 
cycle =   4  energy =   -0.1744645692   ΔE = 0.0007615818 
cycle =   5  energy =   -0.1750628911   ΔE = 0.0005983219 
Removing  S_{5a}
Removing  D_{5a}
Removing  D_{5b}
cycle =   6  energy =   -0.1753395414   ΔE = 0.0002766502 
cycle =   7  energy =   -0.1755271521   ΔE = 0.0001876107 
cycle =   8  energy =   -0.1756242736   ΔE = 0.0000971216 
cycle =   9  energy =   -0.1756854184   ΔE = 0.0000611448 
cycle =  10  energy =   -0.1757189593   ΔE = 0.0000335409 
cycle =  11  energy =   -0.1757392933   ΔE = 0.0000203340 
cycle =  12  energy =   -0.1757507854   ΔE = 0.0000114921 
cycle =  13  energy =   -0.1757576205   ΔE = 0.0000068351 
cycle =  14  energy =   -0.1757615434   ΔE = 0.0000039229 
cycle =  15  energy =   -0.1757638542   ΔE = 0.0000023108 
cycle =  16  energy =   -0.1757651909   ΔE = 0.0000013368 
cycle =  17  energy =   -0.1757659745   ΔE = 0.0000007836 
cycle =  18  energy =   -0.1757664297   ΔE = 0.0000004552 
cycle =  19  energy =   -0.1757666959   ΔE = 0.0000002662 
cycle =  20  energy =   -0.1757668508   ΔE = 0.0000001549 
cycle =  21  energy =   -0.1757669413   ΔE = 0.0000000905 
cycle =  22  energy =   -0.1757669940   ΔE = 0.0000000527 
cycle =  23  energy =   -0.1757670248   ΔE = 0.0000000308 
cycle =  24  energy =   -0.1757670427   ΔE = 0.0000000179 
cycle =  25  energy =   -0.1757670532   ΔE = 0.0000000105 
final CCSD energy correction is  -0.1757670593231534 au
'''
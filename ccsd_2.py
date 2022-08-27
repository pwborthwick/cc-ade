'''
This is an example using cc-ade to perform a CCSD computation without DIIS using
the HF data in a 6-31g basis (r=1.6 bohr).The target level of this compuation is
CCSDT-2 which involves adding diagram S_7 to the singles, D_10a-10c, D_11a-11c 
(T_1 T_3) to the doubles and constructing a triples from T_1a, T_1b, (T_2), T_2a 
and T_2b (T_3) and (T_2^2) T_3a-3e diagrams 
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
user_data['tt'] = np.zeros_like(data['tt'])                 #will need triples to
user_data['o'], user_data['v'] = o, v                       #orbital slices

#fock denominator
eps = np.diag(user_data['fs'])
ds = -eps[v,n] + eps[n,o]
dd = -eps[v,n,n,n] - eps[n,v,n,n] + eps[n,n,o,n] + eps[n,n,n,o]
dt = -eps[v,n,n,n,n,n] - eps[n,v,n,n,n,n] - eps[n,n,v,n,n,n] + eps[n,n,n,o,n,n] + eps[n,n,n,n,o,n] + eps[n,n,n,n,n,o]

#mp2 initialisation
user_data['td'] = user_data['gs'][v,v,o,o] * np.reciprocal(dd)

#run parameters
max_cycle, tol, norm, converged, verbose = 50, 1e-8, 0, False, True

s = expansion('SD', 'S')
s.order_by_reference()

x = expansion(order='S', series='T_3')                       #get S_7 diagram
s.ade += x.ade                                               #append diagram to singles

d = expansion('SD', 'D')
d.order_by_reference()

x = expansion(order='D', series='T_3 + T_1 T_3')              #get D_10a, b and c and D_11a, b, and c
d.ade += x.ade                                                #append D_10a-c and D_11a-c to doubles

t = expansion(order='T', series='T_2 + T_3')                  #get T_1a-b and T_2a-b
t.ade = t.ade[:4]                                             #omit T2c,d,e
x = expansion(order='T', series='T_2 T_2')                    #get T_3a-e 
t.ade += x.ade                                                #append T_2 T_2 diagrams              

vs = validate(s, 'hf')
vd = validate(d, 'hf')
vt = validate(t, 'hf')

e = expansion('SD', 'E')
ve = validate(e, 'hf')

cycle_energy = [0.0]
for cycle in range(cycle_limit):   

    t1 = np.zeros_like(user_data['ts'])
    for diagram in s.ade:
        x =  vs.evaluate(diagram['reference'], user_data)
        t1 += x

    t2 = np.zeros_like(user_data['td'])
    for diagram in d.ade:
        x = vd.evaluate(diagram['reference'], user_data)
        t2 += x

    t3 = np.zeros_like(user_data['tt'])
    for diagram in t.ade:
        x = vt.evaluate(diagram['reference'], user_data)
        t3 += x


    user_data['ts'] = t1 * np.reciprocal(ds) + user_data['ts']
    user_data['td'] = t2 * np.reciprocal(dd) + user_data['td']
    user_data['tt'] = t3 * np.reciprocal(dt) + user_data['tt']

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
cycle =   1  energy =   -0.1723669169   ΔE = 0.0121692966 
cycle =   2  energy =   -0.1734598785   ΔE = 0.0010929616 
cycle =   3  energy =   -0.1764441376   ΔE = 0.0029842592 
cycle =   4  energy =   -0.1770510657   ΔE = 0.0006069280 
cycle =   5  energy =   -0.1779638562   ΔE = 0.0009127905 
cycle =   6  energy =   -0.1782309956   ΔE = 0.0002671394 
cycle =   7  energy =   -0.1785275836   ΔE = 0.0002965880 
cycle =   8  energy =   -0.1786364072   ΔE = 0.0001088236 
cycle =   9  energy =   -0.1787359027   ΔE = 0.0000994955 
cycle =  10  energy =   -0.1787783867   ΔE = 0.0000424840 
cycle =  11  energy =   -0.1788124914   ΔE = 0.0000341048 
cycle =  12  energy =   -0.1788286547   ΔE = 0.0000161633 
cycle =  13  energy =   -0.1788405276   ΔE = 0.0000118729 
cycle =  14  energy =   -0.1788465770   ΔE = 0.0000060494 
cycle =  15  energy =   -0.1788507570   ΔE = 0.0000041800 
cycle =  16  energy =   -0.1788529969   ΔE = 0.0000022399 
cycle =  17  energy =   -0.1788544805   ΔE = 0.0000014836 
cycle =  18  energy =   -0.1788553039   ΔE = 0.0000008234 
cycle =  19  energy =   -0.1788558335   ΔE = 0.0000005296 
cycle =  20  energy =   -0.1788561347   ΔE = 0.0000003012 
cycle =  21  energy =   -0.1788563246   ΔE = 0.0000001898 
cycle =  22  energy =   -0.1788564344   ΔE = 0.0000001098 
cycle =  23  energy =   -0.1788565026   ΔE = 0.0000000682 
cycle =  24  energy =   -0.1788565426   ΔE = 0.0000000400 
cycle =  25  energy =   -0.1788565672   ΔE = 0.0000000246 
cycle =  26  energy =   -0.1788565817   ΔE = 0.0000000145 
final CCSD energy correction is  -0.17885659054031286 au
'''
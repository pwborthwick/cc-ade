import numpy as np
from itertools import permutations, product
import os

class validate(object):
    #class to test coupled-cluster diagrams against Shavitt & Bartlett

    def __init__(self, ade):
        
        current_directory = os.path.dirname(__file__) + '/'

        self.amplitude_file = current_directory + '../validation/hf.npz'
        self.reference_file = current_directory + 'reference.txt'
        self.r = current_directory + 'hf-eom.npz'
        self.has_r = True

        self.occ = 10
        self.ade = ade
        self.amplitudes()
        self.load()

    def load(self):
        #load Shavitt and Bartlett reference data

        self.reference_data = {}

        with open(self.reference_file, 'r') as f:
            for line in f:
                reference_list = line.rstrip('\n').split('&')
                self.reference_data[reference_list[0]] = reference_list[1:]

    def get(self, id, type='sb'):
        #return a value

        if type == 'sb':
            return self.reference_data.get(id, ['no key'])
        if type == 'ade':
            data =  list(filter(lambda id_ref: id_ref['reference'] == id, self.ade.ade))

            #decode permutation information
            p = [i.replace('(','').replace(')','').replace(' ','') for i in data[0]['permutation'].split('\\hat{P}')[1:]][::-1]

            if p == []: p = ['', '']
            if not 'i' in p[0] and p != ['','']: p = [''] + p
            if not 'a' in p[-1] and p != ['','']: p = p + ['']

            return [id, data[0]['code'], p[0], p[1]]

    def _evaluate(self, id, type='sb'):
        #evaluate an diagram with id

        order = {'S':'ai', 'D':'abij', 'T':'abcijk', 'Q':'abcdijkl'}
        level = [i for i in order.keys() if i in id][0]

        #check we have amplitudes in data
        if (level == 'Q') and (self.amplitude_file == 'hf.npz'): 
            print('HF data set has SDT only')
            return 0

        diagram = self.get(id, type)
        if diagram == ['no key']: return -1

        code, i, a = diagram[1:]

        run = 't = ' + code 

        if i.strip() != '':
            run += '\n' + 't = permute(\'' + order[level] + '->' + i + '\', t=t)'
        if a.strip() != '':
            run += '\n' + 't = permute(\'' + order[level] + '->' + a + '\', t=t)'

        #remove explicit numpy reference from ade code
        run = run.replace('np.', '')

        data = {'einsum': np.einsum, 'fs': self.fs, 'gs': self.gs,
                'ts': self.ts, 'td': self.td, 'tt': self.tt,
                'o': slice(None, self.occ), 'v':slice(self.occ, None), 'permute': self.permute }

        if self.has_r: data['rs'], data['rd'] = self.rs, self.rd

        exec(run, {'builtins':None}, data)
        
        return data['t']

    def compare(self, id):

        if id == 'all':
            for data in self.ade.ade:
                if data['labels'] == ('',''): continue

                t = self._evaluate(data['reference'])
                s = self._evaluate(data['reference'], type='ade')
                print('{:8s}  {:12s} {:<18.15f}  {:<18.15f} '.format(data['reference'], data['contraction'], np.linalg.norm(t), np.linalg.norm(s)), end='')
                print(np.isclose(np.linalg.norm(t), np.linalg.norm(s)), np.allclose(t,s))

        elif '_{' in id:
            for data in self.ade.ade:
                if data['reference'] == id: break
            if data['reference'] != id: 
                print('reference ', id, ' not found')
                return
            t = self._evaluate(data['reference'])
            s = self._evaluate(data['reference'], type='ade')
            print('{:8s}  {:12s} {:<18.15f}  {:<18.15f} '.format(data['reference'], data['contraction'], np.linalg.norm(t), np.linalg.norm(s)), end='')
            print(np.isclose(np.linalg.norm(t), np.linalg.norm(s)), np.allclose(t,s))
        else:
            for data in self.ade.ade:
                data_id = data['id']
                hamiltonian = data_id[data_id.find('(')+1:data_id.find(')')]
                if hamiltonian == id:
                    t = self._evaluate(data['reference'])
                    s = self._evaluate(data['reference'], type='ade')
                    print('{:8s}  {:12s} {:<18.15f}  {:<18.15f} '.format(data['reference'], data['contraction'], np.linalg.norm(t), np.linalg.norm(s)), end='')
                    print(np.isclose(np.linalg.norm(t), np.linalg.norm(s)), np.allclose(t,s))

    def evaluate(self, id, data):
        #evaluate diagram using user supplied values

        order = {'E':'','S':'ai', 'D':'abij', 'T':'abcijk', 'Q':'abcdijkl'}
        level = [i for i in order.keys() if i in id][0]

        diagram = self.get(id, 'ade')
        if diagram == ['no key']: return -1

        code, i, a = diagram[1:]

        run = 't = ' + code 
        if i.strip() != '':
            run += '\n' + 't = permute(\'' + order[level] + '->' + i + '\', t=t)'
        if a.strip() != '':
            run += '\n' + 't = permute(\'' + order[level] + '->' + a + '\', t=t)'

        #remove explicit numpy reference from ade code
        run = run.replace('np.', '')

        data['einsum'] = np.einsum; data['permute'] = self.permute

        exec(run, {'builtins':None}, data)
        
        return data['t']


    def amplitudes(self):
        #get saved amplitudes for HF and C

        data = np.load(self.amplitude_file)

        self.ts, self.td, self.tt = data['ts'], data['td'], data['tt']
        self.fs, self.gs = data['f'], data['g']
        if self.has_r: 
            data = np.load(self.r)
            self.rs, self.rd = data['rs'].transpose(1,0), data['r2'].transpose(2,3,0,1)

    def permute(self, p, t='raw'):
        #do permutation p of t - if t not given return raw permutation lists

        def cluster_permutations(p):

            def permutation_parity(p):
                return 1 - 2*(sum(1 for (x,px) in enumerate(p)
                                    for (y,py) in enumerate(p)
                                    if x<y and px>py) % 2)

            def exchange(p, a, b):

                p[b], p[a] = p[a], p[b]

                return p

            #split source indices and required permutation
            source = p[:p.index('-')]
            target = p[p.index('>')+1:]

            #check all target indices exist in source
            valid = not -1 in [source.find(i) for i in target.replace('/', '')]
            if not valid:
                print('alien target index detected')
                return [0], [0]

            #get order of source, index set permutation acting on and permutation type
            order = len(source)//2
            if order > 4:
                print('only coded for orders 2,3 and 4')
                return [0], [0]

            type = 'v' if (target[0] in source[:order]) else '^'
            level = target.count('/')

            #source indices and source index sets
            vector = list(np.arange(2*order))
            lower_id = vector[:order] ; upper_id = vector[order:]

            if level == 0:

                #ensure lexical order of target permutation
                target = ''.join(sorted(target))
                
                #get indices of target labels
                id = [source.index(i) for i in target]

                #permute label indices
                id_perm = list(permutations(id))

                #generate biased permutation list
                cycle = [lower_id + list(p) for p in id_perm] if type == '^' else [list(p) + upper_id for p in id_perm]

                parity = [permutation_parity(p) for p in id_perm]

            if level == 1:

                #ensure lexical order of target
                target = target.split('/')
                target = [sorted(i) for i in target]
                target = sorted(target, key=len)

                assert sum([len(i) for i in target]) == order

                if len(target[0]) == 1:

                    j = source.index(target[0][0])

                    cycle = [exchange(vector[:], j, source.index(i)) for i in target[1]]
                    cycle.insert(0, vector)
                    parity = [permutation_parity(p) for p in cycle]

                else:

                    j = [source.index(i) for i in target[0]]
                    k = [source.index(i) for i in target[1]]

                    cycle = [exchange(vector[:], i[0], i[1]) for i in list(product(j, k))]
                    cycle.insert(0, vector)

                    k = [source.index(i[0]) for i in target[0] + target[1]]
                    cycle.append(exchange(exchange(vector[:], k[0], k[2]), k[1], k[3]))
                    parity = [permutation_parity(p) for p in cycle]

            if (level == 2) and (order == 4):

                #ensure lexical order of target
                target = target.split('/')
                target = [sorted(i) for i in target]
                target = sorted(target, key=len)

                assert sum([len(i) for i in target]) == order

                j, m, n = [source.index(i) for i in target[0]],  [source.index(i) for i in target[1]],  [source.index(i) for i in target[2]]
                cycle = [exchange(vector[:], j[0], i) for i in n] + [exchange(vector[:], m[0], i) for i in n] + [exchange(vector[:], j[0], m[0])]

                cycle += [exchange(exchange(vector[:], m[0], n[1]), j[0], n[0])] + [exchange(exchange(vector[:], m[0], n[0]), j[0], n[1])] 
                cycle += [exchange(i,j[0], m[0]) for i in [exchange(vector[:], j[0], i) for i in n] + [exchange(vector[:], m[0], i) for i in n]]
                cycle.insert(0, vector)     

                parity = [permutation_parity(p) for p in cycle]

            return cycle, parity

        if type(t) == str : return cluster_permutations(p)

        cycle, parity = cluster_permutations(p)
        f = 0
        for i,j  in zip(cycle, parity):
            f += t.transpose(i) * j

        return f

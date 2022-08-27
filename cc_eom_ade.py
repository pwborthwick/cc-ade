from __future__ import division
from sympy import *
from sympy.polys.ring_series import rs_exp
from itertools import permutations
import re
import lablib

class expansion(object):
    #exponential expansion of amplitude operators

    def __init__(self, level='', order='', contraction=''):
        #class initiator 
        #- level (str) the operator types to include 'SD' will include T_1 and T_2
        #- order (int) the maximium degree of powers of operators to include

        self.vertex = {
        'pp': ['+'], 'hh':['-'],'ph':['0'], 'hp':['+', '-'],
        'hhhh':['+', '+'], 'pppp':['-', '-'],'phph':['+', '-'], 'ppph':['+'], 'phpp':['+', '+', '-'],
        'hphh':['-'], 'hhhp':['+', '-', '-'], 'pphh':['0'], 'hhpp':['+', '+', '-', '-']
        }

        self.interaction = {
        'pp':[0],'hh':[0],'ph':[+1],'hp':[-1],'hhhh':[0],'pppp':[0],'phph':[0],'ppph':[+1],
        'phpp':[-1],'hphh':[+1],'hhhp':[-1],'pphh':[+2],'hhpp':[-2]
        }

        #must be either expansion or individual contraction
        assert ((level + order) != '') or (contraction != '')

        if (level + order) != '':
            self.level = level ; self.order = {'S':2, 'D':3, 'T':4, 'Q':5}[order]

            #get the operator set and 
            self.construct()

            #reduce to required order
            self.to_cluster_order()

            #add vertex to operator definition
            self.operator_vertex()

            #make a dictionary of diagram types
            self.diagram_dictionary(order)
        if contraction != '':
            self.contract(contraction)

        diagram_object = diagrams(self)
        self.ade = diagram_object.allocate_labels()

    def construct(self):
        #use Sympy to construct exponential series

        #define symbols for amplitude operators
        R, t1, t2, t3, t4 = ring('T_1, T_2, T_3, T_4', QQ)
        t = [t1, t2, t3, t4]

        #construct raw expansion (Sympy internal)
        expansion_sympy = Rational(1)
        for i in self.level:
            k = 'SDTQ'.find(i)
            expansion_sympy *= rs_exp(t[k], t[k], self.order+1)

        #convert to Python readable form
        raw = srepr(expansion_sympy)

        #parse raw into operator and series blocks
        series_start = raw.find('[')
        operator_block  = raw[:series_start]
        series_block = raw[series_start+1:-2]

        def parse(block, p):
            #parse a sympy block

            end = 0  ; symbols = []
            while True:

                start = block.find(p[0], end)
                if start == -1: break
                start += p[1]
                end = block.find(p[2], start) + p[3]
                symbols.append(block[start:end].replace('PythonRational',''))

            return symbols

        #get operators defined
        self.operators = parse(operator_block, ['Symbol', 8, ')', -1])

        #get term expresions
        self.terms = []
        raw_terms = parse(series_block, ['((', 1, '))', 1])
        for term in raw_terms:
            operator = [int(i) for i in term[1:term.find('),')].split(',')] 
            name = ''.join([str(i) for i in operator]) 
            factor    = [int(i) for i in term[term.find('),')+4:-1].split(',')]

            self.terms.append([name, factor, operator])

        self.len = len(self.terms)

    def to_cluster_order(self):
        #reduce order - Sympy order n is terms containing upto T^(n-1) eg T_1"(n-1) T_2^(n-1)
        # cluster order n means T_m^n contributes m*n so T_1^2 T_2 is order 2+2

        #get valid terms for required order
        valid = [True] * self.len
        for n, term in enumerate(self.terms):

            #running sum of orders over operators in term
            order = sum([i*(m+1) for m, i in enumerate(term[2])])
            n_amp = sum([i for i in term[2]])

            #valid if calculated order is less required order
            valid[n] = (order < self.order + 1)

            #check term has possible contractions
            required_order = self.order - 2 - order
            if valid[n]: valid[n] = self.get('i', required_order) != []
            if valid[n]: (max([len(i) for i in self.get('v', self.get('i', required_order))]) >= n_amp)

            #for connected diagrams no more than 4 terms
            if valid[n]: valid[n] = sum(term[2]) <= 4

        #remove invalid terms
        for n, term in enumerate(self.terms[:]):
            if not valid[n]: self.terms.remove(term)

        #add R to terms
        terms = []
        for series in range(len(self.level)):
            for term in self.terms:

                id = term[0]
                order = sum([i*(m+1) for m, i in enumerate(term[2])])
                if order < self.order - series + 1:
                    terms.append([id + str(series+1)] + term[1:])

        self.terms = terms
        self.len = len(self.terms)

    def operator_vertex(self):
        #append the amplitude operator vertex type 

        enhanced = []
        vertices = [['+', '-'], ['+', '+', '-', '-'], ['+','+','+',  '-', '-', '-'],
                    ['+','+','+','+','-','-','-','-']]
        for n, op in enumerate(self.operators):
            enhanced.append([op, vertices[n]])

        self.operators = enhanced.copy()

    def do_latex(self):
        #construct a Latex string of expansion

        latex = '$'
        for term in self.terms:

            latex += '\\frac{' + str(term[1][0]) + '}{' + str(term[1][1]) + '}'
            for i in range(len(self.operators)):
                if term[2][i] != 0: latex += 'T_' + str(i+1)
                if term[2][i] > 1: latex += '^{' + str(term[2][i]) + '}'

            latex += ' R_{' + str(term[0][4]) + '}'
            latex += ' + '

        latex = latex.replace('\\frac{1}{1}','')[:-2] + '$'

        return latex

    def do_terms(self):
        #evaluate terms for a specific level

        for n_term, term in enumerate(self.terms):

            #vertices   
            vertex = [] ; t = []
            for n, i in enumerate(term[2]):

                if i != 0: 
                    t.append(self.operators[n][1])
                    for j in range(i-1):
                        t.append(self.operators[n][1])

            vertex = t
            vertex.append([['+','-'], ['+','+','-','-'],['+','+','+','-','-','-'],['+','+','+','+','-','-','-','-']][int(term[0][4])-1])
            if vertex == []: vertex = [['0']]

            #interaction number
            l = {'E':0, 'S':1, 'D':2, 'T':3, 'Q':4}[self.type]
            interaction_number = 0

            for v in vertex:
                interaction_number += v.count('+')
            
            interaction_number = l - interaction_number 

            #append vertex and interaction number and R type to terms properties
            self.terms[n_term] =  self.terms[n_term][:3] + [vertex]
            self.terms[n_term] += [interaction_number]

    def do_term(self, id):
        #print information about a term

        for term in self.terms:
            if term[0] == id: break

        print('Level is [', self.level, ']  Order is [', self.order, ']')
        print('Term id  [', term[0], ']   Multiplier is [ ', str(term[1][0]), '/', str(term[1][1]), ' ]')

        latex = ''
        if term[1] != [1,1]:
            if term[1][1] != 1:
                latex += '\\frac{' + str(term[1][0]) + '}{' + str(term[1][1]) + '} '
            else:
                latex += str(term[1][0])

        for n, i in enumerate(term[2]):
            if i == 0: continue
            if i == 1: 
                latex += 'T_' + str(n+1) + ' '
            else :
                latex += 'T_' + str(n+1) + '^{' + str(i) + '} '
            latex += ' R_' + str(term[0][4])

        print('Latex string is [ ', latex,' ]')

        vertex_string = ' | '.join([ ''.join(i) for i in self.get('t', id)[3]])
        print('Vertex string is [ ', vertex_string, ' ]')

        s = '+' if term[4] >= 0 else '-'
        print('Interaction number is [ ', s, abs(term[4]), ' ]')

        diagrams = self.diagrams[id]

        print('\nContractor               Contractions')
        for key in diagrams.keys():
            if diagrams[key] == []: continue
            print('{:<5}                       '.format(key), diagrams[key])

    def contractions(self, id, seed):

        #null vertex
        if id[:4] == '0000' and seed == 0: 
            return seed

        #validate connection counts
        validate = lambda p, m: (p.count('+') < (m+2)) and (p.count('-') < (m+2))
        
        #block into different operator types
        blocks = seed.count('|') + 1

        #operator types 0=T_1 1=T_2 etc
        block_types = [0]*int(id[0]) + [1]*int(id[1]) + [2]*int(id[2]) + [3]*int(id[3]) + [4]*int(id[4])    

        #slice information for each type of operator in p_string
        block_count = [block_types.count(0),
                       block_types.count(0)+block_types.count(1),
                       block_types.count(0)+block_types.count(1)+block_types.count(2),
                       block_types.count(0)+block_types.count(1)+block_types.count(2)+block_types.count(3),
                       block_types.count(0)+block_types.count(1)+block_types.count(2)+block_types.count(3)+block_types.count(4)]

        #raw permutations of seed string
        p_raw = permutations(seed)

        #process raw into ordered and validated list
        p_list = []
        for n, p in enumerate(p_raw):

            #check seperator allocation
            p = ''.join(p)
            if ('||' in p) or p.startswith('|') or p.endswith('|'): continue

            #split into individual operator blocks
            block = p.split('|')
            ordered = []

            #validate each block and order (with block)
            valid = True

            for m, b in enumerate(block):

                #validate each operator in block
                if valid: 
                    k = block_types[m] if block_types[m] != 4 else (int(id[4]) - 1)
                    valid = validate(b, k)

                #sort b and add to the ordered operator list
                q = ''.join(sorted(b))
                ordered.append(q)            
            
            if valid:
                #sort within each operator type block
                ordered = (sorted(ordered[:block_count[0]]) + 
                           sorted(ordered[block_count[0]:block_count[1]]) + 
                           sorted(ordered[block_count[1]:block_count[2]]) +
                           sorted(ordered[block_count[2]:block_count[3]]) +
                           sorted(ordered[block_count[3]:block_count[4]])
                           ) 

                #concatenate operators
                p_list.append('|'.join(ordered))

        #get unique entries
        p_list = list(set(p_list))

        return p_list


    def diagram_dictionary(self, type):
        #make a dictionary out of the diagrams for the terms

        self.diagrams = {}  

        #define variable type
        self.type = type

        #enhance root terms lists
        self.do_terms()

        for term in self.terms:

            id = term[0] ;  interaction_number = term[4]
            interactions = self.get('i', interaction_number)
            vertices = self.get('v', interactions)

            term_diagram = {}
            for n, vertex  in enumerate(vertices):

                #for EOM vertex 0 cannot connect as no 1 vertex
                if vertex == ['0']: 
                    continue

                seed = ''.join(vertex)  + '|'*(len(term[3])-1)

                term_diagram[interactions[n]] = self.contractions(id, seed) 

            self.diagrams[id] = term_diagram


    def get(self, code, source, id = ''):
        #get all Hamiltonians with interaction requested

        if code == 't':
            for term in self.terms:
                if term[0] == source: return term

            return []

        #given an interaction number get a Hamiltonians ids with that number
        if code == 'i':
            interactions = []

            for k in self.interaction.keys():

                if self.interaction[k][0] == source:

                    #connected diagrams only
                    if id != '':

                        term_vertices =  len(self.get('t', id)[3])

                        if len(self.vertex[k]) >= term_vertices: 
                            interactions.append(k)

                            #like vertices - term vertex set
                            vertices = []
                            for i in self.get('t', id)[3]:
                                vertices += i

                            interaction_counts_term    = [vertices.count(s) for s in ['0', '+', '-']]
                            interaction_counts_element = [self.vertex[k].count(s) for s in ['0', '+', '-']]

                            #must be enough vertices in term to accomodate Hamiltonian element
                            valid = True
                            for i in range(3):
                                if valid: valid = interaction_counts_element[i] <= interaction_counts_term[i]

                            if not valid: interactions.pop()
                    else:
                        interactions.append(k)

            return interactions

        #given a list of Hamiltonian ids get list of vertices
        if code == 'v':
            vertices = []
            for id in source:
                vertices.append(self.vertex[id])

            return vertices

        #return from diagram list Shavitt and Bartlett reference record
        if code == 'r':
            valid = True in [i['reference'] == source for i in self.ade]
            value = [i for i in self.ade if i['reference'] == source][0] if valid else {}

            return value

        #return from diagram list s-string record 'Hamiltonian':term eg '+|--:+-|++--'
        if code == 'c':
            h, t = source.split(':')
            amplitudes = t.split('|')
            counts = [str(amplitudes.count(t)) for t in ['+-', '++--', '+++---', '++++----']]
            id = ''.join(counts)  

            #get diagram list entry
            selected = [i for i in self.ade if (i['contraction'] == h) and (i['id'][:4] == id)]

            value = selected[0] if selected != [] else {}

            return value

    def contract(self, contraction):
        #evaluate a single diagram

        #construct self.terms record
        self.terms = []

        h, amplitudes = contraction.split(':')
        amplitude_list = amplitudes.split('|')

        counts = [amplitude_list.count(t) for t in ['+-', '++--', '+++---', '++++----']]

        id = ''.join([str(i) for i in counts])

        t = [re.findall('[+-]', i) for i in  amplitude_list]

        hamiltonian = [k for k in  self.vertex if self.vertex[k] ==  re.findall('[+-]', h.replace('|', ''))][0]
        interaction_number_hamiltonian = self.interaction[hamiltonian][0]

        self.terms.append([id, [1, 1], counts, t, interaction_number_hamiltonian])
        
        #construct self.diagram record
        self.diagrams ={}
        self.diagrams[id] = {hamiltonian:[h]}

        #get level of cluster
        self.order = 2 + amplitudes.count('+') + interaction_number_hamiltonian
        self.type = ['E', 'S', 'D', 'T', 'Q'][self.order-2]

    def order_by_reference(self):
        #impose an order on the diagram definition list

        _nsre = re.compile('([0-9]+)')
        def convert(x):
            return [int(text) if text.isdigit() else text.lower() for text in re.split(_nsre, x['reference'])]

        self.ade.sort(key=convert)

class labels(object):
    #class to manage label assignment

    def __init__(self, type):

        self.type = type
        label_set = {'E':2, 'S':2, 'D':2, 'T':3, 'Q':4}[type]

        cache = {'-': ['i','j','k','l','m','n'], '+': ['a','b','c','d','e','f']}
        self.cache = {'+':{'e': cache['+'][:label_set], 'i':cache['+'][max(label_set, 2):]},
                      '-':{'e': cache['-'][:label_set], 'i':cache['-'][max(label_set, 2):]}}

    def refresh(self):
        #reinstate the cache 

        self.__init__(self.type)

    def get(self, _with, _from, _at=1, nopop=False):
        #get from _with ('p'/'h') _from ('e'/'i') _at position

        if self.cache[_with][_from] == []: return ''

        _chr = ''
        _chr += self.cache[_with][_from][_at-1]
        if not nopop:
            del self.cache[_with][_from][_at-1]

        return _chr

    def replace_all(self, _in, _type):
        #replace all external signs

        for i in _in[:]:
            if (i in _type) and (i == '+'): _in = _in.replace('+', self.get('+', 'e'), 1)
            if (i in _type) and (i == '-'): _in = _in.replace('-', self.get('-', 'e'), 1)

        return _in

    def replace_near(self, op, amp, labels):
        #replace nearest

        amp_list, op_list =  amp, op.split('|')

        a = amp_list[0][::-1].replace(op_list[0], labels[0], 1)[::-1]
        b = amp_list[1].replace(op_list[1], labels[1], 1)

        return a + '|' + b        

class diagrams(object):
    #class to evaluate diagrams

    def __init__(self, expansion):

        self.expansion = expansion
        self.labels = labels(self.expansion.type)

    def diagram_sign(self, diagram_string):
        #return the sign associated with the diagram string: Shavitt & Bartlett Rule 8

        expansion_order =  {'E':0, 'S':1, 'D':2, 'T':3, 'Q':4}[self.expansion.type]
        #down (-) lines
        h = diagram_string.count('-') + expansion_order

        #loops 
        l = expansion_order
        for t in diagram_string.split('|'):
            l += t.count('+-')

        return (-1) ** (h + l), [h, l]

    def equivalent_lines(self, diagram_string):
        #return the number of equivalent lines: Shavitt & Bartlett Rule 6

        equivalent = 0
        for t in diagram_string.split('|'):
            equivalent += t.count('++') + t.count('--')

        return equivalent

    def equivalent_vertices(self, diagram_string, amplitude_string):
        #return the number of equivalent vertices: Shavitt & Bartlett Rule 7

        #if more than one amplitude remve last one so R can't match with a T
        if diagram_string.count('|') != 0:
            diagram_string   = '|'.join(diagram_string.split('|')[:-1])
            amplitude_string = '|'.join(amplitude_string.split('|')[:-1])


        amplitude_type = [len(i)//2 for i in amplitude_string.split('|')]


        counts_p , counts_h, counts_ph  = [0]*max(amplitude_type), [0]*max(amplitude_type), [0]*max(amplitude_type)
        p =  [i == '+'  for i in diagram_string.split('|')]
        h =  [i == '-'  for i in diagram_string.split('|')]
        ph = [i == '+-' for i in diagram_string.split('|')]

        for n, b in enumerate(p):
            if b: counts_p[amplitude_type[n]-1] += 1
        for n, b in enumerate(h):
            if b: counts_h[amplitude_type[n]-1] += 1
        for n, b in enumerate(ph):
            if b: counts_ph[amplitude_type[n]-1] += 1

        equivalent = counts_p.count(2) + counts_h.count(2) + counts_ph.count(2)

        return equivalent

    def diagram_factor(self, diagram_string, amplitude_string):
        #overall factor

        return self.equivalent_lines(diagram_string) + self.equivalent_vertices(diagram_string, '|'.join(amplitude_string))

    def diagram_permutations(self, diagram_labels):
        #determine the permutations factor

        #reduce labels to just external labels
        labels = diagram_labels[0] 
        for label in diagram_labels[1].split('|'):
            labels += ',' + label
        labels = labels.replace('f','').replace('n','')
        if self.expansion.type in ['Q','T', 'D', 'S']: labels = labels.replace('e','').replace('m','')
        if self.expansion.type in ['T','D','S']: labels = labels.replace('d','').replace('l','')
        if self.expansion.type in ['D','S']: labels = labels.replace('c','').replace('k','')

        #copy of labels for particle and hole lines
        p = labels[:].replace('i','').replace('j','').replace('k','').replace('l','').split(',')
        h = labels[:].replace('a','').replace('b','').replace('c','').replace('d','').split(',')

        #remove null strings and sort 
        p = [i for i in p if i != '']
        h = [i for i in h if i != '']
        
        p.sort() ; h.sort()
        permutation = ''
        if self.expansion.type in ['S', 'D']:
            if len(p) > 1: permutation += '\\hat{P}(' + ''.join(p) + ') '
            if len(h) > 1: permutation += '\\hat{P}(' + ''.join(h) + ') '
            return permutation

        elif self.expansion.type in ['T', 'Q']:
            p.sort(key=lambda l: len(l))
            h.sort(key=lambda l: len(l))

            #tidy up a/b/c -> abc
            perm_p, perm_h = '', ''
            if len(p) > 1: 
                perm_p = '\\hat{P}(' + '/'.join(p) + ') '
                if all([len(i) == 1 for i in p]): perm_p = perm_p.replace('/','')
            if len(h) > 1: 
                perm_h += '\\hat{P}(' + '/'.join(h) + ') '
                if all([len(i) == 1 for i in h]): perm_h = perm_h.replace('/','')

            return perm_p + perm_h

        return  ''

    def diagram_sigma(self, diagram_labels):
        #determine indices for summation

        #convert to string and get unique set of labels
        diagram_labels = ''.join(diagram_labels)
        labels = set(diagram_labels)

        sigma = ''
        #test for repeated indices
        for l in labels:
            if (diagram_labels.count(l) > 1): sigma += l

        return sigma.replace('|', '')

    def diagram_latex(self, diagram_definition):
        #generate the Latex expression for diagram

        #factor
        value = diagram_definition['factor'] * 2
        factor = '' if value == 0 else '\\frac{1}{' + str(value) + '}'
        #sign
        sign = '' if diagram_definition['sign'][0] == 1 else '-'
        #permutation
        value = diagram_definition['permutation']
        permutation = '' if value == '' else value
        #summation
        value = diagram_definition['sigma']
        sigma = '' if value == '' else ' \\displaystyle \\sum_{' + value + '}'

        #Hamiltonian operator
        h = diagram_definition['labels'][0]
        operator = ' f_{' + h + '}' if len(h) == 2 else ' \\langle ' + h[:2] +'\\Vert ' + h[2:] + ' \\rangle'

        #amplitude operators
        amplitudes = ''
        a = diagram_definition['labels'][1]
        if a != '':
            for t in a.split('|'):
                type = len(t)//2
                amplitudes += ' t^{' + t[:type] + '}_{' + t[type:] + '}'
        amplitudes = amplitudes[::-1].replace('t', 'r', 1)[::-1]

        return sign + factor + permutation + sigma + operator + amplitudes

    def order_operator(self, type, string):
        #put operator string in correct order

        p = 'cdef' ; h = 'klmn'

        if type == 'pp': return  string[::-1]
        if (type == 'hp') and (string[0] in p): return string[::-1]

        return string

    def diagram_code(self, diagram_definition):
        #generate Python code for diagram

        no_amplitude = diagram_definition['labels'][1] == ''

        #sign and factor
        diagram = '' if diagram_definition['sign'][0] == 1 else '-'
        f = '' if diagram_definition['factor'] == 0 else str(0.5**diagram_definition['factor']) + ' * '
        diagram += f + 'np.einsum(\''

        #einsum string
        diagram += diagram_definition['labels'][0] + ','
        if no_amplitude: diagram = diagram[:-1]
        diagram += diagram_definition['labels'][1].replace('|', ',') + '->'
        diagram += 'abcd'[:self.expansion.order-1] + 'ijkl'[:self.expansion.order-1] + '\', '  #expansion order one less in eom

        #operators
        op = 'fs[' if len(diagram_definition['labels'][0]) == 2 else 'gs['
        for s in diagram_definition['labels'][0]:
            if s in 'ijklmn': op += 'o,'
            if s in 'abcdef': op += 'v,'
        diagram += op[:-1] + '], '

        #amplitudes
        if not no_amplitude:
            diagram += ', '.join(['t' + ['s','d','t','q'][len(i)//2 - 1] for i in diagram_definition['labels'][1].split('|')]) + ')'
            diagram = diagram[::-1].replace('t', 'r', 1)[::-1]
        else:
            diagram = diagram[:-2] + ')'

        return diagram

    def sb_reference(self, type, diagram, amplitudes):
        #return a reference to the diagram from Shavitt & Bartlett

        library_string = diagram + ':' + '|'.join(amplitudes)
        sb_lib = {}

        if type == 'S':
            sb_lib = {'+:+-':'S_{1}', '-|+:+-|+-':'S_{2}', '+-|+:+-|+-':'S_{3}', '+--|+:++--|+-':'S_{4}', '+-|-|+:+-|+-|+-':'S_{5}', '-:+-':'S_{6}',
            '+|-:+-|+-':'S_{7}', '+-|-:+-|+-':'S_{8}', '++-|-:++--|+-':'S_{9}', '+|+-|-:+-|+-|+-':'S_{10}', '+-:+-':'S_{11}', '+|+-:+-|+-':'S_{12}',
            '-|+-:+-|+-':'S_{13}', '+-|+-:++--|+-':'S_{14}', '+|-|+-:+-|+-|+-':'S_{15}', '+-:++--':'S_{16}', '+-|+-:+-|++--':'S_{17}', '++-:++--':'S_{18}',
            '-|++-:+-|++--':'S_{19}', '+--:++--':'S_{20}', '+|+--:+-|++--':'S_{21}'}

        if type == 'D':
            sb_lib = {'+:++--':'D_{1}', '-|+:+-|++--':'D_{2}', '+-|+:+-|++--':'D_{3}', '+--|+:++--|++--':'D_{4}', '+-|-|+:+-|+-|++--':'D_{5}', '-:++--':'D_{6}',
            '+|-:+-|++--':'D_{7}', '+-|-:+-|++--':'D_{8}', '++-|-:++--|++--':'D_{9}', '+|+-|-:+-|+-|++--':'D_{10}', '--:++--':'D_{11}', '+|--:+-|++--':'D_{12}',
            '++|--:++--|++--':'D_{13}', '+|+|--:+-|+-|++--':'D_{14}', '++:++--':'D_{15}', '-|++:+-|++--':'D_{16}', '--|++:++--|++--':'D_{17}', '-|-|++:+-|+-|++--':'D_{18}',
            '+-:++--':'D_{19}', '+|+-:+-|++--':'D_{20}', '-|+-:+-|++--':'D_{21}', '+-|+-:++--|++--':'D_{22}', '+|-|+-:+-|+-|++--':'D_{23}', '+:+-':'D_{24}', 
            '-|+:++--|+-':'D_{25}', '-|+:+-|+-':'D_{26}', '+|+:+-|+-':'D_{27}', '+-|+:++--|+-':'D_{28}', '+|-|+:+-|+-|+-':'D_{29}', '-|+-|+:+-|++--|+-':'D_{30}',
            '+|+|-|-:+-|+-|+-|+-':'D_{31}', '-:+-':'D_{32}', '+|-:++--|+-':'D_{33}', '-|-:+-|+-':'D_{34}', '+|-:+-|+-':'D_{35}', '+-|-:++--|+-':'D_{36}',
            '+|-|-:+-|+-|+-':'D_{37}', '+|+-|-:+-|++--|+-':'D_{38}', '+|-|-|+:+-|+-|+-|+-':'D_{39}', '-|++-:++--|++--':'D_{40}', '+|+--:++--|++--':'D_{41}',
            '-|+-:++--|+-':'D_{42}', '+|+-:++--|+-':'D_{43}', '+|-|+-:+-|++--|+-':'D_{44}', '-|+|+-:+-|++--|+-':'D_{45}', '++|-:++--|+-':'D_{46}', '--|+:++--|+-':'D_{47}',
            '+|+|-:+-|+-|+-':'D_{48}', '-|-|+:+-|+-|+-':'D_{49}', '+|--|+:+-|++--|+-':'D_{50}', '-|++|-:+-|++--|+-':'D_{51}', '+-|-|+:+-|++--|+-':'D_{52}',
            '+-|+|-:+-|++--|+-':'D_{53}'}

        if type == 'T':
            sb_lib = {}

        if type == 'Q':
            sb_lib = {}

        return sb_lib.get(library_string, '')

    def allocate_labels(self):
        #assign the internal and external labels

        ade = []
        #loop over terms
        for term in self.expansion.terms:

            id = term[0] ; term_amplitudes = [''.join(i) for i in term[3]]

            #loop over Hamiltonian types for term
            for hamiltonian in self.expansion.diagrams[id].keys():

                #loop over diagrams
                for diagram in self.expansion.diagrams[id][hamiltonian]:

                    diagram_definition = {'id'    : id + '(' + hamiltonian + ')'}
                    f = {'pp': lablib.pp, 'hh': lablib.hh, 'ph': lablib.ph, 'hp': lablib.hp,
                         'ppph': lablib.ppph, 'hphh': lablib.hphh, 'hhhh': lablib.hhhh, 'pppp': lablib.pppp,
                         'phph': lablib.phph, 'phpp': lablib.phpp, 'hhhp': lablib.hhhp, 'hhpp': lablib.hhpp,
                         'pphh': lablib.pphh}[hamiltonian]

                    #add the labels
                    diagram_definition['labels'] = f(self, hamiltonian, diagram, term_amplitudes)
                    #add sign of diagram
                    diagram_definition['sign']   = self.diagram_sign(diagram)
                    #add multiplicative factor
                    diagram_definition['factor'] = self.diagram_factor(diagram, term_amplitudes)
                    #add summation indices
                    diagram_definition['sigma'] = self.diagram_sigma(diagram_definition['labels'])
                    #add permutation factor to definition
                    diagram_definition['permutation'] = self.diagram_permutations(diagram_definition['labels'])
                    #add latex expression to definition
                    diagram_definition['latex'] = self.diagram_latex(diagram_definition)
                    #add Python code
                    diagram_definition['code'] = self.diagram_code(diagram_definition)
                    #add Shavitt & Bartlett reference
                    diagram_definition['reference'] = self.sb_reference(self.expansion.type, diagram, term_amplitudes)
                    #add operator contraction string
                    diagram_definition['contraction'] = diagram

                    self.labels.refresh()

                    ade.append(diagram_definition)

        return ade

if __name__ == "__main__":

    from eom_validation.validate import validate
    
    s = expansion('SD', 'S')
    s.order_by_reference()
    v = validate(s)
    v.compare('all')

    d = expansion('SD', 'D')
    d.order_by_reference()
    v = validate(d)
    v.compare('all')

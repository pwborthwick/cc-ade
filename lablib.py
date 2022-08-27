def pp(self, id, diagram, amplitudes):
    #pp 1-body Hamiltonian

    internal_label = self.labels.get('+', 'i')
    term_amplitude = amplitudes[0][::-1].replace('+', internal_label, 1)[::-1]
    term_amplitude = self.labels.replace_all(term_amplitude, '+-')

    external_label = self.labels.get('+', 'e')
    operator = external_label + internal_label

    return operator, term_amplitude

def hh(self, id, diagram, amplitudes):
    #hh 1-body Hamiltonian

    internal_label = self.labels.get('-', 'i')
    term_amplitude = amplitudes[0][::-1].replace('-', internal_label, 1)[::-1]
    term_amplitude = self.labels.replace_all(term_amplitude, '+-')

    external_label = self.labels.get('-', 'e')
    operator = internal_label + external_label 

    return operator, term_amplitude

def ph(self, id, diagram, amplitudes):
    #ph 1-body Hamiltonian

    return 'ai', ''

def hp(self, id, diagram, amplitudes):
    #hp 1-body Hamiltonian

    internal_labels = [self.labels.get('+',  'i'),
                       self.labels.get('-', 'i')]

    if diagram == '+-':
        term_amplitude = amplitudes[0][::-1].replace('+', internal_labels[0], 1).replace('-', internal_labels[1], 1)[::-1]
    if diagram == '+|-':
        term_amplitude = self.labels.replace_near(diagram, amplitudes, internal_labels)
    if diagram == '-|+':
        term_amplitude = self.labels.replace_near(diagram, amplitudes, internal_labels[::-1])
 
    term_amplitude = self.labels.replace_all(term_amplitude, '+-')

    operator = internal_labels[1] + internal_labels[0]

    return operator, term_amplitude

def ppph(self, id, diagram, amplitudes):
    #ppph 2-body Hamiltonian

    internal_label =   self.labels.get('+', 'i')

    term_amplitude = amplitudes[0][::-1].replace('+', internal_label, 1)[::-1]
    term_amplitude = self.labels.replace_all(term_amplitude, '+-')

    operator = (self.labels.get('+', 'e') + self.labels.get('+', 'e') + 
                internal_label + self.labels.get('-', 'e'))

    return operator, term_amplitude

def hphh(self, id, diagram, amplitudes):
    #ppph 2-body Hamiltonian

    internal_label = self.labels.get('-', 'i')

    term_amplitude = amplitudes[0][::-1].replace('-', internal_label, 1)[::-1]
    term_amplitude = self.labels.replace_all(term_amplitude, '+-')

    operator = (internal_label + self.labels.get('+', 'e') +
                self.labels.get('-', 'e') + self.labels.get('-', 'e'))

    return operator, term_amplitude

def hhhh(self, id, diagram, amplitudes):
    #hhhh 2-body Hamiltonian

    term_amplitude = '|'.join(amplitudes)

    internal_labels = [self.labels.get('+', 'i'), self.labels.get('+', 'i')]

    if diagram == '++':
        term_amplitude = term_amplitude.replace('++', ''.join(internal_labels[:2]), 1)
    if diagram == '+|+':
        term_amplitude = self.labels.replace_near(diagram, amplitudes, internal_labels)
        amp = term_amplitude.split('|')
        term_amplitude = ''.join(self.labels.replace_all(amp[0], '+') + '|' + amp[1])

    external_labels = [self.labels.get('+', 'e'), self.labels.get('+', 'e')]
    term_amplitude = self.labels.replace_all(term_amplitude, '+-')

    operator = ''.join(external_labels) + ''.join(internal_labels)

    return operator, term_amplitude

def pppp(self, id, diagram, amplitudes):
    #pppp 2-body Hamiltonian

    term_amplitude = '|'.join(amplitudes)

    internal_labels = [self.labels.get('-', 'i'), self.labels.get('-', 'i')]
 
    if diagram == '--':
        term_amplitude = term_amplitude.replace('--',  ''.join(internal_labels), 1)
    if diagram == '-|-':
        term_amplitude = self.labels.replace_near(diagram, amplitudes, internal_labels)
        amp = term_amplitude.split('|')
        term_amplitude = ''.join(self.labels.replace_all(amp[0], '-') + '|' + amp[1])

    external_labels = [self.labels.get('-', 'e'), self.labels.get('-', 'e')]
    term_amplitude = self.labels.replace_all(term_amplitude, '+-')

    operator =   ''.join(internal_labels) + ''.join(external_labels) 

    return operator, term_amplitude

def phph(self, id, diagram, amplitudes):
    #phph 2-body Hamiltonian

    interchange = False
    internal_labels = [self.labels.get('+', 'i'),
                       self.labels.get('-', 'i')]

    if diagram == '+-':
        term_amplitude = amplitudes[0].replace('+', internal_labels[0], 1).replace('-', internal_labels[1], 1)
    if diagram == '+|-':
        term_amplitude = self.labels.replace_near(diagram, amplitudes, internal_labels)
        amp = term_amplitude.split('|')
        term_amplitude = ''.join(self.labels.replace_all(amp[0], '+-')) + '|' + amp[1].replace('+', self.labels.get('+', 'e'), 1)
    if diagram == '-|+':
        term_amplitude = self.labels.replace_near(diagram, amplitudes, internal_labels[::-1])
        amp = term_amplitude.split('|')
        term_amplitude = ''.join(self.labels.replace_all(amp[0], '+-')) + '|' + amp[1]
        interchange = not interchange

    external_labels = [self.labels.get('+', 'e'),
                       self.labels.get('-', 'e')]
    term_amplitude = self.labels.replace_all(term_amplitude, '+-')

    bra = external_labels[0] + internal_labels[1]
    if interchange: bra = bra[::-1]
    operator =  bra + external_labels[1] + internal_labels[0]

    return operator, term_amplitude

def phpp(self, id, diagram, amplitudes):
    #phpp 2-body Hamiltonian

    interchange = False
    internal_labels = [self.labels.get('+', 'i'),
                       self.labels.get('+', 'i'),
                       self.labels.get('-', 'i')]

    if diagram == '++-':
        term_amplitude = amplitudes[0].replace('++', ''.join(internal_labels[:2]), 1).replace('-', internal_labels[2], 1)
        external_label = self.labels.get('+', 'e')
    if diagram == '+-|+':
        term_amplitude = self.labels.replace_near('+|+', amplitudes, internal_labels[:2])
        amp = term_amplitude.split('|') 
        external_label = self.labels.get('+', 'e', _at=len(amp[0])//2)
        term_amplitude = amp[0][::-1].replace('-', internal_labels[2], 1)[::-1] +  '|' + amp[1]
    if diagram == '+|+-':
        term_amplitude = self.labels.replace_near('+|+', amplitudes, internal_labels[:2])
        amp = term_amplitude.split('|')
        term_amplitude = amp[0] + '|' + amp[1].replace('-', internal_labels[2], 1) 
        external_label = self.labels.get('+', 'e', _at=len(amp[0])//2)
        interchange = not interchange
    if diagram == '++|-':
        term_amplitude = self.labels.replace_near(diagram, amplitudes, [internal_labels[1]+internal_labels[0], internal_labels[2]])
        amp = term_amplitude.split('|')  
        external_label = self.labels.get('+', 'e', _at=len(amp[0])//2 - 1)
        interchange = not interchange
    if diagram == '-|++':
        term_amplitude = self.labels.replace_near(diagram, amplitudes, [internal_labels[2], internal_labels[0]+internal_labels[1]])
        amp = term_amplitude.split('|')
        external_label = self.labels.get('+', 'e', _at=len(amp[0])//2  + 1)
    if diagram == '-|+|+':
        term_amplitude = self.labels.replace_near('-|+', amplitudes[:2], [internal_labels[2], internal_labels[0]])
        term_amplitude = term_amplitude + '|' + amplitudes[2].replace('+', internal_labels[1], 1)
        amp = term_amplitude.split('|')
        external_label = self.labels.get('+', 'e', _at=(len(amp[0]) + len(amp[1]))//2)  
    if diagram == '+|-|+':
        term_amplitude = self.labels.replace_near('-|+', amplitudes[1:], internal_labels[1:][::-1])
        term_amplitude = amplitudes[0][::-1].replace('+', internal_labels[0], 1)[::-1] + '|' + term_amplitude
        amp = term_amplitude.split('|')
        external_label = self.labels.get('+', 'e', _at=(len(amp[0]) + len(amp[1]))//2)           
    if diagram == '+|+|-':
        term_amplitude = self.labels.replace_near('+|-', amplitudes[1:], internal_labels[1:])
        term_amplitude = amplitudes[0][::-1].replace('+', internal_labels[0], 1)[::-1] + '|' + term_amplitude
        amp = term_amplitude.split('|')
        external_label = self.labels.get('+', 'e', _at=len(amp[0])//2)           
        interchange = not interchange

    term_amplitude = self.labels.replace_all(term_amplitude, '+-')
    bra = internal_labels[2] + external_label
    if interchange: bra = bra[::-1]
    ket = ''.join(internal_labels[:2])

    return bra + ket, term_amplitude

def hhhp(self, id, diagram, amplitudes):
    #phpp 2-body Hamiltonian

    interchange = False
    internal_labels = [self.labels.get('-', 'i'),
                       self.labels.get('-', 'i'),
                       self.labels.get('+', 'i')]

    if diagram == '+--':
        term_amplitude = amplitudes[0].replace('--', ''.join(internal_labels[:2]), 1).replace('+', internal_labels[2], 1)
        external_label = self.labels.get('-', 'e')
    if diagram == '+-|-':
        term_amplitude = self.labels.replace_near('-|-', amplitudes[:2], internal_labels[:2])
        amp = term_amplitude.split('|')
        term_amplitude = amp[0][::-1].replace('+', internal_labels[2], 1)[::-1] + '|' + amp[1]
        external_label = self.labels.get('-', 'e', _at=len(amp[0])//2)
    if diagram == '-|+-':
        term_amplitude = self.labels.replace_near('-|-', amplitudes[:2], internal_labels[:2])
        amp = term_amplitude.split('|')
        term_amplitude = amp[0] + '|' + amp[1].replace('+', internal_labels[2], 1)
        external_label = self.labels.get('-', 'e', _at=len(amp[0])//2)
        interchange = not interchange
    if diagram == '+|--':
        term_amplitude = self.labels.replace_near(diagram, amplitudes, [internal_labels[2], internal_labels[0]+internal_labels[1]])
        amp = term_amplitude.split('|')
        external_label = self.labels.get('-', 'e', _at=len(amp[0])//2  + 1)
    if diagram == '--|+':
        term_amplitude = self.labels.replace_near(diagram, amplitudes, [internal_labels[1]+internal_labels[0], internal_labels[2]])
        amp = term_amplitude.split('|')
        external_label = self.labels.get('-', 'e', _at=len(amp[0])//2  - 1)
        interchange = not interchange
    if diagram == '+|-|-':
        term_amplitude = self.labels.replace_near('+|-', amplitudes[:2], [internal_labels[2], internal_labels[0]])
        term_amplitude = term_amplitude + '|' + amplitudes[2].replace('-', internal_labels[1], 1)
        amp = term_amplitude.split('|')
        external_label = self.labels.get('-', 'e', _at=(len(amp[0]) + len(amp[1]))//2)  
    if diagram == '-|+|-':
        term_amplitude = self.labels.replace_near('+|-', amplitudes[1:], internal_labels[1:][::-1])
        term_amplitude = amplitudes[0][::-1].replace('-', internal_labels[0], 1)[::-1] + '|' + term_amplitude
        amp = term_amplitude.split('|')
        external_label = self.labels.get('-', 'e', _at=(len(amp[0]) + len(amp[1]))//2)           
    if diagram == '-|-|+':
        term_amplitude = self.labels.replace_near('-|+', amplitudes[1:], internal_labels[1:])
        term_amplitude = amplitudes[0][::-1].replace('-', internal_labels[0], 1)[::-1] + '|' + term_amplitude
        amp = term_amplitude.split('|')
        external_label = self.labels.get('-', 'e', _at=len(amp[0])//2)           
        interchange = not interchange

    term_amplitude = self.labels.replace_all(term_amplitude, '+-')
    ket = internal_labels[2] + external_label
    if interchange: ket = ket[::-1]
    bra = ''.join(internal_labels[:2])

    return bra + ket, term_amplitude

def pphh(self, id, diagram, amplitude):
    #pphh 2-body Hamiltonian

    return 'abij', ''

def hhpp(self, id, diagram, amplitudes):
    #hhpp 2-body Hamiltonian

    interchange = False
    internal_labels = [self.labels.get('+', 'i'),
                       self.labels.get('+', 'i'),
                       self.labels.get('-', 'i'),
                       self.labels.get('-', 'i')]

    if diagram == '++--':
        term_amplitude = amplitudes[0][::-1].replace('++', internal_labels[0]+internal_labels[1], 1).replace('--', internal_labels[3]+internal_labels[2], 1)[::-1]

    if diagram == '+|+--':
        term_amplitude = self.labels.replace_near('+|--', amplitudes, [internal_labels[0], internal_labels[2]+internal_labels[3]])
        amp = term_amplitude.split('|')
        term_amplitude = amp[0] + '|' + amp[1].replace('+', internal_labels[1], 1)
        interchange = not interchange
    if diagram == '+--|+':
        term_amplitude = self.labels.replace_near('--|+', amplitudes, [internal_labels[3]+internal_labels[2], internal_labels[1]])
        amp = term_amplitude.split('|')
        term_amplitude = amp[0][::-1].replace('+', internal_labels[0], 1)[::-1] + '|' + amp[1]
        interchange = not interchange
    if diagram == '++|--':
        term_amplitude = self.labels.replace_near('++|--', amplitudes, [internal_labels[1]+internal_labels[0], 
                                                                        internal_labels[2]+internal_labels[3]])
    if diagram == '--|++':
        term_amplitude = self.labels.replace_near('--|++', amplitudes, [internal_labels[3]+internal_labels[2], 
                                                                        internal_labels[0]+internal_labels[1]])
    if diagram == '-|++-':
        term_amplitude = self.labels.replace_near('-|++', amplitudes, [internal_labels[2], internal_labels[0]+internal_labels[1]])
        amp = term_amplitude.split('|')
        term_amplitude = amp[0] + '|' + amp[1].replace('-', internal_labels[3], 1)
        interchange = not interchange
    if diagram == '++-|-':
        term_amplitude = self.labels.replace_near('++|-', amplitudes, [internal_labels[1]+internal_labels[0], internal_labels[2]])
        amp = term_amplitude.split('|')
        term_amplitude = amp[0][::-1].replace('-', internal_labels[3], 1)[::-1] + '|' + amp[1]  
    if diagram == '+-|+-':
        term_amplitude = self.labels.replace_near('+|-', amplitudes, [internal_labels[0], internal_labels[3]])
        amp = term_amplitude.split('|')
        term_amplitude = amp[0][::-1].replace('-', internal_labels[2], 1)[::-1] + '|' + amp[1].replace('+', internal_labels[1], 1)

    if diagram == '++|-|-':
        term_amplitude =  self.labels.replace_near('-|-', amplitudes[1:], [internal_labels[2], internal_labels[3]])
        term_amplitude = amplitudes[0][::-1].replace('++', internal_labels[1]+internal_labels[0], 1)[::-1]
    if diagram == '--|+|+':
        term_amplitude = self.labels.replace_near('+|+', amplitudes[1:], [internal_labels[0], internal_labels[1]])
        term_amplitude = amplitudes[0][::-1].replace('--', internal_labels[2]+internal_labels[3], 1)[::-1]
    if diagram == '+|+-|-':
        term_amplitude =  self.labels.replace_near('+|-', amplitudes[:2], [internal_labels[0], internal_labels[2]])
        term_amplitude = term_amplitude[::-1].replace('+', internal_labels[1], 1)[::-1] + '|' + amplitudes[2].replace('-', internal_labels[3], 1)
        interchange = not interchange
    if diagram == '-|+-|+':
        term_amplitude =  self.labels.replace_near('-|+', amplitudes[:2], [internal_labels[2], internal_labels[0]])
        term_amplitude = term_amplitude[::-1].replace('-', internal_labels[3], 1)[::-1] + '|' + amplitudes[2].replace('+', internal_labels[1], 1)
        interchange = not interchange
    if diagram == '+|-|+-':
        term_amplitude =  self.labels.replace_near('+|-', amplitudes[:2], [internal_labels[0], internal_labels[2]])
        term_amplitude = term_amplitude + '|' + amplitudes[2].replace('+', internal_labels[1], 1).replace('-', internal_labels[3], 1)
    if diagram == '-|+|+-':
        term_amplitude =  self.labels.replace_near('-|+', amplitudes[:2], [internal_labels[2], internal_labels[0]])
        term_amplitude = term_amplitude + '|' + amplitudes[2].replace('+', internal_labels[1], 1).replace('-', internal_labels[3], 1)
    if diagram == '-|-|++':
        term_amplitude =  self.labels.replace_near('-|-', amplitudes[:2], [internal_labels[2], internal_labels[3]])
        term_amplitude = term_amplitude + '|' + amplitudes[2].replace('++', internal_labels[0]+internal_labels[1], 1)
    if diagram == '-|++|-':
        term_amplitude =  self.labels.replace_near('-|+', amplitudes[:2], [internal_labels[2], internal_labels[0]])
        term_amplitude = term_amplitude[::-1].replace('+', internal_labels[1], 1)[::-1] + '|' + amplitudes[2].replace('-', internal_labels[3], 1)
    if diagram == '+|+|--':
        term_amplitude =  self.labels.replace_near('+|+', amplitudes[:2], [internal_labels[0], internal_labels[1]])
        term_amplitude = term_amplitude + '|' + amplitudes[2].replace('--', internal_labels[2]+internal_labels[3], 1)
    if diagram == '+|--|+':
        term_amplitude = self.labels.replace_near('+|-', amplitudes[:2], [internal_labels[0], internal_labels[2]])
        term_amplitude = term_amplitude[::-1].replace('-', internal_labels[3], 1)[::-1] + '|' + amplitudes[2].replace('+', internal_labels[1], 1)
    if diagram == '+-|-|+':
        term_amplitude = self.labels.replace_near('-|+', amplitudes[1:], [internal_labels[3], internal_labels[1]])
        term_amplitude = amplitudes[0][::-1].replace('-', internal_labels[2], 1).replace('+', internal_labels[0], 1)[::-1] + '|' + term_amplitude
    if diagram == '+-|+|-':
        term_amplitude = self.labels.replace_near('+|-', amplitudes[1:], [internal_labels[1], internal_labels[3]])
        term_amplitude = amplitudes[0][::-1].replace('-', internal_labels[2], 1).replace('+', internal_labels[0], 1)[::-1] + '|' + term_amplitude

    if diagram == '+|+|-|-':
        term_amplitude = self.labels.replace_near('-|-', amplitudes[2:], [internal_labels[2], internal_labels[3]])
        term_amplitude = (amplitudes[0][::-1].replace('+', internal_labels[0], 1)[::-1] + '|' + amplitudes[1][::-1].replace('+', internal_labels[1], 1)[::-1] +
                         '|' + term_amplitude)
    if diagram == '-|-|+|+':
        term_amplitude = self.labels.replace_near('+|+', amplitudes[2:], [internal_labels[0], internal_labels[1]])
        term_amplitude = (amplitudes[0][::-1].replace('-', internal_labels[2], 1)[::-1] + '|' + amplitudes[1][::-1].replace('-', internal_labels[3], 1)[::-1] +
                         '|' + term_amplitude)
    if diagram == '+|-|-|+':
        term_amplitude = self.labels.replace_near('+|-', amplitudes[:2], [internal_labels[0], internal_labels[2]])
        term_amplitude += '|' + self.labels.replace_near('-|+', amplitudes[2:], [internal_labels[3], internal_labels[1]])
    if diagram == '+|-|+|-':
        term_amplitude = self.labels.replace_near('+|-', amplitudes[:2], [internal_labels[0], internal_labels[2]])
        term_amplitude += '|' + self.labels.replace_near('+|-', amplitudes[2:], [internal_labels[1], internal_labels[3]])
    if diagram == '-|+|-|+':
        term_amplitude = self.labels.replace_near('-|+', amplitudes[:2], [internal_labels[2], internal_labels[0]])
        term_amplitude += '|' + self.labels.replace_near('-|+', amplitudes[2:], [internal_labels[3], internal_labels[1]])
    if diagram == '-|+|+|-':
        term_amplitude = self.labels.replace_near('-|+', amplitudes[:2], [internal_labels[2], internal_labels[0]])
        term_amplitude += '|' + self.labels.replace_near('+|-', amplitudes[2:], [internal_labels[1], internal_labels[3]])

    term_amplitude = self.labels.replace_all(term_amplitude, '+-')
    bra = ''.join(internal_labels[2:])
    if interchange: bra = bra[::-1]
    ket = ''.join(internal_labels[:2])

    return bra + ket, term_amplitude    
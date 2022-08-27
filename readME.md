# cc-ade
Algebraic Evaluation of Coupled-Cluster Diagrams

The repo COGUS enabled derivation of various coupled-cluster codes via a second quantization approach. It used Sympy and by way of Wicks and Baker-Cambell-Hausdorff theorems generated coupled-cluster code. This is a purely algebraic approach. In their book 'Many-Body Methods in Chemistry and Physics - MBPT and Coupled-Cluster Theory' Shavitt and Bartlett offer an alternative approach based on a diagrammatic foundation. 

The idea starts with the exponential expansion of cluster operators, for example, if we are looking at singles and doubles level then we consider e<sup>(T_1 + T_2)</sup>. We can again use Sympy to do the work of expanding the exponential power. This gives us a series in powers of products of T<sub>n</sub> operators. Each T<sub>n</sub> is now represented by a sign string (s-string) so T<sub>1</sub> is +-, T<sub>2</sub> is ++-- etc. In diagrammatic terms the + signs are up (particle) lines and the - signs are down (hole) lines. A product of amplitudes are written as +-|++--, which would represent T<sub>1</sub> T<sub>2</sub>. 

One and two body Hamiltonian operators are also represented by s-strings. In diagrammatic terms the s-string is the signs of lines the Hamitonian presents below the interaction vertex. In the figure below we see a T<sub>2</sub> amplitude and an hhhp 2-body Hamiltonian have s-strings of ++-- and +-- respectively. An important quantity is the _interaction number_ - for an amplitude this is the level of excitation or the number of + signs and for a Hamiltonian it is the number of + signs above the vertex minus the number below. For our T<sub>2</sub> amplitude then the interaction number is +2 and for the Hamiltonian it is -1.
![feynman(2)](https://user-images.githubusercontent.com/73105740/184648291-3c7f2627-7560-4e4a-bf60-465dd690699a.png)

Each term in the series then has an interaction number, which is the sum of the interaction numbers of the amplitudes making it up. We also have a target excitation which for singles will be +1, for doubles +2 etc. If a given term has an interaction number of say 3 and we are considering doubles then only Hamiltonian operators with a interaction number of 2-3=-1 can contract with that term to make a connected diagram.

Our task is then to, for each term, get all Hamiltonians that have the correct interaction number. We then must work out all possible valid ways that those Hamiltonians can contract with the terms s-string. For example, if the term is T<sub>1</sub> T<sub>2</sub> and we are considering doubles then our hhhp element with an interaction number of -1 is suitable and can contract as -|+-, +|-- or +-|- with +-|++--. We can use itertools to generate all permutations of the Hamiltonian s-string and prune out the invalid ones eg --|+.

We now have for each term a number of tuples of a Hamiltonian s-string and the terms s-string. Each of these is a connected coupled-cluster diagram. From the Hamiltonian s-string we can determine the sign of the diagram and the multiplying factor. To complete the code generation we need to label the s-strings. Each Hamiltonian s-string has an algorithm for labelling, for +|-- with the hhhp Hamiltonian that algorithm is t<sup>(+)e</sup><sub>(-)</sub> t<sup>(+)</sup><sub>mn[-)</sub> <mn||e\[>. *This is interpreted as (+) contiguously label with external labels, \[+) miss next available label and continue contiguous labelling, mn and ef represent the internal labels available at the level of theory under consideration* so we would first label Hamiltonian with internal labels (for doubles these are {k,l} and {c,d}) as  +<sup>c</sup>-<sup>k</sup>-<sup>l</sup>, then following the labelling algorithm we label the amplitudes, using the external doubles labels {i,j} and {a,b}, as +<sup>c</sup>-<sup>i</sup>|+<sup>a</sup>+<sup>b</sup>-<sup>k</sup>-<sup>l</sup> with <kl||cj>. The permutations are then  **P**(ij) because labels i and j occur on different elements of the expression - there is no (ab) permutation as they appear together.

We have now generated the algebraic equivalent of the coupled-cluster diagram. The cc-ade program has been used to generate all diagrams for S,D,T and Q
levels of the theory although the _lablib_ library of labelling algorithms is applicable for all levels of theory. The generation of code is much much faster than COGUS although less versatile.

The ideas are immediately applicable to EOM_CCSD(TQ) theory and a file *cc-ade_eom* illustrates using *lablib* to generate EOM-CCSD-EE codes. Much of the code in *cc-ade-eom* is the same as *cc-ade* and a better implementation would incorporate plain CC and EOM functionality in one set of code. Validation tests are provided for EOM.
- - -
Code
- - -
+ **cc_ade** - the main class module for expansion and diagram generation
+ **cc_eom_ade** - the main class module for EOM diagrams (a lot of replicated code from cc_ade)
+ **lablib** - a set of labelling algorithms
+ **cc-ade-theory.ipynb** - jupyter notebook explaining the idea behind the program
+ **cc-ade-implementation.ipynb**- jupyter notebook details the use of the program
+ **cc-ade-eom.ipynb** - notes on implementation of EOM
+ **SDTQ-diagrams.ipynb** - reference for all diagrams with Shavitt and Bartlett reference codes
+ **ccsd** - example using cc-ade to compute ccsd energy
+ **qcisd** - example using cc-ade to compute qcisd energy
+ **ccsdt-2** - example of adding custom terms to expansion

**validation folder**
+ **validate** - class to test cc-ade generated code against reference Shavitt and Bartlett 
+ **c.npz** - compressed archive S,D,T and Q amplitudes of carbon atom in 3-21g basis (not converged for testing only!)
+ **hf.npz** - compressed archive S,D and T amplitudes of hydrogen fluoride in 6-21g basis
+ **reference.txt** - einsum and permutation based on Shavitt and Bartlett diagrams
+ **c-3-21g-sdtq.txt** - results of validation test of cc-ade generated diagrams against Shavitt and Bartlett
+ **hf-6-31g-sdt.txt** - results of validation test of cc-ade generated diagrams against Shavitt and Bartlett

**eom_valiation folder**
+ **validate** - class to test cc-eom-ade generated code against reference values 
+ **hf-eom.npz** - r-vectors from HF in 6-31g basis
+ **reference.txt** - einsum and permutation based on reference code diagrams
+ **hf-6-31g-cc-eom.txt** - results of validation test of cc-eom-ade generated diagrams against Shavitt and Bartlett




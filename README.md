# quantum-algorithm-simulation

Barak N. - Independent Work Fall 2016: Simulation of Quantum Algorithms

A quantum algorithm simulation library.

The code that makes up the quantum simulator is divided into separate modules with code for working with vectors, states, operators, and registers. Each module (in its separate file) contains the relevant definitions and functions as well as varying amounts of unit testing code to verify its desired functionality.

## Getting Started

To compile and run tests, use:
```
$ make run_tests
```

This should produce the output

```
Vector tests:           passed.
State tests:            passed.
Operator tests:         passed.
Register tests:         passed.
Input w.
w = 
```

At this point, the tests have passed, and you are prompted for an input `w` with which to initialize the quantum register.

At the prompt, type in a number to see a sample run of the [inversion algorithm](https://arxiv.org/pdf/1511.08253.pdf) using the quantum-algorithm-simulation library.

```
Input w.
w = 24.896

Input number: 24.895508 Actual Inverse:         0.040167889224493
-------------------------
State of 61 qubits in the basis state:
(1.0+0.0i) |0001100011100101010000000000000000000000000000000000000000000 (31ca80000000000)>

Applying CNot(60 ~42 -> 34)
State of 61 qubits in the basis state:
(1.0+0.0i) |0001100011100101010000000000000000000000000000000000000000000 (31ca80000000000)>
[...]
Applying CNot(57 ~42 -> 37)
State of 61 qubits in the basis state:
(1.0+0.0i) |0001100011100101010000010000000000000000000000000000000000000 (31ca82000000000)>
[...]
Applying PauliX(-> 42)
State of 61 qubits in the basis state:
(1.0+0.0i) |0001100011100101010000010000000000000000000000000000000000000 (31ca82000000000)>


Initial approximation: 1/32                      =      0.031250000000000
-----------------------------
State of 61 qubits in the basis state:
(1.0+0.0i) |0001100011100101010000010000000000000000000000000000000000000 (31ca82000000000)>
[...]
x_4: 1/24.97                     =      0.040039062500000
-----------------------------
State of 61 qubits in the basis state:
(1.0+0.0i) |0001100011100101010000010100101000000001101001000110000101001 (31ca82940348c29)>
```

## Built With

* [C        ](https://en.wikipedia.org/wiki/C_(programming_language)) - Written in the C language
* [GMP      ](https://gmplib.org/) -            For support of arbitrary precision integers
* [MakeTutor](http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/) -    Basis of the makefile
* [GDSL     ](http://www.nongnu.org/gdsl/) -    Generic Data Structures Library              


## Authors

* **Barak Nehoran** - *Author* 
* **Iasonas Petras** - *Adviser*

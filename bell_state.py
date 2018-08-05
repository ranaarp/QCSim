from q_operations import *
import numpy as np

a = init_zeros(2)
# a = init_state([1, 0, 0, 0])  # initialize two qubits in |00> state
print(a)
n = 2 # number of qubits

a = np.matmul(np.kron(H(1), I(1)), a) # hadamard on first qubit
print(a)

a = np.matmul(CNOT([1], 2, 2), a) # CNOT from 1st qubit to 2nd
print(a)

measure(a, n)
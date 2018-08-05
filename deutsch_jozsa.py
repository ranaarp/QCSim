from q_operations import *

qubit_num = 2  # basically the size of the function input

a = init_zeros(qubit_num + 1)

print(X(1))
print(np.kron(I(qubit_num), X(1)))
a = np.matmul(np.kron(I(qubit_num), X(1)), a)
print(a)

H_n = H(qubit_num + 1)

a = np.matmul(H_n, a)
print(a)

# oracle = CNOT([1], qubit_num + 1, qubit_num + 1)  # oracle for a balanced function
# print(oracle)
# a = np.matmul(oracle, a)

a = np.matmul(np.kron(I(qubit_num), X(1)), a)   # oracle for a constant function
print(a)


a = np.matmul(np.kron(H(qubit_num), I(1)), a)
print(a)

measure(a, qubit_num + 1)




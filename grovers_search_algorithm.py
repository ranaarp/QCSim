from q_operations import *
import numpy as np

qubit_num = 4  # works for searching in 2**qubit_num of entries


def oracle(num_search , n) :
    # n is the number of qubits you're working on, excluding the ancilla bit
    # remember that the oracle acts on n+1 qubits including the ancilla bit whereas reflection about the mean acts on n qubits

    p = np.array(np.arange(1, n+1))
    print("control bits used in forming the Tofolli gate : ",p)
    Toffoli = CNOT(p, n+1, n+1)
    print("Tofolli gate being used is :\n ",Toffoli)

    binary_number = bin(num_search)[2:]
    if len(binary_number) < n:
        for i in range(n - len(binary_number)):
            binary_number = "0" + binary_number
    # print(binary_number)

    print("binary number is :")
    print(binary_number)

    if binary_number[0] == '1':
        A = I(1)
    else:
        A = X(1)
    print("A now is :")
    print(A)

    for i in binary_number[1:]:
        print("the next bit is :",i)
        if i == '1':
            A = np.kron(A, I(1))
        if i == '0':
            A = np.kron(A, X(1))

    print("A now is (2) :\n", A)
    A = np.kron(A, I(1))
    print("A finally is :\n", A)
    # print(A)
    U = np.matmul(A, Toffoli)
    U = np.matmul(U, A)
    print("\nOracle which gives a \"1\" when it finds",num_search, "is :\n",U)
    return U


def reflection_about_mean(n):
    # reflection about the mean can be done using a matrix that is 2|psi><psi| - I , where |psi> is the state that
    # the hadamard sends the state |000....0>. Multiplying a state vector with the projector matrix |psi><psi| sends
    # sends every entry into the average value.
    # But, here we'll just directly define the matrix operator "reflection about the mean"
    # geometrically vector "v" is sent to "-v + 2a" where a is the average value

    N = 2**n
    Mean = np.full((N, N), 2/N, dtype=complex)
    R2 = -1*I(n) + Mean

    print("\nR2 : \n", R2)
    # print(R2)
    return R2


# print("U : ")
# print(oracle(2, 2))

state = init_zeros(qubit_num + 1)  # the factor +1 is to account for the ancilla bit

U1 = np.kron(I(qubit_num), X(1))

state = np.matmul(U1, state)  # flip the ancilla bit
print("state 1 : \n", state)
# print(state)

U2 = H(qubit_num + 1)

state = np.matmul(U2, state)  # putting everything into a superposition using hadamard
print("state 2 : \n", state)
# print(state)
# print(reflection_about_mean(3))


def grover_subroutine(state, iter):

    # there are two iterative steps that we have to usually apply for sqrt(N) times

    # Step 1 : flip the phase of the required entity that gives f(x) = 1
    # matrix used is : if suppose we are searching for the value 2 in a 3 qubit space, which is "010" in binary.
    # we want our oracle function to give the output "1" only when it gets a 2 as input.
    # to do that we flip all the bits that are zero in the required entity's binary representation and then apply a
    # Toffoli gate that gives 1 only if all the entries are 1. and then we undo our change to the state.
    # in the example we are considering, our oracle will be :(3 data bits, one ancilla bit):
    #       (X(1)*I(1)*X(1)*I(1)).(Toffoli from first three data bits to fourth ancilla bit).(X(1)*I(1)*X(1)*I(1))

    # Step 2 : reflection about the mean


    R1 = oracle(2, qubit_num)  # searching for the number 2

    R2 = np.kron(reflection_about_mean(qubit_num), I(1))

    for i in range(iter):
        state = np.matmul(R1, state)
        print("state after phase inversion. Iteration :", i, "\n", state)
        state = np.matmul(R2, state)
        print("state after inversion about mean. Iteration :", i, "\n", state)
        # measure_magnitudes(state, qubit_num+1, qubit_num, i+1)
        measure(state, qubit_num + 1)


grover_subroutine(state, 24)

#!/usr/bin/python
import math
import sys
import numpy as np
import matplotlib.pyplot as plt

def init_state(state_vector):
    numpy_vector = np.matrix(state_vector)
    numpy_vector = np.transpose(numpy_vector)
    print("Initial state is :")
    print(numpy_vector)
    return numpy_vector


def init_zeros(qubit_num):
    zero_vect = np.zeros(2**qubit_num, dtype=complex)
    zero_vect[0] = 1
    a = np.matrix(zero_vect)
    a = np.transpose(a)
    print(a.shape)
    print("Initial state is :")
    print(a)
    return a


def I(n):
    # identity matrix on n qubits
    I_n = np.identity(2**n, dtype=complex)   # 2**n because it takes 2^n elements to represent n qubits
    return I_n

def i_I(n):
    # identity matrix multiplied with 'i'
    i_I = complex(0, 1)*I(n)
    return i_I

def H(n):
    # hadamard matrix on n qubits
    h_1 = (float(1) / math.sqrt(2)) * np.matrix([[1, 1], [1, -1]], dtype=complex)
    h_n = h_1

    for i in range(n - 1):
        h_n = np.kron(h_1, h_n)

    return h_n


def X(n):
    # bit flip matrix on n qubits (sigma_x)
    N = 2**n
    # temp = [1]
    # for i in range(2**n-1):
    #     temp += [1]
    # x_n = np.diagflat(temp, -0).astype(complex)

    x_n = np.zeros((N, N), dtype = complex)

    i = N-1
    j = 0
    for k in range(N):
        x_n[j, i] = 1
        i -= 1
        j += 1

    return x_n


def Y(n):
    # sigma_y pauli matrix on n qubits
    y_1 = np.matrix([[0, 0-1j], [1j, 0]], dtype=complex)
    y_n = y_1

    for i in range(n - 1):
        y_n = np.kron(y_1, y_n)

    return y_n


def Z(n):
    # phase flip matrix on n qubits (sigma_z)
    z_1 = np.matrix([[1, 0], [0, -1]], dtype=complex)
    z_n = z_1

    for i in range(n - 1):
        z_n = np.kron(z_1, z_n)

    return z_n


def int2vect(num, n):
    # takes in an integer and returns the corresponding vector based on the space size n
    N = 2**n
    vect = []
    for i in range(N):
        vect += [0]
    vect[num] = 1
    vect = np.matrix(vect)
    return vect.T


def flipcheck(control_bits, target_bit, val, n):
    binary = bin(val)[2:]
    # print("binary value is : " + binary)
    if len(binary)<n:
        for i in range(n-len(binary)):
            binary = "0" + binary
    # print("binary value after appending is : " + binary)

    check = 1
    for k in control_bits:
        check *= int(binary[k - 1])

    new_binary = binary
    if check == 1:
        # binary[target_bit-1] = str((int(binary[target_bit-1])+1)%2)
        new_binary = binary[:target_bit-1] + str((int(binary[target_bit-1])+1)%2) + binary[target_bit:]

    return int2vect(int(new_binary, 2), n)


def CNOT(control_bits,target_bit, n):
    # controlled NOT gate applied on n qubits. Function: to flip target_bit when all the control_bits are 1.
    # qubit numbering starts from 1. that is, 1 implies the first qubit, not zero.
    N = 2**n
    cnot_n = np.zeros([N, 1], dtype=float)
    cnot_n[0] = 1
    for i in range(1,N):
        # print("value of i : " + str(i))
        cnot_n = np.concatenate((cnot_n, flipcheck(control_bits, target_bit, i, n)), axis=1)
    # print("CNOT gate being applied is :", cnot_n)
    return cnot_n


def measure(state, n):
    # convert the state from a bra to a ket
    state = np.transpose(state)
    # convert the complex numbers into their absolute values
    state = np.absolute(state)
    state = np.square(state)
    N = 2**n
    x_values = []
    for i in range(N):
        temp = bin(i)[2:]
        if len(temp) < n:
            for i in range(n - len(temp)):
                temp = "0" + temp
                # print(temp)
        x_values.append(temp)
        # print(x_values)

    x = np.arange(N)
    # print("shape of state array right now is : ")
    # print(state.shape)
    # print(state)
    state = np.array(state).reshape((N,))
    # print("shape of state array right now is : ")
    # print(state.shape)
    # print(x_values)
    plt.bar(x, state)
    plt.xlabel('States')
    plt.ylabel('Probabilities')
    plt.yticks(np.arange(0, 1.25, step=0.25))
    plt.xticks(x, x_values)
    plt.show()


# def measure(state, n):
#     # convert the state from a bra to a ket
#     state = np.transpose(state)
#     # convert the complex numbers into their absolute values
#     state = np.absolute(state)
#     state = np.square(state)
#     state = state*2
#     N = 2 ** n
#     x_values = []
#     for i in range(N):
#         temp = bin(i)[2:]
#         if len(temp) < n:
#             for i in range(n - len(temp)):
#                 temp = "0" + temp
#                 # print(temp)
#         x_values.append(temp)
#         # print(x_values)
#
#     x = np.arange(N)
#     # print("shape of state array right now is : ")
#     # print(state.shape)
#     # print(state)
#     state = np.array(state).reshape((N,))
#     # print("shape of state array right now is : ")
#     # print(state.shape)
#     # print(x_values)
#     plt.bar(x, state)
#     plt.xlabel('States')
#     plt.ylabel('Probabilities')
#     plt.yticks(np.arange(0, 1.25, step=0.25))
#     plt.xticks(x, x_values)
#     plt.show()


def measure_magnitudes(state, n, n_m, counter):
    # measures first "n_m" qubits in an "n" qubit register

    diff = 2**(n-n_m)
    N = 2**n                            # total number of of states
    N_m = 2**n_m                        # number of states possible by n_m qubits

    print("\n***||| function measure_magnitudes has been activated |||***")

    state = np.transpose(state)         # changing the state from vertical vector to a horizontal one
    state = np.absolute(state)
    squared_state = np.array(np.square(state))  # array containing the squared values of the state vector (probabilities)
    squared_state = squared_state.reshape((N,))
    print("shape of squared_state is : ", squared_state.shape, squared_state)   # shape is (1, 2**n_m)
    measured_states = []                # array to store all the states of the first n_m qubits
    magnitudes = []             # stores the modulus values of magnitudes of the corresponding states in measured_states
    combined_probability = 0            # this variable will contain the iterated sum of the all the states with first
                                        # n_m qubits same

    j = 1

    x_values =[]                        # array containing the states of the complete qubit register, not just the one's
                                        # we have to measure


    for i in range(N):
        temp = bin(i)[2:]
        if len(temp) < n:
            for i in range(n - len(temp)):
                temp = "0" + temp
                # print(temp)
        x_values.append(temp)

    print(x_values)

    for i in range(N_m):
        # print("type of x_values", type(x_values))
        x_value_now = str(x_values[i])
        # print("type of x_value_now", type(x_values))
        print("\nvalue of x_value : ", x_value_now)
        print("value of j is : ", j)
        if j == 1:
            measured_states.append(x_value_now[0:n_m])
            print("\nmatrix recording measured states is :", measured_states)
            combined_probability = 0

        combined_probability += squared_state[i]                        # increment combined probability

        j += 1

        if j == diff:
            magnitudes.append(np.sqrt(combined_probability))
            print("\nmatrix recording magnitudes is : ", magnitudes)
            j = 1                                                       # set j back to 1
            print("*** changed j to 1 ***")

    print("\n\nfinally measured states are : ", measured_states)        # measured_states will be the x axis
    print("finally measured mod magnitudes are : ", magnitudes)         # magnitudes will be the y axis
    probabilities = np.square(magnitudes)                               # changing magnitudes to probabilities
    print("corresponding probabilities", probabilities)

    plt.bar(measured_states, probabilities)
    plt.xlabel('States')
    plt.ylabel('Probabilities')
    plt.title("iteration : %d"%counter)                               # optional counter for number of iterations
    plt.yticks(np.arange(0, 1.25, step=0.25))
    plt.xticks(np.arange(N_m), measured_states)
    plt.show()


# if __name__ == "__main__":
#     # test program
#
#     # print(CNOT([1], 2, 3))
#     # a = np.array([0.25, 0.1, 0.35, 0.45])
#     # measure(a, 2)
#
#     a = np.matrix([1, 0, 0, 0])
#     n = 2
#     a = np.transpose(a)
#     print(a)
#     a = np.matmul(np.kron(H(1), I(1)), a)
#     print(a)
#     a = np.matmul(CNOT([1], 2, 2), a)
#     print(a)
#     measure(a, n)
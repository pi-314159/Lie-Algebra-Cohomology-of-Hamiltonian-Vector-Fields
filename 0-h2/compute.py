#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sat Jul 5 18:04:28 2021

@author: PI
"""

import itertools
import time

try:
    from sympy import Matrix
except:
    exit("Please install SymPy!")

startTime = time.time()

def GenPair(length, sumTo, start, end):
    # Generates all combinations of size `length` which sum to `sumTo`.
    temp = list(filter(lambda x: sum(x) == sumTo, itertools.product(range(start, end + 1), repeat = length)))
    possiblePairs = []
    # Creates \mathbb Z \times \mathbb Z grading.
    for a in temp:
        for A in temp:
            b = sorted([[1 - x, 1 - y] for x, y in zip(a, A)], reverse = True)
            # (0, 0) does not exist, so we are not going to consider it.
            # We also do not need repeated pairs.
            if ([0, 0] not in b) and (b not in possiblePairs):
                i = 0
                for B in b:
                    # If a term is repeated, then the wedge product is 0;
                    # thus, we will ignore pairs with repeated terms.
                    if b.count(B) > 1:
                        i = 1
                        break
                if i == 0:
                    possiblePairs.append(b)
    return possiblePairs

def LieBracket(a, b):
    return [a[0] + b[0] - 1, a[1] + b[1] - 1, a[0] * b[1] - a[1] * b[0]]

def Cohomology(a):
    if a == 0:
        return "0"
    elif a == 1:
        return "\mathbb C"
    else:
        return "\mathbb C^" + str(a)

print("C^0:\t1")
print("C^1:\t1")

allPairs = [[[[1, 1]]]]
numTestSubsetA = [1]
numTestSubsetB = [1]

n = 2
i = True

#while n < 8:
while i == True:
    myPair = GenPair(n, 0, -2, 1)
    a = len(myPair)
    if a != 0:
        if n % 2 == 0:
            numTestSubsetB.append(a)
        else:
            numTestSubsetA.append(a)
        allPairs.append(myPair)
        print("C^{0}:\t{1}".format(n, a))
        n += 1
    else:
        i = False
        # If we found all basis elements, the Euler characteristic should be 0.
        if sum(numTestSubsetA) - sum(numTestSubsetB) == 0:
            print("Total:", sum(numTestSubsetA) + sum(numTestSubsetB))
        else:
            exit("ERROR")
        del numTestSubsetA, numTestSubsetB

upToDim = len(allPairs)

imageOfD = 0
differentialList = []
kernelList = []

print("\nH^0 = \mathbb C")
for i in range(1, upToDim):
    k = len(allPairs[i - 1])
    # We create an m by n zero matrix to represent differentials.
    differentials = [[0 for j in itertools.repeat(None, k)] for i in itertools.repeat(None, len(allPairs[i]))]
    powerList = list(range(i + 1))
    for j in range(len(allPairs[i])):
        allPairsij = allPairs[i][j]
        # Applying Fuks' differential formula.
        # `ii` and `jj` represent s and t, respectively.
        for ii in range(i):
            for jj in range(ii + 1, i + 1):
                # `plusMinus` not only indicates the sign, but also means the coefficient.
                plusMinus = 1
                workingList = [*powerList]
                workingList.remove(ii)
                workingList.remove(jj)
                lB = LieBracket(allPairsij[ii], allPairsij[jj])
                plusMinus *= lB.pop()
                # If the coefficient is 0 or the bracket gives \alpha_i = \beta_i = 0,
                # then the corresponding matrix element is 0;
                # hence, we do not need to change anything (remember we generated a zero matrix).
                if plusMinus != 0 and lB != [0, 0]:
                    tempA = [lB]
                    for A in workingList:
                        if allPairsij[A] not in tempA:
                            tempA.append(allPairsij[A])
                        else:
                            # Again, repeated means the corresponding matrix element is 0. 
                            plusMinus = 0
                            break
                    if plusMinus != 0:
                        minusOnePower = ii + jj - 1
                        tempB = sorted(tempA, reverse = True)
                        tempBLen = len(tempB) - 1
                        # We reorder the basis with respect to L(...).
                        # We need to change the sign whenever two basis vectors change places.
                        while tempBLen > 0:
                            tempAIndex = tempA.index(tempB[tempBLen])
                            minusOnePower += tempBLen - tempAIndex
                            del tempA[tempAIndex]
                            tempBLen -= 1
                        # The corresponding matrix element is updated with the new coefficient.
                        differentials[j][allPairs[i - 1].index(tempB)] += (-1) ** minusOnePower * plusMinus
    tempBLen = len(differentials[0])
    differentials = Matrix(differentials)
    kernel = differentials.nullspace()
    tempALen = len(kernel)
    print("H^{d_n:d} = {cohomology}".format(d_n = i, cohomology = Cohomology(tempALen - imageOfD)))
    imageOfD = tempBLen - tempALen
    kernelList.append(kernel)
    differentialList.append(differentials)
print("H^{:d} =".format(upToDim), Cohomology(len(allPairs[upToDim - 1]) - imageOfD))

print("\nBasis elements are saved in the variable `allPairs`, \
      \ndifferentials are saved in the variable `differentialList`, \
      \nand kernels are saved in the variable `kernelList`.")
print("\nTime Elapsed: {:.3f} seconds".format(time.time() - startTime))

#! /usr/bin/env python3
import math

#swap all qubits in register
def swap_all(qc,r):
    for i in range(math.floor(len(r)/2)):
        qc.swap(r[i],r[len(r)-1-i])

#apply controlled global phase
def c_global_phase(cntrl,Q,phase):
    Q.u1(phase,cntrl)

#quantum fourier transform
def qft(qc,r):
    nqubits = len(r)
    for i in reversed(range(nqubits)):
        qc.h(r[i])
        for j in reversed(range(i)):
            qc.cu1(2.0*math.pi/2.0**(i+1-j),r[j],r[i])
    swap_all(qc,r)


#inverse quantum fourier transorm 
def iqft(qc,r):
    nqubits=len(r)
    swap_all(qc,r)
    for i in range(nqubits):
        for j in range(0,i):
            #Add controlled Rk gates here as described by Nielsen
            qc.cu1(-2.0*math.pi/2.0**(i+1-j),r[j],r[i])
        qc.h(r[i])


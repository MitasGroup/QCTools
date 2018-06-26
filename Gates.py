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
def qft(q,Q):
    for i in range(len(q)):
        Q.h(q[i])
        for j in range(i+1,len(q)):
            Q.cu1(2.0*math.pi/2**(j+1),q[j],q[i])
    for i in range(int(math.floor(len(q)/2))):
        Q.swap(q[i],q[int(len(q))-1-i])

def iqft(qc,r):
    nqubits=len(r)
    swap_all(qc,r)
    for i in range(nqubits):
        counter=0
        for j in range(0,i):
            #Add controlled Rk gates here as described by Nielsen
            qc.cu1(-2.0*math.pi/2.0**(i+1-j),r[j],r[i])
        qc.h(r[i])


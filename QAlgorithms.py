#! /usr/bin/env python3
import math
from QCTools import Gates

#phase estimation algorithm with controlled unitary function
def phase_estimation(qs,qr,cr,Q,controlled_unitary,*args):
    #apply hadamard to each readout qubit
    Q.h(qr)

    #apply controlled unitary gates
    for r in range(len(qr)):
        for i in range(2**r):
            controlled_unitary(qs,qr[r],Q,*args)
    #apply inverse QFT to readout qubits
    Gates.qft(Q,qr)

    #measure the readout qubits
    for r in range(len(qr)):
        Q.measure(qr[r],cr[r])
    
    # Potentially we output or return the result here

#iterative phase estimation algorithm with controlled unitary function
def iterative_phase_estimation(qs,qr,cr,Q,controlled_unitary,accuracy,*args):
    phase_bits=[]
    phase_factor=0
    rlen=len(qr) 

    n=accuracy
    while n > 0:
        #apply global phase
        Q.u1(-2*math.pi*phase_factor,qr[0])

        #apply hadamard to each readout qubit
        Q.h(qr)

        #apply controlled unitary gates
        for r in range(rlen):
            for i in range(2**(n-1-r)):
                print(i,r)
                controlled_unitary(qs,qr[r],Q,*args)
        #apply inverse QFT to readout qubits
        qift(qr,Q)

        #measure the readout qubits
        Q.measure(qr,cr)
   
        # Need to add the bit values with maximum count to phase_bits

        #add to global phase factor 
        for i in range(n-1,n-1-rlen,-1):
            phase_factor+=phase_bits[i]/(2**(2+i))

        n-=rlen

        Q.initialize(params,qr,qc) # Will break. I am not yet sure what params is refering to here. Read QISkit SDK reference
    # Potentially we output or return the result here

#! /usr/bin/env python3
import math
from QCTools import Gates, QHelperFunctions
import operator

#phase estimation algorithm with controlled unitary function
def phase_estimation(qs,qr,cs,cr,Q,controlled_unitary,*args):
    #apply hadamard to each readout qubit
    Q.h(qr)

    #apply controlled unitary gates
    for r in range(len(qr)):
        for i in range(2**r):
            controlled_unitary(qs,qr[r],Q,*args)
    #apply inverse QFT to readout qubits
    Gates.iqft(Q,qr)

    #x=4
    #measure the readout qubits
    for r in range(len(qr)):
        Q.measure(qr[r],cr[r])
    
    # Potentially we output or return the result here

#iterative phase estimation algorithm with controlled unitary function
def iterative_phase_estimation(qs0,qr0,cs0,cr0,Q0,backend,shots,rsize,controlled_unitary,accuracy,*args):
    from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
    from qiskit import execute, compile
    import re
    import copy

    bit_str=''
    phase_factor=0
    rlen=len(qr0) 
    slen=len(qs0) 

    print("rlen: ",rlen)
    print("acc: ",accuracy)

    n=accuracy
    while n > 0:
        print('Current phase estimate is  0.'+'x'*(n)+str(bit_str))
        
        # Re-initialize Quantum Registers.
        qr = copy.deepcopy(qr0) 
        qs = copy.deepcopy(qs0) 

        # Create a Classical Registers.
        cr = copy.deepcopy(cr0) 
        cs = copy.deepcopy(cs0) 

        # Re-initialize Quantum Circuit
        Q = copy.deepcopy(Q0) 

        #apply hadamard to each readout qubit
        Q.h(qr)

        #apply controlled unitary gates
        for r in range(rlen):
            #apply global phase
            Q.u1(-2.0**(r+1)*math.pi*phase_factor,qr[r])
            #apply controlled powers of U
            for i in range(2**(n-1-r)):
                controlled_unitary(qs,qr[rlen-1-r],Q,*args)
        #apply inverse QFT to readout qubits
        Gates.iqft(Q,qr)

        #measure the readout qubits
        Q.measure(qr,cr)
  
        #execute job
        job = execute(Q, backend, shots=shots)
        QHelperFunctions.watch_job(job,30)
        result = job.result()

        #show the results
        print("simulation: ", result)
        distribution = {k: v / shots for k, v in result.get_counts().items()}

        #concatenating bit string of maximum count to bit_str
        bit_substr=re.compile('[0-1]+').findall(max(distribution.items(), key=operator.itemgetter(1))[0])[1]
        bit_str=bit_substr+bit_str
  
        print(max(distribution.items(), key=operator.itemgetter(1)))

        #add to global phase factor 
        i=0
        phase_factor=0
        for c in bit_str:
            phase_factor+=float(c)/(2**(rlen+1+i))
            i=i+1

        n=n-rlen

    print('Final phase estimate is  0.'+str(bit_str))

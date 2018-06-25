#! /usr/bin/env python3

#quantum fourier transform
def qft(q,Q):
    for i in range(len(q)):
        Q.h(q[i])
        for j in range(i+1,len(q)):
            Q.cu1(2.0*np.pi/2**(j+1),q[j],q[i])
    for i in range(int(math.floor(len(q)/2))):
        Q.swap(q[i],q[int(len(q))-1-i])

#quantum inverse fourier transform
def qift(q,Q):
    for i in range(int(math.floor(len(q)/2))):
        Q.swap(q[i],q[int(len(q))-1-i])
    for i in reversed(range(len(q))):
        for j in reversed(range(i+1,len(q))):
            Q.cu1(-2.0*np.pi/2**(j+1),q[j],q[i])
        Q.h(q[i])

#phase estimation algorithm with controlled unitary function
def phase_estimation(qs,qr,cr,Q,controlled_unitary,*args):
    #apply hadamard to each readout qubit
    Q.h(qr)

    #apply controlled unitary gates
    for r in range(len(qr)):
        for i in range(2**r):
            print(i,r)
            controlled_unitary(qs,qr[r],Q,*args)
    #apply inverse QFT to readout qubits
    qift(qr,Q)

    #measure the readout qubits
    for r in range(len(qr)):
        Q.measure(qr[r],cr[r])
    
    # Potentially we output or return the result here

#iterative phase estimation algorithm with controlled unitary function
def iterative_phase_estimation(qs,qr,cr,Q,controlled_unitary,accuracy,*args):
    phase_bits=[]
    phase_factor=0
    rlen=len(qr) 

    n=accuracy-1
    while n >= 0:
        #apply hadamard to each readout qubit
        Q.h(qr)

        #apply global phase
        Q.u1(-2*pi*phase_factor,qr[0])

        #apply controlled unitary gates
        for r in range(rlen):
            for i in range(2**n):
                print(i,r)
                controlled_unitary(qs,qr[r],Q,*args)
        #apply inverse QFT to readout qubits
        qift(qr,Q)

        #measure the readout qubits
        for r in range(rlen):
            Q.measure(qr[r],cr[r])
   
        # Need to add the bit values with maximum count to phase_bits

        #add to global phase factor 
        for i in range(n,n-rlen,-1):
            phase_factor+=phase_bits[i]/(2**(2+i))

        n-=rlen

        Q.initialize(params,qr,qc) # Will break. I am not yet sure what params is refering to here. Read QISkit SDK reference
    # Potentially we output or return the result here

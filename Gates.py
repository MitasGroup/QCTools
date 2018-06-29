#! /usr/bin/env python3
import math

#swap all qubits in register
def swap_all(qc,r):
    for i in range(math.floor(len(r)/2)):
        qc.swap(r[i],r[len(r)-1-i])

#apply controlled global phase
def c_global_phase(qs,cntrl,Q,phase):
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

def number_operator(tgt,Q,arg):
    Q.u1(-arg,tgt)

def controlled_number_operator(ctl,tgt,Q,arg):
    Q.cu1(-arg,ctl,tgt)

def excitation_operator(p,q,qs,Q,arg):
    assert(p>q)

    Q.h(qs[p])
    Q.h(qs[q])
    for i in reversed(range(q,p)):
        Q.cx(qs[i+1],qs[i])
    Q.rz(arg,qs[q])
    for i in range(q,p):
        Q.cx(qs[i+1],qs[i])
    Q.h(qs[p])
    Q.h(qs[q])

    Q.rx(-pi/2,qs[p]) 
    Q.rx(-pi/2,qs[q]) 
    for i in reversed(range(q,p)):
        Q.cx(qs[i+1],qs[i])
    Q.crz(arg,ctl,qs[q])
    for i in range(q,p):
        Q.cx(qs[i+1],qs[i])
    Q.rx(-pi/2,qs[p]) 
    Q.rx(-pi/2,qs[q]) 

def controlled_excitation_operator(ctl,p,q,qs,Q,arg):
    assert(p>q)

    Q.ch(ctl,qs[p])
    Q.ch(ctl,qs[q])
    for i in reversed(range(q,p)):
        Q.cx(qs[i+1],qs[i])
    Q.crz(arg,ctl,qs[q])
    for i in range(q,p):
        Q.cx(qs[i+1],qs[i])
    Q.ch(ctl,qs[p])
    Q.ch(ctl,qs[q])

    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[p]) #crx of -pi/2
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[q]) #crx of -pi/2
    for i in reversed(range(q,p)):
        Q.cx(qs[i+1],qs[i])
    Q.crz(arg,ctl,qs[q])
    for i in range(q,p):
        Q.cx(qs[i+1],qs[i])
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[p]) #crx of -pi/2
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[q]) #crx of -pi/2

def coulomb_and_exchange_operator(p,q,qs,Q,arg):
    #Supposed to apply global phase, but cannot be measured. 
    #Matters in PEA algorithm, where global phase gets controled by PEA qubit
    Q.rz(arg,qs[p])
    Q.rz(arg,qs[q])
    Q.cx(qs[p],qs[q])
    Q.rz(arg,qs[q])
    Q.cx(qs[p],qs[q])

def controlled_coulomb_and_exchange_operator(ctl,p,q,qs,Q,arg):
    Q.u1(-arg,ctl) #can use control qubit to apply the global phase
    Q.crz(arg,ctl,qs[p])
    Q.crz(arg,ctl,qs[q])
    Q.cx(qs[p],qs[q])
    Q.crz(arg,ctl,qs[q])
    Q.cx(qs[p],qs[q])

def number_excitation_operator(p,q,r,qs,Q,arg):
    assert(p>q)
    assert(q>r)

    # M = H
    Q.h(qs[p])
    Q.h(qs[r])
    for i in reversed(range(q+1,p)):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q+1],qs[q-1])
    for i in reversed(range(r,q-1)):
        Q.cx(qs[i+1],qs[i])
    Q.rz(arg,qs[r])
    for i in range(r,q-1):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q+1],qs[q-1])
    for i in range(q+1,p):
        Q.cx(qs[i+1],qs[i])
    Q.h(qs[p])
    Q.h(qs[r])

    # M = Rx(-pi/2)
    Q.rx(-pi/2,qs[p])
    Q.rx(-pi/2,qs[r])
    for i in reversed(range(q+1,p)):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q+1],qs[q-1])
    for i in reversed(range(r,q-1)):
        Q.cx(qs[i+1],qs[i])
    Q.rz(arg,qs[r])
    for i in range(r,q-1):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q+1],qs[q-1])
    for i in range(q+1,p):
        Q.cx(qs[i+1],qs[i])
    Q.rx(-pi/2,qs[p])
    Q.rx(-pi/2,qs[r])

def controlled_number_excitation_operator(ctl,p,q,r,qs,Q,arg):
    assert(p>q)
    assert(q>r)

    # M = H
    Q.ch(ctl,qs[p])
    Q.ch(ctl,qs[r])
    for i in reversed(range(q+1,p)):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q+1],qs[q-1])
    for i in reversed(range(r,q-1)):
        Q.cx(qs[i+1],qs[i])
    Q.crz(arg,ctl,qs[r])
    for i in range(r,q-1):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q+1],qs[q-1])
    for i in range(q+1,p):
        Q.cx(qs[i+1],qs[i])
    Q.ch(ctl,qs[p])
    Q.ch(ctl,qs[r])

    # M = Rx(-pi/2)
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[p])
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[r])
    for i in reversed(range(q+1,p)):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q+1],qs[q-1])
    for i in reversed(range(r,q-1)):
        Q.cx(qs[i+1],qs[i])
    Q.crz(arg,ctl,qs[r])
    for i in range(r,q-1):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q+1],qs[q-1])
    for i in range(q+1,p):
        Q.cx(qs[i+1],qs[i])
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[p])
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[r])

def double_excitation_operator(p,q,r,s,qs,Q,arg):
    assert(p>q)
    assert(q>r)
    assert(r>s)
    #Note that Y = Rx(-pi/2)

    #(M1,M2,M3,M4) = (H,H,H,H)
    Q.h(qs[p])
    Q.h(qs[q])
    Q.h(qs[r])
    Q.h(qs[s])
    for i in reversed(range(q,p)):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in reversed(range(s,r)):
        Q.cx(qs[i+1],qs[i])
    Q.rz(arg,qs[s])
    for i in range(s,r):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in range(q,p):
        Q.cx(qs[i+1],qs[i])
    Q.h(qs[p])
    Q.h(qs[q])
    Q.h(qs[r])
    Q.h(qs[s])

    #(M1,M2,M3,M4) = (Y,Y,Y,Y)
    Q.rx(-pi/2,qs[p])
    Q.rx(-pi/2,qs[q])
    Q.rx(-pi/2,qs[r])
    Q.rx(-pi/2,qs[s])
    for i in reversed(range(q,p)):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in reversed(range(s,r)):
        Q.cx(qs[i+1],qs[i])
    Q.rz(arg,qs[s])
    for i in range(s,r):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in range(q,p):
        Q.cx(qs[i+1],qs[i])
    Q.rx(pi/2,qs[p])
    Q.rx(pi/2,qs[q])
    Q.rx(pi/2,qs[r])
    Q.rx(pi/2,qs[s])

    #(M1,M2,M3,M4) = (H,Y,H,Y)
    Q.h(qs[p])
    Q.rx(-pi/2,qs[q])
    Q.h(qs[r])
    Q.rx(-pi/2,qs[s])
    for i in reversed(range(q,p)):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in reversed(range(s,r)):
        Q.cx(qs[i+1],qs[i])
    Q.rz(arg,qs[s])
    for i in range(s,r):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in range(q,p):
        Q.cx(qs[i+1],qs[i])
    Q.h(qs[p])
    Q.rx(pi/2,qs[q])
    Q.h(qs[r])
    Q.rx(pi/2,qs[s])

    #(M1,M2,M3,M4) = (Y,H,Y,H)
    Q.rx(-pi/2,qs[p])
    Q.h(qs[q])
    Q.rx(-pi/2,qs[r])
    Q.h(qs[s])
    for i in reversed(range(q,p)):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in reversed(range(s,r)):
        Q.cx(qs[i+1],qs[i])
    Q.rz(arg,qs[s])
    for i in range(s,r):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in range(q,p):
        Q.cx(qs[i+1],qs[i])
    Q.rx(pi/2,qs[p])
    Q.h(qs[q])
    Q.rx(pi/2,qs[r])
    Q.h(qs[s])

    #(M1,M2,M3,M4) = (Y,Y,H,H)
    Q.rx(-pi/2,qs[p])
    Q.rx(-pi/2,qs[q])
    Q.h(qs[r])
    Q.h(qs[s])
    for i in reversed(range(q,p)):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in reversed(range(s,r)):
        Q.cx(qs[i+1],qs[i])
    Q.rz(arg,qs[s])
    for i in range(s,r):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in range(q,p):
        Q.cx(qs[i+1],qs[i])
    Q.rx(pi/2,qs[p])
    Q.rx(pi/2,qs[q])
    Q.h(qs[r])
    Q.h(qs[s])

    #(M1,M2,M3,M4) = (H,H,Y,Y)
    Q.h(qs[p])
    Q.h(qs[q])
    Q.rx(-pi/2,qs[r])
    Q.rx(-pi/2,qs[s])
    for i in reversed(range(q,p)):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in reversed(range(s,r)):
        Q.cx(qs[i+1],qs[i])
    Q.rz(arg,qs[s])
    for i in range(s,r):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in range(q,p):
        Q.cx(qs[i+1],qs[i])
    Q.h(qs[p])
    Q.h(qs[q])
    Q.rx(pi/2,qs[r])
    Q.rx(pi/2,qs[s])

    #(M1,M2,M3,M4) = (Y,H,H,Y)
    Q.rx(-pi/2,qs[p])
    Q.h(qs[q])
    Q.h(qs[r])
    Q.rx(-pi/2,qs[s])
    for i in reversed(range(q,p)):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in reversed(range(s,r)):
        Q.cx(qs[i+1],qs[i])
    Q.rz(arg,qs[s])
    for i in range(s,r):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in range(q,p):
        Q.cx(qs[i+1],qs[i])
    Q.rx(pi/2,qs[p])
    Q.h(qs[q])
    Q.h(qs[r])
    Q.rx(pi/2,qs[s])

    #(M1,M2,M3,M4) = (H,Y,Y,H)
    Q.h(qs[p])
    Q.rx(-pi/2,qs[q])
    Q.rx(-pi/2,qs[r])
    Q.h(qs[s])
    for i in reversed(range(q,p)):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in reversed(range(s,r)):
        Q.cx(qs[i+1],qs[i])
    Q.rz(arg,qs[s])
    for i in range(s,r):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in range(q,p):
        Q.cx(qs[i+1],qs[i])
    Q.h(qs[p])
    Q.rx(pi/2,qs[q])
    Q.rx(pi/2,qs[r])
    Q.h(qs[s])

def controlled_double_excitation_operator(ctl,p,q,r,s,qs,Q,arg):
    assert(p>q)
    assert(q>r)
    assert(r>s)

    #Note that Y = Rx(-pi/2)

    #(M1,M2,M3,M4) = (H,H,H,H)
    Q.ch(ctl,qs[p])
    Q.ch(ctl,qs[q])
    Q.ch(ctl,qs[r])
    Q.ch(ctl,qs[s])
    for i in reversed(range(q,p)):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in reversed(range(s,r)):
        Q.cx(qs[i+1],qs[i])
    Q.crz(arg,ctl,qs[s])
    for i in range(s,r):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in range(q,p):
        Q.cx(qs[i+1],qs[i])
    Q.ch(ctl,qs[p])
    Q.ch(ctl,qs[q])
    Q.ch(ctl,qs[r])
    Q.ch(ctl,qs[s])

    #(M1,M2,M3,M4) = (Y,Y,Y,Y)
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[p])
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[q])
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[r])
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[s])
    for i in reversed(range(q,p)):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in reversed(range(s,r)):
        Q.cx(qs[i+1],qs[i])
    Q.crz(arg,ctl,qs[s])
    for i in range(s,r):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in range(q,p):
        Q.cx(qs[i+1],qs[i])
    Q.cu3(pi/2,-pi/2,pi/2,ctl,qs[p])
    Q.cu3(pi/2,-pi/2,pi/2,ctl,qs[q])
    Q.cu3(pi/2,-pi/2,pi/2,ctl,qs[r])
    Q.cu3(pi/2,-pi/2,pi/2,ctl,qs[s])

    #(M1,M2,M3,M4) = (H,Y,H,Y)
    Q.ch(ctl,qs[p])
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[q])
    Q.ch(ctl,qs[r])
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[s])
    for i in reversed(range(q,p)):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in reversed(range(s,r)):
        Q.cx(qs[i+1],qs[i])
    Q.crz(arg,ctl,qs[s])
    for i in range(s,r):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in range(q,p):
        Q.cx(qs[i+1],qs[i])
    Q.ch(ctl,qs[p])
    Q.cu3(pi/2,-pi/2,pi/2,ctl,qs[q])
    Q.ch(ctl,qs[r])
    Q.cu3(pi/2,-pi/2,pi/2,ctl,qs[s])

    #(M1,M2,M3,M4) = (Y,H,Y,H)
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[p])
    Q.ch(ctl,qs[q])
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[r])
    Q.ch(ctl,qs[s])
    for i in reversed(range(q,p)):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in reversed(range(s,r)):
        Q.cx(qs[i+1],qs[i])
    Q.crz(arg,ctl,qs[s])
    for i in range(s,r):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in range(q,p):
        Q.cx(qs[i+1],qs[i])
    Q.cu3(pi/2,-pi/2,pi/2,ctl,qs[p])
    Q.ch(ctl,qs[q])
    Q.cu3(pi/2,-pi/2,pi/2,ctl,qs[r])
    Q.ch(ctl,qs[s])

    #(M1,M2,M3,M4) = (Y,Y,H,H)
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[p])
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[q])
    Q.ch(ctl,qs[r])
    Q.ch(ctl,qs[s])
    for i in reversed(range(q,p)):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in reversed(range(s,r)):
        Q.cx(qs[i+1],qs[i])
    Q.crz(arg,ctl,qs[s])
    for i in range(s,r):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in range(q,p):
        Q.cx(qs[i+1],qs[i])
    Q.cu3(pi/2,-pi/2,pi/2,ctl,qs[p])
    Q.cu3(pi/2,-pi/2,pi/2,ctl,qs[q])
    Q.ch(ctl,qs[r])
    Q.ch(ctl,qs[s])

    #(M1,M2,M3,M4) = (H,H,Y,Y)
    Q.ch(ctl,qs[p])
    Q.ch(ctl,qs[q])
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[r])
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[s])
    for i in reversed(range(q,p)):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in reversed(range(s,r)):
        Q.cx(qs[i+1],qs[i])
    Q.crz(arg,ctl,qs[s])
    for i in range(s,r):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in range(q,p):
        Q.cx(qs[i+1],qs[i])
    Q.ch(ctl,qs[p])
    Q.ch(ctl,qs[q])
    Q.cu3(pi/2,-pi/2,pi/2,ctl,qs[r])
    Q.cu3(pi/2,-pi/2,pi/2,ctl,qs[s])

    #(M1,M2,M3,M4) = (Y,H,H,Y)
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[p])
    Q.ch(ctl,qs[q])
    Q.ch(ctl,qs[r])
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[s])
    for i in reversed(range(q,p)):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in reversed(range(s,r)):
        Q.cx(qs[i+1],qs[i])
    Q.crz(arg,ctl,qs[s])
    for i in range(s,r):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in range(q,p):
        Q.cx(qs[i+1],qs[i])
    Q.cu3(pi/2,-pi/2,pi/2,ctl,qs[p])
    Q.ch(ctl,qs[q])
    Q.ch(ctl,qs[r])
    Q.cu3(pi/2,-pi/2,pi/2,ctl,qs[s])

    #(M1,M2,M3,M4) = (H,Y,Y,H)
    Q.ch(ctl,qs[p])
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[q])
    Q.cu3(-pi/2,-pi/2,pi/2,ctl,qs[r])
    Q.ch(ctl,qs[s])
    for i in reversed(range(q,p)):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in reversed(range(s,r)):
        Q.cx(qs[i+1],qs[i])
    Q.crz(arg,ctl,qs[s])
    for i in range(s,r):
        Q.cx(qs[i+1],qs[i])
    Q.cx(qs[q],qs[r])
    for i in range(q,p):
        Q.cx(qs[i+1],qs[i])
    Q.ch(ctl,qs[p])
    Q.cu3(pi/2,-pi/2,pi/2,ctl,qs[q])
    Q.cu3(pi/2,-pi/2,pi/2,ctl,qs[r])
    Q.ch(ctl,qs[s])

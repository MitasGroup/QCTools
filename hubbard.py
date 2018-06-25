#!/usr/bin/env python3

def global_H(qc,q):

    for i in range(len(q)):
        qc.h(q[i])

def adiabatic_evolve(qc,q0,q1,delta,tau,steps):

    for m in range(1,steps+1):

        # Apply Rx(-2*(delta) gate on qubits 0 and 1 
        qc.u3(-2.0*delta,-pi/2.0,pi/2.0,q0)
        qc.u3(-2.0*delta,-pi/2.0,pi/2.0,q1)

        # Apply CNOT from qubit 0 to qubit 1
        qc.cx(q0,q1)
        
        # Apply Rz(-2*(delta) gate on qubit 1
        qc.u1(m*delta**2.0/tau,q1)
        
        # Apply CNOT from qubit 0 to qubit 1
        qc.cx(q0,q1)

def measure(qc,q,c,axis):

    for i in range(len(q)):
        if(axis[i]=='z'):
            pass
        elif(axis[i]=='x'):
            qc.h(q[i])
        elif(axis[i]=='y'):
            qc.sdg(q[i])
            qc.h(q[i])
        else:
            sys.exit(axis[i]+" direction not implemented")
        print("Measuring qubit "+str(i)+" along "+axis[i]+" axis")
    qc.measure(q,c)

def E(probability,qubits):
   
    val = 0.0

    if(qubits[0] and not qubits[1]):
        for k, v in probability.items():
            val = val+v*(-1)**int(k[0])
    elif(qubits[1] and not qubits[0]):
        for k, v in probability.items():
            val = val+v*(-1)**int(k[1])
    else:
        for k, v in probability.items():
            val = val+v*(-1)**int(k[0])*(-1)**int(k[1])

    return val

def lowest_pending_jobs():
    """Returns the backend with lowest pending jobs."""
    list_of_backends = available_backends(
        {'local': False, 'simulator': False})
    device_status = [get_backend(backend).status
                     for backend in list_of_backends]

    best = min([x for x in device_status if x['available'] is True],
               key=lambda x: x['pending_jobs'])
    return best['name']

def watch_job(job):
    lapse = 0
    interval = 10
    while not job_sim.done:
        print('Status @ {} seconds'.format(interval * lapse))
        print(job_sim.status)
        time.sleep(interval)
        lapse += 1

def create_circuit(delta,tau,steps):
    # Create a Quantum Register with 2 qubits.
    q = QuantumRegister(2)
    # Create a Classical Register with 2 bits.
    c = ClassicalRegister(2)
    # Create a Quantum Circuit
    qc = QuantumCircuit(q, c)
    #global_H(qc,q)
    #qc.initialize(q)
    qc.h(q)
    adiabatic_evolve(qc,q[0],q[1],delta,tau,steps)
    return qc,q,c

if __name__ == "__main__":

    import sys, time, getpass

    try:
        # Add path to Qconfig.py if not already present
        # sys.path.append("../../") # go to parent dir
        import Qconfig_qexperience as Qconfig
        qx_config = {
           "APItoken": Qconfig.APItoken,
           "url": Qconfig.config['url'],
           "hub": Qconfig.config['hub'],
           "group": Qconfig.config['group'],
           "project": Qconfig.config['project']}
        print('Qconfig loaded from %s.' % Qconfig.__file__)
    except:
        print('Qconfig.py not found; Qconfig loaded using user input.')
        APItoken = getpass.getpass('Please input your token and hit enter: ')
        qx_config = {
            "APItoken": APItoken,
            "url":"https://quantumexperience.ng.bluemix.net/api"}


    # Import the QISKit SDK
    from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
    from qiskit import available_backends, execute, register, get_backend
    from qiskit.tools.visualization import circuit_drawer
    from math import pi
    import numpy as np

    register(qx_config['APItoken'], qx_config['url'],qx_config['hub'],qx_config['group'],qx_config['project'])#IBM Q-Console

    nshots = 4096
    for i in range(1,11):

        # Method 1
        # U = i
        # steps = i
        # tau = 0.1
        # delta = 0.1

        # Method 2
        U = i*0.5
        steps = 5
        tau = 1.25/U
        delta = 0.25

        # x1
        print("\nPreparing x1 circuit")
        qc,q,c = create_circuit(delta,tau,steps)

        # Measure along X1 and Z2 axis
        measure(qc,q,c,['x','z'])

        # Set backend
        #backend = lowest_pending_jobs()
        #print("the best backend is " + backend)
        backend = 'ibmqx4'

        #circuit_drawer(qc).save("qc.pdf")

        # Compile and run the Quantum circuit on a simulator backend
        job_sim = execute(qc, backend, shots=nshots, max_credits=10)
        watch_job(job_sim)
        sim_result = job_sim.result()

        # Show the results
        print("simulation: ", sim_result)
        distribution = {k: v / nshots for k, v in sim_result.get_counts().items()}
        print(distribution)
        x1_exp = E(distribution,[False,True])

        # x2
        print("\nPreparing x2 circuit")
        qc,q,c = create_circuit(delta,tau,steps)
        # Measure along Z1 and X2 axis
        measure(qc,q,c,['z','x'])

        # Set backend
        backend = lowest_pending_jobs()
        print("the best backend is " + backend)
        backend = 'ibmqx4'

        # Compile and run the Quantum circuit on a simulator backend
        job_sim = execute(qc, backend, shots=nshots, max_credits=10)
        watch_job(job_sim)
        print(job_sim.status)
        sim_result = job_sim.result()
        
        # Show the results
        print("simulation: ", sim_result)
        distribution = {k: v / nshots for k, v in sim_result.get_counts().items()}
        print(distribution)
        x2_exp = E(distribution,[True,False])

        # z1z2
        print("\nPreparing z1z2 circuit")
        qc,q,c = create_circuit(delta,tau,steps)

        # Measure along Z1 and Z2 axis
        measure(qc,q,c,['z','z'])
        backend = lowest_pending_jobs()
        print("the best backend is " + backend)
        backend = 'ibmqx4'

        # Compile and run the Quantum circuit on a simulator backend
        job_sim = execute(qc, backend, shots=nshots, max_credits=10)
        watch_job(job_sim)
        sim_result = job_sim.result()

        # Show the results
        print("simulation: ", sim_result)
        distribution = {k: v / nshots for k, v in sim_result.get_counts().items()}
        print(distribution)
        z1z2_exp = E(distribution,[True,True])

        h_exp = -(x1_exp+x2_exp)+0.5*steps*delta/tau*z1z2_exp

        print("\n <X1> = %12.8f" % (x1_exp))
        print("   <X2> = %12.8f" % (x2_exp))
        print(" <Z1Z2> = %12.8f" % (z1z2_exp))

        print("      U = ",U)
        print("\n  <H> = %12.8f" % (h_exp))
 
    print("\n EXPERIMENT COMPLETE. \n")

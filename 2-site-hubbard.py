#! This should simulate the results from arXiv:1712.08581

from qiskit import QuantumProgram, QISKitError,QuantumJob,QuantumRegister,ClassicalRegister,QuantumCircuit
from qiskit import available_backends, execute, register, get_backend
from qiskit.tools.visualization import circuit_drawer
import numpy as np
import sys,time,getpass

shots = 4096

try:
    sys.path.append('/home/camelto2/Research/quantum_computing/')
    import Qconfig_ibmq_experience as Qconfig
    qx_config = {
            "APItoken": Qconfig.APItoken,
            "url"     : Qconfig.config['url'],
            "hub"     : Qconfig.config['hub'],
            "group"   : Qconfig.config['group'],
            "project" : Qconfig.config['project']
            }
    print('Qconfig loaded from %s.' % Qconfig.__file__ )
    register(qx_config['APItoken'],qx_config['url'],qx_config['hub'],qx_config['group'],qx_config['project'])
except:
    APItoken = getpass.getpass('Please input your token and hit enter: ')
    qx_config = {
            "APItoken": APItoken,
            "url": "https://quantumexperience.ng.bluemix.net/api"
            }
    print("Qconfig not found. Loaded from user input")

def lowest_pending_jobs():
    """Returns the backend with lowest pending jobs."""
    list_of_backends = available_backends(
        {'local': False, 'simulator': False})
    device_status = [get_backend(backend).status
                     for backend in list_of_backends]

    best = min([x for x in device_status if x['available'] is True],
                key=lambda x: x['pending_jobs'])
    return best['name']

def exact(U):
    return 0.5*(U-np.sqrt(16+U**2))-U/2

def setup():
    q = QuantumRegister(2)
    c = ClassicalRegister(2)
    Q = QuantumCircuit(q,c)
    return q,c,Q

def ground_state(q,Q):
    Q.h(q) #applies to each one

def adiabatic_evolution(q,Q,delta,tau,Nsteps):
    #Adiabatic evolution
    for m in range(1,Nsteps+1):
        #apply exp(i delta (X1+X2) ) 
        for i in range(len(q)):
            Q.rx(-2*delta,q[i])
    
        #apply exp(-i m delta^2/ (2tau) Z1 Z2) 
        Q.cx(q[0],q[1])
        Q.rz(m*delta**2/tau,q[1])
        Q.cx(q[0],q[1])

def measure(q,c,Q,measurement):
    if measurement == 'X1':
        Q.h(q[0])
    elif measurement == 'X2':
        Q.h(q[1])
    Q.measure(q,c)

    #backend = lowest_pending_jobs()
    #print("the best backend is " + backend)

    #backend = 'local_qasm_simulator'
    backend = 'ibmqx4'
    job = execute(Q,backend=backend,shots=shots)

    lapse = 0
    interval = 10
    while not job.done:
        print('Status @ {} seconds'.format(interval*lapse))
        print(job.status)
        time.sleep(interval)
        lapse += 1
    print(job.status)

    result = job.result()
    return result

def analyze(H1,H2,H3,U):
    states=['00','10','01','11']
    scaledH1 = { k: v/shots for k,v in H1.get_counts().items()}
    scaledH2 = { k: v/shots for k,v in H2.get_counts().items()}
    scaledH3 = { k: v/shots for k,v in H3.get_counts().items()}

    for state in states:
        if state not in scaledH1:
            scaledH1[state] = 0.0
        if state not in scaledH2:
            scaledH2[state] = 0.0
        if state not in scaledH3:
            scaledH3[state] = 0.0

    X1 = scaledH1['00']+scaledH1['10'] - (scaledH1['01']+scaledH1['11'])
    X2 = scaledH2['00']+scaledH2['01'] - (scaledH2['10']+scaledH2['11'])
    Z1Z2 = scaledH3['00'] + scaledH3['11'] - scaledH3['01'] - scaledH3['10']

    return -(X1+X2)+0.5*U*Z1Z2


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import pickle
    data = {'M1': [], 'M2': []}
    #Us = [0,1,2,3,4,5,6]
    Us = [0,2,4,6]
    data['U'] = Us
    for U in Us:
        q,c,Q = setup() #for measuring <X1>
        ground_state(q,Q)
        adiabatic_evolution(q,Q,0.1,0.1,U)
        res1 = measure(q,c,Q,'X1')
        q,c,Q = setup() #for measuring <X1>
        ground_state(q,Q)
        adiabatic_evolution(q,Q,0.1,0.1,U)
        res2 = measure(q,c,Q,'X2')
        q,c,Q = setup() #for measuring <X1>
        ground_state(q,Q)
        adiabatic_evolution(q,Q,0.1,0.1,U)
        res3 = measure(q,c,Q,'Z1Z2')
        E = analyze(res1,res2,res3,U)
        data['M1'].append(E)

    for U in Us:
        q,c,Q = setup() #for measuring <X1>
        ground_state(q,Q)
        if U != 0:
            adiabatic_evolution(q,Q,0.25,1.25/U,5)
        res1 = measure(q,c,Q,'X1')
        q,c,Q = setup() #for measuring <X1>
        ground_state(q,Q)
        if U != 0:
            adiabatic_evolution(q,Q,0.25,1.25/U,5)
        res2 = measure(q,c,Q,'X2')
        q,c,Q = setup() #for measuring <X1>
        ground_state(q,Q)
        if U != 0:
            adiabatic_evolution(q,Q,0.25,1.25/U,5)
        res3 = measure(q,c,Q,'Z1Z2')
        E = analyze(res1,res2,res3,U)
        data['M2'].append(E)

    plt.scatter(data['U'],data['M1'],label='M1')
    plt.scatter(data['U'],data['M2'],label='M2')
    plt.plot(np.linspace(0,6,100),exact(np.linspace(0,6,100)),label='exact')
    plt.legend()
    plt.show()
    with open('data.p','wb') as f:
        pickle.dump(data,f)

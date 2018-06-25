
def lowest_pending_jobs():
    from qiskit import available_backends, get_backend
    list_of_backends = available_backends(
            {'local': False, 'simulator': False})
    device_status = [get_backend(backend).status
            for backend in list_of_backends]

    best = min([ x for x in device_status if x['available'] is True],
            key = lambda x: x['pending_jobs'])
    return best['name']

def watch_job(job,interval=10):
    lapse = 0
    while not job.done:
        print('Status @ {} seconds'.format(interval * lapse))
        print(job.status)
        time.sleep(interval)
        lapse += 1

def configure(Qconfig_path,Qconfig_name):
    import time,sys,getpass
    from qiskit import register
    try:
        sys.path.append(Qconfig_path)
        Qconfig = __import__(Qconfig_name)
        qx_config = {
                "APItoken": Qconfig.APItoken,
                "url"     : Qconfig.config['url'],
                "hub"     : Qconfig.config['hub'],
                "group"   : Qconfig.config['group'],
                "project" : Qconfig.config['project']
                }
        print('Qconfig loaded from %s. ' % Qconfig.__file__)
    except:
        APItoken = getpass.getpass('Please input your token and hit enter: ')
        qx_config = {
                "APItoken": APItoken,
                "url"     : "https://quantumexperience.ng.bluemix.net/api",
                "hub"     : None,
                "group"   : None,
                "project" : None,
                }
        print("Qconfig not found. Loaded from user input")

    register(qx_config['APItoken'],qx_config['url'],qx_config['hub'],qx_config['group'],qx_config['project'])


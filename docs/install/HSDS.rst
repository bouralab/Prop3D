===================================
Highly Scalable Data Service (HSDS)
===================================

**Note: If you only want to use the dataset, this is unnecessary. See instructions on how to use the** :doc:`../precomputed_datasets`

*This documentation has been modified from `HOBO-request <>` documentation by LF Murillo.*

*Accompanying files can be found in the* ``docs/k3s-conf`` *directory after cloning Prop3D* (``git clone https://github.com/bouralab/Prop3D``)


You can set up HSDS on any cloud platform or a single machine using Docker or on a cluster using Kubernetes (or AKS on Microsoft Azure). If you are not using a single machine, please follow the official `Kubernetes instruction instructions <https://kubernetes.io/docs/setup/>`_.

Single Machine Setup
--------------------

For single machine setup, please clone the `HSDS <https://github.com/HDFGroup/hsds>`_ repository.

Running HSDS (Highly Scalable Data Service) on K3s
++++++++++++++++++++++++++++++++++++++++++++++++++

This document describes the installation of HSDS on K3s, 
a certified Kubernetes distribution. In our example, we
will install HSDS on a single host running K3s for testing
purposes, but you can run it on as many hosts as you need.

Requirements
++++++++++++

K3s is pre-packaged for various Linux distributions. Before
getting started, `make sure you satisfy the requirements for 
k3s with your distro <https://rancher.com/docs/k3s/latest/en/advanced>`_
which include cgroups, legacy iptables, and more.

Our example is based on k3s version `1.21.5+k3s2`.

Installing K3s 
++++++++++++++

If your distro does not have K3s packaged, you can run the 
installation script yourself with the following command:

.. code-block:: bash

    curl -sfL https://get.k3s.io | sh -


After it finishes the script, you should have k3s up and running
as a systemd service:

.. code-block:: bash

    k3s.service - Lightweight Kubernetes
    Loaded: loaded (/etc/systemd/system/k3s.service; enabled; vendor preset: enab
    Active: active (running) since Wed 2021-11-03 23:41:55 EDT; 2min 21s ago
        Docs: https://k3s.io
    Process: 2446 ExecStartPre=/bin/sh -xc ! /usr/bin/systemctl is-enabled --quiet
    Process: 2448 ExecStartPre=/sbin/modprobe br_netfilter (code=exited, status=0/
    Process: 2449 ExecStartPre=/sbin/modprobe overlay (code=exited, status=0/SUCCE
    Main PID: 2451 (k3s-server)
        Tasks: 99
    Memory: 1.3G
    CGroup: /system.slice/k3s.service
    ...


You can check if everything is ready for hsds installation with kubectl, as one
would for the official distribution of Kubernetes:

.. code-block:: bash

    sudo kubectl get nodes

In our case, we will only get one node listed as we are just running a single-node
install for testing purposes:

.. code-block:: bash

    NAME             STATUS   ROLES                  AGE     VERSION
    $YOUR_HOSTNAME   Ready    control-plane,master   5m55s   v1.21.5+k3s2

Looking good! Time to continue on to deploy HSDS.

Deploying HSDS (Highly Scalable Data Service) 
+++++++++++++++++++++++++++++++++++++++++++++

For our test, we will deploy HSDS with 8 data nodes (pods) on a regular
filesystem (described as a POSIX storage install in the HSDS documentation).

Instead of using the yml distributed with HSDS, you will need use to
deployment files that were modified to work with K3s. You will find these
files in the directory ``doc/k3s-conf`` in the root of the ``hobo-connect``
repository on Gitlab.

HSDS has been 'containerized' for deployment with Kubernetes, so we can skip
the step of preparing an application to be deployed.

Now, we can proceed to the first step, which is to apply the configuration
for 'role-based access control' (RBAC): 

.. code-block:: bash    

    sudo kubectl apply -f k8s_rbac.yml

This will enable HSDS pods to "find each other" so that requests can 
be accelerated by distributing reads and writes across multiple pods.

Next, edit the file ``override.yml`` with the parameters that you need and
proceed to run a command that will create a 'configuration map' that can 
be used by HSDS. This is the approach that Kubernetes uses to separate 
configurations that are specific to your use-case and environment from
your container images (which have the standalone application only):

.. code-block:: bash

    sudo kubectl create configmap hsds-config --from-file=config.yml --from-file=override.yml


HSDS can use a password file to authenticate users using the HTTP Basic Auth protocol 
(authentication using OpenID or Azure Active Directory is also supported).
To construct a password file, create a text file like the following (hopefully using more
secure passwords!):

.. code-block:: 

    admin:admin
    test_user1:test
    test_user2:test

Each line in the file is in the format &lt;username&gt;:&lt;password&gt;.  You'll need the 
admin user for performing certain tasks like setting up top-level folders.

Once you have the password file, run:

.. code-block:: bash

    kubectl create secret generic user-password --from-file=<passwd_file>


You can always check to see if everything has been loaded properly, but you really do 
not have to. Using ``kubectl`` you can ``get`` (request) info about configmap:

.. code-block:: bash

    sudo kubectl get configmap

The output should include the configmap we just loaded:

.. code-block:: bash

    NAME               DATA   AGE
    kube-root-ca.crt   1      61m
    hsds-config        2      21s

Awesome, just a few more commands and we are done!
 
We need to configure storage by 'claiming a persistent volume' (PVC).
Before running the command, edit the file `hsds_pvc.yml` and set the disk
space that you want to use. Then, proceed with the following command:

.. code-block:: bash

    sudo kubectl apply -f hsds_pvc.yml

K3s will output ``persistentvolumeclaim/hsds-pvc created``, so you know you
are onto something good.

Now, let's run two more commands for deploying the HSDS container and
expose the HSDS service in the cluster, respectively:

.. code-block:: bash

    sudo kubectl apply -f hsds_deployment.yml
    sudo kubectl apply -f hsds_service.yml

You may want to check if everything is good on the K3s side:

You will see that it will change its ``STATUS`` from ``ContainerCreating``
to ``Running``:

.. code-block:: bash

    NAME                    READY   STATUS    RESTARTS   AGE
    hsds-857754bf58-p2n8b   2/2     Running   0          103s


You may want to look into the service that has been provisioned as well:

.. code-block:: bash

    sudo kubectl get services

And see something along the lines of this output:

.. code-block:: bash

    NAME         TYPE        CLUSTER-IP      EXTERNAL-IP   PORT(S)          AGE
    kubernetes   ClusterIP   10.43.0.1       <none>        443/TCP          78m
    hsds         NodePort    10.43.173.154   <none>        5101:32613/TCP   3m

You can try sending a request to the service like so:

.. code-block:: bash

    $ curl http://10.43.173.154:5101/about


Replace the IP above with the cluster IP value for the HSDS service.

If you see a JSON response with a "status" key of "READY",
congratulations, you have HSDS up and running on K3s! 

If not, you can review the HSDS logs to see what the problem is.  Each pod
has two containers named "sn" and "dn" which support frontend and backend aspects of the HSDS service.

To display the log for the sn container you can run:

.. code-block:: bash

    kubectl logs hsds-857754bf58-p2n8b sn

Where hsds-8577... is the pod id.  Similarly use "dn" to see the dn logs.

Since HSDS can be rather chatty, it can be useful to filter by just ERROR or WARN
entries:

.. code-block:: bash

    kubectl logs hsds-857754bf58-p2n8b sn | grep ERROR

You can tweak the number of pods (instances) of HSDS, scaling it up
or down as needed with the following command. For our tests, we will
run 8 replicas:

.. code-block:: bash

    sudo kubectl scale --replicas=8 deployment/hsds

The number of HSDS pods you will be able to create depends on the amount 
of memory available on the machine.  If you see one or more pods that stay
in "Pending" status, it's likely there's not sufficient system resources to 
support that number of pods and you'll need to scale down a bit.


But... before we say good-bye, make sure you configure the users of your 
test instance. Below, we will create the directory for the 'admin' user:

.. code-block:: bash

    pip install h5pyd
    hstouch -u admin -p admin_password -e http://<ip>:5101 /home/
    # run the following for each user who will need a "home" folder:
    hstouch -u admin -p admin_password -e http://<ip>:5101 -o <username> /home/<username>

Now each user who will be interacting with the system can run: hsconfigure.
They will be prompted for server endpoint, username, and password.  Information will
be stored in a file ``.hscfg`` in their home directory.  This will be used to authenticate 
with the server when using tools like: ``hsinfo``, ``hsls``, ``hsload``, etc., and also when 
using h5pyd in Python scripts.


Happy K3s + HSDS testing!

--
Sign-off: LF Murilo, 11-03-2021

Recommended HSDS changes
++++++++++++++++++++++++

.. code-block:: yaml

    dn_ram: 6g
    sn_ram: 6g
    max_tcp_connections: 1000
    max_task_count: 1000
    aio_max_pool_connections: 264
    metadata_mem_cache_size: 1g
    chunk_mem_cache_size: 1g
    data_cache_size: 1g
    timeout: 120


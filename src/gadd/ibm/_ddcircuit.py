#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import numpy as np
from qiskit import (QuantumCircuit, IBMQ, execute, transpile,
                    schedule as build_schedule)

class DdCircuit(QuantumCircuit):
    """
    DdCircuit is class for Dynamical Decoupling Circuits.
    Inherits all mehods/ valid data from qiskit's QuantumCircuit
    """

#TODO: See why this causes "more than one backend matches criteria bug
#    def load_backend(name='ibmq_armonk'):
#        # obtain backend and get useful properties
#        provider = IBMQ.load_account()
#        self.backend = provider.get_backend(name)
#        self.inst_map = self.backend.defaults().instruction_schedule_map
#        # list of backend methods for possible future use
#        #backend.properties()
#        #backend.configuration()
#        #backend_configuration.dt
#        #backend_defaults()
#        #backend.job_limit()
#        #assert backend_configuration().open_pulse

    def set_backend(self, name='ibmq_armonk'):
        """
        Sets backend member data to (name) by actually loading 
        """
        provider = IBMQ.load_account()
        self.backend = provider.get_backend(name)
        return

    def get_backend(self):
        return self.backend

    def gate_frev(self, n=1):
        """
        Appends free evolution to circuit

        Inputs:
        -------------------------
        n -- int, number of identity gates to freely evolve

        Retuns:
        -------------------------
        ng -- int, number of gates appended

        Result:
        -------------------------
        Appends n identity gates to self
        """
        ng = 0
        for i in range(n):
            self.barrier(0)
            self.id(0)
            ng += 1

        self.barrier(0)

        return ng

    def gate_xy4(self, ncyc=1, ni=0):
        """
        Appends xy4 sequence to circuit

        Inputs:
        -------------------------
        ncyc -- int, number of times to cycle through xy4 sequence
        ni -- int, number of identity gates between pauli pulses

        Returns
        -------------------------
        ng -- int, number of gates in appended

        Result:
        -------------------------
        Applies Y I...I X I...I Y I...I X I...I sequence ncyc times with ni identities between pauli pulses
        """
        ng = 0
        for i in range(ncyc):
            self.barrier(0)
            self.y(0)
            ng += 1
            for j in range(ni):
                self.barrier(0)
                self.id(0)
                ng += 1

            self.barrier(0)
            self.x(0)
            ng += 1
            for j in range(ni):
                self.barrier(0)
                self.id(0)
                ng += 1
            self.barrier(0)
            self.y(0)
            ng += 1
            for j in range(ni):
                self.barrier(0)
                self.id(0)
                ng += 1
            self.barrier(0)
            self.x(0)
            ng += 1
            for j in range(ni):
                self.barrier(0)
                self.id(0)
                ng += 1

        self.barrier(0)
        return ng

    def gate_xy4test(self, ncyc=1, ni=0):
        """
        Appends xy4 sequence to circuit the same as gate_xy4 BUT uses 
        native gates directlry. 

        NOTE: I just wanted to see if for some reason QC gave better results
        while using native gates... the answer is no which makes sense with
        transpiling doing the work anyway.
        """
        ng = 0
        for i in range(ncyc):
            self.barrier(0)
            self.add_y(0)
            ng += 1
            for j in range(ni):
                self.barrier(0)
                self.id(0)
                ng += 1
            self.barrier(0)
            self.add_x(0)
            ng += 1
            for j in range(ni):
                self.barrier(0)
                self.id(0)
                ng += 1
            self.barrier(0)
            self.add_y(0)
            ng += 1
            for j in range(ni):
                self.barrier(0)
                self.id(0)
                ng += 1
            self.barrier(0)
            self.add_x(0)
            ng += 1
            for j in range(ni):
                self.barrier(0)
                self.id(0)
                ng += 1

        self.barrier(0)
        return ng

    def gate_xz4(self, ncyc=1, ni=0):
        """
        Appends xz4 sequence. Z is applied much faster, so I thought I would
        try it instead of y and see if results improve.
        """
        ng = 0
        for i in range(ncyc):
            self.barrier(0)
            self.z(0)
            ng += 1
            for j in range(ni):
                self.barrier(0)
                self.id(0)
                ng += 1
            self.barrier(0)
            self.x(0)
            ng += 1
            for j in range(ni):
                self.barrier(0)
                self.id(0)
                ng += 1
            self.barrier(0)
            self.z(0)
            ng += 1
            for j in range(ni):
                self.barrier(0)
                self.id(0)
                ng += 1
            self.barrier(0)
            self.x(0)
            ng += 1
            for j in range(ni):
                self.barrier(0)
                self.id(0)
                ng += 1

        self.barrier(0)
        return ng

    # Adding functions which implement all the gates we need in terms of the
    # IBM native gates. These were done with armonk as of 1 May 2020 using
    # transpile(circ, backend) for some basic gates.
    def add_x(self, qidx):
        """
        Appends x gate to circuit on qubit qidx in terms of native u3 gate.
        """
        self.u3(np.pi, 0, np.pi, qidx)
        return

    def add_y(self, qidx):
        """
        Appends y gate to circuit on qubit qidx in terms of native u3 gate.
        """
        self.u3(np.pi, np.pi/2, np.pi/2, qidx)
        return

    def add_z(self, qidx):
        """
        Appends z gate to circuit on qubit qidx in terms of native u1 gate.
        """
        #NOTE:  u1(pi) = u3(0, pi, 0)
        self.u1(np.pi, qidx)
        return

    def add_xb(self, qidx):
        """
        Appends x bar to circ on qubit qidx in terms of native gates.
        xbar comes from RGAn sequences and is pi rotation of x around
        all axes.
        """
        #NOTE: u1(pi) = u3(0, pi, 0) = z
        self.u1(np.pi, qidx)
        return

    def add_yb(self, qidx):
        """
        Appends y bar to circ on qubit qidx in terms of native gates.
        ybar comes from RGAn sequences and is pi rotation of y around
        all axes.
        """
        #NOTE: u1(3pi) = u1(pi) = z
        self.u1(3*np.pi, qidx)
        return 

    def add_zb(self, qidx):
        """
        Appends z bar to circ on qubit qidx in terms of native gates.
        zbar comes from RGAx sequences and is pi rotation of z around
        all axes.
        """
        #NOTE: u3(pi, 0, pi) = x
        self.u3(np.pi, 0, np.pi, qidx)
        return

    def add_H(self, qidx):
        """
        Appends hadamard to circ on qubit qidx in terms of native u2 gate.
        """
        #NOTE: u2(0, pi) = u3(pi/2, 0, pi)
        self.u2(0, np.pi, qidx)
        return

    

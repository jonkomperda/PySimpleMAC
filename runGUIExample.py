#! /usr/bin/env python
"""
This script behaves as an executable.
It's intent is to to be run from the terminal, and it will launch the exampleGUI script
which runs an interactive demonstration of PySimpleMAC.
.. warning::
This script as well as exampleGUI require MPI4PY as well as a working build of OpenMPI 
or MPICH on the users system. You may need to edit the path in this script to the correct
location of 'mpirun' on your system.
"""
import os

os.system('/usr/local/mpifort/bin/mpirun -n 2 python exampleGUI.py')
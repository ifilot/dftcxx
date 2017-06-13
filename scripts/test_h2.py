#!/usr/bin/env python2

from PyQuante import Molecule
from PyQuante.dft import dft

h2 = Molecule('h2', [(1,(0,0,0.143885)),(1,(0,0,0.856115))], units="Angstrom")
en,orbe,orbs = dft(h2, functional="LDA", basis="sto-3g")

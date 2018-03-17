#!/usr/bin/env python

from PyQuante import Molecule
from PyQuante.dft import dft

h2 = Molecule('h2', [(1,(0,0,-.36700000000000000000)),(1,(0,0,0.36700000000000000000))], units="Angstrom")
en,orbe,orbs = dft(h2, functional="LDA", basis="sto-3g")

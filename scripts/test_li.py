#!/usr/bin/env python

from PyQuante import Molecule
from PyQuante.dft import dft

he = Molecule('lih', [(3,(0,0,0.403635)), (1,(0,0,-1.210905))], units="Bohr")
en,orbe,orbs = dft(he, functional="LDA", basis="sto-3g")

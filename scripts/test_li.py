#!/usr/bin/env python

from PyQuante import Molecule
from PyQuante.dft import dft

lih = Molecule('LiH', [(1,(0,0,-1.210905)), (3,(0,0,0.403635))], units="Bohr")
en,orbe,orbs = dft(lih, functional="LDA", basis="sto-3g")

#!/usr/bin/env python

from PyQuante import Molecule
from PyQuante.dft import dft

n2 = Molecule('n2', [(7,(0,0,0)),(7,(0,0,1.097600))], units="Angstrom")
en,orbe,orbs = dft(n2, functional="LDA", basis="sto-3g")

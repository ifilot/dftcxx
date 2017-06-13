#!/usr/bin/env python2

from PyQuante import Molecule
from PyQuante.dft import dft

he = Molecule('he', [(2,(0,0,0))], units="Angstrom")
en,orbe,orbs = dft(he, functional="LDA", basis="sto-3g")

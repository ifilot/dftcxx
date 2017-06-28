#!/usr/bin/env python2

from PyQuante import Molecule
from PyQuante.dft import dft

be = Molecule('be', [(4,(0,0,0))], units="Angstrom")
en,orbe,orbs = dft(be, functional="LDA", basis="sto-3g")

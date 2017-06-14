#!/usr/bin/env python

from PyQuante import Molecule
from PyQuante.dft import dft

h2o = Molecule('h2o', [(8,(0,0,0)),(1,(0.7570,0.5860,0.0)),(1,(-0.757,0.5860,0.0))], units="Bohr")
en,orbe,orbs = dft(h2o, functional="LDA", basis="sto-3g")

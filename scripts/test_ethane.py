#!/usr/bin/env python

from PyQuante import Molecule
from PyQuante.dft import dft

coord = [
(1,(0.0998334, 0.995004, -0.6)),
(1,(0.911616, 0.411044, 1.6)),
(1,(0.811782, -0.58396, -0.6)),
(1,(-0.0998334, -0.995004, 1.6)),
(1,(-0.911616, -0.411044, -0.6)),
(1,(-0.811782, 0.583961, 1.6)),
(6, (0, 0, 0)),
(6, (0, 0, 1))
]

ethane = Molecule('ethane', coord, units="Angstrom")
en,orbe,orbs = dft(ethane, functional="LDA", basis="sto-3g")

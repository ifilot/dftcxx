/**************************************************************************
 *   dftcxx.cpp  --  This file is part of DFTCXX.                         *
 *                                                                        *
 *   Copyright (C) 2016, Ivo Filot                                        *
 *                                                                        *
 *   DFTCXX is free software:                                             *
 *   you can redistribute it and/or modify it under the terms of the      *
 *   GNU General Public License as published by the Free Software         *
 *   Foundation, either version 3 of the License, or (at your option)     *
 *   any later version.                                                   *
 *                                                                        *
 *   DFTCXX is distributed in the hope that it will be useful,            *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty          *
 *   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 *   See the GNU General Public License for more details.                 *
 *                                                                        *
 *   You should have received a copy of the GNU General Public License    *
 *   along with this program.  If not, see http://www.gnu.org/licenses/.  *
 *                                                                        *
 **************************************************************************/

#include <chrono>
#include <boost/format.hpp>

#include "molecule.h"
#include "dft.h"

int main(int argc, char** argv) {

    if(argc != 2) {
        std::cerr << "Please specify an input file" << std::endl;
        exit(-1);
    }

    auto start = std::chrono::system_clock::now(); //toc

    auto mol = std::make_shared<Molecule>(argv[1]);

    DFT dft;
    dft.add_molecule(mol);
    dft.scf();

    auto end = std::chrono::system_clock::now(); //toc
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << boost::format("Total elapsed time: %f ms\n") % elapsed.count();

    return 0;
}

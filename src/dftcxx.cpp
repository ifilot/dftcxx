/**************************************************************************
 *   This file is part of DFTCXX.                                         *
 *                                                                        *
 *   Author: Ivo Filot <ivo@ivofilot.nl>                                  *
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
#include <tclap/CmdLine.h>

#include "molecule.h"
#include "dft.h"
#include "config.h"

int main(int argc, char** argv) {

    try {
        TCLAP::CmdLine cmd("Perform DFT calculation.", ' ', PROGRAM_VERSION);

        // input filename
        TCLAP::ValueArg<std::string> arg_input_filename("i","input","Input file (i.e. h2.in)",true,"h2.in","filename");
        cmd.add(arg_input_filename);

        cmd.parse(argc, argv);

        const std::string input_filename = arg_input_filename.getValue();

        std::cout << "--------------------------------------------------------------" << std::endl;
        std::cout << std::endl;
        std::cout << "Executing DFTCXX v." << PROGRAM_VERSION << std::endl;
        std::cout << "Author: Ivo Filot <ivo@ivofilot.nl>" << std::endl;
        std::cout << std::endl;
        std::cout << "--------------------------------------------------------------" << std::endl;
        std::cout << std::endl;

        auto start = std::chrono::system_clock::now();
        auto mol = std::make_shared<Molecule>(input_filename);

        DFT dft;
        dft.add_molecule(mol);
        dft.set_hartree_evaluation(DFT::BECKE_GRID);
        dft.scf();

        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << boost::format("Total elapsed time: %f ms\n") % elapsed.count();

        return 0;

    } catch (TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() <<
                     " for arg " << e.argId() << std::endl;
        return -1;
    }
}

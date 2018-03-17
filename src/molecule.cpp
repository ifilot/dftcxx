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

#include "molecule.h"

Atom::Atom() :
    atnr(0),
    position(vec3(0,0,0))
{
}

Atom::Atom(unsigned int _atnr, const vec3& _position):
    atnr(_atnr),
    position(_position) {
}

Molecule::Molecule(const std::string& filename) {
    this->read_molecule_from_file(filename);
}

void Molecule::read_molecule_from_file(const std::string& filename) {
    std::cout << "           Reading input file           " << std::endl;
    std::cout << "========================================" << std::endl;

    static const double angstrom_to_bohr = 1.889725989;

    std::cout << "Reading file:\t\t" << filename << std::endl;
    std::cout << std::endl;

    std::ifstream infile(filename);
    std::string line;

    std::getline(infile, line); // grab name
    std::cout << "Calculating system:\t" << line << std::endl;

    std::getline(infile, line); // grab basis set
    std::string basis_set = line;
    std::cout << "Using basis set:\t" << basis_set << std::endl;

    bool unit_angstrom = false;
    std::getline(infile, line); // length unit
    if(line.compare("angstrom") == 0) {
        unit_angstrom = true;
    }

    std::getline(infile, line); // grab number of atoms
    unsigned int nratoms = boost::lexical_cast<unsigned int>(line);
    std::cout << "Atoms in system:\t" << nratoms << std::endl;

    std::vector<std::string> pieces;
    for(unsigned int i=0; i<nratoms; i++) {
        std::getline(infile, line);
        boost::split(pieces, line, boost::is_any_of(" \t"), boost::token_compress_on);
        unsigned int atnr = this->get_atom_number_from_string(pieces[0]);
        double x = boost::lexical_cast<double>(pieces[1]);
        double y = boost::lexical_cast<double>(pieces[2]);
        double z = boost::lexical_cast<double>(pieces[3]);

        if(unit_angstrom) {
            x *= angstrom_to_bohr;
            y *= angstrom_to_bohr;
            z *= angstrom_to_bohr;
        }

        this->add_atom(Atom(atnr, vec3(x, y, z) ) );
    }

    infile.close();

    this->set_basis_set(basis_set);
    std::cout << "========================================" << std::endl;

    for(unsigned int i=0; i<nratoms; i++) {
        std::cout << boost::format("%i  %12.6f  %12.6f  %12.6f")
                     % this->atoms[i]->get_charge()
                     % this->atoms[i]->get_position()[0]
                     % this->atoms[i]->get_position()[1]
                     % this->atoms[i]->get_position()[2] << std::endl;
    }

    std::cout << "========================================" << std::endl;
    std::cout << std::endl;
}

void Molecule::set_basis_set(const std::string& basis_set) {
    // define origin vector
    static const vec3 origin(0,0,0);

    // collect for which atom types we need to fetch basis
    // set information
    unsigned int highest_atom = 0;
    for(unsigned int i=0; i<this->atoms.size(); i++) {
        highest_atom = std::max(highest_atom, this->atoms[i]->get_charge());
    }

    // open the file
    std::ifstream infile("../" + basis_set);
    std::string line;
    std::vector<CGF> bs_cgfs;
    while(std::getline(infile, line)) {

        // skip comment lines
        if(line[0] == '#') {
            continue;
        }

        std::vector<std::string> pieces;
        boost::split(pieces, line, boost::is_any_of(" \t"), boost::token_compress_on);
        unsigned int atnr =  boost::lexical_cast<unsigned int>(pieces[0]);
        unsigned int ncgfs = boost::lexical_cast<unsigned int>(pieces[1]);
        bs_cgfs.clear();

        // collect cgfs
        for(unsigned int i=0; i<ncgfs; i++) {
            std::getline(infile, line);
            boost::split(pieces, line, boost::is_any_of(" \t"), boost::token_compress_on);
            char type = pieces[0][0];
            unsigned int ngtos = boost::lexical_cast<unsigned int>(pieces[1]);

            // allocate cgfs in vector
            switch(type) {
                case 'S':
                    bs_cgfs.resize(bs_cgfs.size() + 1);
                break;

                case 'P':
                    bs_cgfs.resize(bs_cgfs.size() + 3);
                break;

                case 'D':
                    bs_cgfs.resize(bs_cgfs.size() + 6);
                break;
            }

            // collect gtos for cgf
            unsigned int cgfs_size = bs_cgfs.size();
            for(unsigned int j=0; j<ngtos; j++) {
                std::getline(infile, line);
                boost::split(pieces, line, boost::is_any_of(" \t"), boost::token_compress_on);

                double exponent = boost::lexical_cast<double>(pieces[1]);
                double coefficient = boost::lexical_cast<double>(pieces[2]);

                switch(type) {
                    case 'S':
                        bs_cgfs[cgfs_size - 1].add_gto(CGF::GTO_S, exponent, coefficient, origin);
                    break;

                    case 'P':
                        bs_cgfs[cgfs_size-3].add_gto(CGF::GTO_PX, exponent, coefficient, origin);
                        bs_cgfs[cgfs_size-2].add_gto(CGF::GTO_PY, exponent, coefficient, origin);
                        bs_cgfs[cgfs_size-1].add_gto(CGF::GTO_PZ, exponent, coefficient, origin);
                    break;

                    case 'D':
                        bs_cgfs[cgfs_size-6].add_gto(CGF::GTO_DX2, exponent, coefficient, origin);
                        bs_cgfs[cgfs_size-5].add_gto(CGF::GTO_DXY, exponent, coefficient, origin);
                        bs_cgfs[cgfs_size-4].add_gto(CGF::GTO_DXZ, exponent, coefficient, origin);
                        bs_cgfs[cgfs_size-3].add_gto(CGF::GTO_DY2, exponent, coefficient, origin);
                        bs_cgfs[cgfs_size-2].add_gto(CGF::GTO_DYZ, exponent, coefficient, origin);
                        bs_cgfs[cgfs_size-1].add_gto(CGF::GTO_DZ2, exponent, coefficient, origin);
                    break;
                }
            }
        }

        // add basis functions to the atoms
        for(unsigned int i=0; i<this->atoms.size(); i++) {
            if(this->atoms[i]->get_charge() == atnr) {
                for(unsigned int j=0; j<bs_cgfs.size(); j++) {
                    bs_cgfs[j].set_position(this->atoms[i]->get_position());
                    this->add_cgf(i, bs_cgfs[j]);
                }
            }
        }

        if(atnr == highest_atom) { // highest atom reached, break reading file
            break;
        }
    }

    infile.close();
}

unsigned int Molecule::get_nr_elec() const {
    unsigned int nr_elec = 0;

    for(unsigned int i=0; i<this->atoms.size(); i++) {
        nr_elec += this->atoms[i]->get_charge();
    }

    return nr_elec;
}

unsigned int Molecule::get_atom_number_from_string(std::string el) {
    if(el.compare("H") == 0) {
        return 1;
    }
    if(el.compare("He") == 0) {
        return 2;
    }
    if(el.compare("Li") == 0) {
        return 3;
    }
    if(el.compare("Be") == 0) {
        return 4;
    }
    if(el.compare("B") == 0) {
        return 5;
    }
    if(el.compare("C") == 0) {
        return 6;
    }
    if(el.compare("N") == 0) {
        return 7;
    }
    if(el.compare("O") == 0) {
        return 8;
    }
    if(el.compare("F") == 0) {
        return 9;
    }
    if(el.compare("Ne") == 0) {
        return 10;
    }
    if(el.compare("Ar") == 0) {
        return 18;
    }

    return 0;
}

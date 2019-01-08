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

/**
 * @brief      atom constructor
 */
Atom::Atom() :
    atnr(0),
    position(vec3(0,0,0))
{
}

/**
 * @brief      atom constructor
 *
 * @param[in]  atnr       atom number (i.e. number of protons)
 * @param[in]  _position  atom position
 */
Atom::Atom(unsigned int _atnr, const vec3& _position):
    atnr(_atnr),
    position(_position) {
}

/**
 * @brief      molecule constructor
 *
 * @param[in]  filename   input file
 * @param[in]  _settings  pointer to settings object
 */
Molecule::Molecule(const std::string& filename, const std::shared_ptr<Settings>& _settings) :
settings(_settings)
{
    this->read_molecule_from_file(filename);
}

/**
 * @brief      Reads a molecule from file.
 *
 * @param[in]  filename  The filename
 */
void Molecule::read_molecule_from_file(const std::string& filename) {
    std::cout << "             Reading input file              " << std::endl;
    std::cout << "=============================================" << std::endl;

    static const double angstrom_to_bohr = 1.889725989;
    bool unit_angstrom;
    try {
        unit_angstrom = this->settings->get_value("units").compare("angstrom") == 0 ? true : false;
    } catch(const std::exception& e) {
        unit_angstrom = false;
    }

    std::cout << "Reading file:\t\t" << filename << std::endl;
    std::cout << std::endl;

    std::cout << "System name: " << this->settings->get_value("name") << std::endl;
    std::string basis_set = "basis/" + this->settings->get_value("basis") + ".dat";
    std::cout << "Basis set: " << basis_set << std::endl;
    std::cout << std::endl;

    // ignore all settings information
    std::ifstream infile(filename);
    std::string line;
    while(std::getline(infile, line)) {
        if(line.substr(0,7).compare("system:") == 0) {
            break;
        }
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
    this->print_geometry();
}

/**
 * @brief      perturb atoms of molecule
 *
 * @param[in]  p     perturbation vector
 */
void Molecule::perturb_atoms(const VectorXd& p) {
    for(unsigned int i=0; i<this->get_nr_atoms(); i++) {
        vec3 newpos = this->get_atom(i)->get_position();
        for(unsigned int j=0; j<3; j++) {
            newpos[j] += p[i*3 + j];
        }
        this->get_atom(i)->set_position(newpos);
    }

    this->update_cgfs();
}

/**
 * @brief      print the current geometry
 */
void Molecule::print_geometry() const {
    std::cout << "=============================================" << std::endl;

    for(unsigned int i=0; i<this->get_nr_atoms(); i++) {
        std::cout << boost::format("%i  %12.6f  %12.6f  %12.6f")
                     % this->atoms[i]->get_charge()
                     % this->atoms[i]->get_position()[0]
                     % this->atoms[i]->get_position()[1]
                     % this->atoms[i]->get_position()[2] << std::endl;
    }

    std::cout << "=============================================" << std::endl;
    std::cout << "Total number of GTOs: " << this->get_nr_gtos() << std::endl;
    std::cout << std::endl;
}

/**
 * @brief      Sets the basis set.
 *
 * @param[in]  basis_set  The basis set
 */
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
    if (!boost::filesystem::exists("../" + basis_set)) {
        std::cerr << "Please make sure you are running dftcxx from the build directory..." << std::endl;
        throw std::runtime_error("Cannot open ../" + basis_set + "!");
    }

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

/**
 * @brief      Gets the number of elec.
 *
 * @return     The number of elec.
 */
unsigned int Molecule::get_nr_elec() const {
    unsigned int nr_elec = 0;

    for(unsigned int i=0; i<this->atoms.size(); i++) {
        nr_elec += this->atoms[i]->get_charge();
    }

    return nr_elec;
}

/**
 * @brief      the total number of gtos in this molecule
 *
 * @return     total number of gtos
 */
unsigned int Molecule::get_nr_gtos() const {
    unsigned int sz = 0;

    for(unsigned int i=0; i<this->cgfs.size(); i++) {
        sz += this->cgfs[i].size();
    }

    return sz;
}

/**
 * @brief      Gets the atom number from string.
 *
 * @param[in]  el    element name
 *
 * @return     atom number
 */
unsigned int Molecule::get_atom_number_from_string(std::string el) {
    static const std::vector<std::string> names = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na",
                                                   "Mg", "Al", "Si", "P", "S", "Cl", "Ar"};

    for(unsigned int i=0; i<names.size(); i++) {
        if(names[i].compare(el) == 0) {
            return (i+1);
        }
    }

    throw std::runtime_error("Unknown element: " + el);
}

/**
 * @brief      update the cgfs after a geometry change
 */
void Molecule::update_cgfs() {
    this->cgfs.clear();
    this->set_basis_set("basis/" + this->settings->get_value("basis") + ".dat");
}

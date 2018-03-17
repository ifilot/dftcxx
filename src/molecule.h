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

#ifndef _MOLECULE_H
#define _MOLECULE_H

#include <Eigen/Dense>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <string>
#include <bitset>
#include <fstream>

#include "cgf.h"

typedef Eigen::Vector3d vec3;

class Atom {
public:
    Atom();

    Atom(unsigned int atnr, const vec3& _position);

    inline const vec3& get_position() const {
        return this->position;
    }

    inline const unsigned int get_charge() const {
        return this->atnr;
    }

private:
    unsigned int atnr;      // atomic number
    vec3 position;          // atom position
};

class Molecule {
public:
    Molecule(const std::string& filename);

    void read_molecule_from_file(const std::string& filename);

    inline unsigned int get_nr_atoms() const {
        return this->atoms.size();
    }

    inline unsigned int get_nr_bfs() const {
        return this->cgfs.size();
    }

    inline const CGF& get_cgf(unsigned int i) const {
        return this->cgfs[i];
    }

    inline const std::vector<CGF>* get_cgfs() const {
        return &this->cgfs;
    }

    inline const std::shared_ptr<Atom>& get_atom(unsigned int i) const {
        return this->atoms[i];
    }

    inline const vec3& get_atomic_position(unsigned int i) const {
        return this->atoms[i]->get_position();
    }

    inline const unsigned int get_atomic_charge(unsigned int i) const {
        return this->atoms[i]->get_charge();
    }

    inline void add_atom(const Atom& atom) {
        this->atoms.emplace_back(std::make_shared<Atom>(atom));
    }

    inline void add_cgf(unsigned int atid, const CGF cgf) {
        this->cgfs.push_back(cgf);
    }

    void set_basis_set(const std::string& basis_set);

    unsigned int get_nr_elec() const;

private:
    std::vector<std::shared_ptr<Atom>> atoms;
    std::vector<CGF> cgfs;  // vector of cgfs of this molecule

    unsigned int get_atom_number_from_string(std::string el);
};

#endif //_MOLECULE_H

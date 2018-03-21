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
#include "settings.h"

typedef Eigen::Vector3d vec3;

/**
 * @brief      Class for atom.
 */
class Atom {
public:

    /**
     * @brief      atom constructor
     */
    Atom();

    /**
     * @brief      atom constructor
     *
     * @param[in]  atnr       atom number (i.e. number of protons)
     * @param[in]  _position  atom position
     */
    Atom(unsigned int atnr, const vec3& _position);

    /**
     * @brief      Gets the position.
     *
     * @return     The position.
     */
    inline const vec3& get_position() const {
        return this->position;
    }

    /**
     * @brief      Gets the charge.
     *
     * @return     The charge.
     */
    inline const unsigned int get_charge() const {
        return this->atnr;
    }

private:
    unsigned int atnr;      // atomic number
    vec3 position;          // atom position
};

/**
 * @brief      Class for molecule.
 */
class Molecule {
private:
    std::vector<std::shared_ptr<Atom>> atoms;   // atoms in the molecule
    std::vector<CGF> cgfs;                      // vector of cgfs of this molecule
    std::shared_ptr<Settings> settings;         // pointer to settings file

public:

    /**
     * @brief      molecule constructor
     *
     * @param[in]  filename   input file
     * @param[in]  _settings  pointer to settings object
     */
    Molecule(const std::string& filename, const std::shared_ptr<Settings>& _settings);

    /**
     * @brief      Reads a molecule from file.
     *
     * @param[in]  filename  The filename
     */
    void read_molecule_from_file(const std::string& filename);

    /**
     * @brief      Gets the number of atoms.
     *
     * @return     The number of atoms.
     */
    inline unsigned int get_nr_atoms() const {
        return this->atoms.size();
    }

    /**
     * @brief      Gets the number of bfs.
     *
     * @return     The number of bfs.
     */
    inline unsigned int get_nr_bfs() const {
        return this->cgfs.size();
    }

    /**
     * @brief      get contracted Gaussian function (CGF) i
     *
     * @param[in]  i     cgf identifier
     *
     * @return     The cgf.
     */
    inline const CGF& get_cgf(unsigned int i) const {
        return this->cgfs[i];
    }

    /**
     * @brief      get the vector of all contracted Gaussian functions (cgf)
     *
     * @return     The cgfs.
     */
    inline const std::vector<CGF>* get_cgfs() const {
        return &this->cgfs;
    }

    /**
     * @brief      get pointer to atom
     *
     * @param[in]  i     atom id
     *
     * @return     pointer to atom class
     */
    inline const std::shared_ptr<Atom>& get_atom(unsigned int i) const {
        return this->atoms[i];
    }

    /**
     * @brief      Gets the atomic position.
     *
     * @param[in]  i     atom id
     *
     * @return     The atomic position.
     */
    inline const vec3& get_atomic_position(unsigned int i) const {
        return this->atoms[i]->get_position();
    }

    /**
     * @brief      Gets the atomic charge.
     *
     * @param[in]  i     atom id
     *
     * @return     The atomic charge.
     */
    inline const unsigned int get_atomic_charge(unsigned int i) const {
        return this->atoms[i]->get_charge();
    }

    /**
     * @brief      Adds an atom.
     *
     * @param[in]  atom  The atom
     */
    inline void add_atom(const Atom& atom) {
        this->atoms.emplace_back(std::make_shared<Atom>(atom));
    }

    /**
     * @brief      Adds a contracted Gaussian functional
     *
     * @param[in]  atid  atom id
     * @param[in]  cgf   The cgf
     */
    inline void add_cgf(unsigned int atid, const CGF cgf) {
        this->cgfs.push_back(cgf);
    }

    /**
     * @brief      Sets the basis set.
     *
     * @param[in]  basis_set  The basis set
     */
    void set_basis_set(const std::string& basis_set);

    /**
     * @brief      Gets the number of elec.
     *
     * @return     The number of elec.
     */
    unsigned int get_nr_elec() const;

    /**
     * @brief      the total number of gtos in this molecule
     *
     * @return     total number of gtos
     */
    unsigned int get_nr_gtos() const;

private:

    /**
     * @brief      Gets the atom number from string.
     *
     * @param[in]  el    element name
     *
     * @return     atom number
     */
    unsigned int get_atom_number_from_string(std::string el);
};

#endif //_MOLECULE_H

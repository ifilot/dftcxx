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

#ifndef _SETTINGS_H
#define _SETTINGS_H

#include <string>
#include <unordered_map>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem.hpp>

class Settings {
private:
    std::string name;                   // name of the system
    std::string basis_set;              // basis set (i.e. sto-3g, sto-6g)
    unsigned int grid_level;            // fineness of the grid
    unsigned int radial_points;         // number of radial points
    unsigned int lebedev_order;         // lebedev order
    unsigned int lmax;                  // lmax
    unsigned int hartree_evaluation;    // Hartree evaluation method

    std::unordered_map<std::string, std::string> key_values;    // map to store settings

public:

    enum {  // defines the type of evaluation for the Hartree potential
        BECKE_GRID,
        TWO_ELECTRON_INTEGRALS,

        NUM_HARTEE_EVALUATION
    };

    enum {  // defines fineness of the numerical integration
        GRID_COARSE,
        GRID_MEDIUM,
        GRID_FINE,
        GRID_ULTRAFINE,

        NR_GRID_RESOLUTIONS
    };

    enum {  // lebedev orders
        LEBEDEV_6,
        LEBEDEV_14,
        LEBEDEV_26,
        LEBEDEV_38,
        LEBEDEV_50,
        LEBEDEV_74,
        LEBEDEV_86,
        LEBEDEV_110,
        LEBEDEV_146,
        LEBEDEV_170,
        LEBEDEV_194,

        NUM_LEBEDEV_POINTS
    };

    /**
     * @brief      constructor
     *
     * @param[in]  filename  input file
     */
    Settings(const std::string& filename);

    /**
     * @brief      obtain a value by a key from settings
     *
     * @param[in]  key   key
     *
     * @return     value
     */
    const std::string& get_value(const std::string& key) const;

    /**
     * @brief      Gets the hartree evaluation method
     *
     * @return     The hartree evaluation method
     */
    inline unsigned int get_hartree_evaluation_method() const {
        return this->hartree_evaluation;
    }

    /**
     * @brief      get the number of radial points
     *
     * @return     number of radial points
     */
    inline unsigned int get_radial_points() const {
        return this->radial_points;
    }

    /**
     * @brief      get the Lebedev order
     *
     * @return     The lebedev order.
     */
    inline unsigned int get_lebedev_order() const {
        return this->lebedev_order;
    }

    /**
     * @brief      return maximum angular value
     *
     * @return     The lmax.
     */
    inline unsigned int get_lmax() const {
        return this->lmax;
    }

private:

    /**
     * @brief      Sets the default settings.
     */
    void set_default_settings();

    /**
     * @brief      Reads a settings file.
     *
     * @param[in]  filename  input file
     */
    void read_settings_file(const std::string& filename);

    /**
     * @brief      Sets the grid fineness.
     *
     * @param[in]  fineness  fineness constant
     */
    void set_grid_fineness(unsigned int fineness);
};

#endif // _SETTINGS_H

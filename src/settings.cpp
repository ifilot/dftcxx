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

#include "settings.h"

/**
 * @brief      constructor
 *
 * @param[in]  filename  input file
 */
Settings::Settings(const std::string& filename) {
    this->read_settings_file(filename);
    this->set_default_settings();
}

/**
 * @brief      Reads a settings file.
 *
 * @param[in]  filename  input file
 */
void Settings::read_settings_file(const std::string& filename) {
    // check if filename exists
    if (!boost::filesystem::exists(filename)) {
        throw std::runtime_error("Cannot open " + filename + "!");
    }

    // static regex patterns
    const boost::regex regex_system("^system:\\s*$");
    boost::smatch what;

    // read from file
    std::ifstream infile(filename);
    std::string line;
    std::vector<std::string> pieces;

    while(std::getline(infile, line)) {
        boost::split(pieces, line, boost::is_any_of("="), boost::token_compress_on);
        if(pieces.size() == 2) {
            std::string keyword = pieces[0];
            std::string value = pieces[1];

            boost::trim(keyword);
            boost::trim(value);

            this->key_values.emplace(keyword, value);
        }

        if(boost::regex_match(line, what, regex_system)) {
            break;
        }
    }
}

/**
 * @brief      obtain a value by a key from settings
 *
 * @param[in]  key   key
 *
 * @return     value
 */
const std::string& Settings::get_value(const std::string& key) const {
    auto got = this->key_values.find(key);
    if(got != this->key_values.end()) {
        return got->second;
    } else {
        throw std::logic_error("Could not find " + key);
    }
}

/**
 * @brief      Sets the default settings.
 */
void Settings::set_default_settings() {
    // *****************************
    // set hartree evaluation method
    // *****************************
    std::string he;
    try {
        he = this->get_value("hartree_evaluation");
    } catch(const std::exception& e) {
        he = "becke_grid";
    }

    if(he.compare("two_electron_integrals") == 0) {
        this->hartree_evaluation = Settings::TWO_ELECTRON_INTEGRALS;
    } else if(he.compare("becke_grid")) {
        this->hartree_evaluation = Settings::BECKE_GRID;
    } else {
        Settings::BECKE_GRID;
    }

    // *****************************
    // set fineness of the grid
    // *****************************
    std::string grid_fineness;
    try {
        grid_fineness = this->get_value("grid");
    } catch(const std::exception& e) {
        grid_fineness = "medium";
    }

    if(grid_fineness.compare("coarse") == 0) {
        this->set_grid_fineness(GRID_COARSE);
    } else if(grid_fineness.compare("medium") == 0) {
        this->set_grid_fineness(GRID_MEDIUM);
    } else if(grid_fineness.compare("fine") == 0) {
        this->set_grid_fineness(GRID_FINE);
    } else if(grid_fineness.compare("ultrafine") == 0) {
        this->set_grid_fineness(GRID_ULTRAFINE);
    }

    // ************************************
    // overrule grid parameters if present
    // ************************************

    try {
        this->radial_points = boost::lexical_cast<unsigned int>(this->get_value("radial_points"));
    } catch(const std::exception& e) {
        // do nothing
    }

    try {
        this->lebedev_order = boost::lexical_cast<unsigned int>(this->get_value("lebedev_order"));
    } catch(const std::exception& e) {
        // do nothing
    }

    try {
        this->lmax = boost::lexical_cast<unsigned int>(this->get_value("lmax"));
    } catch(const std::exception& e) {
        // do nothing
    }
}

/**
 * @brief      Sets the grid fineness.
 *
 * @param[in]  fineness  fineness constant
 */
void Settings::set_grid_fineness(unsigned int fineness) {
    // set the resolution of the grid
    switch(fineness) {
        case GRID_COARSE:
            this->radial_points = 10;
            this->lebedev_order = LEBEDEV_50;
            this->lmax = 5;
        break;
        case GRID_MEDIUM:
            this->radial_points = 15;
            this->lebedev_order = LEBEDEV_110;
            this->lmax = 8;
        break;
        case GRID_FINE:
            this->radial_points = 20;
            this->lebedev_order = LEBEDEV_146;
            this->lmax = 10;
        break;
        case GRID_ULTRAFINE:
            this->radial_points = 30;
            this->lebedev_order = LEBEDEV_194;
            this->lmax = 11;
        break;
        default: // medium settings
            this->radial_points = 15;
            this->lebedev_order = LEBEDEV_110;
            this->lmax = 8;
        break;
    }
}

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

#ifndef _UTIL_PNG_H
#define _UTIL_PNG_H

#include <fstream>
#include <iostream>
#include <vector>
#include <png.h>

namespace PNG {

    void write_image_buffer_to_png(const std::string& filename, const std::vector<uint8_t>& buffer, unsigned int width, unsigned int height, unsigned int col);
    void write_file_callback(png_structp png_ptr, png_bytep data, png_size_t count);
}

#endif //_UTIL_PNG_H

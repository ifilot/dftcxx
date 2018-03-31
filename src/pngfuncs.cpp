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

#include "pngfuncs.h"

void PNG::write_image_buffer_to_png(const std::string& filename, const std::vector<uint8_t>& buffer, unsigned int width, unsigned int height, unsigned int col) {
    png_structp png_ptr;
    png_infop info_ptr;

    /* create file */
    std::ofstream ofile(filename.c_str(), std::ios::binary);

    if (!ofile.is_open() ) {
        std::cerr << "[write_png_file] File " << filename.c_str() << " could not be opened for reading" << std::endl;
        exit(-1);
    }

    /* initialize stuff */
    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

    if (!png_ptr) {
        std::cerr << "[write_png_file] png_create_write_struct failed" << std::endl;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        std::cerr << "[write_png_file] png_create_info_struct failed" << std::endl;
    }

    if (setjmp(png_jmpbuf(png_ptr))) {
        std::cerr << "[write_png_file] Error during init_io";
    }

    png_set_write_fn(png_ptr, (void *)&ofile, write_file_callback, NULL);

    /* write header */
    if (setjmp(png_jmpbuf(png_ptr))) {
        std::cerr << "[write_png_file] Error during writing header" << std::endl;
    }

    png_set_IHDR(png_ptr, info_ptr, width, height,
                 8, col, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    png_write_info(png_ptr, info_ptr);


    /* write bytes */
    if (setjmp(png_jmpbuf(png_ptr))) {
        std::cerr << "[write_png_file] Error during writing bytes" << std::endl;
    }

    png_bytep *row_pointers;
    if(col == PNG_COLOR_TYPE_GRAY) {
        row_pointers = new png_bytep[height];
        for(unsigned int i=0; i<height; i++) {
            row_pointers[i] = new unsigned char[width];
            for(unsigned int j=0; j<width; j++) {
                row_pointers[i][j] = buffer[i * width + j];
            }
        }
    } else {
        static const unsigned int coldepth = 4;
        row_pointers = new png_bytep[height];
        for(unsigned int i=0; i<height; i++) {
            row_pointers[i] = new unsigned char[width * coldepth];
            for(unsigned int j=0; j<width; j++) {
                for(unsigned int p=0; p<coldepth; p++) {
                    // note that height needs to inverted for correct image capture
                    row_pointers[i][j * coldepth + p] = buffer[((height - i - 1) * width + j) * coldepth + p];
                }
            }
        }
    }
    png_write_image(png_ptr, row_pointers);

    /* end write */
    if (setjmp(png_jmpbuf(png_ptr))) {
        std::cerr << "[write_png_file] Error during end of write" << std::endl;
    }

    png_write_end(png_ptr, NULL);

    for(unsigned int i=0; i<height; i++) {
        delete[] row_pointers[i];
    }
    delete[] row_pointers;

    png_destroy_write_struct(&png_ptr, &info_ptr);

    ofile.close();
}

void PNG::write_file_callback(png_structp png_ptr, png_bytep data, png_size_t count) {
    std::ofstream *outfile = (std::ofstream*)png_get_io_ptr(png_ptr);
    outfile->write((char*)data, count);
}

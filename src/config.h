#/**************************************************************************
 #  config.h.in  --  This file is part of DFTCXX.                        *
 #                                                                       *
 #  Copyright (C) 2016, Ivo Filot                                        *
 #                                                                       *
 #  DFTCXX is free software:                                             *
 #  you can redistribute it and/or modify it under the terms of the      *
 #  GNU General Public License as published by the Free Software         *
 #  Foundation, either version 3 of the License, or (at your option)     *
 #  any later version.                                                   *
 #                                                                       *
 #  DFTCXX is distributed in the hope that it will be useful,            *
 #  but WITHOUT ANY WARRANTY; without even the implied warranty          *
 #  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 #  See the GNU General Public License for more details.                 *
 #                                                                       *
 #  You should have received a copy of the GNU General Public License    *
 #  along with this program.  If not, see http://www.gnu.org/licenses/.  *
 #                                                                       *
#**************************************************************************/

#ifndef _CONFIG_H
#define _CONFIG_H

#define VERSION_MAJOR 1
#define VERSION_MINOR 1
#define VERSION_MICRO 2
#define VERSION "1.1.2"

static const std::string PROGRAM_VERSION(VERSION);
static const unsigned int PROGRAM_VERSION_MAJOR = VERSION_MAJOR;
static const unsigned int PROGRAM_VERSION_MINOR = VERSION_MINOR;
static const unsigned int PROGRAM_VERSION_MICRO = VERSION_MICRO;

#endif // _CONFIG_H

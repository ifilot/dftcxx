 #*************************************************************************
 #   CMakeLists.txt  --  This file is part of DFTCXX.                     *
 #                                                                        *
 #   Author: Ivo Filot <i.a.w.filot@tue.nl>                               *
 #                                                                        *
 #   DFTCXX is free software: you can redistribute it and/or modify       *
 #   it under the terms of the GNU General Public License as published    *
 #   by the Free Software Foundation, either version 3 of the License,    *
 #   or (at your option) any later version.                               *
 #                                                                        *
 #   DFTCXX is distributed in the hope that it will be useful,            *
 #   but WITHOUT ANY WARRANTY; without even the implied warranty          *
 #   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 #   See the GNU General Public License for more details.                 *
 #                                                                        *
 #   You should have received a copy of the GNU General Public License    *
 #   along with this program.  If not, see http://www.gnu.org/licenses/.  *
 #                                                                        *
 #*************************************************************************/

add_library(unittest STATIC unittest.cpp)

# set executables
SET(EXECUTABLES TestMolecules)

#######################################################
# Add executables
#######################################################
add_executable(TestMolecules test_molecules.cpp)

#######################################################
# Link mkmsources and other dependencies
#######################################################
if(CMAKE_BUILD_TYPE MATCHES Debug)
    SET(CPPUNIT_LIB ${CPPUNIT_DEBUG_LIBRARY})
else()
    SET(CPPUNIT_LIB ${CPPUNIT_LIBRARY})
endif()

# configure common settings for test executables
foreach(testexec ${EXECUTABLES})
    target_link_libraries(${testexec}
                          unittest
                          dftcxxsources
                          ${Boost_LIBRARIES}
                          ${CPPUNIT_LIB})
    add_test(NAME ${testexec} COMMAND ${testexec})
    set_tests_properties(${testexec} PROPERTIES FIXTURES_REQUIRED Dataset)
endforeach()

//---------------------------------------------------------------------------------------
//
//	CAVVIS - A Field-of-View Geometric Algorithm for Protein Cavity Detection
//
//  Copyright (C) 2019 Instituto de Telecomunicações (www.it.pt)
//  Copyright (C) 2019 Universidade da Beira Interior (www.ubi.pt)
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//  Contacts:
//      Tiago Simões <tiago.simoes@it.ubi.pt>
//      Abel Gomes <agomes@di.ubi.pt>
//---------------------------------------------------------------------------------------
#ifndef __CavVis__Utils__
#define __CavVis__Utils__
#include <stdio.h>
#include <stdio.h>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include <complex.h>
#include <unistd.h>
#include "DataStructures.h"
using namespace std;
namespace CavVis{
    namespace Utils {
        /* String to bool */
        bool string2bool (const std::string & v);
        /* Split string */
        vector<string> split(string line);
        /* Remove spaces of string */
        string removeSpaces(string input);
        /* Remove char of string */
        string removeChar(string input, char x);
        /* Bool var to string */
        string to_string(bool b);
        /* If necessary adds a slash to the end of a directory name */
        bool directoryNameHasEndSlash(string str);
        /* Check if a directory string as a slash at the end */
        string directoryEndSlashFix(string str);
        /* Get max point as a float. Set to true which coordinate you need */
        float getMax(vector<Point> * points, bool x, bool y, bool z);
        /* Get min point as a float. Set to true which coordinate you need */
        float getMin(vector<Point> * points, bool x, bool y, bool z);
        /* Get max point as a Point. Set to true which coordinate you need */
        Point getMaxPoint(vector<Point> * points, bool x, bool y, bool z);
        /* Get min point as a Point. Set to true which coordinate you need */
        Point getMinPoint(vector<Point> * points, bool x, bool y, bool z);
        /* Check if float a and b are equal using an epsilon=1000000*/
        bool isEqual(float f1, float f2);
        /* Return current date and time */
        const std::string currentDateTime();
    }
}


#endif /* defined(__CavVis__Utils__) */

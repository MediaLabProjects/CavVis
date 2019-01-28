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
#ifndef __CavVis__Parsers__
#define __CavVis__Parsers__
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <iomanip>
#include "DataStructures.h"
#include "Utils.h"
using namespace std;
namespace CavVis{
    class Parser{
    public:
        Parser(){};
        /* Creates .xyz files for each cavity detected by the CavVis algorithm. */
        void createXYZ(vector<Protein> * proteins, string outputDirectory);
        /* Split parameters of a PDB file */
        vector<string> parametersSplit(string line);
        /* Read contents of a PDB File */
        bool read(string file, string PDB_id, Protein * p);
    private:
        //
    };
}

#endif /* defined(__CavVis__Parsers__) */

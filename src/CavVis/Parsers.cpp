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
#include "Parsers.h"
namespace CavVis{
    
    /* Creates .xyz files for each cavity detected by the CavVis algorithm. */
    void Parser::createXYZ(vector<Protein> * proteins, string outputDirectory){
        outputDirectory = Utils::directoryEndSlashFix(outputDirectory);
        for(int i=0; i < proteins->size(); i++){
            Protein * protein = &proteins->at(i);
            for(int j=0; j < protein->cavities.size(); j++){
                string xyz = outputDirectory + protein->id + "_" + to_string(j) + ".xyz";
                std::ofstream o(xyz.c_str());
                Cavity * c = &protein->cavities[j];
                int natoms      = (int) c->atoms.size();
                o << natoms-1 << endl;
                for(int k=0; k < c->atoms.size(); k++)
                o << "H " << c->atoms[k].coord.x << " " << c->atoms[k].coord.y << " " << c->atoms[k].coord.z << endl;
                o.close();
            }
        }
    }
    
    /* Split parameters of a PDB file */
    vector<string> Parser::parametersSplit(string line){
        bool debug=false;
        vector<string> lines;
        
        // if its an HETATM and not an ATOM line replace to ajust columns positions equal to the ATOM line
        if (line[0] == 'H' && line[1] == 'E' && line[2] == 'T' && line[3] == 'A' && line[4] == 'T' && line[5] == 'M'){
            string r("HETATM");
            line.replace(line.find(r),r.length(),"ATOM  ");
        }
        
        // ATOM
        //Example ATOM Column 1 to Column 4
        string ATOM = std::string() + line[0] + line[1] + line[2] + line[3];
        ATOM = Utils::removeSpaces(ATOM);
        lines.push_back(ATOM); // 0
        if (debug) cout << ATOM << "\n";
        
        // ATOM serial Number
        string serialNumber = string() + line[6] +line[7] + line[8] +line[9] + line[10];
        serialNumber = Utils::removeSpaces(serialNumber);
        lines.push_back(serialNumber); // 1
        if (debug) cout << serialNumber << "\n";
        
        // ATOM Name
        string name = string() + line[12] + line[13] + line[14] + line[15];
        name = Utils::removeSpaces(name);
        lines.push_back(name); // 2
        if (debug) cout << name << "\n";
        
        // Alternate location indicator
        string loc = string() + line[16];
        loc = Utils::removeSpaces(loc);
        lines.push_back(loc); // 3
        if (debug) cout << loc << "\n";
        
        // Residue Name
        string residueName = string() + line[17] + line[18] + line[19];
        residueName = Utils::removeSpaces(residueName);
        lines.push_back(residueName); // 4
        if (debug) cout << residueName << "\n";
        
        // Chain identifier
        string chainIdentifier = string() + line[21];
        chainIdentifier = Utils::removeSpaces(chainIdentifier);
        lines.push_back(chainIdentifier); // 5
        if (debug) cout << chainIdentifier << "\n";
        
        // Residue sequence number
        string residueSequenceNumber = string()+line[22]+line[23]+line[24]+line[25];
        residueSequenceNumber = Utils::removeSpaces(residueSequenceNumber);
        lines.push_back(residueSequenceNumber); // 6
        if (debug) cout << residueSequenceNumber << "\n";
        
        // Code for insertions of residues
        string insertionCode = string() + line[26];
        insertionCode = Utils::removeSpaces(insertionCode);
        lines.push_back(insertionCode); // 7
        if (debug) cout << insertionCode << "\n";
        
        // X coordinate
        string x = string() + line[30] + line[31] + line[32] + line[33] + line[34] + line[35] + line[36] + line[37];
        x = Utils::removeSpaces(x);
        lines.push_back(x); // 8
        if (debug) cout << x << "\n";
        
        // Y coordinate
        string y = string() + line[38] + line[39] + line[40] + line[41] + line[42] + line[43] + line[44] + line[45];
        y = Utils::removeSpaces(y);
        lines.push_back(y); // 9
        if (debug) cout << y << "\n";
        
        // Z coordinate
        string z = string() + line[46] + line[47] + line[48] + line[49] + line[50] + line[51] + line[52]  + line[53];
        z = Utils::removeSpaces(z);
        lines.push_back(z); // 10
        if (debug) cout << z << "\n";
        
        // Occupancy
        string occupancy = string() + line[54] + line[55] + line[56] + line[57] + line[58] + line[59];
        occupancy = Utils::removeSpaces(occupancy);
        lines.push_back(occupancy); // 11
        if (debug) cout << occupancy << "\n";
        
        // Temperature factor
        string temperature = string() + line[60] + line[61] + line[62] + line[63] + line[64] + line[65];
        temperature = Utils::removeSpaces(temperature);
        lines.push_back(temperature); // 12
        if (debug) cout << temperature << "\n";
        
        // Segment identifier
        string segmentIdentifier = string() + line[72] + line[73] + line[74] + line[75];
        segmentIdentifier = Utils::removeSpaces(segmentIdentifier);
        lines.push_back(segmentIdentifier); // 13
        if (debug) cout << segmentIdentifier << "\n";
        
        // Element symbol
        string elementSymbol = string() + line[76] + line[77];
        elementSymbol = Utils::removeSpaces(elementSymbol);
        lines.push_back(elementSymbol); // 14
        if (debug) cout << elementSymbol << "\n";
        
        // Charge
        string charge = string() + line[78] + line[79];
        charge = Utils::removeSpaces(charge);
        lines.push_back(charge); // 15
        if (debug) cout << charge << "\n";
        
        if (debug) cout << line << "\n";
        return lines;
    }
    
    /* Read contents of a PDB File */
    bool Parser::read(string file, string PDB_id, Protein * p){
        string line;
        float minx=50,miny=50,minz=50;
        float maxx=-50,maxy=-50,maxz=-50;
        int k=0; // count number of atoms;
        p->id = PDB_id; // Set id of the protein
        ifstream myfile (file.c_str());
        if (myfile.is_open())
        {
            while (getline (myfile,line)){
                if (line.length() == 0) continue;
                if ((line[0] == 'A' && line[1] == 'T' && line[2] == 'O' && line[3] == 'M')){
                    vector<string> pars = parametersSplit(line);
                    Atom a;
                    // Set x, y, z of the current atom
                    std::stringstream(pars[1] )  >> a.id;
                    
                    std::stringstream(pars[8] )  >> a.coord.x;
                    std::stringstream(pars[9] )  >> a.coord.y;
                    std::stringstream(pars[10])  >> a.coord.z;
                    
                    minx=a.coord.x<minx?a.coord.x:minx;
                    maxx=a.coord.x>maxx?a.coord.x:maxx;
                    miny=a.coord.y<miny?a.coord.y:miny;
                    maxy=a.coord.y>maxy?a.coord.y:maxy;
                    minz=a.coord.z<minz?a.coord.z:minz;
                    maxz=a.coord.z>maxz?a.coord.z:maxz;
                    // Add min, max values
                    p->min_x = minx;
                    p->max_x = maxx;
                    
                    p->min_y = miny;
                    p->max_y = maxy;
                    
                    p->min_z = minz;
                    p->max_z = maxz;
                    
                    // Set type (element symbol) of the current atom
                    a.type = pars[14];
                    if (a.type == "H"){a.color.R = 227; a.color.G = 232; a.color.B = 225; a.radius = hydrogen_H;}
                    if (a.type == "C"){a.color.R = 92; a.color.G = 92; a.color.B = 92;a.radius = carbon_C;}
                    if (a.type == "N"){a.color.R = 92; a.color.G = 92; a.color.B = 92;a.radius = nitrogen_N;}
                    if (a.type == "O"){a.color.R = 168; a.color.G = 0; a.color.B = 0;a.radius = oxygen_O;}
                    if (a.type == "P"){a.color.R = 168; a.color.G = 13; a.color.B = 92;a.radius = phosphorous_P;}
                    if (a.type == "S"){a.color.R = 229; a.color.G = 182; a.color.B = 0;a.radius = sulfur_S;}
                    if (a.type == "F"){a.color.R = 0; a.color.G = 97; a.color.B = 0;a.radius = fluorine_F;}
                    if (a.type == "Ci" || a.type == "CI"){a.color.R = 0; a.color.G = 184; a.color.B = 0;a.radius = chlorine_Ci;}
                    if (a.type == "Cu" || a.type == "CU"){a.color.R = 191; a.color.G = 134; a.color.B = 0;a.radius = copper_Cu;}
                    p->atoms.push_back(a);
                    
                    // Save maximum radius for this protein
                    if (a.radius > p->max_radius) p->max_radius = a.radius;
                    p->numberOfatoms++; // Increment number of available - atoms -
                    k++; // Increment number of available atoms
                }
                
                if (line[0] == 'H' && line[1] == 'E' && line[2] == 'T' && line[3] == 'A' && line[4] == 'T' && line[5] == 'M')
                    p->numberOfHeteroAtoms++; // Increment number of available - hetero atoms -
                
                
            }
            
            myfile.close();
        }else{
            cout << "Unable to open file [" + file + "].\n\n";
            return false;
        }
        return true;
        
    }
}
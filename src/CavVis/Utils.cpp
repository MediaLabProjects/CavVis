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
#include "Utils.h"
namespace CavVis{
    namespace Utils{
        
        /* Check if string is a number */
        bool is_number(const std::string &s) {
            return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
        }
        
        /* String to bool */
        bool string2bool (const std::string & v){
            return !v.empty () &&
            (strcasecmp (v.c_str (), "true") == 0 ||
             atoi (v.c_str ()) != 0);
        }
        
        /* Split string */
        vector<string> splitString(string line){
            string buf;
            stringstream ss(line);
            vector<string> tokens;
            while (ss >> buf)
                tokens.push_back(buf);
            return(tokens);
        }
        
        /* Remove spaces of string */
        string removeSpaces(string input){
            input.erase(std::remove(input.begin(),input.end(),' '),input.end());
            return input;
        }
        
        /* Remove char of string */
        string removeChar(string input, char x){
            input.erase(std::remove(input.begin(),input.end(),x),input.end());
            return input;
        }
        
        /* Bool var to string */
        string to_string(bool b){
            if (b) return "true"; else return "false";
        }
        
        /* Check if a directory string as a slash at the end */
        bool directoryNameHasEndSlash(string str){
            if (str[str.length()-1] == '/') return true; return false;
        }
        
        /* If necessary adds a slash to the end of a directory name */
        string directoryEndSlashFix(string str){
            if ( directoryNameHasEndSlash(str)) return str; else
                return(str + "/");
        }
        
        /* Get max point as a float. Set to true which coordinate you need */
        float getMax(vector<Point> * points, bool x, bool y, bool z){
            if (x){
                auto max = std::max_element(points->begin(), points->end(), [] (Point const& lhs, Point const& rhs) {return lhs.coord.x < rhs.coord.x;});
                return(max->coord.x);
            }
            if (y){
                auto max = std::max_element(points->begin(), points->end(), [] (Point const& lhs, Point const& rhs) {return lhs.coord.y < rhs.coord.y;});
                return(max->coord.y);
            }
            if (z){
                auto max = std::max_element(points->begin(), points->end(), [] (Point const& lhs, Point const& rhs) {return lhs.coord.z < rhs.coord.z;});
                return(max->coord.z);
            }
            return(0);
        }
        
        /* Get min point as a float. Set to true which coordinate you need */
        float getMin(vector<Point> * points, bool x, bool y, bool z){
            if (x){
                auto min = std::max_element(points->begin(), points->end(), [] (Point const& lhs, Point const& rhs) {return lhs.coord.x > rhs.coord.x;});
                return(min->coord.x);
            }
            if (y){
                auto min = std::max_element(points->begin(), points->end(), [] (Point const& lhs, Point const& rhs) {return lhs.coord.y > rhs.coord.y;});
                return(min->coord.y);
            }
            if (z){
                auto min = std::max_element(points->begin(), points->end(), [] (Point const& lhs, Point const& rhs) {return lhs.coord.z > rhs.coord.z;});
                return(min->coord.z);
            }
            return(0);
        }
        
        /* Get max point as a Point. Set to true which coordinate you need */
        Point getMaxPoint(vector<Point> * points, bool x, bool y, bool z){
            float max=0;
            Point Amax;
            Point atm;
            int n= (int) points->size();
            atm = points->at(0);
            if (x) max = atm.coord.x;
            if (y) max = atm.coord.y;
            if (z) max = atm.coord.z;
            Amax = atm;
            for(int i=1; i < n;i++){
                atm = points->at(i);
                if (x) if (max < atm.coord.x){ max = atm.coord.x; Amax = atm;}
                if (y) if (max < atm.coord.y){ max = atm.coord.y; Amax = atm;}
                if (z) if (max < atm.coord.z){ max = atm.coord.z; Amax = atm;}
            }
            return(Amax);
        }
        
        /* Get min point as a Point. Set to true which coordinate you need */
        Point getMinPoint(vector<Point> * points, bool x, bool y, bool z){
            float min=0;
            Point atm;
            Point Amax;
            int n=(int) points->size();
            atm = points->at(0);
            if (x) min = atm.coord.x;
            if (y) min = atm.coord.y;
            if (z) min = atm.coord.z;
            Amax = atm;
            for(int i=1; i < n;i++){
                atm = points->at(i);
                if (x) if (min > atm.coord.x){ min = atm.coord.x; Amax = atm; }
                if (y) if (min > atm.coord.y){ min = atm.coord.y; Amax = atm; }
                if (z) if (min > atm.coord.z){ min = atm.coord.z; Amax = atm; }
            }
            return(Amax);
        }
        
        /* Check if float a and b are equal using an epsilon=1000000*/
        bool isEqual(float f1, float f2){
            return ( (int)(f1 *1000000)) == ((int)(f2 * 1000000) );
        }
        
        /* Return current date and time */
        const std::string currentDateTime() {
            time_t     now = time(0);
            struct tm  tstruct;
            char       buf[80];
            tstruct = *localtime(&now);
            // format time
            strftime(buf, sizeof(buf), "%d-%m-%Y -> %H:%M:%S", &tstruct);
            return buf;
        }
        
    }
}

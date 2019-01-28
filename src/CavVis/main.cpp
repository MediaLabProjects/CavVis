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
#include "main.h"

int main(int argc, char * argv[]) {
    
    // --------------------
    // ! Read command line parameters
    // --------------------
    if (!processArgs(argc, argv)) return 1;
    protein.type = Protein::type::UNSPECIFIED; // default
    
    // ! Get the id of the protein from the filename
    PDB_id = input.substr(input.find_last_of("/\\")+1, input.length());
    PDB_id = PDB_id.substr(0, PDB_id.find_last_of("."));
    
    printf("-----------------------------------------\n");
    cout <<"              CAVVIS v" << version <<"  \n" ;
    printf("-----------------------------------------\n");

    // --------------------
    // ! Read PDB Protein into object
    // --------------------
    Parser pdbParser;
    if (!pdbParser.read(input, PDB_id, &protein)) return 1;
    cout << "Processing protein [" << PDB_id << "] with "<< protein.numberOfatoms+protein.numberOfHeteroAtoms << " atoms. \n";
    cout << "Please wait, this process can take some time.\n\n";
    
    // --------------------
    // ! Start CavVis
    // --------------------
    Algorithm CavVis(&points, &protein, gridSpacing, top, tdistance, dareas, ffilter, false, false);

    // --------------------
    // ! Prepare results and create output files
    // --------------------
    vector<Protein> proteins;
    proteins.push_back(protein); // for now CavSeeker processes only one protein at a time
    pdbParser.createXYZ(&proteins, output);
    
    cout << "- Results successfully generated for protein " << protein.id << ". Thank you for using CavVis.\n\n";

    // --------------------
    // ! OpenGL molecular Viewer (if requested)
    // --------------------
    if (!opengl_viewer) return 0;
    cout << "- Molecular Viewer: Generating Surface for the Viewer, Please Wait...\n";
    
    // Generate Protein Surface
    Grid vgrid; // Main protein grid
    vgrid.gridSpacing = v_gridSpacing;
    vgrid.iso = 1.0;
    vgrid.C   = 0.33; // Smooth parameter
    vgrid.d   = 2.35; // Decay rate parameter of the gaussian surface
    Surface MC;
    MC.MarchingCubes(true, &protein, &vgrid, false);
    mblob = &MC.mblob;
    
    // Add colors to cavity regions
    for(int i=0; i < protein.cavities.size(); i++){
        protein.cavities[i].regionColor = new Color();
        protein.cavities[i].regionColor->random();
    }
    
    string title = "CavVis " + version + " ["+ PDB_id +"] "+ "(" + to_string(protein.numberOfatoms+protein.numberOfHeteroAtoms) + " Atoms)";
    title += ":: To exit this window press [Q] or check terminal for more commands.";
    cout << "  Use in the OpenGL Window [W] to zoom in, [S] to zoom out,\n";
    cout << "  [R] to rotate, [J] to show cavity surfaces, [C] change colors, or [Q] to quit.\n";
    glut(argc, argv, title, 800, 900);
    glutMainLoop();
    
    return 0;
}

/* Generate Cavity Surfaces */
void cavitySurfaces(){
    printf("\n- Generating Cavity surfaces, please wait this may take some time.\n");
    bool debug          = false;
    bool output_mc      = false;    // Output marching cubes information
    float radius        = 1.0;      // Filling/dummy atom radius
    float cgridspacing  = v_gridSpacing;   // Grid spacing used to generate each cavity surface
    
    float minx=FLT_MAX,miny=FLT_MAX,minz=FLT_MAX;
    float maxx=FLT_MIN,maxy=FLT_MIN,maxz=FLT_MIN;
    vector<Protein> cavities;
    
    // Generate data for the marching cubes algorithm using the Protein class
    for(int i=0; i < protein.cavities.size(); i++){
        Protein cavity; // store cavity data as a protein object to the use the marching cubes
        cavity.color = protein.cavities[i].regionColor;
        if (debug){
            printf("Cavity [%d]:\n", i);
            printf("%d\n", (int) protein.cavities[i].atoms.size()-1);
        }
        for(int j=0; j < protein.cavities[i].atoms.size(); j++){
            Atom a = protein.cavities[i].atoms[j];
            a.id = j;
            a.radius = radius;
            minx=a.coord.x<minx?a.coord.x:minx;
            maxx=a.coord.x>maxx?a.coord.x:maxx;
            miny=a.coord.y<miny?a.coord.y:miny;
            maxy=a.coord.y>maxy?a.coord.y:maxy;
            minz=a.coord.z<minz?a.coord.z:minz;
            maxz=a.coord.z>maxz?a.coord.z:maxz;
            cavity.atoms.push_back(a);
            if (debug) printf("H   %f  %f  %f\n", a.coord.x, a.coord.y, a.coord.z);
        }
        // Store min, max values
        cavity.min_x = minx;
        cavity.max_x = maxx;
        cavity.min_y = miny;
        cavity.max_y = maxy;
        cavity.min_z = minz;
        cavity.max_z = maxz;
        cavity.max_radius = 2.0;
        cavity.numberOfatoms = (int) cavity.atoms.size();
        // Store cavity
        cavities.push_back(cavity);
        // reset
        minx=FLT_MAX;miny=FLT_MAX;minz=FLT_MAX;
        maxx=FLT_MIN;maxy=FLT_MIN;maxz=FLT_MIN;
        if (debug){
            printf(" - Number of atoms = %d\n", (int) cavity.numberOfatoms);
            printf(" - min_x=%f, max_x=%f\n", cavity.min_x, cavity.max_x);
            printf(" - min_y=%f, max_y=%f\n", cavity.min_y, cavity.max_y);
            printf(" - min_z=%f, max_z=%f\n", cavity.min_z, cavity.max_z);
            printf("\n");
        }
        
    }
    
    // Generate surface for each cavity previously populated
    for(int i=0; i < cavities.size(); i++){
        Protein * cavity = &cavities[i];
        vector<mTriangle> cblob;
        Grid cgrid; // Main protein grid
        cgrid.gridSpacing = cgridspacing;
        cgrid.iso = 1.0;
        cgrid.C   = 0.33; // Smooth parameter
        cgrid.d   = 2.35; // Decay rate parameter of the gaussian surface
        Surface MC;
        MC.MarchingCubes(true, cavity, &cgrid, output_mc);
        MC.mblob[0].color = cavity->color;
        surfaces.push_back(MC.mblob);
    }
}

/* OpenGL Draw objects in embedded viewer */
void drawObjects(){
    // Change colors in real time
    if (getC_show() == 1){
        for(int i=0; i < protein.cavities.size(); i++){
            protein.cavities[i].regionColor = new Color();
            protein.cavities[i].regionColor->random();
        }
        if (getJ_show() == 1){
            surfaces.clear();
            cavitySurfaces();
        }
        setC_show(0);
    }
    
    // Draw protein surface
    Color surfaceColor(0.827451, 0.827451, 0.827451);
    Viewer::Graphics::drawSurface(mblob, surfaceColor, false, false);
    
    // Draw protein cavity filling spheres
    for(int i=0; i < protein.cavities.size() && getJ_show() == 0; i++){
        Color * regionColor = protein.cavities[i].regionColor;
        // Draw cavity filling spheres
        for(int j=0; j < protein.cavities[i].atoms.size() && getJ_show() == 0; j++){
            Atom * fAtom = &protein.cavities[i].atoms[j];   // filling or dummy atom
            Viewer::Graphics::drawSphere(fAtom, 0.8, regionColor->R, regionColor->G, regionColor->B, false);
        }
    }
    
    // Draw surfaces of each cavity
    if (surfaces.size() == 0 && getJ_show() == 1) cavitySurfaces();
    
    for(int i=0; i < surfaces.size() && getJ_show() == 1; i++){
      vector<mTriangle> * blob = &surfaces[i];
      Color * regionColor = blob->at(0).color;
      Viewer::Graphics::drawSurface(blob, *regionColor, false, false);
    }
}

/* OpenGL core functions */
void displayFunc(void){
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // Draw custom objects
    drawObjects();
    glFlush();
    glutSwapBuffers();
}

/* OpenGL core functions */
void init() {
    glClearColor( 1, 1, 1, 1);
    glEnable(GL_NORMALIZE);
    glEnable(GL_MULTISAMPLE);
    glShadeModel(GL_SMOOTH);
    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, high_shininess);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);
}

/* OpenGL core functions */
void glut(int argc, char * argv[], string title, int width, int height){
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowPosition(50,25);
    glutInitWindowSize(width, height);
    glutCreateWindow(title.c_str());
    glutDisplayFunc(displayFunc);
    glScalef(0.04,0.04,0.04); // default
    zprInit(protein.max_x,protein.min_x, protein.max_y, protein.min_y, protein.max_z, protein.min_z);
    zprSelectionFunc(drawObjects);
    init();
}

/* Handle terminal related stuff */
void printHelp(){
    cout << "\nUsage: CavVis -f <PDB filename> -o <output directory> [options]" << endl;
    cout << "\n";
    cout << "[options]:" << endl;
    cout << "  -s <grid spacing value> or --spacing <grid spacing value>  : Defines the grid spacing\n";
    cout << "                                                               Typically 0.6 is used for best accuracy\n";
    cout << "                                                               and 0.8 for best performance." << endl;
    cout << "  -d <cluster distance> or --distance <cluster distance>     : Defines the cluster distance threshold." << endl;
    cout << "  -a or --areas                                              : Displays cavity areas." << endl;
    cout << "  -n or --nofilter                                           : Do not filter filling spheres (best performance)" << endl;
    cout << "  -v <viewer grid spacing> or --viewer <viewer grid spacing> : Displays an OpenGL viewer." << endl;
    cout << "\n";
    cout << "Use examples:" << endl;
    cout << "  CavVis -f 1mzl.pdb -o output/" << endl;
    cout << "  CavVis -f 1mzl.pdb -o output/ -s 0.6" << endl;
    cout << "  CavVis -f 1mzl.pdb -o output/ -s 0.6 -d 1.4" << endl;
    cout << "  CavVis -f 1mzl.pdb -o output/ -v 0.6 \n\n";
}

/* Handle terminal related stuff */
bool processArgs(int argc, char** argv){
    bool fail=false;
    const char* const short_opts = ":f:o:s:d:t:v:nah";
    const option long_opts[] = {
        {"file", required_argument, nullptr, 'f'},
        {"output", required_argument, nullptr, 'o'},
        
        {"spacing", required_argument, nullptr, 's'},
        {"distance", required_argument, nullptr, 'd'},
        {"viewer", required_argument, nullptr, 'v'},
        {"top", required_argument, nullptr, 't'},
        
        {"nofilter", no_argument, nullptr, 'n'},
        {"areas", no_argument, nullptr, 'a'},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, no_argument, nullptr, 0}
    };
    
    while (true){
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);
        if (-1 == opt) break;
        switch (opt){
            // --------------------------------------------------------------------------------
            // Input PDB File
            case 'f':
                input = std::string(optarg);
                // Filename validation
                if (input.find(".pdb")==std::string::npos){
                    cout << "CavVis: invalid file specified (use a *.pdb file instead)\n";
                    return false;
                }
                break;
                // Filter filling spheres
            case 'n':
                ffilter=false;
                break;
                // Areas
            case 'a':
                dareas=true;
                break;
                // Output directory
            case 'o':
                output = std::string(optarg);
                break;
                // CavVis Grid Spacing value
            case 's':
                try{ gridSpacing = stof(optarg);
                }catch (const std::invalid_argument& ex){ cout << "CavVis: invalid grid spacing specified\n"; return false;}
                cout << "Grid spacing set to " << gridSpacing << endl;
                break;
            case 'd':
                try{ tdistance = stof(optarg);
                }catch (const std::invalid_argument& ex){ cout << "CavVis: invalid distance specified\n"; return false;}
                cout << "Cluster distance merging set to " << tdistance << endl;
                break;
                break;
                // CavVis tops outputted
            case 't':
                try{ top = stoi(optarg);
                }catch (const std::invalid_argument& ex){ cout << "CavVis: invalid top specified\n"; return false;}
                cout << "Top set to " << top << endl;
                break;
                // OpenGL viewer grid spacing
            case 'v':
                try{ v_gridSpacing = stof(optarg);
                }catch (const std::invalid_argument& ex){ cout << "CavVis: invalid viewer grid spacing specified\n"; return false;}
                cout << "Viewer Grid spacing set to " << v_gridSpacing << endl;
                opengl_viewer = true;
                cout << "Using OpenGL viewer" << endl;
                break;
                // --------------------------------------------------------------------------------
            case 'h':
                fail=true;
                break;
            case '?':
                printf("CavVis: unknown option: %c\n", optopt);
                fail=true;
                break;
            default:
                printf("CavVis: -%c needs an argument\n", optopt);
                fail=true;
                break;
        }
    }
    
    if (input.size() == 0){
        cout << "CavVis: input file not specified\n";
        fail=true;
    }
    
    if (output.size() == 0){
        fail=true;
        cout << "CavVis: output directory not specified\n";
    }
    
    if (fail) printHelp();
    return(!fail);
}


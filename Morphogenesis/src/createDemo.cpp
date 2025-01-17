#define __CL_ENABLE_EXCEPTIONS
#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_TARGET_OPENCL_VERSION 300
#define SDK_SUCCESS 0
#define SDK_FAILURE 1

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <inttypes.h>
#include <errno.h>
#include <math.h>
#include <regex>
#include <jsoncpp/json/json.h>
#define CHECK_ERROR(status) if (status != CL_SUCCESS) { printf("Error: %d\n", status); exit(1); }
#include "fluid.h"
#include "fluid_system.h"
#include "host_CL.cpp"
#include <chrono>


int main ( int argc, const char** argv )
{
     /*
     * Arguments: num_particles, spacing, x_dim, y_dim, z_dim, demo_space, simSpace
     * Inputs:       125            1       6       6     6        0          8
     *
     * Additional paths for old version that merged make_demo and make_demo2:
     * /home/goldi/Documents/KDevelop\ Projects/Morphogenesis/Morphogenesis/src/demo
     * /home/goldi/Documents/KDevelop\ Projects/Morphogenesis/Morphogenesis/src/demo/out
     */

    //Initialize Variables-----------------
    uint num_particles, demoType, simSpace;
    char input_folder[256];
    char output_folder[256];
    float spacing, x_dim, y_dim, z_dim;
    //-------------------------------------

    //Find working directory-----------------------------
    fs::path Path = fs::current_path().parent_path();
    const char* directory = Path.c_str();
    std::cout << "Directory: " << directory << std::endl;
    //---------------------------------------------------

    //Initialize JSON---------------------------------------------------
    char jsonPath[512]; sprintf(jsonPath, "%s/%s", directory, argv[8]);
    cout << "Concatenated path: " << jsonPath << std::endl;
    ifstream ifs(jsonPath);
    Json::Reader reader;
    Json::Value obj_;
    Json::Value obj;
    obj["verbosity"] = 1;
    obj["opencl_platform"] = 0;
    obj["opencl_device"] = 0;
    bool b = reader.parse(ifs, obj);
    if (!b) {cout << "Error: " << reader.getFormattedErrorMessages();}
    else    {cout << "NB lists .json file entries alphabetically: \n" << obj ;}
    cout << "\n\n\n" << endl;
    const char* path = obj["demo_path"].asCString();

    //------------------------------------------------------------------

    //Initialize Program------------------------
    FluidSystem fluid(obj);
    fluid.Initialize();     // Clear all buffers
    uint verbosity = obj["verbosity"].asUInt();
    //------------------------------------------

    //Initialize specifications---------------------------------------------------------------
    if ( argc != 9 && argc != 1 ) {
        printf ( "usage: create_demo num_particles spacing x_dim y_dim z_dim \n \
        demoType(0:free falling, 1: remodelling & actuation, 2: diffusion & epigenetics.) \n \
        simSpace(0:regression test, 1:tower, 2:wavepool, 3:small dam break, 4:dual-wavepool, 5: microgravity, \n \
            6:Morphogenesis small demo  7:use SpecificationFile.txt  8:parameter sweep default )\n" );
        return 0;
    } else if (argc == 9) {
        num_particles = atoi(argv[1]);          printf ( "num_particles = %u\n", num_particles );
        spacing =       atof(argv[2]);          printf ( "spacing = %f\n",             spacing );
        x_dim =         atof(argv[3]);          printf ( "x_dim = %f\n",                 x_dim );
        y_dim =         atof(argv[4]);          printf ( "y_dim = %f\n",                 y_dim );
        z_dim =         atof(argv[5]);          printf ( "z_dim = %f\n",                 z_dim );
        demoType =      atof(argv[6]);          printf ( "demoType = %u",             demoType );               printf (" (0:free falling, 1: remodelling & actuation, 2: diffusion & epigenetics.)\n");
        simSpace =      atof(argv[7]);          printf ( "simSpace = %u",              simSpace);               printf (" (0:regression test, 1:tower, 2:wavepool, 3:small dam break, 4:dual-wavepool, 5: microgravity, \n 6:Morphogenesis small demo  7:use SpecificationFile.txt  8:parameter sweep default )\n\n");

        sprintf ( input_folder, "%s/%s", directory, path);
        sprintf ( output_folder, "%s/%s", directory, path);

        // Check if output_folder and input_folder (if provided) exist
        if (!fs::exists(input_folder)) {
            if (!fs::create_directories(input_folder)) {std::cerr << "Error: Failed to create output folder: " << input_folder << std::endl; return 1;}}
        if (!fs::exists(output_folder)) {
            if (!fs::create_directories(output_folder)) {std::cerr << "Error: Failed to create output folder: " << output_folder << std::endl; return 1;}}

        printf ( "input_folder = %s ,\noutput_folder = %s\n", input_folder, output_folder );

        }  else {
            num_particles = 4000;          printf ( "num_particles = %u\n",    num_particles );
            spacing =        1.0;          printf ( "spacing = %f\n",                spacing );
            x_dim =         10.0;          printf ( "x_dim = %f\n",                    x_dim );
            y_dim =         10.0;          printf ( "y_dim = %f\n",                    y_dim );
            z_dim =            3;          printf ( "z_dim = %f\n",                    z_dim );
            demoType =         0;          printf ( "demoType = %u",                demoType );               printf (" (0:free falling, 1: remodelling & actuation, 2: diffusion & epigenetics.)\n");
        simSpace =             8;          printf ( "simSpace = %u",                 simSpace);               printf (" (0:regression test, 1:tower, 2:wavepool, 3:small dam break, 4:dual-wavepool, 5: microgravity, \n 6:Morphogenesis small demo  7:use SpecificationFile.txt  8:parameter sweep default )\n\n");
    }
    //----------------------------------------------------------------------------------------

    //initialize simulation--------------------------------------------------------------------------
    fluid.WriteDemoSimParams(output_folder, GPU_OFF, CPU_YES , num_particles, spacing, x_dim, y_dim, z_dim, demoType, simSpace, verbosity);/*const char * relativePath*/

    if(argc !=1){                                           // i.e not relying on defaults in simspace 8
    fluid.launchParams.num_particles    = num_particles;    // Write default values to fluid.launchParams...
    fluid.launchParams.demoType         = demoType;
    fluid.launchParams.simSpace         = 7;                // i.e. use the Specfile.txt generated.
    fluid.launchParams.x_dim            = x_dim;
    fluid.launchParams.y_dim            = y_dim;
    fluid.launchParams.z_dim            = z_dim;

    fluid.launchParams.num_files        = 400;
    fluid.launchParams.steps_per_InnerPhysicalLoop = 3;
    fluid.launchParams.steps_per_file   = 6;
    fluid.launchParams.freeze_steps     = 1;
    fluid.launchParams.verbosity            = 0;
    fluid.launchParams.file_num         = 0;

    fluid.launchParams.save_ply         = 'n';
    fluid.launchParams.save_csv         = 'n';
    fluid.launchParams.save_vtp         = 'y';
    fluid.launchParams.gene_activity    = 'n';
    fluid.launchParams.remodelling      = 'n';
    }

    string paramsPath(  string(path)  + "/SimParams.txt");    // Set file paths relative to data/ , where SpecfileBatchGenerator will be run.
    string genomePath(  string(path)  +    "/genome.csv");
    string outPath(     string(path)  +           "/out");
    string pointsPath(path);

    for(int i=0;i<paramsPath.length();i++)fluid.launchParams.paramsPath[i] = paramsPath [i];
    for(int i=0;i<pointsPath.length();i++)fluid.launchParams.pointsPath[i] = pointsPath [i];
    for(int i=0;i<genomePath.length();i++)fluid.launchParams.genomePath[i] = genomePath [i];
    for(int i=0;i<outPath.length(); i++)  fluid.launchParams.outPath[i]    =    outPath [i];

    printf("paramsPath: %s\n", fluid.launchParams.paramsPath);
    printf("pointsPath: %s\n", fluid.launchParams.pointsPath);
    printf("genomePath: %s\n", fluid.launchParams.genomePath);
    printf("outPath:    %s\n", fluid.launchParams.outPath);

    fluid.WriteExampleSpecificationFile(input_folder);
    printf("\ncreateDemo finished.\n");
    //fluid.Exit_no_CL ();
    return 0;
}

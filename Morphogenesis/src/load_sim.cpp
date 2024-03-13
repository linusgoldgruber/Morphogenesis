#define __CL_ENABLE_EXCEPTIONS
#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_TARGET_OPENCL_VERSION 300
#define SDK_SUCCESS 0
#define SDK_FAILURE 1

#include <stdio.h>
#include <inttypes.h>
#include <errno.h>
#include <string.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <regex>
#include <filesystem>
#include <jsoncpp/json/json.h>
#include "host_CL.cpp"
#include "fluid_system.h"
#include <CL/cl.h>

#define CHECK_ERROR(err) if (err != CL_SUCCESS) { printf("Error: %d\n", err); exit(1); }


int main ( int argc, const char** argv )
{
    char paramsPath[256];
    char pointsPath[256];
    char genomePath[256];
    char outPath[256];
    uint num_files, steps_per_file, freeze_steps, debug;
    int file_num=0;
    char save_ply, save_csv, save_vtp,  gene_activity, remodelling;
    if ( argc != 13 ) {
        printf ( "usage: load_sim  simulation_data_folder  output_folder  num_files  steps_per_file  freeze_steps save_ply(y/n)  save_csv(y/n)  save_vtp(y/n)  debug(0-5) gene_activity(y/n)  remodelling(y/n) jason_file \n\ndebug: 0=full speed, 1=current special output,  2=host cout, 3=device printf, 4=SaveUintArray(), 5=save .csv after each kernel.\n\n" );
        return 0;
    } else {
        sprintf ( paramsPath, "%s/SimParams.txt", argv[1] );
        printf ( "simulation parameters file = %s\n", paramsPath );

        sprintf ( pointsPath, "%s/particles_pos_vel_color100001.csv", argv[1] );
        printf ( "simulation points file = %s\n", pointsPath );

        sprintf ( genomePath, "%s/genome.csv", argv[1] );
        printf ( "simulation genome file = %s\n", genomePath );

        sprintf ( outPath, "%s/", argv[2] );
        printf ( "output folder = %s\n", outPath );
        
        num_files = atoi(argv[3]);
        printf ( "num_files = %u\n", num_files );
        
        steps_per_file = atoi(argv[4]);
        printf ( "steps_per_file = %u\n", steps_per_file );
        
        freeze_steps = atoi(argv[5]);
        
        save_ply = *argv[6];
        printf ( "save_ply = %c\n", save_ply );
        
        save_csv = *argv[7];
        printf ( "save_csv = %c\n", save_csv );
        
        save_vtp = *argv[8];
        printf ( "save_vtp = %c\n", save_vtp );
        
        debug = atoi(argv[9]);
        printf ("debug = %u\n", debug );
        
        gene_activity = *argv[10];
        printf ("gene_activity = %c\n", gene_activity );
        
        remodelling = *argv[11];
        printf ("remodelling = %c\n", remodelling );

        ifstream ifs(argv[12]);
        printf ("JSON file = %s\n", argv[12]);


    }

std::cout <<"\nchk load_sim_0.1\n"<<std::flush;  

    // Initialize
    ifstream ifs(argv[12]);
    Json::Reader reader;
    Json::Value obj_;
    Json::Value obj;
    obj["verbosity"] = 1;
    obj["opencl_platform"] = 0;
    obj["opencl_device"] = 0;

    bool b = reader.parse(ifs, obj);
    //if (!b) { cout << "Error: " << reader.getFormattedErrorMessages();}   else {cout << "NB lists .json file entries alphabetically: \n" << obj ;}
    cout << "\n\n\n" << endl;
    //obj["kernel_filepath"] = "/home/goldi/Documents/KDevelop Projects/Morphogenesis/Morphogenesis/src/kernelTest.cl";
    FluidSystem fluid(obj);
    fluid.Initialize();
    fluid.InitializeOpenCL();

std::cout <<"\nchk load_sim_0.2\n"<<std::flush;
    
    fluid.ReadSimParams ( paramsPath );
std::cout <<"\nchk load_sim_0.2.1\n"<<std::flush;

    fluid.ReadGenome ( genomePath );    
std::cout <<"\nchk load_sim_0.2.2\n"<<std::flush;

    fluid.ReadPointsCSV2 ( pointsPath, GPU_SINGLE, CPU_YES );    // NB currently GPU allocation is by Allocate particles, called by ReadPointsCSV.
std::cout <<"\nchk load_sim_0.3\n"<<std::flush;

    fluid.Init_CLRand();
std::cout <<"\nchk load_sim_1.0\n"<<std::flush;

    auto old_begin = std::chrono::steady_clock::now();
    
    fluid.TransferFromCL ();
    fluid.SavePointsCSV2 ( outPath, file_num );
    if(save_vtp=='y') fluid.SavePointsVTP2( outPath, file_num);
    file_num++;
    
    fluid.TransferPosVelVeval ();

std::cout <<"\nchk load_sim_2.0\n"<<std::flush;
    fluid.setFreeze(true);
    for (int k=0; k<freeze_steps; k++){
        if(debug>1) std::cout<<"\n\nFreeze()"<<k<<"\n"<<std::flush;
         /*
        fluid.Freeze (outPath, file_num);                   // save csv after each kernel - to investigate bugs
        file_num+=10;
         */
        //fluid.Freeze (outPath, file_num, (debug=='y'), (gene_activity=='y'), (remodelling=='y')  );       // creates the bonds // fluid.Freeze(outPath, file_num) saves file after each kernel,, fluid.Freeze() does not.         // fluid.Freeze() creates fixed bond pattern, triangulated cubic here. 
        
        fluid.Run2Simulation();
        fluid.TransferPosVelVeval (); // Freeze movement until heal() has formed bonds, over 1st n timesteps.
        if(save_csv=='y'||save_vtp=='y') fluid.TransferFromCL ();
        if(save_csv=='y') fluid.SavePointsCSV2 ( outPath, file_num+90);
        if(save_vtp=='y') fluid.SavePointsVTP2 ( outPath, file_num+90);
        file_num+=100;
    }
    fluid.setFreeze(false);

    if(debug>1) printf("\n\nFreeze finished, starting normal Run ##############################################\n\n");
    
    for ( ; file_num<num_files; file_num+=100 ) {
        
        for ( int j=0; j<steps_per_file; j++ ) {//, bool gene_activity, bool remodelling 
            
            fluid.Run2Simulation();  // run the simulation  // Run(outPath, file_num) saves file after each kernel,, Run() does not.
        }// 0:start, 1:InsertParticles, 2:PrefixSumCellsCL, 3:CountingSortFull, 4:ComputePressure, 5:ComputeForce, 6:Advance, 7:AdvanceTime

        //fluid.SavePoints (i);                         // alternate file formats to write
        // TODO flip mutex
        auto begin = std::chrono::steady_clock::now();
        if(save_csv=='y'||save_vtp=='y') fluid.TransferFromCL ();
        if(save_csv=='y') fluid.SavePointsCSV2 ( outPath, file_num+90);
        if(save_vtp=='y') fluid.SavePointsVTP2 ( outPath, file_num+90);
        cout << "\n File# " << file_num << ". " << std::flush;
        
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> time = end - begin;
        std::chrono::duration<double> begin_dbl = begin - old_begin;
        if(debug>0) std::cout << "\nLoop duration : "
                    << begin_dbl.count() <<" seconds. Time taken to write files for "
                    << fluid.NumPoints() <<" particles : " 
                    << time.count() << " seconds\n" << std::endl;
        old_begin = begin;
        
        //fluid.WriteParticlesToHDF5File(i);
        //if(debug>1) printf ( "\nsaved file_num=%u, frame number =%i \n",file_num,  file_num*steps_per_file );
    }
    
    file_num++;
    fluid.WriteSimParams ( outPath ); 
    fluid.WriteGenome( outPath );
  //  fluid.SavePointsCSV2 ( outPath, file_num );                   //fluid.SavePointsCSV ( outPath, 1 );
  //  fluid.SavePointsVTP2 ( outPath, file_num );

    printf ( "\nClosing load_sim.\n" );
    fluid.Exit ();                                                  // Clean up and close
    
    /*
    size_t   free1, free2, total;
    cudaMemGetInfo(&free1, &total);
<<<<<<< HEAD:Morphogenesis/src/load_sim.cpp
    if(debug>0) printf("\nCuda Memory, before cuCtxDestroy(clContext): free=%lu, total=%lu.\t",free1,total);
   // clCheck(clFinish(), "load_sim.cpp ", "clFinish", "before cuCtxDestroy(clContext)", 1/_*mbDebug*_/);
    
    cl_int clResult = cuCtxDestroy ( clContext ) ;
=======
    if(debug>0) printf("\nCuda Memory, before cuCtxDestroy(cuContext): free=%lu, total=%lu.\t",free1,total);
   // clCheck(clFinish(), "load_sim.cpp ", "clFinish", "before cuCtxDestroy(cuContext)", 1/_*mbDebug*_/);  
    
    cl_int clResult = cuCtxDestroy ( cuContext ) ;
>>>>>>> 75eada6585054e07bf9262a150f34af03aa68428:src/load_sim.cpp
    if ( clResult!=0 ) {printf ( "error closing, clResult = %i \n",clResult );}
    
   // clCheck(clFinish(), "load_sim.cpp ", "clFinish", "after cudaDeviceReset()", 1/_*mbDebug*_/); 
    cudaMemGetInfo(&free2, &total);
    if(debug>0) printf("\nAfter cuCtxDestroy(clContext): free=%lu, total=%lu, released=%lu.\n",free2,total,(free2-free1) );
    if(debug>0) printf ( "\nClosing load_sim.\n" );
    */
    return 0;
}

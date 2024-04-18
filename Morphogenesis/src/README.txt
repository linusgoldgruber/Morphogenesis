
| 1. SETUP PIELINE: set the name of your "demo" folder within the Morphogenesis/src/ directory in the config.json file (line 4)


| 2. EXECUTION PIPELINE: createDemo -> checkDemo -> createDemo2

    NB: Creating multiple executables in CMakeLists.txt runs into issues with VTK atm,
        thus example CMakeLists.txt files are provided (see Morphogenesis/src/ direcotry).


| 3. LAUNCH SETTINGS:


createDemo:
-----------------------------------------------
Launch-Configurations:

    Executable:
        build/createDemo

    Arguments:          num_particles, spacing, x_dim, y_dim, z_dim, demoType, simSpace, JSON file

        for example:         20            1       4      4      4       0       6       config.json      (Linus' basic demo)
                 or:        125            1       6      6      6       0       5       config.json      ("free falling" from NH89's github)

----------------------------------------------

checkDemo:
----------------------------------------------
Launch-Configurations:

    Executable:
        build/checkDemo

    Arguments:          JSON file

        for example:   config.json

----------------------------------------------


createDemo2:
----------------------------------------------
Launch-Configurations:

    Executable:
        build/createDemo2

    Arguments:          JSON file

        for example:   config.json

----------------------------------------------


| 4. Environment variables:
----------------------------------------------
COLLECT_GCC_OPTIONS=            '-E' '-v' '-o' '/dev/null' '-mtune=generic' '-march=x86-64'

COMPILER_PATH=                  /usr/lib/gcc/x86_64-linux-gnu/11/:/usr/lib/gcc/x86_64-linux-gnu/11/:/usr/lib/gcc/x86_64-linux-gnu/:/usr/lib/gcc/x86_64-linux-gnu/11/:/usr/lib/gcc/x86_64-linux-gnu/

CPATH=                          /usr/local/cuda/include:/usr/local/include:/usr/include/x86_64-linux-gnu/c++/11:/usr/include/c++/11:

C_INCLUDE_PATH=                 /usr/local/cuda/include:/usr/local/include:/usr/include:

C_PLUS_INCLUDE_PATH=            /usr/local/cuda/include:/usr/local/include:/usr/include/x86_64-linux-gnu/c++/11:/usr/include/c++/11/tr1:/usr/include/c++/11:

LD_LIBRARY_PATH=$LIBRARY_PATH:  /opt/rocm/opencl/lib:/opt/intel/oneapi/lib

LIBRARY_PATH=                   /usr/lib/gcc/x86_64-linux-gnu/11/:/usr/lib/gcc/x86_64-linux-gnu/11/../../../x86_64-linux-gnu/:/usr/lib/gcc/x86_64-linux-gnu/11/../../../../lib/:/lib/x86_64-linux-gnu/:/lib/../lib/:/usr/lib/x86_64-linux-gnu/:/usr/lib/../lib/:/usr/lib/gcc/x86_64-linux-gnu/11/../../../:/lib/:/usr/lib/

-----------------------------------------------











NICK (old)
0) Allocate buffers (might be automatic, check)
    (i)  particles   FELASTIDX, FNERVEIDX, FCONC, FEPIGEN
    (ii) genome
    

1) Initalize correct UID for each particle

2) Initialize Buffer[FELASTIDX],   FNERVEIDX   , FCONC , FEPIGEN

3) UpdateGenome(); //  sends genome to device. // NB need to initialize genome from file, or something.

4) Need script to generate simulations + functions to read them from file. ? what are the available i/o functions ? currently FluidSystem::SavePoints* () 
Also see how files are read in OpenCL-SPH, and how I planned to use Json, Yaml, hdf5 etc...
Probably use a folder with named csv files.

5) 


-----------------------------------------------------------------------------------------------------------------------------------------------------
CMAKE:

Where is the source code: 
/home/goldi/Documents/KDevelop Projects/Morphogenesis/Morphogenesis/src

Where to build the binaries:
/home/goldi/Documents/KDevelop Projects/Morphogenesis/Morphogenesis/src/build

----
load_sim.cpp

Executable:
/home/goldi/Documents/KDevelop Projects/Morphogenesis/Morphogenesis/src/build/load_sim

Arguments:
/home/goldi/Documents/KDevelop\ Projects/Morphogenesis/Morphogenesis/src/demo /home/goldi/Documents/KDevelop\ Projects/Morphogenesis/Morphogenesis/src/demo/out 3 10 n n y y 3 n n /home/goldi/Documents/KDevelop\ Projects/Morphogenesis/Morphogenesis/src/opencl_config.json


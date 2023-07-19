#include <allocator_OpenCL.h>
#include <iostream>
#include <cstdlib>
#include <assert.h>

bool openCLCheck(cl_int launch_stat, const char* obj, const char* method, const char* apicall, const char* arg, bool bDebug)
{
    cl_int kern_stat = CL_SUCCESS;

    if (bDebug) {
        #   since the StartOpenCL function does not have access to a queue variable, it would not be able to call the original version of the openCLCheck function with bDebug set to true. clFinish also blocks until the queue is empty, but has different error codes that it outputs (like CL_SUCCESS)
        kern_stat = clFinish(queue);
    }
    if (kern_stat != CL_SUCCESS || launch_stat != CL_SUCCESS) {
        
        # Assumes getErrorString to be useable by OpenCL
        const char* launch_statmsg = getErrorString(launch_stat);
        const char* kern_statmsg = getErrorString(kern_stat);
        std::cout << "GVDB OpenCL ERROR:\n";
        std::cout << "  Launch status: " << getErrorString(launch_stat) << "\n";
        std::cout << "  Kernel status: " << getErrorString(kern_stat) << "\n";
        std::cout << "  Caller: " << obj << "::" << method << "\n";
        std::cout << "  Call:   " << apicall << "\n";
        std::cout << "  Args:   " << arg << "\n";

        if (bDebug) {
            std::cout << "  Generating assert so you can examine call stack.\n";
            assert(0); // debug - trigger break (see call stack)
        } else {
            std::cout << "Error. Application will exit.\n"; // exit - return 0
        }
        return false;
    }
    return true;
}

void StartOpenCL(int devsel, cl_context ctxsel, cl_device_id& dev, cl_context& ctx, cl_command_queue* strm, bool verbose)
{

    int version = 0;
    char name[128];

    cl_uint cnt = 0;
    cl_device_id dev_id;
    # No need to explicitly initialize in OpenCL (Cuda: cuInit(0);)
        
    //--- List devices
    clGetDeviceIDs(NULL, CL_DEVICE_TYPE_ALL, 0, NULL, &cnt);
    
    if (cnt == 0) {
        std::cout << "ERROR: No OpenCL devices found.\n";
        dev = (int)NULL;
        ctx = NULL;
        std::cout << "Error. Application will exit.\n";
        return;
    }
    if (verbose)
        std::cout << "  Device List:\n";
    std::vector<cl_device_id> devices(cnt);
    clGetDeviceIDs(NULL, CL_DEVICE_TYPE_ALL, cnt, devices.data(), NULL);
    for (int n = 0; n < cnt; n++) {
        dev_id = devices[n];
        clGetDeviceInfo(dev_id, CL_DEVICE_NAME, sizeof(name), name, NULL);

        size_t pi;
        clGetDeviceInfo(dev_id, CL_DEVICE_IMAGE3D_MAX_WIDTH, sizeof(pi), &pi, NULL);
        if (verbose)
            std::cout << "Max. image3D width: " << pi << "\n";

        clGetDeviceInfo(dev_id, CL_DEVICE_IMAGE3D_MAX_HEIGHT, sizeof(pi), &pi, NULL);
        if (verbose)
            std::cout << "Max. image3D height: " << pi << "\n";

        clGetDeviceInfo(dev_id, CL_DEVICE_IMAGE3D_MAX_DEPTH, sizeof(pi), &pi, NULL);
        if (verbose)
            std::cout << "Max. image3D depth: " << pi << "\n";
        if (verbose)
            std::cout << "   " << n << ". " << name << "\n";
    }

    //--- Create new context with Driver API
    dev = devices[devsel];
    ctx = clCreateContext(NULL, 1, &dev, NULL, NULL, NULL);

    clGetDeviceInfo(dev, CL_DEVICE_NAME, sizeof(name), name, NULL);
    if (verbose)
        std::cout << "   Using Device: " << (int)dev << ", " << name << ", Context: " << (void*)ctx << "\n";
}

Vector3DF openCLGetMemUsage()
{
    Vector3DF mem;
    cl_ulong free, total;
    clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(total), &total, NULL);
    clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_AVAILABLE, sizeof(free), &free, NULL);
    free /= (1024.0f * 1024.0f); // MB
    total /= (1024.0f * 1024.0f);
    mem.x = total - free; // used
    mem.y = free;
    mem.z = total;
    return mem;
}

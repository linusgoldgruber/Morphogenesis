#ifndef DEF_ALLOCATOR
#define DEF_ALLOCATOR

#include <vector>
#include <CL/cl.h>
#include "vector.h"

// Maximum number of GVDB Pool levels
#define MAX_POOL 10

// Global OpenCL helpers
#define MIN_RUNTIME_VERSION 4010
#define MIN_COMPUTE_VERSION 0x20
extern void StartOpenCL(int devsel, cl_context ctxsel, cl_device_id &dev, cl_context &ctx, cl_command_queue *strm, bool verbose);
extern bool openclCheck(cl_int e, const char *obj, const char *method, const char *apicall, const char *arg, bool bDebug);
extern Vector3DF openclGetMemUsage();

#endif //DEF_ALLOCATOR
//#include <helper_math.h>

extern "C" {

__global__ void kernel_vexp10(double *v)
{
    int i = blockIdx.x *blockDim.x + threadIdx.x;
    v[i]=exp10(v[i]);
}

}

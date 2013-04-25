Description : Measure the floating point performance of CPU
              with SSE and AVX instruction set

Compilation :

  icc -o measure measure_GFLOPS.cpp -fast -fopenmp -lstdc++ -lrt 
 
  or with PAPI 
  
  icc -o measure measure_GFLOPS.cpp -fast -fopenmp -lstdc++ -DUSE_PAPI -lrt


Note :
 Use directly OpenMP to measure time


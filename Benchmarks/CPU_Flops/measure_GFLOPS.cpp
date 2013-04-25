/**
*   Description : Measure the floating point performance of CPU
*                 with SSE and AVX instruction set
*   Note        :
*                 Found at  http://stackoverflow.com/questions/8389648/ 
*
*   License     : Creative Commons SA (2.5)
*                 http://creativecommons.org/licenses/by-sa/2.5/
*
*   Author      : Alexander J. Yee (2012)
*                 Paul P. Hilscher (2013) - added PAPI interface comparison
*
*
**/ 


#include "xmmintrin.h"
#include "immintrin.h"
#include "stdint.h"
#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include <stdlib.h>
#include <omp.h>
#include <papi.h> 
  
struct Measure {
   float rtime;  ///< Real time
   float ptime;  ///< Processor time
   long long flpops; ///< Total floating point operations
   float mflops;   ///< Million Floating Operations Per Second
} M;

double test_dp_mac_AVX(double x,double y,size_t iterations){

    register __m256d r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,rA,rB,rC,rD,rE,rF;

    //  Generate starting data.
    r0 = _mm256_set1_pd(x);
    r1 = _mm256_set1_pd(y);

    r8 = _mm256_set1_pd(-0.0);

    r2 = _mm256_xor_pd(r0,r8);
    r3 = _mm256_or_pd(r0,r8);
    r4 = _mm256_andnot_pd(r8,r0);
    r5 = _mm256_mul_pd(r1,_mm256_set1_pd(0.37796447300922722721));
    r6 = _mm256_mul_pd(r1,_mm256_set1_pd(0.24253562503633297352));
    r7 = _mm256_mul_pd(r1,_mm256_set1_pd(4.1231056256176605498));
    r8 = _mm256_add_pd(r0,_mm256_set1_pd(0.37796447300922722721));
    r9 = _mm256_add_pd(r1,_mm256_set1_pd(0.24253562503633297352));
    rA = _mm256_sub_pd(r0,_mm256_set1_pd(4.1231056256176605498));
    rB = _mm256_sub_pd(r1,_mm256_set1_pd(4.1231056256176605498));

    rC = _mm256_set1_pd(1.4142135623730950488);
    rD = _mm256_set1_pd(1.7320508075688772935);
    rE = _mm256_set1_pd(0.57735026918962576451);
    rF = _mm256_set1_pd(0.70710678118654752440);

    uint64_t iMASK = 0x800fffffffffffffull;
    __m256d MASK = _mm256_set1_pd(*(double*)&iMASK);
    __m256d vONE = _mm256_set1_pd(1.0);

    // Only inner loop is counted
    size_t c = 0;
    while (c < iterations){
        size_t i = 0;
        while (i < 1000){
            //  Here's the meat - the part that really matters.

            // 48 instructions with 4 flops / instruction
            r0 = _mm256_mul_pd(r0,rC);
            r1 = _mm256_add_pd(r1,rD);
            r2 = _mm256_mul_pd(r2,rE);
            r3 = _mm256_sub_pd(r3,rF);
            r4 = _mm256_mul_pd(r4,rC);
            r5 = _mm256_add_pd(r5,rD);
            r6 = _mm256_mul_pd(r6,rE);
            r7 = _mm256_sub_pd(r7,rF);
            r8 = _mm256_mul_pd(r8,rC);
            r9 = _mm256_add_pd(r9,rD);
            rA = _mm256_mul_pd(rA,rE);
            rB = _mm256_sub_pd(rB,rF);

            r0 = _mm256_add_pd(r0,rF);
            r1 = _mm256_mul_pd(r1,rE);
            r2 = _mm256_sub_pd(r2,rD);
            r3 = _mm256_mul_pd(r3,rC);
            r4 = _mm256_add_pd(r4,rF);
            r5 = _mm256_mul_pd(r5,rE);
            r6 = _mm256_sub_pd(r6,rD);
            r7 = _mm256_mul_pd(r7,rC);
            r8 = _mm256_add_pd(r8,rF);
            r9 = _mm256_mul_pd(r9,rE);
            rA = _mm256_sub_pd(rA,rD);
            rB = _mm256_mul_pd(rB,rC);

            r0 = _mm256_mul_pd(r0,rC);
            r1 = _mm256_add_pd(r1,rD);
            r2 = _mm256_mul_pd(r2,rE);
            r3 = _mm256_sub_pd(r3,rF);
            r4 = _mm256_mul_pd(r4,rC);
            r5 = _mm256_add_pd(r5,rD);
            r6 = _mm256_mul_pd(r6,rE);
            r7 = _mm256_sub_pd(r7,rF);
            r8 = _mm256_mul_pd(r8,rC);
            r9 = _mm256_add_pd(r9,rD);
            rA = _mm256_mul_pd(rA,rE);
            rB = _mm256_sub_pd(rB,rF);

            r0 = _mm256_add_pd(r0,rF);
            r1 = _mm256_mul_pd(r1,rE);
            r2 = _mm256_sub_pd(r2,rD);
            r3 = _mm256_mul_pd(r3,rC);
            r4 = _mm256_add_pd(r4,rF);
            r5 = _mm256_mul_pd(r5,rE);
            r6 = _mm256_sub_pd(r6,rD);
            r7 = _mm256_mul_pd(r7,rC);
            r8 = _mm256_add_pd(r8,rF);
            r9 = _mm256_mul_pd(r9,rE);
            rA = _mm256_sub_pd(rA,rD);
            rB = _mm256_mul_pd(rB,rC);

            i++;
        }

        //  Need to renormalize to prevent denormal/overflow.
        r0 = _mm256_and_pd(r0,MASK);
        r1 = _mm256_and_pd(r1,MASK);
        r2 = _mm256_and_pd(r2,MASK);
        r3 = _mm256_and_pd(r3,MASK);
        r4 = _mm256_and_pd(r4,MASK);
        r5 = _mm256_and_pd(r5,MASK);
        r6 = _mm256_and_pd(r6,MASK);
        r7 = _mm256_and_pd(r7,MASK);
        r8 = _mm256_and_pd(r8,MASK);
        r9 = _mm256_and_pd(r9,MASK);
        rA = _mm256_and_pd(rA,MASK);
        rB = _mm256_and_pd(rB,MASK);
        r0 = _mm256_or_pd(r0,vONE);
        r1 = _mm256_or_pd(r1,vONE);
        r2 = _mm256_or_pd(r2,vONE);
        r3 = _mm256_or_pd(r3,vONE);
        r4 = _mm256_or_pd(r4,vONE);
        r5 = _mm256_or_pd(r5,vONE);
        r6 = _mm256_or_pd(r6,vONE);
        r7 = _mm256_or_pd(r7,vONE);
        r8 = _mm256_or_pd(r8,vONE);
        r9 = _mm256_or_pd(r9,vONE);
        rA = _mm256_or_pd(rA,vONE);
        rB = _mm256_or_pd(rB,vONE);

        c++;
    }

    r0 = _mm256_add_pd(r0,r1);
    r2 = _mm256_add_pd(r2,r3);
    r4 = _mm256_add_pd(r4,r5);
    r6 = _mm256_add_pd(r6,r7);
    r8 = _mm256_add_pd(r8,r9);
    rA = _mm256_add_pd(rA,rB);

    r0 = _mm256_add_pd(r0,r2);
    r4 = _mm256_add_pd(r4,r6);
    r8 = _mm256_add_pd(r8,rA);

    r0 = _mm256_add_pd(r0,r4);
    r0 = _mm256_add_pd(r0,r8);

    //  Prevent Dead Code Elimination
    double out = 0;
    out += ((double*)&r0)[0];
    out += ((double*)&r0)[1];
    out += ((double*)&r0)[2];
    out += ((double*)&r0)[3];

    return out;
}


double test_dp_mac_SSE(double x,double y,size_t iterations){
    register __m128d r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,rA,rB,rC,rD,rE,rF;

    //  Generate starting data.
    r0 = _mm_set1_pd(x);
    r1 = _mm_set1_pd(y);

    r8 = _mm_set1_pd(-0.0);

    r2 = _mm_xor_pd(r0,r8);
    r3 = _mm_or_pd(r0,r8);
    r4 = _mm_andnot_pd(r8,r0);
    r5 = _mm_mul_pd(r1,_mm_set1_pd(0.37796447300922722721));
    r6 = _mm_mul_pd(r1,_mm_set1_pd(0.24253562503633297352));
    r7 = _mm_mul_pd(r1,_mm_set1_pd(4.1231056256176605498));
    r8 = _mm_add_pd(r0,_mm_set1_pd(0.37796447300922722721));
    r9 = _mm_add_pd(r1,_mm_set1_pd(0.24253562503633297352));
    rA = _mm_sub_pd(r0,_mm_set1_pd(4.1231056256176605498));
    rB = _mm_sub_pd(r1,_mm_set1_pd(4.1231056256176605498));

    rC = _mm_set1_pd(1.4142135623730950488);
    rD = _mm_set1_pd(1.7320508075688772935);
    rE = _mm_set1_pd(0.57735026918962576451);
    rF = _mm_set1_pd(0.70710678118654752440);

    uint64_t iMASK = 0x800fffffffffffffull;
    __m128d MASK = _mm_set1_pd(*(double*)&iMASK);
    __m128d vONE = _mm_set1_pd(1.0);

    size_t c = 0;
    while (c < iterations){
        size_t i = 0;
        while (i < 1000){
            //  Here's the meat - the part that really matters.

            // 48 instructions with 2 flops/instruction
            r0 = _mm_mul_pd(r0,rC);
            r1 = _mm_add_pd(r1,rD);
            r2 = _mm_mul_pd(r2,rE);
            r3 = _mm_sub_pd(r3,rF);
            r4 = _mm_mul_pd(r4,rC);
            r5 = _mm_add_pd(r5,rD);
            r6 = _mm_mul_pd(r6,rE);
            r7 = _mm_sub_pd(r7,rF);
            r8 = _mm_mul_pd(r8,rC);
            r9 = _mm_add_pd(r9,rD);
            rA = _mm_mul_pd(rA,rE);
            rB = _mm_sub_pd(rB,rF);

            r0 = _mm_add_pd(r0,rF);
            r1 = _mm_mul_pd(r1,rE);
            r2 = _mm_sub_pd(r2,rD);
            r3 = _mm_mul_pd(r3,rC);
            r4 = _mm_add_pd(r4,rF);
            r5 = _mm_mul_pd(r5,rE);
            r6 = _mm_sub_pd(r6,rD);
            r7 = _mm_mul_pd(r7,rC);
            r8 = _mm_add_pd(r8,rF);
            r9 = _mm_mul_pd(r9,rE);
            rA = _mm_sub_pd(rA,rD);
            rB = _mm_mul_pd(rB,rC);

            r0 = _mm_mul_pd(r0,rC);
            r1 = _mm_add_pd(r1,rD);
            r2 = _mm_mul_pd(r2,rE);
            r3 = _mm_sub_pd(r3,rF);
            r4 = _mm_mul_pd(r4,rC);
            r5 = _mm_add_pd(r5,rD);
            r6 = _mm_mul_pd(r6,rE);
            r7 = _mm_sub_pd(r7,rF);
            r8 = _mm_mul_pd(r8,rC);
            r9 = _mm_add_pd(r9,rD);
            rA = _mm_mul_pd(rA,rE);
            rB = _mm_sub_pd(rB,rF);

            r0 = _mm_add_pd(r0,rF);
            r1 = _mm_mul_pd(r1,rE);
            r2 = _mm_sub_pd(r2,rD);
            r3 = _mm_mul_pd(r3,rC);
            r4 = _mm_add_pd(r4,rF);
            r5 = _mm_mul_pd(r5,rE);
            r6 = _mm_sub_pd(r6,rD);
            r7 = _mm_mul_pd(r7,rC);
            r8 = _mm_add_pd(r8,rF);
            r9 = _mm_mul_pd(r9,rE);
            rA = _mm_sub_pd(rA,rD);
            rB = _mm_mul_pd(rB,rC);

            i++;
        }

        //  Need to renormalize to prevent denormal/overflow.
        r0 = _mm_and_pd(r0,MASK);
        r1 = _mm_and_pd(r1,MASK);
        r2 = _mm_and_pd(r2,MASK);
        r3 = _mm_and_pd(r3,MASK);
        r4 = _mm_and_pd(r4,MASK);
        r5 = _mm_and_pd(r5,MASK);
        r6 = _mm_and_pd(r6,MASK);
        r7 = _mm_and_pd(r7,MASK);
        r8 = _mm_and_pd(r8,MASK);
        r9 = _mm_and_pd(r9,MASK);
        rA = _mm_and_pd(rA,MASK);
        rB = _mm_and_pd(rB,MASK);
        r0 = _mm_or_pd(r0,vONE);
        r1 = _mm_or_pd(r1,vONE);
        r2 = _mm_or_pd(r2,vONE);
        r3 = _mm_or_pd(r3,vONE);
        r4 = _mm_or_pd(r4,vONE);
        r5 = _mm_or_pd(r5,vONE);
        r6 = _mm_or_pd(r6,vONE);
        r7 = _mm_or_pd(r7,vONE);
        r8 = _mm_or_pd(r8,vONE);
        r9 = _mm_or_pd(r9,vONE);
        rA = _mm_or_pd(rA,vONE);
        rB = _mm_or_pd(rB,vONE);

        c++;
    }

    r0 = _mm_add_pd(r0,r1);
    r2 = _mm_add_pd(r2,r3);
    r4 = _mm_add_pd(r4,r5);
    r6 = _mm_add_pd(r6,r7);
    r8 = _mm_add_pd(r8,r9);
    rA = _mm_add_pd(rA,rB);

    r0 = _mm_add_pd(r0,r2);
    r4 = _mm_add_pd(r4,r6);
    r8 = _mm_add_pd(r8,rA);

    r0 = _mm_add_pd(r0,r4);
    r0 = _mm_add_pd(r0,r8);

    //  Prevent Dead Code Elimination
    double out = 0;
    out += ((double*)&r0)[0];
    out += ((double*)&r0)[1];

    return out;
}

void test_dp(int tds,size_t iterations, int pc){

    double *sum = (double*) malloc(tds * sizeof(double));

    // get time 
    timespec ts_start, ts_end;
      
    clock_gettime(CLOCK_REALTIME, &ts_start); 
   

#ifdef USE_PAPI
    if(PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT)    std::cout << "ERROR" << std::endl;

    int num_hwcntrs = 0;
    if (( num_hwcntrs = PAPI_num_counters()) < PAPI_OK) std::cout << "ERROR" << std::endl;
    std::cout << "HW counters : " << num_hwcntrs << std::endl;

    int events[] = { PAPI_DP_OPS, PAPI_VEC_DP, PAPI_FP_OPS };
    int num      = sizeof(events)/sizeof(int);

    PAPI_start_counters(events, num);
    long long start_usec( PAPI_get_real_usec() );
#endif // USE_PAPI

    // Main benchmark part
#pragma omp parallel num_threads(tds)
    {
        double ret;
        
        if     (pc == 2) ret = test_dp_mac_SSE(1.1,2.1,iterations);
        else if(pc == 4) ret = test_dp_mac_AVX(1.1,2.1,iterations);
        sum[omp_get_thread_num()] = ret;
    }

#ifdef USE_PAPI
    long long end_usec( PAPI_get_real_usec() );
    long long values[num];
    if(PAPI_read_counters(values, num) != PAPI_OK) std::cout << "Error" << std::endl;
   
    double dtime = (end_usec - start_usec) * 1.e-6;
    for(int n = 0; n < num; n++)  std::cout << "Values : " << values[n]/1.e9 << " " << values[n]*1.e-9/dtime  << std::endl;
    if(PAPI_stop_counters(values, num) != PAPI_OK) std::cout << "Error" << std::endl;       
#endif // USE_PAPI

    //PAPI_flops( &M.rtime, &M.ptime, &M.flpops, &M.mflops);
    clock_gettime(CLOCK_REALTIME, &ts_end); 
    
    double secs  = (ts_end.tv_sec - ts_start.tv_sec);
    double nsecs = (ts_end.tv_nsec - ts_start.tv_nsec);
	 
    double time = secs + nsecs * 1.e-9;
    
    double gflop = 48 * 1000 * iterations * tds * pc * 1.e-9;
    std::cout << gflop << " GFLOP in " << secs << " seconds. GFLOPS = "  
              << std::setprecision(5) << std::setw(5) << gflop / time << std::endl;
    
     // as measured by PAPI
    //std::cout << " PAPI : " << M.flpops/1.e9 << " GFLOP." << "GFLOPS(1) : " << M.mflops/1.e3 << "GFLOPS. GFLOPS = " << (M.flpops/1.e9/time)  << std::endl;

    double out = 0;
    int c = 0;
    while (c < tds) out += sum[c++];
   
    std::cout << "sum = " << out << std::endl;
    std::cout << std::endl;

    free(sum);
}

int main()
{
  std::cout << "Benchmark SSE" << std::endl; test_dp(1,1000000, 2);
  std::cout << "Benchmark AVX" << std::endl; test_dp(1,1000000, 4);

}

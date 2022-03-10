//
// Fast Cross-Correlation algorithm
//
// Author: Marc-Antoine Maheux
// Email: marc-antoine.maheux@usherbrooke.ca
//
// Refer to the following paper for details:
//
// Grondin, F., Maheux, M.-A., Lauzon, J.-S. and Michaud, F., "Fast Cross-Correlation
// for TDoA Estimation on Small Aperture Microphone Arrays", arXiV
//

#ifndef __FCCPHAT_SIMD
#define __FCCPHAT_SIMD

#ifdef FCCPHAT_USE_SIMD

    #ifdef __cplusplus
    extern "C" {
    #endif


    #define FLOAT_SIMD_ALIGNMENT 32
    #define FLOAT_SIMD_ALIGNMENT_ATTRIBUTE __attribute__((aligned(FLOAT_SIMD_ALIGNMENT)))
    #define FLOAT8_SIZE 8

    #ifdef __clang__

        typedef float float8_t __attribute__ ((ext_vector_type (8)));
        typedef float float4_t __attribute__((vector_size (4 * sizeof(float))));

        #define SHUFFLE_FLOAT8_T(dest, src, i0, i1, i2, i3, i4, i5, i6, i7) \
            dest = __builtin_shufflevector((src), (src), (i0), (i1), (i2), (i3), (i4), (i5), (i6), (i7))

    #elif __GNUC__

        typedef float float8_t __attribute__((vector_size (8 * sizeof(float))));
        typedef float float4_t __attribute__((vector_size (4 * sizeof(float))));
        typedef int32_t int32x8_t __attribute__((vector_size (8 * sizeof(int32_t))));

        #define SHUFFLE_FLOAT8_T(dest, src, i0, i1, i2, i3, i4, i5, i6, i7) do {   \
                int32x8_t mask = {(i0), (i1), (i2), (i3), (i4), (i5), (i6), (i7)}; \
                (dest) = __builtin_shuffle((src), mask);                           \
            } while(0)

    #else

        #error "Vector extension is not supported by the compiler."

    #endif

    typedef union float_vector_union {
        float8_t f8;
        float4_t f4[2];
    } float_vector_union;

    #define SUM_8(vv) ((vv)[0] + (vv)[1] + (vv)[2] + (vv)[3] + (vv)[4] + (vv)[5] + (vv)[6] + (vv)[7]);

    #ifdef __cplusplus
    }
    #endif



#endif



#endif

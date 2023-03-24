/**
* @file dft.h
* @author Melih Altun @2015
**/

#include <string.h>
#include <stdlib.h>

#ifndef DFT_INCLUDE
#define DFT_INCLUDE

#define PI 3.141592653589793

//#define PHASE
//#define MAGNITUDE

#define MAX_DFT_SIZE 65535

typedef enum
{
	DFT_NO_ERROR,
	DFT_NULL_PTR_ERROR,
	DFT_COMPUTATION_ERROR
} dftReturn;

typedef struct _dft_instance{
	size_t dftSize;
	double *Re;
	double *Im;
#if defined (MAGNITUDE)
	double *abs;
#endif
#if defined (PHASE)
	double *angle;
#endif
}dft_instance;

/* Memory allocations */
dftReturn set_dft_instance(dft_instance *dft, int size);

/* DFT with real input */
dftReturn dft(dft_instance *inst, double x1[], int clk, size_t n, size_t m);

/* DFT with complex input */
dftReturn dft_complex(dft_instance *inst, double xr[], double xi[], int clk, size_t n, size_t m);

/* DFT of certain selected frequency bins for real inputs */
dftReturn dftSelectedFreq(dft_instance *inst, double x[], int f[], int clk, size_t N, size_t F);

/* DFT of certain selected frequency bins for complex inputs */
dftReturn dftSelectedFreq_complex(dft_instance *inst, double xr[], double xi[], int f[], int clk, size_t N, size_t F);

/* Inverse DFT with a real output */
dftReturn idft_real(double x[], dft_instance *inst, int clk, size_t N, size_t M);

/* Inverse DFT with complex outputs */
dftReturn idft_complex(double xr[], double xi[], dft_instance *inst, int clk, size_t N, size_t M);

/* Frequency bin of the fundemental frequency */
dftReturn dft_w0(int *maxInd, dft_instance *inst);

/* Memory deallocations*/
dftReturn delete_dft_instance (dft_instance *dft);
#endif

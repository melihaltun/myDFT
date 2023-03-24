/**
* @file dft.c
* @author Melih Altun @2015
**/

#include <math.h>
#include "dft.h"

// auxilary function: replacement for % operator, since % operator can be problematic when inputs are negative
static int mod(int a, int b) {
	if (b == 0)
		return -1;
	if (b < 0)
		return mod(-a, -b);
	int c = a % b;
	return (c < 0) ? c + b : c;
}

/* allocates memory for DFT operations */
dftReturn set_dft_instance(dft_instance *inst, int size)
{
	if (inst == NULL)
		return DFT_NULL_PTR_ERROR;
	if (size < 1 || size > MAX_DFT_SIZE)
		return DFT_COMPUTATION_ERROR;
	inst->dftSize = size;
	inst->Re = calloc(size, sizeof(double));
	inst->Im = calloc(size, sizeof(double));
#if defined (MAGNITUDE)
	inst->abs = calloc(size, sizeof(double));
#endif
#if defined (PHASE)
	inst->angle = calloc(size, sizeof(double));
#endif
	return DFT_NO_ERROR;
}

/* Compute the first M fft coeficients of real series x[] of lenght N. x[] can be a circular buffer with clk pointing to the buffer index. Set clk to 0 if x[] is not a circular buffer.
 For real valued inputs, calculating first N/2 coefficients is sufficient for obtaining a full dft:
 if x = real{x} -> real{fft{x}} is even symmetric and imag{fft{x}} is odd symmetric. */
dftReturn dft_real(dft_instance *inst, double x[], int clk, size_t N, size_t M)
{
	size_t k, n;
	double arg;
	int i;

	if (inst == NULL || x == NULL)
		return DFT_NULL_PTR_ERROR;
	if (N < 1 || N > MAX_DFT_SIZE || M < 1)
		return DFT_COMPUTATION_ERROR;

	if (M > inst->dftSize)
		M = inst->dftSize;  //M cannot be larger than dft size

	//X[k] = Sum (n = 0 -> N-1) x[n] * e^-(j 2 pi k n / N)
	//Re{X[k]} = Sum (n = 0 -> N-1) x[n] * cos(2 pi k n / N)
	//Im{X[k]} = Sum (n = 0 -> N-1) x[n] * -sin(2 pi k n / N)
	for (k = 0; k < M; k++) {
		inst->Re[k] = 0;
		inst->Im[k] = 0;
		arg = 2.0f * (double)PI * (double)k / (double)N;
		for (n = 0; n < N; n++) {
			i = mod(clk - ((int)N - n), (int)N); // % operator is problematic when inputs are negative
			inst->Re[k] += (x[i] * cosf(n * arg));
			inst->Im[k] -= (x[i] * sinf(n * arg));
		}
#if defined (MAGNITUDE)
		inst->abs[k] = sqrtf(inst->Re[k] * inst->Re[k] + inst->Im[k] * inst->Im[k]);
#endif
#if defined (PHASE)
		inst->angle[k] = atan2(inst->Im[k], inst->Re[k]);
#endif
	}
	return DFT_NO_ERROR;
}

/* Compute the first M fft coeficients of  complex series  x = xr[] + j xi[] of lenght N. xr[] and xi[] can be a circular buffers with clk pointing to the buffer index
Set clk to 0 if xr[] and xi[] are not a circular buffers. */
dftReturn dft_complex(dft_instance *inst, double xr[], double xi[], int clk, size_t N, size_t M)
{
	size_t k, n;
	double arg;
	int i;

	if (inst == NULL || xr == NULL || xi == NULL)
		return DFT_NULL_PTR_ERROR;
	if (N < 1 || N > MAX_DFT_SIZE || M < 1)
		return DFT_COMPUTATION_ERROR;

	if (M > inst->dftSize)
		M = inst->dftSize;  //M cannot be larger than dft size

	//X[k] = Sum (n = 0 -> N-1) x[n] * e^-(j 2 pi k n / N)
	//Re{X[k]} = Sum (n = 0 -> N-1) Re{x[n]} * cos(2 pi k n / N) + j^2 * Im{x[n]} * -sin(2 pi k n / N)
	//Im{X[k]} = Sum (n = 0 -> N-1) Im{x[n]} * cos(2 pi k n / N) + Re{x[n]} * -sin(2 pi k n / N)
	for (k = 0; k < M; k++) {
		inst->Re[k] = 0;
		inst->Im[k] = 0;
		arg = (double)(2.0 * PI * (double)k / (double)N);
		for (n = 0; n < N; n++) {
			i = mod(clk - ((int)N - n), (int)N);
			inst->Re[k] += (xr[i] * cosf(n * arg)) + (xi[i] * sinf(n * arg));
			inst->Im[k] += (xi[i] * cosf(n * arg)) - (xr[i] * sinf(n * arg));
		}
#if defined (MAGNITUDE)
		inst->abs[k] = (double)sqrt(inst->Re[k] * inst->Re[k] + inst->Im[k] * inst->Im[k]);
#endif
#if defined (PHASE)
		inst->angle[k] = atan2(inst->Im[k], inst->Re[k]);
#endif
	}
	return DFT_NO_ERROR;
}

/* Compute dft coeficients belonging to specficied frequencies in f[] of length F from series
 x[] of lenght N. x1[] can be a circular buffer with clk being the first buffer index */
dftReturn dftSelectedFreq(dft_instance *inst, double x[], int f[], int clk, size_t N, size_t F)
{
	size_t k, n;
	double arg;
	int i;

	if (inst == NULL || x == NULL || f == NULL)
		return DFT_NULL_PTR_ERROR;
	if (N < 1 || N > MAX_DFT_SIZE || F < 1 || F > MAX_DFT_SIZE )
		return DFT_COMPUTATION_ERROR;

	for (k = 0; k < F; k++) {
		inst->Re[k] = 0;
		inst->Im[k] = 0;
		if (f[k] > (int)inst->dftSize)
			return DFT_COMPUTATION_ERROR;  // frequency index cannot be larger than dft size.
		arg = (double)(2.0 * PI * (double)f[k] / (double)N);
		for (n = 0; n < N; n++) {
			i = mod(clk - ((int)N - n), (int)N);
			inst->Re[k] += (x[i] * cosf(n * arg));
			inst->Im[k] -= (x[i] * sinf(n * arg));
		}
#if defined (MAGNITUDE)
		inst->abs[k] = (double)sqrt(inst->Re[k] * inst->Re[k] + inst->Im[k] * inst->Im[k]);
#endif
#if defined (PHASE)
		inst->angle[k] = atan2(inst->Im[k], inst->Re[k]);
#endif
	}
	return DFT_NO_ERROR;
}

/* Compute dft coeficients belonging to specficied frequencies in f[] of length F from series
 x = xr[] + xi[] of lenght N. xr[] and xi[] can be a circular buffers with clk pointing to the buffer index 
 Set clk to 0 if xr[] and xi[] are not a circular buffers.*/
dftReturn dftSelectedFreq_complex(dft_instance *inst, double xr[], double xi[], int f[], int clk, size_t N, size_t F)
{
	size_t k, n;
	double arg;
	int i;

	if (inst == NULL || xr == NULL || xi == NULL || f == NULL)
		return DFT_NULL_PTR_ERROR;
	if (N < 1 || N > MAX_DFT_SIZE || F < 1 || F > MAX_DFT_SIZE)
		return DFT_COMPUTATION_ERROR;

	for (k = 0; k < F; k++) {
		inst->Re[k] = 0;
		inst->Im[k] = 0;
		if (f[k] > (int)inst->dftSize)
			return DFT_COMPUTATION_ERROR;  // frequency index cannot be larger than dft size.
		arg = (double)(2.0 * PI * (double)f[k] / (double)N);
		for (n = 0; n < N; n++) {
			i = mod(clk - ((int)N - n), (int)N);
			inst->Re[k] += (xr[i] * cosf(n * arg)) + (xi[i] * sinf(n * arg));
			inst->Im[k] += (xi[i] * cosf(n * arg)) - (xr[i] * sinf(n * arg));
		}
#if defined (MAGNITUDE)
		inst->abs[k] = (double)sqrt(inst->Re[k] * inst->Re[k] + inst->Im[k] * inst->Im[k]);
#endif
#if defined (PHASE)
		inst->angle[k] = atan2(inst->Im[k], inst->Re[k]);
#endif
	}
	return DFT_NO_ERROR;
}

/* inverse dft function to compute abs(x) from FT{X}
 Function uses first M FT coefficients from a series of size N to reconstruct x */
dftReturn idft_real(double x[], dft_instance *inst, int clk, size_t N, size_t M)
{
	size_t k, n;
	double arg, imgPart;
	int i;

	if (inst == NULL || x == NULL )
		return DFT_NULL_PTR_ERROR;
	if (N < 1 || N > MAX_DFT_SIZE || M < 1 )
		return DFT_COMPUTATION_ERROR;

	if (M > inst->dftSize)
		M = inst->dftSize;  //M cannot be larger than dft size

	//x[n] = 1/N Sum (k = 0 -> N-1) X[k] e^(j 2 pi k n / N)
	//x[n] = 1/N Sum (n = 0 -> N-1) X[k] (cos(2 pi k n / N) + j sin(2 pi k n / N))
	//x[n] = 1/N Sum (n = 0 -> N-1) re{X[k]} cos(2 pi k n / N) - im{X[k]} sin(2 pi k n / N)) + j re{X[k]} sin(2 pi k n / N) + j im{X[k]} cos(2 pi k n / N)

	for (n = 0; n < N; n++) {
		i = mod(clk - ((int)N - n), (int)N); // % operator is problematic when inputs are negative
		x[i] = 0;
		imgPart = 0;
		arg = (double)(2.0 * PI * (double)n / (double)N);
		for (k = 0; k < M; k++) {
			x[i] += 1 / (double)N * (inst->Re[k] * cosf(k * arg) - inst->Im[k] * sinf(k * arg));
			imgPart += 1 / (double)N * (inst->Re[k] * sinf(k * arg) + inst->Im[k] * cosf(k * arg));
		}
	}
	return DFT_NO_ERROR;
}


/* inverse dft function to compute x = re{x} + j im{x} from FT{X}
Function uses first M FT coefficients from a series of size N to reconstruct x */
dftReturn idft_complex(double xr[], double xi[], dft_instance *inst, int clk, size_t N, size_t M)
{
	size_t k, n;
	double arg;
	int i;

	if (inst == NULL || xr == NULL || xi == NULL)
		return DFT_NULL_PTR_ERROR;
	if (N < 1 || N > MAX_DFT_SIZE || M < 1 )
		return DFT_COMPUTATION_ERROR;

	if (M > inst->dftSize)
		M = inst->dftSize;  //M cannot be larger than dft size

	//x[n] = 1/N Sum (k = 0 -> N-1) X[k] e^(j 2 pi k n / N)
	//x[n] = 1/N Sum (n = 0 -> N-1) X[k] (cos(2 pi k n / N) + j sin(2 pi k n / N))
	//x[n] = 1/N Sum (n = 0 -> N-1) re{X[k]} cos(2 pi k n / N) - im{X[k]} sin(2 pi k n / N)) + j re{X[k]} sin(2 pi k n / N) + j im{X[k]} cos(2 pi k n / N)

	for (n = 0; n < N; n++) {
		i = mod(clk - ((int)N - n), (int)N);
		xr[i] = 0;
		xi[i] = 0;
		arg = (double)(2.0 * PI * (double)n / (double)N);
		for (k = 0; k < M; k++) {
			xr[i] += 1 / (double)N * (inst->Re[k] * cosf(k * arg) - inst->Im[k] * sinf(k * arg));
			xi[i] += 1 / (double)N * (inst->Re[k] * sinf(k * arg) + inst->Im[k] * cosf(k * arg));
		}
	}
	return DFT_NO_ERROR;
}

/*returns the frequency index with max magnitute value - i.e. fundemental frequency */
dftReturn dft_w0(int *maxInd, dft_instance *inst)
{
	size_t imax = 1, i;  //imax=1 -> ignore DC

	if (inst == NULL || maxInd == NULL)
		return DFT_NULL_PTR_ERROR;

	for (i = 2; i < inst->dftSize; i++)
#if defined (MAGNITUDE)
		imax = (inst->abs[i] > inst->abs[imax]) ? i : imax;
#else
		imax = ((inst->Re[i]*inst->Re[i]) + (inst->Im[i] * inst->Im[i]) > (inst->Re[imax] * inst->Re[imax]) + (inst->Im[imax] * inst->Im[imax])) ? i : imax;
#endif
	*maxInd = imax;
	return DFT_NO_ERROR;
}

/* releases memory */
dftReturn delete_dft_instance(dft_instance* dft)
{
	if (dft == NULL)
		return DFT_NULL_PTR_ERROR;
	if (dft->Re != NULL) {
		free(dft->Re);
		dft->Re = NULL;
	}
	if (dft->Im != NULL) {
		free(dft->Im);
		dft->Im = NULL;
	}
#if defined (MAGNITUDE)
	if (dft->abs != NULL) {
		free(dft->abs);
		dft->abs = NULL;
	}
#endif
#if defined (PHASE)
	if (dft->angle != NULL) {
		free(dft->angle);
		dft->angle = NULL;
	}
#endif
	dft->dftSize = 0;
	return DFT_NO_ERROR;
}

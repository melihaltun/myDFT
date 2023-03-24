/**
* @file dft_test.c
* @author Melih Altun @2023
**/

#include "dft.h"

void main() {
	dft_instance dft1;
	double a[] = { 3.0, -6.0, 4.0, 8.0, 2.0, 5.0, 4.0, -2.0, 1.0, -4.0 };
	int selected_f[] = { 1, 3 };
	double f1r, f1i, f3r, f3i, Ar[10], Ai[10], b[10];

	set_dft_instance(&dft1, 10);

	dftSelectedFreq(&dft1, a, selected_f, 0, 10, 2);

	f1r = dft1.Re[0];   // should be -15.2532889043741
	f1i = dft1.Im[0];   // should be -10.0125937026671j
	f3r = dft1.Re[1];   // should be 3.75328890437411 
	f3i = dft1.Im[1];   // should be 11.4454343449828j

	dft_real(&dft1, a, 0, 10, 10);
	memcpy(Ar, dft1.Re, 10 * sizeof(double));
	memcpy(Ai, dft1.Im, 10 * sizeof(double));

	//Expected results
	//Ar = { 15.0, -15.2532889043741, -2.13525491562421, 3.75328890437411, 14.6352549156242, 13.0, 14.6352549156242, 3.75328890437411, -2.13525491562421, -15.2532889043741 }
	//Ai = { 0.0, -10.0125937026671, 4.11449676604731, 11.4454343449828, -6.65739561406607, 0.0, 6.65739561406607, -11.4454343449828, -4.11449676604731, 10.0125937026671 }

	idft_real(b, &dft1, 0, 10, 10);

	delete_dft_instance(&dft1);
}


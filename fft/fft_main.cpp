#include <stdio.h>
#include <stdlib.h>

#include "fft_test.h"
#include "fft_core.h"

int main(void)
{
	setvbuf(stdout, NULL, _IONBF, 0);

	testForwardAndBackward2D();
	/*complex<double> test_fn_1[] =
	{
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
	};
	complex<double> test_fn_1_spectrum[32];

	discreteFourierFast(test_fn_1, 32, test_fn_1_spectrum, ftdFunctionToSpectrum);
	discreteFourierFast(test_fn_1_spectrum, 32, test_fn_1, ftdSpectrumToFunction);

	for (int i = 0; i < 32; i++)
	{
		printf("F%d = (r:%.3f, i:%.3f);\n", i, test_fn_1[i].real(), test_fn_1[i].imag());
	}*/


	printf("       N\t time_triv\t  time_fft\tvel diff by\t precision\n");
	compareVelocityDiscreteFourier();
	testVelocityDiscreteFourierFast();

	return EXIT_SUCCESS;
}

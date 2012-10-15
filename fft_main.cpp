#include <stdio.h>
#include <stdlib.h>

#include "fft_test.h"
#include "fft_core.h"

int main(void)
{
	setvbuf(stdout, NULL, _IONBF, 0);

	printf("       N\t time_triv\t  time_fft\tvel diff by\t precision\n");
	compareVelocityDiscreteFourier();
	testVelocityDiscreteFourierFast();
	printf("\n");
	testForwardAndBackward2D();

	return EXIT_SUCCESS;
}

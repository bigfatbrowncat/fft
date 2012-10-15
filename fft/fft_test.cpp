#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "fft_core.h"

#include <complex>

using namespace std;

const complex<float> I(0, 1);
const float PI = 3.1415926535;

void discreteFourierTrivial(const complex<float>* fn, int N, complex<float>* Fn, fourier_transform_direction ftd)
{
	if (N <= 0) throw INCORRECT_SPECTRUM_SIZE_FOR_FFT;

	float norm, exp_dir;
	switch (ftd)
	{
	case ftdFunctionToSpectrum:
		norm = 1;
		exp_dir = -1;
		break;
	case ftdSpectrumToFunction:
		norm = 1.0 / N;
		exp_dir = 1;
		break;
	default:
		throw UNSUPPORTED_FTD;
	}

	for (int n = 0; n < N; n++)
	{
		Fn[n] = 0;
		for (int k = 0; k < N; k++)
		{
			complex<float> t = (exp_dir * 2*PI/N) * k*n * I;
			Fn[n] += fn[k] * exp(t);
		}
		Fn[n] *= norm;
	}
}

void compareVelocityDiscreteFourierCycle(int k)
{
	int N = 1 << k;

	complex<float>* test_fn = new complex<float>[N];
	complex<float>* test_fn_trivial_spectrum = new complex<float>[N];
	complex<float>* test_fn_fast_spectrum = new complex<float>[N];

	for (int i = 0; i < N / 2; i++)
	{
		test_fn[i] = 1;
		test_fn[N / 2 + i] = -1;
	}

	float t_triv = 0, t_fft = 0;
	clock_t before, after_fast, tmp;

	int rep;

	before = clock();
	rep = 0;
	do
	{
		discreteFourierFast(test_fn, N, test_fn_fast_spectrum, ftdFunctionToSpectrum);
		tmp = clock();
		rep ++;
	}
	while ((tmp - before) / CLOCKS_PER_SEC < 0.1);
	after_fast = tmp;
	t_fft = ((float)after_fast - before) / rep / CLOCKS_PER_SEC * 1000;

	before = clock();
	rep = 0;
	do
	{
		discreteFourierTrivial(test_fn, N, test_fn_trivial_spectrum, ftdFunctionToSpectrum);
		tmp = clock();
		rep ++;
	}
	while ((tmp - before) / CLOCKS_PER_SEC < 0.1);
	after_fast = tmp;
	t_triv = ((float)after_fast - before) / rep / CLOCKS_PER_SEC * 1000;

	// Calculating dispersion
	float Dr = 0, Di = 0;
	for (int i = 0; i < N; i++)
	{
		float ddr = real(test_fn_fast_spectrum[i]) - real(test_fn_trivial_spectrum[i]);
		float ddi = imag(test_fn_fast_spectrum[i]) - imag(test_fn_trivial_spectrum[i]);
		ddr *= ddr; ddi *= ddi;
		Dr += ddr; Di += ddi;
	}

	float D = sqrt(max(Dr, Di) / N);
	float prec = -round(log10(D)) - 1;

	if (t_triv == 0 && t_fft == 0)
	{
		printf("%8d\t%10.3f\t%10.3f\t<indefinite>", N, t_triv, t_fft);
	}
	else if (t_triv != 0 && t_fft != 0)
	{
		printf("%8d\t%10.3f\t%10.3f\t%11.2f", N, t_triv, t_fft, (float)t_triv / t_fft);
	}
	else
	{
		printf("%8d\t%10.3f\t%10.3f\t<too fast>", N, t_triv, t_fft);
	}

	printf("\t%*.0f\n", 10, prec);

	delete [] test_fn;
	delete [] test_fn_trivial_spectrum;
}

void testVelocityDiscreteFourierFastCycle(int k)
{
	int N = 1 << k;

	complex<float>* test_fn = new complex<float>[N];
	complex<float>* test_fn_spectrum = new complex<float>[N];

	for (int i = 0; i < N / 2; i++)
	{
		test_fn[i] = 1;
		test_fn[N / 2 + i] = -1;
	}

	float t_fft = 0;
	clock_t before, after_fast;

	before = clock();
	discreteFourierFast(test_fn, N, test_fn_spectrum, ftdFunctionToSpectrum);
	after_fast = clock();

	t_fft = after_fast - before;

	printf("%8d\t<too slow>\t%10.3f\t      <n/a>\t     <n/a>\n", N, t_fft);


	delete [] test_fn;
	delete [] test_fn_spectrum;
}

void compareVelocityDiscreteFourier()
{
	for (int k = 3; k < 13; k++)
	{
		compareVelocityDiscreteFourierCycle(k);
	}
}

void testVelocityDiscreteFourierFast()
{
	for (int k = 13; k < 23; k++)
	{
		testVelocityDiscreteFourierFastCycle(k);
	}
}

void test_discreteFourierTrivial()
{
	complex<float> test_fn_1[] =
	{
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
	};
	complex<float> test_fn_1_spectrum[32];
	discreteFourierTrivial(test_fn_1, 32, test_fn_1_spectrum, ftdFunctionToSpectrum);

	for (int i = 0; i < 32; i++)
	{
		printf("F%d = (r:%.3f, i:%.3f);\n", i, test_fn_1_spectrum[i].real(), test_fn_1_spectrum[i].imag());
	}

}
void test_discreteFourierFast()
{
	complex<float> test_fn_1[] =
	{
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
	};
	complex<float> test_fn_1_spectrum[32];
	discreteFourierFast(test_fn_1, 32, test_fn_1_spectrum, ftdFunctionToSpectrum);

	for (int i = 0; i < 32; i++)
	{
		printf("F%d = (r:%.3f, i:%.3f);\n", i, test_fn_1_spectrum[i].real(), test_fn_1_spectrum[i].imag());
	}

}

void testForwardAndBackwardCycle(int k)
{
	int N = 1 << k;

	complex<float>* test_fn = new complex<float>[N];
	complex<float>* test_fn_spectrum = new complex<float>[N];

	for (int i = 0; i < N; i++)
	{
		test_fn[i] = (float)(i - N/2) / (N/2);
	}

	discreteFourierFast(test_fn, N, test_fn_spectrum, ftdFunctionToSpectrum);
	discreteFourierFast(test_fn_spectrum, N, test_fn, ftdSpectrumToFunction);

	// Calculating dispersion
	float D = 0;
	for (int i = 0; i < N; i++)
	{
		float ideal = (float)(i - N/2) / (N/2);
		float dd = ideal - real(test_fn[i]);
		dd *= dd;
		D += dd;
	}
	D = sqrt(D / N);

	float prec = -round(log10(D)) - 1;

	printf("N = %d, precision = %.0f\n", N, prec);

	delete [] test_fn;
	delete [] test_fn_spectrum;
}

void testForwardAndBackward()
{
	for (int k = 3; k < 23; k++)
	{
		testForwardAndBackwardCycle(k);
	}
}

void testForwardAndBackward2DCycle(int k, int m)
{
	int i_max = 1 << k;
	int j_max = 1 << m;

	complex<float>* test_fn = new complex<float>[i_max * j_max];
	complex<float>* test_fn_spectrum = new complex<float>[i_max * j_max];

	float ij_max = (i_max + j_max) / 2;
	for (int i = 0; i < i_max; i++)
	for (int j = 0; j < j_max; j++)
	{
		test_fn[j * i_max + i] = ((float)i + j - ij_max) / ij_max;
	}

	clock_t before, after;

	before = clock();

	discreteFourierFast2D(test_fn, i_max, j_max, test_fn_spectrum, ftdFunctionToSpectrum);
	discreteFourierFast2D(test_fn_spectrum, i_max, j_max, test_fn, ftdSpectrumToFunction);

	after = clock();

	// Calculating dispersion
	float D = 0;
	for (int i = 0; i < i_max; i++)
	for (int j = 0; j < j_max; j++)
	{
		float ideal = ((float)i + j - ij_max) / ij_max;
		float dd = ideal - real(test_fn[j * i_max + i]);
		dd *= dd;
		D += dd;
	}
	D = sqrt(D / (i_max * j_max));

	float prec = -round(log10(D)) - 1;

	printf("%6d\t%6d\t%10.0f\t%6ld\n", i_max, j_max, prec, after - before);

	delete [] test_fn;
	delete [] test_fn_spectrum;
}

void testForwardAndBackward2D()
{
	printf(" i_max\t j_max\t precision\t  time\n");
	for (int k = 3; k < 13; k++)
	{
		testForwardAndBackward2DCycle(k, k - 1);
	}
}

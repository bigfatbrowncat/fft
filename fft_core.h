#ifndef FFT_CORE_H_
#define FFT_CORE_H_

#include <complex>

using namespace std;

enum fourier_transform_direction { ftdFunctionToSpectrum, ftdSpectrumToFunction };

#define INCORRECT_SPECTRUM_SIZE_FOR_FFT 			1
#define UNSUPPORTED_FTD								2

void discreteFourierFast(const complex<double>* f, int i_max, complex<double>* F, fourier_transform_direction ftd);
void discreteFourierFast2D(const complex<double>* f, int i_max, int j_max, complex<double>* F, fourier_transform_direction ftd);


#endif /* FFT_CORE_H_ */

# Cross Spectrum Example

A simple example showing cross power spectrum computation between two complex valued input signals `s0` and `s1`. Cross power spectrum is computed using Welch' method with Hamming winnow function.

Implementation of cross power spectrum estimation is done using C++ and SIMD intrinsics, so it is assumed that computer hardware supports AVX intrinsics. For both cases elapsed times are measured and difference is shown.

FFT computations are done using [`fftw3`](http://www.fftw.org//) library, so that must be installed before actually using the code.

**NB!** The code in this repo is just for illustration purposes of how cross power spectrum can be estimated and for showing some optimizations of the critical parts. it is not meant for using in real life production code since some preprocessing parts are left out.

## Building

Builing can be done using `cmake` in the root directory of the project as follows:
```
mkdir build
cd build
cmake .. && make
```

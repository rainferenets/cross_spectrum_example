#pragma once

#include <string>
#include <cmath>
#include <random>
#include <iostream>
#include <fftw3.h>
#include <immintrin.h>
#include "definitions.hpp"

class DSP {
private:
    static constexpr size_t win_len = 512;
    static constexpr size_t stp_len = 256;
    static constexpr size_t steps = (signal_length - win_len)/stp_len + 1;

    static constexpr float f0 = 2000.0f;
    static constexpr float df = -0.2f;
    static constexpr float pi = static_cast<float>(M_PI);
    static constexpr float pi2 = static_cast<float>(M_PI*2.0);

    float* win_func = nullptr;

    fftwf_complex* buf_s0 = nullptr;
    fftwf_complex* buf_s1 = nullptr;
    fftwf_complex* buf_fft_s0 = nullptr;
    fftwf_complex* buf_fft_s1 = nullptr;

    fftwf_complex* cross_spectrum = nullptr;

    fftwf_plan plan_s0 = nullptr;
    fftwf_plan plan_s1 = nullptr;

    float norm_coef;

    void reset_data(fftwf_complex *data);
    void create_hamming_win();
    void accumulate_cross_spectrum(const fftwf_complex *fft1, const fftwf_complex *fft2, fftwf_complex *cross_spect);
    void accumulate_cross_spectrum_avx(const float* fft1, const float* fft2, float* cross_spect);

public:
    DSP();
    ~DSP();

    void get_memory(float** data);
    void free_memory(float** data);
    void generate_data(float** data);
    void get_cross_spectrum_mag(float** data_in, float* data_out);
    void get_cross_spectrum_mag_avx(float** data_in, float* data_out);
};


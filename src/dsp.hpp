#pragma once

#include <string>
#include <cmath>
#include <random>
#include <fftw3.h>
#include "definitions.hpp"

class DSP {
private:
    static constexpr size_t win_len = 256;
    static constexpr size_t fft_len = 1024;

    static constexpr float f0 = 2000.0f;
    static constexpr float df = -0.2f;
    static constexpr float pi = static_cast<float>(M_PI);
    static constexpr float pi2 = static_cast<float>(M_PI*2.0);

    float* win_func = nullptr;
    fftwf_complex* fft_s0 = nullptr;
    fftwf_complex* fft_s1 = nullptr;

    void create_hamming_win();

public:
    DSP();
    ~DSP();

    void get_memory(float** data);
    void free_memory(float** data);
    void generate_data(float** data);
};


#include <iostream>
#include <chrono>
#include "dsp.hpp"

int main(int argc, char* argv[])
{
    constexpr size_t N = 100;

    DSP dsp;

    float *data_in[4];
    float data_out_0[fft_len];
    float data_out_1[fft_len];
    double time_orig, time_avx;
    double diff = 0.0;

    dsp.get_memory(data_in);
    dsp.generate_data(data_in);

    auto t0 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < N; i++) {
        dsp.get_cross_spectrum_mag(data_in, data_out_0);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
    time_orig = static_cast<double>(duration)/1000.0;

    std::printf("\nCPP f32 time   = %6.3f ms\n", time_orig);

    t0 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < N; i++) {
        dsp.get_cross_spectrum_mag_avx(data_in, data_out_1);
    }
    t1 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
    time_avx = static_cast<double>(duration)/1000.0;

    std::printf("AVX f32 time   = %6.3f ms\n", time_avx);

    for (size_t i = 0; i < fft_len; i++) {
        double delta = (double)data_out_0[i] - (double)data_out_1[i];
        diff += delta*delta;
    }
    diff /= fft_len;

    std::printf("Change in time = %4.3f%%\n", (time_orig-time_avx)/time_orig*100.0);
    std::printf("MSE, CPP<->AVX = %e\n\n", diff);

    dsp.free_memory(data_in);

    return 0;
}
#include <iostream>
#include <fstream>
#include <chrono>
#include "dsp.hpp"

int main(int argc, char* argv[])
{
    DSP dsp;

    float *data_in[4];

    dsp.get_memory(data_in);
    dsp.generate_data(data_in);

    for (size_t i = 0; i < signal_length; i++) {
        std::printf("%8.5f %8.5f %8.5f %8.5f\n", data_in[0][i], data_in[1][i], data_in[2][i], data_in[3][i]);
    }

    dsp.free_memory(data_in);

    return 0;
}
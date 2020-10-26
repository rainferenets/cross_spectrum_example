#include "dsp.hpp"

DSP::DSP()
{
}

DSP::~DSP()
{
}


void DSP::get_memory(float** data)
{
    for (size_t i = 0; i < 4; i++) {
       data[i] = static_cast<float*>(std::aligned_alloc(32, signal_length*sizeof(float)));
    }
}


void DSP::free_memory(float** data)
{
    for (size_t i = 0; i < 4; i++) {
       std::free(data[i]);
    }
}


void DSP::generate_data(float** data)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> noise(0.0f, 0.2f);

    for (size_t i = 0; i < signal_length; i++) {
        data[0][i] = std::cos(pi2*f0*i/fs) + noise(gen);
        data[1][i] = std::sin(pi2*f0*i/fs) + noise(gen);
        data[2][i] = std::cos(pi2*f0*i/fs + 0.2f) + noise(gen);
        data[3][i] = std::sin(pi2*f0*i/fs + 0.2f) + noise(gen);
    }
}
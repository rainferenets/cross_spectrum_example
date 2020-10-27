#include "dsp.hpp"

DSP::DSP()
{
    win_func = static_cast<float*>(std::aligned_alloc(32, win_len*sizeof(float)));
    fft_s0 = fftwf_alloc_complex(fft_len);
    fft_s1 = fftwf_alloc_complex(fft_len);

    create_hamming_win();
}

DSP::~DSP()
{
    std::free(win_func);
    if (fft_s0) {
        fftwf_free(fft_s0);
    }
    if (fft_s1) {
        fftwf_free(fft_s1);
    }
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
    std::normal_distribution<> noise(0.0f, 0.1f);

    for (size_t i = 0; i < signal_length; i++) {
        data[0][i] = 0.20f*std::cos(pi2*f0*i/fs) + noise(gen);
        data[1][i] = 0.25f*std::sin(pi2*f0*i/fs) + noise(gen);
        data[2][i] = 0.30f*std::cos(pi2*f0*i/fs + df) + noise(gen);
        data[3][i] = 0.25f*std::sin(pi2*f0*i/fs + df) + noise(gen);
    }
}


void DSP::create_hamming_win()
{
    for (size_t i = 0; i < win_len; i++) {
        win_func[i] = 0.54f - 0.46f*std::cos(2.0f*pi*(float)i/(float)(win_len-1));
    }
}
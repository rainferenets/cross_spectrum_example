#include "dsp.hpp"


DSP::DSP()
{
    std::printf("\n");
    std::printf("signal length  = %ld\n", signal_length);
    std::printf("window length  = %ld\n", win_len);
    std::printf("step length    = %ld\n", stp_len);
    std::printf("Welch's steps  = %ld\n", steps);

    win_func = static_cast<float*>(std::aligned_alloc(32, win_len*sizeof(float)));

    buf_s0 = fftwf_alloc_complex(fft_len);
    buf_s1 = fftwf_alloc_complex(fft_len);
    buf_fft_s0 = fftwf_alloc_complex(fft_len);
    buf_fft_s1 = fftwf_alloc_complex(fft_len);
    cross_spectrum = fftwf_alloc_complex(fft_len);

    create_hamming_win();

    norm_coef = 0.0f;
    for (size_t i = 0; i < win_len; i++) {
        norm_coef += win_func[i]*win_func[i];
    }

    norm_coef *= fft_len*steps;
    norm_coef = 1.0f/norm_coef;

    std::printf("FFTW wisdom    = ");
    int r = fftwf_import_wisdom_from_filename("fft_wisdom.dat");
    if (r != 0) {
        std::printf("yes\n");
    }
    else {
        std::printf("no\n");
    }

    plan_s0 = fftwf_plan_dft_1d(fft_len, buf_s0, buf_fft_s0, FFTW_FORWARD, FFTW_MEASURE);
    plan_s1 = fftwf_plan_dft_1d(fft_len, buf_s1, buf_fft_s1, FFTW_FORWARD, FFTW_MEASURE);

    fftwf_export_wisdom_to_filename("fft_wisdom.dat");

    reset_data(buf_s0);
    reset_data(buf_s1);
}


DSP::~DSP()
{
    std::free(win_func);
    if (buf_s0) fftwf_free(buf_s0);
    if (buf_s1) fftwf_free(buf_s1);
    if (buf_fft_s0) fftwf_free(buf_fft_s0);
    if (buf_fft_s1) fftwf_free(buf_fft_s1);
    if (cross_spectrum) fftwf_free(cross_spectrum);
    if (plan_s0) fftwf_destroy_plan(plan_s0);
    if (plan_s1) fftwf_destroy_plan(plan_s1);
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


void DSP::reset_data(fftwf_complex *data)
{
    for (size_t i = 0; i < fft_len; i++) {
        data[i][0] = data[i][1] = 0.0f;
    }
}


void DSP::create_hamming_win()
{
    for (size_t i = 0; i < win_len; i++) {
        win_func[i] = 0.54f - 0.46f*std::cos(2.0f*pi*(float)i/(float)(win_len-1));
    }
}


void DSP::accumulate_cross_spectrum(const fftwf_complex *fft1, const fftwf_complex *fft2, fftwf_complex *cross_spect)
{
#pragma GCC ivdep 
    for (size_t i = 0; i < fft_len; i++) { 
        cross_spect[i][0] += (fft2[i][0]*fft1[i][0] + fft2[i][1]*fft1[i][1]);
        cross_spect[i][1] += (-fft2[i][0]*fft1[i][1] + fft2[i][1]*fft1[i][0]);
    }
}


void DSP::accumulate_cross_spectrum_avx(const float* fft1, const float* fft2, float* cross_spect)
{
    const __m256 odd_signbits  = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f);

    for (size_t i = 0; i < fft_len*2; i += 8) {
        __m256 res = _mm256_load_ps(cross_spect+i);
        __m256 vec1 = _mm256_load_ps(fft1+i);
        __m256 vec2 = _mm256_load_ps(fft2+i);

        __m256 vec3 = _mm256_mul_ps(vec1, vec2);

        vec2 = _mm256_permute_ps(vec2, 0xB1);
        vec1 = _mm256_xor_ps(vec1, odd_signbits);
        __m256 vec4 = _mm256_mul_ps(vec1, vec2);

        vec1 = _mm256_hadd_ps(vec3, vec4);
        vec1 = _mm256_permute_ps(vec1, 0xD8);
        res = _mm256_add_ps(res, vec1);

        _mm256_store_ps(cross_spect+i, res); 
    }
}


void DSP::get_cross_spectrum_mag(float** data_in, float* data_out)
{
    reset_data(cross_spectrum);

    for (size_t i = 0; i < steps; i++) {
        size_t i0 = i*stp_len;

        for (size_t j = 0; j < win_len; j++) {
            buf_s0[j][0] = data_in[0][i0+j]*win_func[j];
            buf_s0[j][1] = data_in[1][i0+j]*win_func[j];
            buf_s1[j][0] = data_in[2][i0+j]*win_func[j];
            buf_s1[j][1] = data_in[3][i0+j]*win_func[j];
        }

        fftwf_execute(plan_s0);
        fftwf_execute(plan_s1);

        accumulate_cross_spectrum(buf_fft_s0, buf_fft_s1, cross_spectrum);
    }

    for (size_t i = 0; i < fft_len; i++) {
        data_out[i] = 10.0f*std::log10(std::sqrt(cross_spectrum[i][0]*cross_spectrum[i][0] + cross_spectrum[i][1]*cross_spectrum[i][1])*norm_coef);
    }
}


void DSP::get_cross_spectrum_mag_avx(float** data_in, float* data_out)
{
    reset_data(cross_spectrum);

    for (size_t i = 0; i < steps; i++) {
        size_t i0 = i*stp_len;

        for (size_t j = 0; j < win_len; j++) {
            buf_s0[j][0] = data_in[0][i0+j]*win_func[j];
            buf_s0[j][1] = data_in[1][i0+j]*win_func[j];
            buf_s1[j][0] = data_in[2][i0+j]*win_func[j];
            buf_s1[j][1] = data_in[3][i0+j]*win_func[j];
        }

        fftwf_execute(plan_s0);
        fftwf_execute(plan_s1);

        accumulate_cross_spectrum_avx(reinterpret_cast<float*>(buf_fft_s0), reinterpret_cast<float*>(buf_fft_s1), reinterpret_cast<float*>(cross_spectrum));
    }

    for (size_t i = 0; i < fft_len; i++) {
        data_out[i] = 10.0f*std::log10(std::sqrt(cross_spectrum[i][0]*cross_spectrum[i][0] + cross_spectrum[i][1]*cross_spectrum[i][1])*norm_coef);
    }
}
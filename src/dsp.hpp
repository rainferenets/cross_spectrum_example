#pragma once

#include <string>
#include <cmath>
#include <random>
#include "definitions.hpp"

class DSP {
private:
    static constexpr float f0 = 2000.0f;
    static constexpr float f1 = 3000.0f;
    static constexpr float pi = static_cast<float>(M_PI);
    static constexpr float pi2 = static_cast<float>(M_PI*2.0);

    
public:
    DSP();
    ~DSP();

    void get_memory(float** data);
    void free_memory(float** data);
    void generate_data(float** data);
};


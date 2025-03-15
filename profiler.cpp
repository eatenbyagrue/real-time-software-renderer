#include "profiler.hpp"
#include <chrono>
#include <iostream>
#include <ostream>

void Profiler::register_profile(Profile profile, std::string label) {
    this->labels[profile] = label;
    this->timings[profile] = 0;
}

void Profiler::start(Profile profile) {
    this->current_profile = profile;
    this->t0 = std::chrono::high_resolution_clock::now();
}
void Profiler::stop() {
    this->t1 = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
    this->timings[current_profile] += elapsed_time;
}
void Profiler::draw_profiler_card() {
    for (const auto &[key, value] : this->labels) {
        std::cout << this->labels[key] << (double)this->timings[key] / 1000.0 << "\n";
    }
}

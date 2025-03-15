#ifndef PROFILER_H
#define PROFILER_H

#include <chrono>
#include <map>
#include <string>

enum class Profile { TRANSFORMING, CLIPPING, DRAWING };

/*Basic-ass profiler. Register profiles with a label to be printed, then take the time with start
 * and stop. Only one profile can be timed at any one moment. A profile is 'what computing time is
 * spent on'.*/

class Profiler {
  private:
    Profile current_profile;
    std::map<Profile, std::string> labels;
    std::map<Profile, long int> timings;
    std::chrono::time_point<std::chrono::high_resolution_clock> t0;
    std::chrono::time_point<std::chrono::high_resolution_clock> t1;

  public:
    void register_profile(Profile profile, std::string label);
    void start(Profile profile);
    void stop();
    void draw_profiler_card();
};

#endif

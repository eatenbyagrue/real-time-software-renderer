#include <algorithm>
#include <array>
#include <iostream>
int main() {
    std::array<double, 3> d = {1.23, -23, 91};
    std::cout << d.at(0) << " " << d.at(1) << " " << d.at(2) << std::endl;
    std::sort(d.begin(), d.end());
    std::cout << d.at(0) << " " << d.at(1) << " " << d.at(2) << std::endl;

    int count_positive = std::count_if(d.begin(), d.end(), [](int d) { return d >= 0; });
    std::cout << count_positive << std::endl;
    return 0;
}

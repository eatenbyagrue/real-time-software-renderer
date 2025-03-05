
#include <iostream>
struct Test {
    int a;
    double b;
};

int main() {
    Test test{1, 2.2};
    Test test2 = test;
    Test &test3 = test;
    std::cout << "original: " << test.a << " " << test.b << " copy: " << test2.a << " " << test2.b
              << " reference" << test3.a << " " << test3.b << std::endl;
    test.a = 5;
    std::cout << "original: " << test.a << " " << test.b << " copy: " << test2.a << " " << test2.b
              << " reference" << test3.a << " " << test3.b << std::endl;
    return 0;
}

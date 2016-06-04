//#include "checkpoint.h"
#include "simpson.h"
#include <array>            // for std::array
#include <cmath>            // for std::sqrt
#include <iomanip>
#include <iostream>
#include <string>

namespace {
    static auto constexpr DIGIT = 15U;
    static auto constexpr LOOPMAX = 1000000000UL;
    static auto constexpr N = 10001U;
}

int main()
{
    auto const func = myfunctional::make_functional([](double x) { return 1.0 / (2.0 * std::sqrt(x)); });
	
    simpson::Simpson s(100);
    s.simpson<decltype(func), 0>(func, 1.0, 4.0);

    //std::cout << "正確な値：\t" << std::setprecision(DIGIT) << exact  << '\n';
	//std::cout << "AVX無効：\t" << std::setprecision(DIGIT) << res[0] << '\n';
    //std::cout << "AVX有効：\t" << std::setprecision(DIGIT) << res[1] << '\n';

	return 0;
}

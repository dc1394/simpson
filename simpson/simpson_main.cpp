#include "checkpoint.h"
#include "simpson_integral.h"
#include <array>        // for std::array
#include <cmath>        // for std::sqrt
#include <iomanip>      // for std::setprecision
#include <iostream>     // for std::ios_base::fixed, std::ios_base::floatfield, std::cout

namespace {
    static auto constexpr DIGIT = 15U;
    static auto constexpr N = 100000000;
}

int main()
{
    checkpoint::CheckPoint cp;

    auto const func = myfunctional::make_functional([](double x) { return 1.0 / (2.0 * std::sqrt(x)); });
	simpson::Simpson<decltype(func)> s(func, N, 1.0, 4.0);
    std::array<double, 7> res{};

    cp.checkpoint("処理開始", __LINE__);

    res[0] = s.simpson<simpson::ParallelType::NoParallel>();

    cp.checkpoint("並列化無し", __LINE__);

    res[1] = s.simpson<simpson::ParallelType::StdAsync>();

    cp.checkpoint("std::asyncで並列化", __LINE__);
    
    res[2] = s.simpson<simpson::ParallelType::Tbb>();

    cp.checkpoint("TBBで並列化", __LINE__);

    res[3] = s.simpson<simpson::ParallelType::Tbb2>();

    cp.checkpoint("TBBで並列化2", __LINE__);

    res[4] = s.simpson<simpson::ParallelType::Ppl>();

    cp.checkpoint("PPLで並列化", __LINE__);
    
    res[5] = s.simpson<simpson::ParallelType::Cilk>();

    cp.checkpoint("Cilkで並列化", __LINE__);

    res[6] = s.simpson<simpson::ParallelType::OpenMp>();

    cp.checkpoint("OpenMPで並列化", __LINE__);
    
    std::cout << "積分点: " << N << '\n';

    cp.checkpoint_print();
    checkpoint::usedmem();

    std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
    std::cout << "並列化無し：\t" << std::setprecision(DIGIT) << res[0]  << '\n';
    std::cout << "std::asyncで並列化：" << std::setprecision(DIGIT) << res[1] << '\n';
    std::cout << "TBBで並列化：\t" << std::setprecision(DIGIT) << res[2] << '\n';
    std::cout << "TBBで並列化2：\t" << std::setprecision(DIGIT) << res[3] << '\n';
    std::cout << "PPLで並列化：\t" << std::setprecision(DIGIT) << res[4] << '\n';
    std::cout << "Cilkで並列化：\t" << std::setprecision(DIGIT) << res[5] << '\n';
    std::cout << "OpenMPで並列化：" << std::setprecision(DIGIT) << res[6] << '\n';

	return 0;
}

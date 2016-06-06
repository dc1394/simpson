#include "checkpoint.h"
#include "simpson.h"
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
	simpson::Simpson s(N);
    std::array<double, 6> res{};

    cp.checkpoint("処理開始", __LINE__);

    res[0] = s.simpson<decltype(func), simpson::ParallelType::NoParallel>(func, 1.0, 4.0);

    cp.checkpoint("並列化無し", __LINE__);

    res[1] = s.simpson<decltype(func), simpson::ParallelType::Tbb>(func, 1.0, 4.0);

    cp.checkpoint("TBBで並列化", __LINE__);

    res[2] = s.simpson<decltype(func), simpson::ParallelType::Tbb2>(func, 1.0, 4.0);

    cp.checkpoint("TBBで並列化2", __LINE__);

    res[3] = s.simpson<decltype(func), simpson::ParallelType::Ppl>(func, 1.0, 4.0);

    cp.checkpoint("PPLで並列化", __LINE__);
    
    res[4] = s.simpson<decltype(func), simpson::ParallelType::Cilk>(func, 1.0, 4.0);

    cp.checkpoint("Cilkで並列化", __LINE__);

    res[5] = s.simpson<decltype(func), simpson::ParallelType::OpenMp>(func, 1.0, 4.0);

    cp.checkpoint("OpenMPで並列化", __LINE__);
        
    cp.checkpoint_print();
    checkpoint::usedmem();

    std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
    std::cout << "並列化無し：\t" << std::setprecision(DIGIT) << res[0]  << '\n';
	std::cout << "TBBで並列化：\t" << std::setprecision(DIGIT) << res[1] << '\n';
    std::cout << "TBBで並列化2：\t" << std::setprecision(DIGIT) << res[2] << '\n';
    std::cout << "PPLで並列化：\t" << std::setprecision(DIGIT) << res[3] << '\n';
    std::cout << "Cilkで並列化：\t" << std::setprecision(DIGIT) << res[4] << '\n';
    std::cout << "OpenMPで並列化：" << std::setprecision(DIGIT) << res[5] << '\n';

	return 0;
}

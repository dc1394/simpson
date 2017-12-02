/*! \file paralleltype.h
    \brief 並列化の手法を表す列挙型の宣言

    Copyright ©  2016 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#ifndef _PARALLELTYPE_H_
#define _PARALLELTYPE_H_

#pragma once

#include <cstdint>  // for std::int32_t

namespace simpson {
    enum class ParallelType : std::int32_t {
        Cilk = 1,
		CPP17 = 2,
        NoParallel = 3,
        OpenMp = 4,
        Ppl = 5,
        StdAsync = 6,
        Tbb = 7,
        Tbb2 = 8
    };
}

#endif  // _PARALLELTYPE_H_

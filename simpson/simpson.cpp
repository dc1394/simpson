/*! \file simpson.cpp
    \brief simpsonの公式で数値積分を行うクラスの実装

    Copyright ©  2014 @dc1394 All Rights Reserved.
*/
#include "simpson.h"

namespace simpson {
    Simpson::Simpson(std::int32_t n)
        : avxSupported(availableAVX()), n_(n)
    {
    }
    
    bool Simpson::availableAVX() const
    {
#if (_MSC_FULL_VER >= 160040219)
        std::array<std::int32_t, 4> cpuInfo = { 0 };
        ::__cpuid(cpuInfo.data(), 1);

        auto const osUsesXSAVE_XRSTORE = cpuInfo[2] & (1 << 27) || false;
        auto const cpuAVXSuport = cpuInfo[2] & (1 << 28) || false;

        if (osUsesXSAVE_XRSTORE && cpuAVXSuport)
        {
            // Check if the OS will save the YMM registers
            auto const xcrFeatureMask = _xgetbv(_XCR_XFEATURE_ENABLED_MASK);
            return (xcrFeatureMask & 0x6) || false;
        }
#endif
        return false;
    }
}

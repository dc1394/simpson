/*! \file simpson.h
    \brief simpsonの公式で数値積分を行うクラスの宣言

    Copyright ©  2016 @dc1394 All Rights Reserved.
*/
#ifndef _SIMPSON_H_
#define _SIMPSON_H_

#pragma once

#include "functional.h"
#include <array>                            // for std::array
#include <cstdint>                          // for std::int32_t
#include <vector>                           // for std::vector
#include <intrin.h>                         // for _xgetbv
#include <boost/align/aligned_allocator.hpp>// for boost::alignment::aligned_allocator
#include <boost/mpl/int.hpp>
#include <dvec.h>                           // for F64vec4, F64vec2

namespace simpson {
    //! A class.
    /*!
        simpsonの公式で数値積分を行うクラス
    */
	class Simpson final
	{
    public:
        // #region コンストラクタ

        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param n simpsonの公式の分点
        */
        Simpson(std::int32_t n);

        // #endregion コンストラクタ

        // #region メンバ関数

        template <typename FUNCTYPE, int V>
        double simpson(myfunctional::Functional<FUNCTYPE> const & func, double x1, double x2) {
            return operator()(func, x1, x2, boost::mpl::int_<V>());
        }

    private:
        //! A private member function (template function).
        /*!
            Simpsonの公式によって数値積分を実行する
            \param func 被積分関数
            \param x1 積分の下端
            \param x2 積分の上端
            \return 積分値
        */
        template <typename FUNCTYPE>
        double operator()(myfunctional::Functional<FUNCTYPE> const & func, double x1, double x2, boost::mpl::int_<0>) const;

        //! A private member function.
        /*!
            AVX命令が使用可能かどうかをチェックする
            \return AVX命令が使用可能ならtrue、使用不可能ならfalse
        */
        bool availableAVX() const;

        // #endregion メンバ関数

        // #region メンバ変数

        //! A private member variable (constant).
        /*!
            AVX命令が使用可能かどうか
        */
        bool const avxSupported;

        //! A private member variable (constant).
        /*!
            simpsonの公式の分点
        */
        std::uint32_t const n_;

        // #endregion メンバ変数

        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        Simpson() = delete;

        //! A private copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
        */
		Simpson(Simpson const &) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト
            \return コピー元のオブジェクト
        */
		Simpson & operator=(Simpson const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
	};
    
    template <typename FUNCTYPE>
    inline double Simpson::operator()(myfunctional::Functional<FUNCTYPE> const & func, double x1, double x2, boost::mpl::int_<0>) const
    {
        return 0.0;
    }
}

#endif  // _SIMPSON_H_


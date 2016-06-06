/*! \file simpson.h
    \brief simpsonの公式で数値積分を行うクラスの宣言と実装

    Copyright ©  2016 @dc1394 All Rights Reserved.
*/
#ifndef _SIMPSON_H_
#define _SIMPSON_H_

#pragma once

#include "functional.h"
#include "paralleltype.h"
#include <cstdint>                  // for std::int32_t
#include <functional>               // for std::plus            
#include <tuple>                    // for std::tuple
#include <omp.h>                    // for pragma omp parallel for
#include <ppl.h>                    // for concurrency::parallel_for
#include <boost/mpl/int.hpp>        // for boost::mpl::int_
#include <cilk/cilk.h>              // for cilk_for
#include <cilk/reducer_opadd.h>     // for cilk::reducer_opadd
#include <tbb/blocked_range.h>      // for tbb:blocked_range
#include <tbb/combinable.h>         // for tbb::combinable
#include <tbb/parallel_reduce.h>    // for tbb:parallel_reduce
#include <tbb/parallel_for.h>       // for tbb::parallel_for

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
        Simpson::Simpson(std::int32_t n) : n_(n) {}

        // #endregion コンストラクタ

        // #region publicメンバ関数

        //! A public member function (template function).
        /*!
            Simpsonの公式によって数値積分を実行する
            \param func 被積分関数
            \param x1 積分の下端
            \param x2 積分の上端
            \return 積分値
        */
        template <typename FUNCTYPE, ParallelType N>
        double simpson(myfunctional::Functional<FUNCTYPE> const & func, double x1, double x2) const
        {
            return (*this)(func, x1, x2, boost::mpl::int_<static_cast<std::int32_t>(N)>());
        }

        // #endregion publicメンバ関数

    private:

        // #region privateメンバ関数

        //! A private member function (template function).
        /*!
            Simpsonの公式によって数値積分を実行する（Cilkで並列化）
            \param func 被積分関数
            \param x1 積分の下端
            \param x2 積分の上端
            \param テンプレート部分特殊化用の引数
            \return 積分値
        */
        template <typename FUNCTYPE>
        double operator()(myfunctional::Functional<FUNCTYPE> const & func, double x1, double x2, boost::mpl::int_<static_cast<std::int32_t>(ParallelType::Cilk)>) const;

        //! A private member function (template function).
        /*!
            Simpsonの公式によって数値積分を実行する（並列化なし）
            \param func 被積分関数
            \param x1 積分の下端
            \param x2 積分の上端
            \param テンプレート部分特殊化用の引数
            \return 積分値
        */
        template <typename FUNCTYPE>
        double operator()(myfunctional::Functional<FUNCTYPE> const & func, double x1, double x2, boost::mpl::int_<static_cast<std::int32_t>(ParallelType::NoParallel)>) const;

        //! A private member function (template function).
        /*!
            Simpsonの公式によって数値積分を実行する（OpenMPで並列化）
            \param func 被積分関数
            \param x1 積分の下端
            \param x2 積分の上端
            \param テンプレート部分特殊化用の引数
            \return 積分値
        */
        template <typename FUNCTYPE>
        double operator()(myfunctional::Functional<FUNCTYPE> const & func, double x1, double x2, boost::mpl::int_<static_cast<std::int32_t>(ParallelType::OpenMp)>) const;

        //! A private member function (template function).
        /*!
            Simpsonの公式によって数値積分を実行する（PPLで並列化）
            \param func 被積分関数
            \param x1 積分の下端
            \param x2 積分の上端
            \param テンプレート部分特殊化用の引数
            \return 積分値
        */
        template <typename FUNCTYPE>
        double operator()(myfunctional::Functional<FUNCTYPE> const & func, double x1, double x2, boost::mpl::int_<static_cast<std::int32_t>(ParallelType::Ppl)>) const;

        //! A private member function (template function).
        /*!
            Simpsonの公式によって数値積分を実行する（TBBで並列化）
            \param func 被積分関数
            \param x1 積分の下端
            \param x2 積分の上端
            \param テンプレート部分特殊化用の引数
            \return 積分値
        */
        template <typename FUNCTYPE>
        double operator()(myfunctional::Functional<FUNCTYPE> const & func, double x1, double x2, boost::mpl::int_<static_cast<std::int32_t>(ParallelType::Tbb)>) const;
        
        //! A private member function (template function).
        /*!
            Simpsonの公式によって数値積分を実行する（TBBで並列化その2）
            \param func 被積分関数
            \param x1 積分の下端
            \param x2 積分の上端
            \param テンプレート部分特殊化用の引数
            \return 積分値
        */
        template <typename FUNCTYPE>
        double operator()(myfunctional::Functional<FUNCTYPE> const & func, double x1, double x2, boost::mpl::int_<static_cast<std::int32_t>(ParallelType::Tbb2)>) const;
        
        // #endregion privateメンバ関数

        // #region メンバ変数

        //! A private member variable (constant).
        /*!
            simpsonの公式の分点
        */
        std::int32_t const n_;

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

    // #region templateメンバ関数の実装

    template <typename FUNCTYPE>
    double Simpson::operator()(myfunctional::Functional<FUNCTYPE> const & func, double x1, double x2, boost::mpl::int_<static_cast<std::int32_t>(ParallelType::Cilk)>) const
    {
        cilk::reducer_opadd<double> sum;
        sum.set_value(0.0);

        auto const dh = (x2 - x1) / static_cast<double>(n_);
        cilk_for (auto i = 0; i < n_; i += 2) {
            auto const f0 = func(x1 + static_cast<double>(i) * dh);
            auto const f1 = func(x1 + static_cast<double>(i + 1) * dh);
            auto const f2 = func(x1 + static_cast<double>(i + 2) * dh);
            sum += (f0 + 4.0 * f1 + f2);
        }

        return sum.get_value() * dh / 3.0;
    }
    
    template <typename FUNCTYPE>
    double Simpson::operator()(myfunctional::Functional<FUNCTYPE> const & func, double x1, double x2, boost::mpl::int_<static_cast<std::int32_t>(ParallelType::NoParallel)>) const
    {
        auto sum = 0.0;
        auto const dh = (x2 - x1) / static_cast<double>(n_);

        for (auto i = 0; i < n_; i += 2) {
            auto const f0 = func(x1 + static_cast<double>(i) * dh);
            auto const f1 = func(x1 + static_cast<double>(i + 1) * dh);
            auto const f2 = func(x1 + static_cast<double>(i + 2) * dh);
            sum += (f0 + 4.0 * f1 + f2);
        }

        return sum * dh / 3.0;
    }

    template <typename FUNCTYPE>
    double Simpson::operator()(myfunctional::Functional<FUNCTYPE> const & func, double x1, double x2, boost::mpl::int_<static_cast<std::int32_t>(ParallelType::OpenMp)>) const
    {
        auto sum = 0.0;
        auto const dh = (x2 - x1) / static_cast<double>(n_);

        #pragma omp parallel for reduction(+:sum)
        for (auto i = 0; i < n_; i += 2) {
            auto const f0 = func(x1 + static_cast<double>(i) * dh);
            auto const f1 = func(x1 + static_cast<double>(i + 1) * dh);
            auto const f2 = func(x1 + static_cast<double>(i + 2) * dh);
            sum += (f0 + 4.0 * f1 + f2);
        }

        return sum * dh / 3.0;
    }

    template <typename FUNCTYPE>
    double Simpson::operator()(myfunctional::Functional<FUNCTYPE> const & func, double x1, double x2, boost::mpl::int_<static_cast<std::int32_t>(ParallelType::Ppl)>) const
    {
        auto const dh = (x2 - x1) / static_cast<double>(n_);
        concurrency::combinable<double> sum;

        concurrency::parallel_for<std::int32_t>(
            0,
            n_ / 2,
            [&](auto i) {
                auto const f0 = func(x1 + static_cast<double>(i * 2) * dh);
                auto const f1 = func(x1 + static_cast<double>(i * 2 + 1) * dh);
                auto const f2 = func(x1 + static_cast<double>(i * 2 + 2) * dh);
                sum.local() += (f0 + 4.0 * f1 + f2);
        });

        return sum.combine(std::plus<double>()) * dh / 3.0;
    }

    template <typename FUNCTYPE>
    double Simpson::operator()(myfunctional::Functional<FUNCTYPE> const & func, double x1, double x2, boost::mpl::int_<static_cast<std::int32_t>(ParallelType::Tbb)>) const
    {
        auto const dh = (x2 - x1) / static_cast<double>(n_);

        auto const sum = tbb::parallel_reduce(
            tbb::blocked_range<std::int32_t>(0, n_ / 2),
            0.0,
            [&](auto const & range, auto sumlocal) {
                for (auto && i = range.begin(); i != range.end(); ++i) {
                    auto const f0 = func(x1 + static_cast<double>(i * 2) * dh);
                    auto const f1 = func(x1 + static_cast<double>(i * 2 + 1) * dh);
                    auto const f2 = func(x1 + static_cast<double>(i * 2 + 2) * dh);
                    sumlocal += (f0 + 4.0 * f1 + f2);
                }
                return sumlocal;
            },
            std::plus<double>());

        return sum * dh / 3.0;
    }

    template <typename FUNCTYPE>
    double Simpson::operator()(myfunctional::Functional<FUNCTYPE> const & func, double x1, double x2, boost::mpl::int_<static_cast<std::int32_t>(ParallelType::Tbb2)>) const
    {
        auto const dh = (x2 - x1) / static_cast<double>(n_);
        tbb::combinable<double> sum;

        tbb::parallel_for(
            tbb::blocked_range<std::int32_t>(0, n_ / 2),
            [&](auto const & range) {
                auto sumlocal = 0.0;
                for (auto && i = range.begin(); i != range.end(); ++i) {
                    auto const f0 = func(x1 + static_cast<double>(i * 2) * dh);
                    auto const f1 = func(x1 + static_cast<double>(i * 2 + 1) * dh);
                    auto const f2 = func(x1 + static_cast<double>(i * 2 + 2) * dh);
                    sumlocal += (f0 + 4.0 * f1 + f2);
                }
                sum.local() += sumlocal;
        });

        return sum.combine(std::plus<double>()) * dh / 3.0;
    }
    
    // #endregion templateメンバ関数の実装
}

#endif  // _SIMPSON_H_


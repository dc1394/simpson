/*! \file simpson_integral.h
    \brief simpsonの公式で数値積分を行うクラスの宣言と実装

    Copyright ©  2016 @dc1394 All Rights Reserved.
*/
#ifndef _SIMPSON_INTEGRAL_H_
#define _SIMPSON_INTEGRAL_H_

#pragma once

#include "functional.h"
#include "paralleltype.h"
#include <algorithm>                // for std::max
#include <cstdint>                  // for std::int32_t
#include <functional>               // for std::plus
#include <future>                   // for std::async, std::future
#include <thread>                   // for std::thread::hardware_concurrency
#include <vector>                   // for std::vector
#include <omp.h>                    // for pragma omp parallel for
#include <ppl.h>                    // for concurrency::parallel_for
#include <boost/mpl/int.hpp>        // for boost::mpl::int_
#include <boost/range/numeric.hpp>  // for boost::accumulate
#include <cilk/cilk.h>              // for cilk_for
#include <cilk/reducer_opadd.h>     // for cilk::reducer_opadd
#include <tbb/blocked_range.h>      // for tbb:blocked_range
#include <tbb/combinable.h>         // for tbb::combinable
#include <tbb/parallel_for.h>       // for tbb::parallel_for
#include <tbb/parallel_reduce.h>    // for tbb:parallel_reduce

namespace simpson {
    //! A template class.
    /*!
        simpsonの公式で数値積分を行うクラス
    */
    template <typename FUNCTYPE>
	class Simpson final {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param func_ 被積分関数
            \param n simpsonの公式の分割数
            \param x1 積分の下限
            \param x2 積分の上限
        */
        Simpson(myfunctional::Functional<FUNCTYPE> const & func_, std::int32_t n, double x1, double x2)
            : func_(func_), n_(n), x1_(x1), dh_((x2 - x1) / static_cast<double>(n_)) {}

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~Simpson() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region publicメンバ関数

        //! A public member function (template function).
        /*!
            Simpsonの公式によって数値積分を実行する
            \return 積分値
        */
        template <ParallelType N>
        double simpson() const
        {
            return (*this)(boost::mpl::int_<static_cast<std::int32_t>(N)>());
        }

        // #endregion publicメンバ関数

        // #region privateメンバ関数

    private:
        //! A private member function (template function).
        /*!
            Simpsonの公式によって数値積分を実行する（Cilkで並列化）
            \param テンプレート部分特殊化用の引数
            \return 積分値
        */
        double operator()(boost::mpl::int_<static_cast<std::int32_t>(ParallelType::Cilk)>) const;

        //! A private member function (template function).
        /*!
            Simpsonの公式によって数値積分を実行する（並列化なし）
            \param テンプレート部分特殊化用の引数
            \return 積分値
        */
        double operator()(boost::mpl::int_<static_cast<std::int32_t>(ParallelType::NoParallel)>) const;

        //! A private member function (template function).
        /*!
            Simpsonの公式によって数値積分を実行する（OpenMPで並列化）
            \param テンプレート部分特殊化用の引数
            \return 積分値
        */
        double operator()(boost::mpl::int_<static_cast<std::int32_t>(ParallelType::OpenMp)>) const;

        //! A private member function (template function).
        /*!
            Simpsonの公式によって数値積分を実行する（PPLで並列化）
            \param テンプレート部分特殊化用の引数
            \return 積分値
        */
        double operator()(boost::mpl::int_<static_cast<std::int32_t>(ParallelType::Ppl)>) const;

        //! A private member function (template function).
        /*!
            Simpsonの公式によって数値積分を実行する（std::asyncで並列化）
            \param テンプレート部分特殊化用の引数
            \return 積分値
        */
        double operator()(boost::mpl::int_<static_cast<std::int32_t>(ParallelType::StdAsync)>) const;

        //! A private member function (template function).
        /*!
            Simpsonの公式によって数値積分を実行する（TBBで並列化）
            \param テンプレート部分特殊化用の引数
            \return 積分値
        */
        double operator()(boost::mpl::int_<static_cast<std::int32_t>(ParallelType::Tbb)>) const;
        
        //! A private member function (template function).
        /*!
            Simpsonの公式によって数値積分を実行する（TBBで並列化その2）
            \param テンプレート部分特殊化用の引数
            \return 積分値
        */
        double operator()(boost::mpl::int_<static_cast<std::int32_t>(ParallelType::Tbb2)>) const;
        
        // #endregion privateメンバ関数

        // #region メンバ変数

        //! A private member variable (constant).
        /*!
            被積分関数
        */
        myfunctional::Functional<FUNCTYPE> const func_;

        //! A private member variable (constant).
        /*!
            simpsonの公式の積分点
        */
        std::int32_t const n_;

        //! A private member variable (constant).
        /*!
            積分の下限
        */
        double const x1_;

        //! A private member variable (constant).
        /*!
            積分の微小区間
        */
        double const dh_;

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
    double Simpson<FUNCTYPE>::operator()(boost::mpl::int_<static_cast<std::int32_t>(ParallelType::Cilk)>) const
    {
        cilk::reducer_opadd<double> sum;
        sum.set_value(0.0);

        cilk_for (auto i = 0; i < n_; i += 2) {
            auto const f0 = func_(x1_ + static_cast<double>(i) * dh_);
            auto const f1 = func_(x1_ + static_cast<double>(i + 1) * dh_);
            auto const f2 = func_(x1_ + static_cast<double>(i + 2) * dh_);
            sum += (f0 + 4.0 * f1 + f2);
        }

        return sum.get_value() * dh_ / 3.0;
    }
    
    template <typename FUNCTYPE>
    double Simpson<FUNCTYPE>::operator()(boost::mpl::int_<static_cast<std::int32_t>(ParallelType::NoParallel)>) const
    {
        auto sum = 0.0;

        for (auto i = 0; i < n_; i += 2) {
            auto const f0 = func_(x1_ + static_cast<double>(i) * dh_);
            auto const f1 = func_(x1_ + static_cast<double>(i + 1) * dh_);
            auto const f2 = func_(x1_ + static_cast<double>(i + 2) * dh_);
            sum += (f0 + 4.0 * f1 + f2);
        }

        return sum * dh_ / 3.0;
    }

    template <typename FUNCTYPE>
    double Simpson<FUNCTYPE>::operator()(boost::mpl::int_<static_cast<std::int32_t>(ParallelType::OpenMp)>) const
    {
        auto sum = 0.0;

        #pragma omp parallel for reduction(+:sum)
        for (auto i = 0; i < n_; i += 2) {
            auto const f0 = func_(x1_ + static_cast<double>(i) * dh_);
            auto const f1 = func_(x1_ + static_cast<double>(i + 1) * dh_);
            auto const f2 = func_(x1_ + static_cast<double>(i + 2) * dh_);
            sum += (f0 + 4.0 * f1 + f2);
        }

        return sum * dh_ / 3.0;
    }

    template <typename FUNCTYPE>
    double Simpson<FUNCTYPE>::operator()(boost::mpl::int_<static_cast<std::int32_t>(ParallelType::Ppl)>) const
    {
        concurrency::combinable<double> sum;

        concurrency::parallel_for<std::int32_t>(
            0,
            n_ / 2,
            [&](auto i) {
                auto const f0 = func_(x1_ + static_cast<double>(i * 2) * dh_);
                auto const f1 = func_(x1_ + static_cast<double>(i * 2 + 1) * dh_);
                auto const f2 = func_(x1_ + static_cast<double>(i * 2 + 2) * dh_);
                sum.local() += (f0 + 4.0 * f1 + f2);
        });

        return sum.combine(std::plus<double>()) * dh_ / 3.0;
    }

    template <typename FUNCTYPE>
    double Simpson<FUNCTYPE>::operator()(boost::mpl::int_<static_cast<std::int32_t>(ParallelType::StdAsync)>) const
    {
        auto const numthreads = static_cast<std::int32_t>((std::max)(std::thread::hardware_concurrency(), 1u));

        std::vector<std::future<double>> future(numthreads);

        auto const nmax = n_ / 2;
        for (auto i = 0; i < numthreads; i++) {
            std::int32_t localnmax;
            if (i == numthreads - 1) {
                localnmax = nmax;
            }
            else {
                localnmax = nmax / numthreads * (i + 1);
            }

            auto const localnmin = nmax / numthreads * i;
            
            future[i] = std::async(
                std::launch::async,
                [this](std::int32_t nmin, std::int32_t nmax) {
                    auto sum = 0.0;

                    for (auto i = nmin; i < nmax; i++) {
                        auto const f0 = func_(x1_ + static_cast<double>(2 * i) * dh_);
                        auto const f1 = func_(x1_ + static_cast<double>(2 * i + 1) * dh_);
                        auto const f2 = func_(x1_ + static_cast<double>(2 * i + 2) * dh_);
                        sum += (f0 + 4.0 * f1 + f2);
                    }

                    return sum;
                },
                localnmin,
                localnmax);
        }

        std::vector<double> result;
        result.reserve(numthreads);

        for (auto && f : future) {
            result.push_back(f.get());
        }

        return boost::accumulate(result, 0.0) * dh_ / 3.0;
    }

    template <typename FUNCTYPE>
    double Simpson<FUNCTYPE>::operator()(boost::mpl::int_<static_cast<std::int32_t>(ParallelType::Tbb)>) const
    {
        auto const sum = tbb::parallel_reduce(
            tbb::blocked_range<std::int32_t>(0, n_ / 2),
            0.0,
            [&](auto const & range, auto sumlocal) {
                for (auto && i = range.begin(); i != range.end(); ++i) {
                    auto const f0 = func_(x1_ + static_cast<double>(i * 2) * dh_);
                    auto const f1 = func_(x1_ + static_cast<double>(i * 2 + 1) * dh_);
                    auto const f2 = func_(x1_ + static_cast<double>(i * 2 + 2) * dh_);
                    sumlocal += (f0 + 4.0 * f1 + f2);
                }
                return sumlocal;
            },
            std::plus<double>());

        return sum * dh_ / 3.0;
    }

    template <typename FUNCTYPE>
    double Simpson<FUNCTYPE>::operator()(boost::mpl::int_<static_cast<std::int32_t>(ParallelType::Tbb2)>) const
    {
        tbb::combinable<double> sum;

        tbb::parallel_for(
            tbb::blocked_range<std::int32_t>(0, n_ / 2),
            [&](auto const & range) {
                auto sumlocal = 0.0;
                for (auto && i = range.begin(); i != range.end(); ++i) {
                    auto const f0 = func_(x1_ + static_cast<double>(i * 2) * dh_);
                    auto const f1 = func_(x1_ + static_cast<double>(i * 2 + 1) * dh_);
                    auto const f2 = func_(x1_ + static_cast<double>(i * 2 + 2) * dh_);
                    sumlocal += (f0 + 4.0 * f1 + f2);
                }
                sum.local() += sumlocal;
        });

        return sum.combine(std::plus<double>()) * dh_ / 3.0;
    }
    
    // #endregion templateメンバ関数の実装
}

#endif  // _SIMPSON_INTEGRAL_H_

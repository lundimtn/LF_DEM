//
//  TensorBase.h
//  LF_DEM
//
//  Copyright (c) 2017 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class TensorBase
 \brief Base tensor object
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef LF_DEM_TensorBase__
#define LF_DEM_TensorBase__
#include <array>
#include <iterator>


template<std::size_t size_, typename DerivedTensor>
class TensorBase {
public:
	std::array<double, size_> elm;
	/*
	 * (12, 13, 23)
	 */
	double & operator [](unsigned i) {
		return elm[i];
	}

	const double & operator [](unsigned i) const {
		return elm[i];
	}

	auto begin() noexcept ->decltype(std::begin(elm)) {
		return elm.begin();
	}

	auto end() noexcept ->decltype(std::end(elm)) {
		return elm.end();
	}

	constexpr std::size_t size() const {
		return size_;
	}

	DerivedTensor& operator += (const DerivedTensor& s)
	{
		for (unsigned i=0; i<elm.size(); i++) {
			elm[i] += s.elm[i];
		}
		return *this;
	}

	DerivedTensor& operator -= (const DerivedTensor& s)
	{
		for (unsigned i=0; i<elm.size(); i++) {
			elm[i] -= s.elm[i];
		}
		return *this;
	}

	template <typename T>
	DerivedTensor& operator *= (const T& d)
	{
		for (auto &e: elm) {
			e *= d;
		}
		return *this;
	}

	template <typename T>
	DerivedTensor& operator /= (const T& d)
	{
		for (auto &e: elm) {
			e /= d;
		}
		return *this;
	}

	void reset() {
		elm.fill(0);
	}
};

template<std::size_t size_, typename DerivedTensor>
inline DerivedTensor operator + (const TensorBase<size_, DerivedTensor> &s)
{
	return s;
}

template<std::size_t size_, typename DerivedTensor>
inline DerivedTensor operator + (const TensorBase<size_, DerivedTensor> &a1,
                                 const TensorBase<size_, DerivedTensor> &a2)
{
	DerivedTensor a3;
	for (unsigned i=0; i<size_; i++) {
		a3[i] = a1[i] + a2[i];
	}
	return a3;
}

/* subtraction */
template<std::size_t size_, typename DerivedTensor>
inline DerivedTensor operator - (const TensorBase<size_, DerivedTensor> &s)
{
	DerivedTensor ms;
	for (unsigned i=0; i<size_; i++) {
		ms[i] = -s[i];
	}
	return ms;
}

template<std::size_t size_, typename DerivedTensor>
inline DerivedTensor operator - (const TensorBase<size_, DerivedTensor> &a1,
                                 const TensorBase<size_, DerivedTensor> &a2)
{
	DerivedTensor a3;
	for (unsigned i=0; i<size_; i++) {
		a3[i] = a1[i] - a2[i];
	}
	return a3;
}

// output stream operator
template<std::size_t size_, typename DerivedTensor>
inline std::ostream& operator << (std::ostream& out,
                                  const TensorBase<size_, DerivedTensor> &s)
{
	for (auto e: s) {
		out << e << ' ';
	}
	return out;
}

template<typename T, std::size_t size_, typename DerivedTensor>
inline DerivedTensor operator * (const T& a,
                                 const TensorBase<size_, DerivedTensor> &s)
{
	DerivedTensor ms;
	for (unsigned i=0; i<size_; i++) {
		ms[i] = a*s[i];
	}
	return ms;
}

template<typename T, std::size_t size_, typename DerivedTensor>
inline DerivedTensor operator * (const TensorBase<size_, DerivedTensor> &s,
                                 const T &a)
{
	return a*s;
}

template<typename T, std::size_t size_, typename DerivedTensor>
inline DerivedTensor operator / (const TensorBase<size_, DerivedTensor> &s,
                                 const T &a)
{
	DerivedTensor ms;
	for (unsigned i=0; i<size_; i++) {
		ms[i] = s[i]/a;
	}
	return ms;
}

#endif // #ifndef LF_DEM_TensorBase__

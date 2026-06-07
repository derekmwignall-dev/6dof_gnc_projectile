#pragma once

// Generic Taylor expansion
template <int N> // number of iterations/terms
class TaylorExpansion
{
public:
	constexpr double expTaylor(double x = 0.0 ) // evalulated as e^x = sum^N_(k=0){x^k/k!}
	{
		double sum{ 1 }; // total expansion sum, initiated at 1 for x^0/0!
		double term{ 1 }; // value of current term

		for (int k{ 1 }; k <= N; ++k){ // loops k over N terms
			term *= (x / k);
			sum += term;
		}

		return sum;
	}
};

namespace Conversion {
	constexpr double degToRad(double degrees) // converts degrees to radians
	{
		return degrees * (Constants::pi / 180.0);
	}
	constexpr double radToDeg(double radians) // converts radians to degrees
	{
		return radians * (180.0 / Constants::pi);
	}
}
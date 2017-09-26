// ConsoleApplication2.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
//


#include"stdafx.h"
#include<iostream>
#include "Eigen/Dense"
#include <cmath>


Eigen::MatrixXd S(const int step, const double h) {
	Eigen::MatrixXd ret = Eigen::MatrixXd::Zero(step - 2, step - 2);
	for (int i = 0; i < ret.cols(); ++i) {
		ret(i, i) = 2.0 / h;
		if (i < ret.cols() - 1) ret(i+1, i) = -1.0 / h;
		if (i>0) ret(i-1, i) = -1.0 / h;
	}
	return ret;
}

Eigen::MatrixXd B(const int step, const double h) {
	Eigen::MatrixXd ret = Eigen::MatrixXd::Zero(step - 2, step - 2);
	for (int i = 0; i < ret.cols(); ++i) {
		ret(i, i) = 0;
		if (i < ret.cols() -1 ) ret(i , i + 1) = 1.0 / 2.0;
		if (i > 0) ret(i , i -1 ) = -1.0 / 2.0;
	}
	return ret;
}

Eigen::MatrixXd M(const int step, const double h) {
	Eigen::MatrixXd ret = Eigen::MatrixXd::Zero(step - 2, step - 2);
	for (int i = 0; i < ret.cols(); ++i) {
		ret(i, i) = h/6.0 * 4;
		if( i < ret.cols() - 1) ret(i+1 , i) = h /6.0;
		if (i > 0) ret(i-1, i) = h/6.0;
	}
	return ret;
}

Eigen::MatrixXd A(const int step, const double h , const double sigma , const double rate) {
	return sigma * sigma / 2.0 * S(step, h)
		+ (sigma * sigma / 2.0 - rate) * B(step, h) + rate * M(step, h);
}



Eigen::VectorXd u(const int step, const double h , double K,const double smin) {
	Eigen::VectorXd ret = Eigen::VectorXd::Zero(step - 2);
	for (int i = 0; i < ret.rows(); ++i) {
		ret(i) = std::max(std::exp(h * (i + 1) + smin) - K, 0.0);
	}
	return ret;
}


int main()
{
	
	const int step = 500;
	const double smin = -std::log(3);
	const double smax = std::log(110.0);
	const double h = (smax - smin) / step;
	
	
	const int tstep = step;
	const double tmin = 0.0;
	const double tmax = 1.0;
	const double k = (tmax - tmin) / tstep;

	const double sigma = 0.3;
	const double rate = 0.05;
	const double K = 100;

	std::cout << u(step, h, K, smin);

	const Eigen::MatrixXd lhs = M(step, h) + k / 2.0 * A(step, h, sigma, rate);
	const Eigen::MatrixXd rhs = M(step, h) - k / 2.0 * A(step, h, sigma, rate);

	const Eigen::FullPivLU<Eigen::MatrixXd> lu = lhs.fullPivLu();
	Eigen::MatrixXd x = u(step, h, K, smin);

	for (int i = 0; i < tstep; ++i) {
		x = lu.solve(rhs * x);
	}

	std::cout << x << std::endl << std::endl;

    return 0;
}


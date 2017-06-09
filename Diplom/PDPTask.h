#pragma once

#include "stdafx.h"

using namespace std;

class PDPTask
{
public:
	PDPTask(const double a, const double b, const vector<double> u0, const vector<double> leftBorder, const vector<double> rightBorder);
	const double getA() const;
	const double getB() const;
	const vector<double> getU0() const;
	const vector<double> getLeftBorder() const;
	const vector<double> getRightBorder() const;

private:
	PDPTask();
	const double a;
	const double b;
	const vector<double> u0;
	const vector<double> leftBorder;
	const vector<double> rightBorder;
};
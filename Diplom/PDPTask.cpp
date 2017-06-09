#include "stdafx.h"
#include "PDPTask.h"

PDPTask::PDPTask(const double a, const double b, const vector<double> u0, const vector<double> leftBorder, const vector<double> rightBorder) : a(a), b(b), u0(u0), leftBorder(leftBorder), rightBorder(rightBorder)
{
}

const double PDPTask::getA() const
{
	return a;
}

const double PDPTask::getB() const
{
	return b;
}

const vector<double> PDPTask::getU0() const
{
	return u0;
}

const vector<double> PDPTask::getLeftBorder() const
{
	return leftBorder;
}

const vector<double> PDPTask::getRightBorder() const
{
	return rightBorder;
}

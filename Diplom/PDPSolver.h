#pragma once

#include "stdafx.h"
#include "PDPTask.h"

using namespace std;

 class PDPSolver
{
public:
	static vector<vector<double>> getSerialPolynomialSolution(const PDPTask task, const int n, const int m, const double tau);
	static vector<vector<double>> getParallelPolynomialSolution1(const PDPTask task, const int n, const int m, const double tau);
	static vector<vector<double>> getParallelPolynomialSolution2(const PDPTask task, const int n, const int m, const double tau, const int numberOfThreads);
	static vector<vector<double>> getSerialTrigonometricSolution(const PDPTask task, const int n, const int m, const double tau);
	static vector<vector<double>> getParallelTrigonometricSolution1(const PDPTask task, const int n, const int m, const double tau);
	static vector<vector<double>> getParallelTrigonometricSolution2(const PDPTask task, const int n, const int m, const double tau, const int numberOfThreads);

private:
	PDPSolver();
};
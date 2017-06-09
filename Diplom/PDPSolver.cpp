#include "stdafx.h"
#include "PDPSolver.h"

PDPSolver::PDPSolver()
{
}

vector<vector<double>> PDPSolver::getSerialPolynomialSolution(const PDPTask task, const int n, const int m, const double tau)
{
	const double h = (task.getB() - task.getA()) / n;
	vector<vector<double>> solution(n + 1, vector<double>(m));
	for (int i = 0; i <= n; i++)
		solution[i][0] = task.getU0()[i];
	double tmp = tau / h / h;
	double koef1 = tmp;
	double koef2 = 1 - 2 * tmp;
	for (int k = 0; k < m - 1; k++)
	{
		solution[0][k + 1] = task.getLeftBorder()[k + 1];
		solution[n][k + 1] = task.getRightBorder()[k + 1];
		for (int i = 1; i < n; i++)
			solution[i][k + 1] = koef1 * (solution[i + 1][k] + solution[i - 1][k]) + koef2 * solution[i][k];
	}
	return solution;
}

vector<vector<double>> PDPSolver::getParallelPolynomialSolution1(const PDPTask task, const int n, const int m, const double tau)
{
	const double h = (task.getB() - task.getA()) / n;
	vector<vector<double>> solution(n + 1, vector<double>(m));
	for (int i = 0; i <= n; i++)
		solution[i][0] = task.getU0()[i];
	double tmp = tau / h / h;
	double koef1 = tmp;
	double koef2 = 1 - 2 * tmp;
	for (int k = 0; k < m - 1; k++)
	{
		solution[0][k + 1] = task.getLeftBorder()[k + 1];
		solution[n][k + 1] = task.getRightBorder()[k + 1];
#pragma omp parallel for
		for (int i = 1; i < n; i++)
			solution[i][k + 1] = koef1 * (solution[i + 1][k] + solution[i - 1][k]) + koef2 * solution[i][k];
	}
	return solution;
}

vector<vector<double>> PDPSolver::getParallelPolynomialSolution2(const PDPTask task, const int n, const int m, const double tau, const int numberOfThreads)
{
	const double h = (task.getB() - task.getA()) / n;
	vector<vector<double>> solution(n + 1, vector<double>(m));
	for (int i = 0; i <= n; i++)
		solution[i][0] = task.getU0()[i];
	double tmp = tau / h / h;
	double koef1 = tmp;
	double koef2 = 1 - 2 * tmp;
#pragma omp parallel for
	for (int i = 0; i < numberOfThreads; i++)
	{
		int th = omp_get_thread_num();
		int p = n / numberOfThreads;
		for (int k = 0; k < m - 1; k++)
		{
			solution[0][k + 1] = task.getLeftBorder()[k + 1];
			solution[n][k + 1] = task.getRightBorder()[k + 1];
			for (int i = th * p + 1; i < n && i < (th + 1) * p; i++)
				solution[i][k + 1] = koef1 * (solution[i + 1][k] + solution[i - 1][k]) + koef2 * solution[i][k];
		}
	}
	return solution;
}

vector<vector<double>> PDPSolver::getSerialTrigonometricSolution(const PDPTask task, const int n, const int m, const double tau)
{
	const double h = (task.getB() - task.getA()) / n;
	vector<vector<double>> solution(n + 1, vector<double>(m));
	for (int i = 0; i <= n; i++)
		solution[i][0] = task.getU0()[i];
	double tmp = tau / (1 - cos(h));
	double koef1 = tmp / 2;
	double koef2 = 1 - tmp;
	for (int k = 0; k < m - 1; k++)
	{
		solution[0][k + 1] = task.getLeftBorder()[k + 1];
		solution[n][k + 1] = task.getRightBorder()[k + 1];
		for (int i = 1; i < n; i++)
			solution[i][k + 1] = koef1 * (solution[i + 1][k] + solution[i - 1][k]) + koef2 * solution[i][k];
	}
	return solution;
}

vector<vector<double>> PDPSolver::getParallelTrigonometricSolution1(const PDPTask task, const int n, const int m, double tau)
{
	const double h = (task.getB() - task.getA()) / n;
	vector<vector<double>> solution(n + 1, vector<double>(m));
	for (int i = 0; i <= n; i++)
		solution[i][0] = task.getU0()[i];
	double tmp = tau / (1 - cos(h));
	double koef1 = tmp / 2;
	double koef2 = 1 - tmp;
	for (int k = 0; k < m - 1; k++)
	{
		solution[0][k + 1] = task.getLeftBorder()[k + 1];
		solution[n][k + 1] = task.getRightBorder()[k + 1];
#pragma omp parallel for
		for (int i = 1; i < n; i++)
			solution[i][k + 1] = koef1 * (solution[i + 1][k] + solution[i - 1][k]) + koef2 * solution[i][k];
	}
	return solution;
}

vector<vector<double>> PDPSolver::getParallelTrigonometricSolution2(const PDPTask task, const int n, const int m, const double tau, const int numberOfThreads)
{
	const double h = (task.getB() - task.getA()) / n;
	vector<vector<double>> solution(n + 1, vector<double>(m));
	for (int i = 0; i <= n; i++)
		solution[i][0] = task.getU0()[i];
	double tmp = tau / (1 - cos(h));
	double koef1 = tmp / 2;
	double koef2 = 1 - tmp;
#pragma omp parallel for
	for (int i = 0; i < numberOfThreads; i++)
	{
		int th = omp_get_thread_num();
		int p = n / numberOfThreads;
		for (int k = 0; k < m - 1; k++)
		{
			solution[0][k + 1] = task.getLeftBorder()[k + 1];
			solution[n][k + 1] = task.getRightBorder()[k + 1];
			for (int i = th * p + 1; i < n && i < (th + 1) * p; i++)
				solution[i][k + 1] = koef1 * (solution[i + 1][k] + solution[i - 1][k]) + koef2 * solution[i][k];
		}
	}
	return solution;
}

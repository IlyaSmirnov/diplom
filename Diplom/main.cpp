#include "stdafx.h"
#include "PDPTask.h"
#include "PDPSolver.h"

using namespace std;

const double pi = 3.141592653589793238463;

void Task1();
void Task2();
void Task3();

int _tmain(int argc, char*argv[])
{
	cout << "Task 1" << endl;
	Task1();
	cout << "Task 2" << endl;
	Task2();
	cout << "Task 3" << endl;
	Task3();

    return 0;
}

void Task1()
{
	const int numberOfThreads = 4;
	omp_set_num_threads(numberOfThreads);
	const double a = 0;
	const double b = 1;
	const int n = 100;
	const int m = 500;
	const double h = (b - a) / n;
	const double tau = h * h / 2;
	vector<double> u0(n + 1);
	for (int i = 0; i <= n; i++)
		u0[i] = cos(0.5 * (a + i * h)) + (1 - (a + i * h)) * (a + i * h);
	vector<double> leftBorder(m);
	for (int i = 0; i < m; i++)
		leftBorder[i] = exp(-0.25 * i * tau);
	vector<double> rightBorder(m);
	for (int i = 0; i < m; i++)
		rightBorder[i] = exp(-0.25 * i * tau) * cos(0.5);
	double start_time = omp_get_wtime();
	PDPSolver::getSerialPolynomialSolution(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau);
	cout << "time of serial polynomial algorithm " << omp_get_wtime() - start_time << endl;
	start_time = omp_get_wtime();
	PDPSolver::getParallelPolynomialSolution1(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau);
	cout << "time of parallel polynomial algorithm 1 " << omp_get_wtime() - start_time << endl;
	start_time = omp_get_wtime();
	PDPSolver::getParallelPolynomialSolution2(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau, numberOfThreads);
	cout << "time of parallel polynomial algorithm 2 " << omp_get_wtime() - start_time << endl;
	start_time = omp_get_wtime();
	PDPSolver::getSerialTrigonometricSolution(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau);
	cout << "time of serial trigonometric algorithm " << omp_get_wtime() - start_time << endl;
	start_time = omp_get_wtime();
	PDPSolver::getParallelTrigonometricSolution1(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau);
	cout << "time of parallel trigonometric algorithm 1 " << omp_get_wtime() - start_time << endl;
	start_time = omp_get_wtime();
	PDPSolver::getParallelTrigonometricSolution2(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau, numberOfThreads);
	cout << "time of parallel trigonometric algorithm 2 " << omp_get_wtime() - start_time << endl;
}

void Task2()
{
	const int numberOfThreads = 4;
	omp_set_num_threads(numberOfThreads);
	const double a = 0;
	const double b = 2;
	const int n = 100;
	const int m = 500;
	const int mSmall = m * 10;
	const double h = (b - a) / n;
	const double tau = h * h / 2;
	const double tauSmall = h * h / 2 / 10;
	vector<double> u0(n + 1);
	for (int i = 0; i <= n; i++)
		if ((a + i * h) >= 0 && (a + i * h) <= 1)
			u0[i] = a + i * h;
		else
			u0[i] = 2 - (a + i * h);
	vector<double> leftBorder(m);
	for (int i = 0; i < m; i++)
		leftBorder[i] = 0;
	vector<double> rightBorder(m);
	for (int i = 0; i < m; i++)
		rightBorder[i] = 0;
	vector<double> leftBorderSmall(mSmall);
	for (int i = 0; i < m; i++)
		leftBorderSmall[i] = 0;
	vector<double> rightBorderSmall(mSmall);
	for (int i = 0; i < m; i++)
		rightBorderSmall[i] = 0;
	double start_time = omp_get_wtime();
	PDPSolver::getSerialPolynomialSolution(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau);
	cout << "time of serial polynomial algorithm " << omp_get_wtime() - start_time << endl;
	start_time = omp_get_wtime();
	PDPSolver::getParallelPolynomialSolution1(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau);
	cout << "time of parallel polynomial algorithm 1 " << omp_get_wtime() - start_time << endl;
	start_time = omp_get_wtime();
	PDPSolver::getParallelPolynomialSolution2(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau, numberOfThreads);
	cout << "time of parallel polynomial algorithm 2 " << omp_get_wtime() - start_time << endl;
	start_time = omp_get_wtime();
	PDPSolver::getSerialTrigonometricSolution(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau);
	cout << "time of serial trigonometric algorithm " << omp_get_wtime() - start_time << endl;
	start_time = omp_get_wtime();
	PDPSolver::getParallelTrigonometricSolution1(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau);
	cout << "time of parallel trigonometric algorithm 1 " << omp_get_wtime() - start_time << endl;
	start_time = omp_get_wtime();
	PDPSolver::getParallelTrigonometricSolution2(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau, numberOfThreads);
	cout << "time of parallel trigonometric algorithm 2 " << omp_get_wtime() - start_time << endl;
	PDPSolver::getSerialPolynomialSolution(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau);
	start_time = omp_get_wtime();
	PDPSolver::getSerialPolynomialSolution(PDPTask(a, b, u0, leftBorderSmall, rightBorderSmall), n, mSmall, tauSmall);
	cout << "time of serial polynomial algorithm on small grid " << omp_get_wtime() - start_time << endl;
	start_time = omp_get_wtime();
	PDPSolver::getSerialTrigonometricSolution(PDPTask(a, b, u0, leftBorderSmall, rightBorderSmall), n, mSmall, tauSmall);
	cout << "time of serial trigonometric algorithm on small grid " << omp_get_wtime() - start_time << endl;
	ofstream out;
	out.open("Task2Pol.txt");
	vector<vector<double>> solution(n + 1, vector<double>(m));
	solution = PDPSolver::getSerialPolynomialSolution(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau);
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{
			out << solution[i][j] << ' ';

		}
		out << solution[n][j] << endl;
	}
	out.close();
	out.open("Task2Trig.txt");
	solution = PDPSolver::getSerialTrigonometricSolution(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau);
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{
			out << solution[i][j] << ' ';

		}
		out << solution[n][j] << endl;
	}
	out.close();
	vector<vector<double>> solutionSmall(n + 1, vector<double>(mSmall));
	out.open("Task2PolSmall.txt");
	solutionSmall = PDPSolver::getSerialPolynomialSolution(PDPTask(a, b, u0, leftBorderSmall, rightBorderSmall), n, mSmall, tauSmall);
	for (int j = 0; j < mSmall; j += 10)
	{
		for (int i = 0; i < n; i++)
		{
			out << solutionSmall[i][j] << ' ';

		}
		out << solutionSmall[n][j] << endl;
	}
	out.close();
	out.open("Task2TrigSmall.txt");
	solutionSmall = PDPSolver::getSerialTrigonometricSolution(PDPTask(a, b, u0, leftBorderSmall, rightBorderSmall), n, mSmall, tauSmall);;
	for (int j = 0; j < mSmall; j += 10)
	{
		for (int i = 0; i < n; i++)
		{
			out << solutionSmall[i][j] << ' ';

		}
		out << solutionSmall[n][j] << endl;
	}
	out.close();
}

void Task3()
{
	const int numberOfThreads = 4;
	omp_set_num_threads(numberOfThreads);
	const double a = 0;
	const double b = 1;
	const int n = 100;
	const int m = 500;
	const int mSmall = m * 10;
	const double h = (b - a) / n;
	const double tau = h * h / 2;
	const double tauSmall = h * h / 2 / 10;
	vector<double> u0(n + 1);
	for (int i = 0; i <= n; i++)
			u0[i] = 2 * sin(pi * (a + i * h));
	vector<double> leftBorder(m);
	for (int i = 0; i < m; i++)
		leftBorder[i] = 0;
	vector<double> rightBorder(m);
	for (int i = 0; i < m; i++)
		rightBorder[i] = 0;
	vector<double> leftBorderSmall(mSmall);
	for (int i = 0; i < m; i++)
		leftBorderSmall[i] = 0;
	vector<double> rightBorderSmall(mSmall);
	for (int i = 0; i < m; i++)
		rightBorderSmall[i] = 0;
	double start_time = omp_get_wtime();
	PDPSolver::getSerialPolynomialSolution(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau);
	cout << "time of serial polynomial algorithm " << omp_get_wtime() - start_time << endl;
	start_time = omp_get_wtime();
	PDPSolver::getParallelPolynomialSolution1(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau);
	cout << "time of parallel polynomial algorithm 1 " << omp_get_wtime() - start_time << endl;
	start_time = omp_get_wtime();
	PDPSolver::getParallelPolynomialSolution2(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau, numberOfThreads);
	cout << "time of parallel polynomial algorithm 2 " << omp_get_wtime() - start_time << endl;
	start_time = omp_get_wtime();
	PDPSolver::getSerialTrigonometricSolution(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau);
	cout << "time of serial trigonometric algorithm " << omp_get_wtime() - start_time << endl;
	start_time = omp_get_wtime();
	PDPSolver::getParallelTrigonometricSolution1(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau);
	cout << "time of parallel trigonometric algorithm 1 " << omp_get_wtime() - start_time << endl;
	start_time = omp_get_wtime();
	PDPSolver::getParallelTrigonometricSolution2(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau, numberOfThreads);
	cout << "time of parallel trigonometric algorithm 2 " << omp_get_wtime() - start_time << endl;
	PDPSolver::getSerialPolynomialSolution(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau);
	start_time = omp_get_wtime();
	PDPSolver::getSerialPolynomialSolution(PDPTask(a, b, u0, leftBorderSmall, rightBorderSmall), n, mSmall, tauSmall);
	cout << "time of serial polynomial algorithm on small grid " << omp_get_wtime() - start_time << endl;
	start_time = omp_get_wtime();
	PDPSolver::getSerialTrigonometricSolution(PDPTask(a, b, u0, leftBorderSmall, rightBorderSmall), n, mSmall, tauSmall);
	cout << "time of serial trigonometric algorithm on small grid " << omp_get_wtime() - start_time << endl;
	ofstream out;
	out.open("Task3Pol.txt");
	vector<vector<double>> solution(n + 1, vector<double>(m));
	solution = PDPSolver::getSerialPolynomialSolution(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau);
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{
			out << solution[i][j] << ' ';

		}
		out << solution[n][j] << endl;
	}
	out.close();
	out.open("Task3Trig.txt");
	solution = PDPSolver::getSerialTrigonometricSolution(PDPTask(a, b, u0, leftBorder, rightBorder), n, m, tau);
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{
			out << solution[i][j] << ' ';

		}
		out << solution[n][j] << endl;
	}
	out.close();
	vector<vector<double>> solutionSmall(n + 1, vector<double>(mSmall));
	out.open("Task3PolSmall.txt");
	solutionSmall = PDPSolver::getSerialPolynomialSolution(PDPTask(a, b, u0, leftBorderSmall, rightBorderSmall), n, mSmall, tauSmall);
	for (int j = 0; j < mSmall; j += 10)
	{
		for (int i = 0; i < n; i++)
		{
			out << solutionSmall[i][j] << ' ';

		}
		out << solutionSmall[n][j] << endl;
	}
	out.close();
	out.open("Task3TrigSmall.txt");
	solutionSmall = PDPSolver::getSerialTrigonometricSolution(PDPTask(a, b, u0, leftBorderSmall, rightBorderSmall), n, mSmall, tauSmall);;
	for (int j = 0; j < mSmall; j += 10)
	{
		for (int i = 0; i < n; i++)
		{
			out << solutionSmall[i][j] << ' ';

		}
		out << solutionSmall[n][j] << endl;
	}
	out.close();
}
//
// Created by sin_een on 1/15/19.
//

#include <iostream>
#include "Matrix.hpp"
#include <eigen3/Eigen/Dense>
#include <stack>
#include <chrono>

using Eigen::MatrixXd;

using namespace std;

//std::stack<clock_t> tictoc_stack;
stack<chrono::time_point<chrono::system_clock>> tictoc_stack;

/**
 *
 */
void tic()
{
	//tictoc_stack.push(clock());
	tictoc_stack.push(chrono::system_clock::now());
}

/**
 *
 */
void toc()
{
	// std::cout << "Time elapsed: "
	// 		<< ((double)(clock() - tictoc_stack.top())) / CLOCKS_PER_SEC
	// 		<< std::endl;
	// tictoc_stack.pop();
	chrono::duration<double> elapsed_seconds = chrono::system_clock::now() - tictoc_stack.top();
	cout << elapsed_seconds.count() << "\n";
	tictoc_stack.pop();
}

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		throw "wrong number of parameters";
	}
	unsigned int size = (unsigned int)stoi(argv[1]);
	if (size > 500)
	{
		throw "wrong argument";
	}
	MatrixXd mE1 = MatrixXd::Random(size, size);
	MatrixXd mE2 = MatrixXd::Random(size, size);

	Matrix m = Matrix<int>(size, size);
	for ( unsigned int r = 0 ; r < size; r++ )
	{
		for(unsigned int c = 0; c < size; c++)
		{
			m(r, c) = 1;
		}
	}
	// Eigen

	cout << "size " << size << endl ;



	tic();
	mE1 * mE2;
	cout << "eigen mult ";
	toc();

//	flush;

	tic();
	mE1 + mE2;
	cout << "eigen add ";
	toc();

	//mine

//	flush;

	tic();
	m* m;
	cout << "matlib mult ";
	toc();

//	flush;

	tic();
	m + m;
	cout << "matlib add ";
	toc();

	return 0;
}
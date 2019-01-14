// ============================ Includes ============================

#include "Matrix.hpp"
#include <iostream>
#include <sstream>

// ================== Change this variable if you want ==================

// "true" if you want to debug one error at a time
bool oneErrorAtATime = false;

// ======================================================================

int numOfErrors = 0;

void checkPrint(std::stringstream &buffer, const std::string &answer, const std::string &tested)
{
	if (buffer.str() != answer)
	{
		numOfErrors++;
		std::cout << "ERROR " << numOfErrors << ":" << std::endl;
		std::cout << tested << std::endl;
		std::cout << "Your answer: " << std::endl << "\"" << buffer.str() << "\"" << std::endl;
		std::cout << "Right answer: " << std::endl << "\"" << answer << "\"" << std::endl << std::endl;
		std::cout.flush();
		if (oneErrorAtATime)
		{
			std::exit(1);
		}
	}
	buffer.str(std::string());
	buffer.clear();
}

void checkBool(bool statement, std::string &message)
{
	if (statement)
	{
		numOfErrors++;
		std::cout << "ERROR " << numOfErrors << ":" << std::endl;
		std::cout << message << std::endl << std::endl;
		if (oneErrorAtATime)
		{
			std::exit(1);
		}
	}
}

void testMatrix()
{

	std::stringstream buffer;
	Matrix<int> a(3, 1, {4, 5, 3});
	Matrix<int> b(1, 3, {4, 5, 3});
	Matrix<int> ab(3, 1, {4, 5, 3});
	Matrix<double> c(3, 2, {4, 5, 6, 7, 8, 9});
	Matrix<double> d(2, 3, {9, 8, 7, 6, 5, 4});
	Matrix<Complex> e(2, 2, {{1, 3}, {4, 4}, {2, 3}, {2, 1}});
	Matrix<Complex> f(2, 2, {{7, 7}, {6, 6}, {5, 5}, {4, 4}});
	bool test;
	std::string message;

	// Test zero constructor
	Matrix<int> m;
	test = m.rows() != 1 || m.cols() != 1 || m(0, 0) != 0;
	message = "Default constructor should init a matrix of size 1 x 1 with zero value\n";
	message += "Yours initiated a matrix of " + std::to_string(m.rows()) + " x " + std::to_string(m.cols());
	message += " with M[0][0]=" + std::to_string(m(0, 0));
	checkBool(test, message);

	// Test (int,int) constructor
	Matrix<int> n(3, 3);
	buffer << n;
	checkPrint(buffer, "0\t0\t0\t\n0\t0\t0\t\n0\t0\t0\t\n", "Test (3,3) constructor");

	// Test (int,int,vector) constructor
	Matrix<int> l(3, 3, {1, 2, 3, 4, 5, 6, 7, 8, 9});
	buffer << l;
	checkPrint(buffer, "1\t2\t3\t\n4\t5\t6\t\n7\t8\t9\t\n", "Test (3,3,vector) constructor");

	// Test copy constructor
	Matrix<int> p(l);
	buffer << p;
	checkPrint(buffer, "1\t2\t3\t\n4\t5\t6\t\n7\t8\t9\t\n", "Test copy constructor");

	// Test operator=
	Matrix<int> z(2, 2);
	l = z;
	buffer << l;
	checkPrint(buffer, "0\t0\t\n0\t0\t\n", "Test operator=");
	test = l.cols() != 2 || l.rows() != 2;
	message = "Operator= did not change the sizes of row and col variables after copy.";
	checkBool(test, message);

	// Make sure copy constructor made a "deep copy"
	buffer << p;
	checkPrint(buffer, "1\t2\t3\t\n4\t5\t6\t\n7\t8\t9\t\n", "You did not make a deep copy in copy constructor");

	// Test operator+
	Matrix<Complex> g = e + f;
	buffer << g;
	checkPrint(buffer, "8 + 10i\t10 + 10i\t\n7 + 8i\t6 + 5i\t\n", "Test operator+");

	// Test operator-
	g = e - f;
	buffer << g;
	checkPrint(buffer, "-6 - 4i\t-2 - 2i\t\n-3 - 2i\t-2 - 3i\t\n", "Test operator-");

	// Test operator*
	Matrix<double> t = c * d;
	buffer << t;
	checkPrint(buffer, "66\t57\t48\t\n96\t83\t70\t\n126\t109\t92\t\n", "Test operator*");

	// Test operator== and operator!=
	test = a == b || !(a != b) || !(a == ab) || a != ab;
	message = "Operator== and/or Operator!= are wrong\n";
	message += "Matrix 3 x 2 and matrix 2 x 3 with same values are not equal";
	checkBool(test, message);

	// Test operator() (const and non-const)
	const Matrix<int> checkConst(3, 1, {4, 5, 3});
	int i = checkConst(1, 0); // If this line has an error, your const operator() is wrong
	test = d(1, 2) != 4 || i != 5;
	message = "Operator() is wrong";
	checkBool(test, message);

	// Test iterators
	auto begin = d.begin();
	auto end = d.end();
	for (; begin != end; begin++)
	{
		buffer << *begin << "\t";
	}
	checkPrint(buffer, "9\t8\t7\t6\t5\t4\t", "Test iterator");
	// Make sure iterator is a const iterator
	begin = d.begin();

	// =======================================================
	// Uncomment next line and make sure it has a compile error
	// =======================================================
	// *begin = 5;

	// Check isSquareMatrix
	test = !e.isSquareMatrix() || d.isSquareMatrix();
	message = "isSquareMatrix method is wrong";
	checkBool(test, message);

	// Check rows and cols
	test = d.rows() != 2 || d.cols() != 3;
	message = "rows method or cols method are wrong";
	checkBool(test, message);

	// Check template transpose method
	Matrix<double> mat(4, 4, {7, 12, 1, 14, 2, 13, 8, 11, 16, 3, 10, 5, 9, 6, 15, 4});
	Matrix<double> trans = mat.trans();
	buffer << trans;
	checkPrint(buffer, "7\t2\t16\t9\t\n12\t13\t3\t6\t\n1\t8\t10\t15\t\n14\t11\t5\t4\t\n",
			   "Test template trans()");

	// Check complex transpose method
	Matrix<Complex> com(2, 2, {{1, -3}, {4, 4}, {2, -3}, {2, -1}});
	Matrix<Complex> comTrans = com.trans();
	buffer << comTrans;
	checkPrint(buffer, "1 + 3i\t2 + 3i\t\n4 - 4i\t2 + 1i\t\n", "Test complex trans()");
}

int main()
{
	testMatrix();
	if (numOfErrors == 0)
	{
		std::cout << "All tests passed";
	}
	else
	{
		std::cout << "Num of errors: " << numOfErrors << std::endl;
	}
	return 0;
}
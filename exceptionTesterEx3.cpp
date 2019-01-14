// Evyatar.
// Ver 1.1

#include <iostream>
#include "Matrix.hpp"
#include <vector>


using namespace std;

int main()
{

	vector<int> missed;

	cout << endl;
	cout << "This tester only test the existence of exceptions, and their 'what()' functions."
		 << endl;
	cout << "It is recommended to also use Idan Gabay's tester." << endl << endl;

//
//	cout << "Exception num " << 1 << ":" << endl;
//	try
//	{
//		Matrix<int> m1(5, 5, {1, 2, 3});
//		cout << "You missed this exception" << endl;
//		missed.push_back(1);
//	}
//	catch (exception &e)
//	{
//		cout << "Problem:\t\t\t" << "vector is size is not rows*cols" << endl;
//		cout << "Your explanation:\t" << e.what() << endl;
//	}
//	cout << endl;
//
//	cout << "Exception num " << 2 << ":" << endl;
//	try
//	{
//		Matrix<int> m1(2, 2, {1, 2, 3, 4, 5});
//		cout << "You missed this exception" << endl;
//		missed.push_back(2);
//
//	}
//	catch (exception &e)
//	{
//		cout << "Problem:\t\t\t" << "vector is size is not rows*cols" << endl;
//		cout << "Your explanation:\t" << e.what() << endl;
//	}
//	cout << endl;

	cout << "Exception num " << 3 << ":" << endl;
	try
	{
		Matrix<int> m1(2, 2), m2(2, 3);
		m1 + m2;
		cout << "You missed this exception" << endl;
		missed.push_back(3);

	}
	catch (exception &e)
	{
		cout << "Problem:\t\t\t" << "Matrices not the same size +" << endl;
		cout << "Your explanation:\t" << e.what() << endl;
	}
	cout << endl;

	cout << "Exception num " << 4 << ":" << endl;
	try
	{
		Matrix<int> m1(3, 3), m2(2, 3);
		m1 - m2;
		cout << "You missed this exception" << endl;
		missed.push_back(4);

	}
	catch (exception &e)
	{
		cout << "Problem:\t\t\t" << "Matrices not the same size -" << endl;
		cout << "Your explanation:\t" << e.what() << endl;
	}
	cout << endl;

	cout << "Exception num " << 5 << ":" << endl;
	try
	{
		Matrix<int> m1(3, 4), m2(3, 4);
		m1 * m2;
		cout << "You missed this exception" << endl;
		missed.push_back(5);

	}
	catch (exception &e)
	{
		cout << "Problem:\t\t\t" << "Matrices not fitting size *" << endl;
		cout << "Your explanation:\t" << e.what() << endl;
	}
	cout << endl;

	cout << "Exception num " << 6 << ":" << endl;
	try
	{
		Matrix<int> m1(3, 4);
		m1.trans();
		cout << "You missed this exception" << endl;
		missed.push_back(6);

	}
	catch (exception &e)
	{
		cout << "Problem:\t\t\t" << "matrix is not squared" << endl;
		cout << "Your explanation:\t" << e.what() << endl;
	}
	cout << endl;

	cout << "Exception num " << 7 << ":" << endl;
	try
	{
		Matrix<Complex> m1(3, 4);
		m1.trans();
		cout << "You missed this exception" << endl;
		missed.push_back(7);

	}
	catch (exception &e)
	{
		cout << "Problem:\t\t\t" << "matrix is not squared (Complex)" << endl;
		cout << "Your explanation:\t" << e.what() << endl;
	}
	cout << endl;

	cout << "Exception num " << 8 << ":" << endl;
	try
	{
		Matrix<Complex> m1(5, 7);
		m1(-1, 3);
		missed.push_back(8);

		cout << "You missed this exception" << endl;
	}
	catch (exception &e)
	{
		cout << "Problem:\t\t\t" << "array index out of bounds" << endl;
		cout << "Your explanation:\t" << e.what() << endl;
	}
	cout << endl;

	cout << "Exception num " << 9 << ":" << endl;
	try
	{
		Matrix<Complex> m1(5, 7);
		m1(4, -4);
		cout << "You missed this exception" << endl;
		missed.push_back(9);

	}
	catch (exception &e)
	{
		cout << "Problem:\t\t\t" << "array index out of bounds" << endl;
		cout << "Your explanation:\t" << e.what() << endl;
	}
	cout << endl;

	cout << "Exception num " << 10 << ":" << endl;
	try
	{
		Matrix<Complex> m1(5, 7);
		m1(0, 7);
		cout << "You missed this exception" << endl;
		missed.push_back(10);

	}
	catch (exception &e)
	{
		cout << "Problem:\t\t\t" << "array index out of bounds" << endl;
		cout << "Your explanation:\t" << e.what() << endl;
	}
	cout << endl;

	cout << "Exception num " << 11 << ":" << endl;
	try
	{
		Matrix<Complex> m1(5, 7);
		m1(5, 2);
		cout << "You missed this exception" << endl;
		missed.push_back(11);

	}
	catch (exception &e)
	{
		cout << "Problem:\t\t\t" << "array index out of bounds" << endl;
		cout << "Your explanation:\t" << e.what() << endl;
	}
	cout << endl;

	cout << "Exception num " << 12 << ":" << endl;
	try
	{
		Matrix<Complex> m1(-2, 7);
		cout << "You missed this exception" << endl;
		missed.push_back(12);

	}
	catch (exception &e)
	{
		cout << "Problem:\t\t\t" << "allocation problem" << endl;
		cout << "Your explanation:\t" << e.what() << endl;
	}
	cout << endl;

	cout << "Exception num " << 13 << ":" << endl;
	try
	{
		Matrix<Complex> m1(-2, -2, {1, 2, 3, 4});
		cout << "You missed this exception" << endl;
		missed.push_back(13);

	}
	catch (exception &e)
	{
		cout << "Problem:\t\t\t" << "allocation problem" << endl;
		cout << "Your explanation:\t" << e.what() << endl;
	}
	cout << endl;

//	cout << "Exception num " << 14 << ":" << endl;
//	try
//	{
//		Matrix<Complex> m1(5, 5);
//		auto begin = m1.begin();
//		--begin;
//		missed.push_back(14);
//
//	}
//	catch (exception &e)
//	{
//		cout << "Problem:\t\t\t" << "reached end of iterator" << endl;
//		cout << "Your explanation:\t" << e.what() << endl;
//	}
//	cout << endl;
//
//	cout << "Exception num " << 15 << ":" << endl;
//	try
//	{
//		Matrix<Complex> m1(5, 5);
//		auto end = m1.end();
//		++end;
//		missed.push_back(15);
//
//	}
//	catch (exception &e)
//	{
//		cout << "Problem:\t\t\t" << "reached end of iterator" << endl;
//		cout << "Your explanation:\t" << e.what() << endl;
//	}
//	cout << endl;
//==================================================================================
	if (missed.empty())
	{
		cout << "You found all of the faults! Now go up and make sure that you print the relevant"
				" text." << endl;
	} else
	{
		cout << "You missed the exceptions:\t";
		for (int &i : missed)
		{
			cout << i << ' ';
		}
		cout << endl;
	}

	return 0;
}

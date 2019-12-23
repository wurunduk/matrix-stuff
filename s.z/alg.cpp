#include <math.h>

#include "alg.h"
#include "matrix.h"

const double epsilon = 1e-100;

void TridiagonalRotation(double* a, int n)
{
	double a_ii;
	double a_ij;
	double a_ji;
	double a_jj;

	for (int i = 1; i < n - 1; ++i)
	{
		for (int j = i + 1; j < n; ++j)
		{
			//get (x y) vector whcih is gonna be transformed to (1 0)*N
			double x = a[i * n + i - 1];
			double y = a[j * n + i - 1];

			if (fabs(y) < epsilon)
				continue;

			double len = sqrt(x * x + y * y);

			if (len < epsilon)
				continue;

			double cos = x / len;
			double sin = -y / len;

			//as a aresult we get (1 0)*N
			a[i * n + i - 1] = len;
			a[j * n + i - 1] = .0;

			//matrix is symmetrical, we will get the same result if mirrored against xy axis
			a[(i - 1) * n + i] = len;
			a[(i - 1) * n + j] = .0;

			for (int k = i + 1; k < n; k++)
			{
				if (k == j)
					continue;

				x = a[i * n + k];
				y = a[j * n + k];

				a[i * n + k] = x * cos - y * sin;
				a[j * n + k] = x * sin + y * cos;

				a[k * n + i] = x * cos - y * sin;
				a[k * n + j] = x * sin + y * cos;
			}

			x = a[i * n + i];
			y = a[j * n + j];

			double z = a[i * n + j];
			double v = a[j * n + i];

			a_ii = x * cos - v * sin;
			a_jj = z * sin + y * cos;

			a_ij = z * cos - y * sin;
			a_ji = x * sin + v * cos;

			a[i * n + i] = a_ii * cos - a_ij * sin;
			a[j * n + j] = a_ji * sin + a_jj * cos;

			a[i * n + j] = a_ii * sin + a_ij * cos;
			a[j * n + i] = a_ii * sin + a_ij * cos;
		}
	}
}

int sign_count(const double* a, int n, double lambda)
{
	int signs = 0;

	double b = a[0] - lambda;
	if(b < .0) signs++;

	for (int i = 1; i < n; ++i)
	{
		if (fabs(b) < epsilon)
			b = 1e-10;

		b = a[i * n + i] - lambda - a[i * n + i - 1] * a[(i - 1) * n + i] / b;

		if (b < 0)
			signs++;
	}

	return signs;
}

int FindValues(double* a, double* values, int n, double eps)
{
	int count = 0;
	int iterations = 0;

    double norm = Length(a, n, n);
    
	double right_bound = norm + eps;
	double left_bound = -right_bound;

	double current_left = left_bound;
	double current_right = right_bound;
    
	//Transfrom matrix to tridiagonal form
	TridiagonalRotation(a, n);

	for (int i = 0; i < n;)
	{
		while (current_right - current_left > eps)
		{
			double d = (current_left + current_right)/2.0;
			if (sign_count(a, n, d) < i + 1)
				current_left = d;
			else
				current_right = d;
            
			iterations++;
		}

		count = sign_count(a, n, current_right) - sign_count(a, n, current_left);

		for (int j = 0; j < count; j++)
			values[i + j] = (current_left + current_right)/2.0;

		i += count;

		current_left = (current_left + current_right)/2.0;
		current_right = right_bound;
	}

	return iterations;
}

#include "halton.h"
static int _halton_i = 40;
static int _halton_b = 2;

double Halton(int b, int i)
{
	_halton_i = i;
	_halton_b = b;
	
	double f = 1;
	double r = 0;
	
	while(i > 0)
	{
		f = f / b;
		r = r + f * (i % b);
		i = i / b;
	}
	return r;
}



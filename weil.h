#pragma once
#include "types.h"

two_int weil(int enc1, int enc2, Curve * group1, Curve * group2)
{
	int real = 1;
	int imag = 2;
	two_int temp = {real, imag};
	return temp;
}

two_int weil_one(int enc1, Point generator, Curve * group1, Curve * group2)
{
	int real = 1;
	int imag = 2;
	two_int temp = { real, imag };
	return temp;
}
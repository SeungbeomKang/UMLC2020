#pragma once
#include "types.h"

Curve* EC_extend(struct Curve * point)
{
	struct Curve* empty;
	static struct Curve p1;
	p1.a = 0;
	p1.b = 0;
	p1.x = 0;
	p1.y = 0;
	p1.order = 0;
	empty = &p1;

	return empty;
}
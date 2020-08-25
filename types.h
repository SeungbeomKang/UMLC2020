#pragma once
typedef struct Point
{
	int a;
	int b;
	int x;
	int y;
	int order;
}Point;

typedef struct Curve
{
	int a;
	int b;
	int x;
	int y;
	int order;
	struct Point gen;
}Curve;

typedef struct two_int
{
	int first;
	int second;
}two_int;
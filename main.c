#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "types.h"
#include "EC_points.h"
#include "EC_extend.h"
#include "max_order.h"
#include "rat.h"
#include "merge.h"
#include "lcm.h"
#include "rmod.h"
#include "poly_calculate.h"
#include "ismember.h"
#include "EC_pmult.h"
#include "weil.h"
#include "cmod.h"

/* System Params */
int main(void)
{
	int p = 43;
	int q = p * p - 1;
	int k = 2;
	int a = 1;
	int b = 0;

	Point init_point = {0};

	Curve curve;
	curve.a = a;
	curve.b = b;
	curve.x = 0;
	curve.y = 0;
	curve.order = 0;
	curve.gen = init_point;

	Curve * Group1, * Group2;
	Group1 = (Curve *)malloc(sizeof(Curve)*p);
	Group2 = (Curve *)malloc(sizeof(Curve)*p);
	Group1 = EC_points(curve);
	Group2 = EC_extend(Group1);
	
	int num_points = sizeof(Group1) / sizeof(Curve);
	int* Orders = (int *)malloc(sizeof(int)*num_points);

	for (int i = 0; i < num_points; i++)
	{
		Orders[i] = Group1[i].order;
	}

	// maximum Value and its Index
	two_int vi = max_order(Orders);

	Point g = Group1[vi.second].gen;
	Point h = Group2[vi.second].gen;

	for (int i = 0; i < num_points; i++)
	{
		Group1[i].gen = g;
		Group2[i].gen = h;
	}



/* From function to R1CS */
	int A[6][4] =
	{
		{0,0,0,5},
		{1,0,1,0},
		{0,0,0,0},
		{0,1,0,0},
		{0,0,1,0},
		{0,0,0,1}
	};

	int B[6][4] =
	{
		{0,0,1,1},
		{1,1,0,0},
		{0,0,0,0},
		{0,0,0,0},
		{0,0,0,0},
		{0,0,0,0}
	};

	int C[6][4] =
	{
		{0,0,0,0},
		{0,0,0,0},
		{0,0,0,1},
		{1,0,0,0},
		{0,1,0,0},
		{0,0,1,0}
	};

	int Z[5] = {24, -50, 35, -10, 1};
	int R[6] = {1, 3, 35, 9, 27, 30};

	int Ind_mid[3] = { 4, 5, 6 };
	int Ind_IO[3] = { 1, 2, 3 };

	int NumWires = sizeof(A)/sizeof(A[0]);
	int NumGates = sizeof(A[0])/sizeof(int);



/* From R1CS to QAP */
	double Ap[6][4] =
	{
		{-5.0, 9.166666666666666, -5.0, 0.8333333333333334},
		{8.0, -11.333333333333332, 5.0, -0.6666666666666666},
		{0.0, 0.0, 0.0, 0.0},
		{-6.0, 9.5, -4.0, 0.5},
		{4.0, -7.0, 3.5, -0.5},
		{-1.0, 1.8333333333333333, -1.0, 0.16666666666666666},
	};

	double Bp[6][4] =
	{
		{3.0, -5.166666666666667, 2.5, -0.33333333333333337},
		{-2.0, 5.166666666666667, -2.5, 0.33333333333333337},
		{0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0}
	};

	double Cp[6][4] =
	{
		{0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0 },
		{-1.0, 1.8333333333333333, -1.0, 0.16666666666666666},
		{4.0, -4.333333333333333, 1.5, -0.16666666666666666},
		{-6.0, 9.5, -4.0, 0.5},
		{4.0, -7.0, 3.5, -0.5}
	};

	double ABp[12][4] = { 0 };
	double ABCp[18][4] = { 0 };
	merge(Ap, Bp, ABp);
	merge(ABp, Cp, ABCp);

	int NUM[18][4] = { 0 };
	int DEN[18][4] = { 0 };
	
	rat(ABCp, NUM);
	rat(ABCp, DEN);

	int mul = lcm(DEN);

	int mul2 = mul * mul;
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			Ap[i][j] *= mul;
			Bp[i][j] *= mul;
			Cp[i][j] *= mul2;
		}
	}

	rmod_2(Ap, q);
	rmod_2(Bp, q);
	rmod_2(Cp, q);

	rmod_1(Z, q);
	rmod_1(R, q);

	double linear_R_Ap[6][4] = { 0 };
	double linear_R_Bp[6][4] = { 0 };
	double linear_R_Cp[6][4] = { 0 };
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			linear_R_Ap[i][j] = R[i] * Ap[i][j];
			linear_R_Bp[i][j] = R[i] * Bp[i][j];
			linear_R_Cp[i][j] = R[i] * Cp[i][j];
		}
	}
	double multiplied_poly[4] = { 0 };
	polymult(linear_R_Ap, linear_R_Bp, multiplied_poly);
	double P_x[4] = { 0 };
	polysubtract(multiplied_poly, linear_R_Cp, P_x);

	rmod_1(P_x, q);

	double H_x[4] = { 0 };
	int rem = 0;
	polydiv(P_x, Z, q, H_x, &rem);
	rmod_1(H_x, q);



/* Setup ( CRS creation ) */
	srand((unsigned int)time(NULL));
	int alpha = ((rand() << 15 | rand()) % (p*p - 4)) + 2;
	int beta = ((rand() << 15 | rand()) % (p*p - 4)) + 2;
	int gamma = ((rand() << 15 | rand()) % (p*p - 4)) + 2;
	int delta = ((rand() << 15 | rand()) % (p*p - 4)) + 2;

	int x_val = (rand() % (q - 1)) + 1;
	int check[4] = { 1,2,3,4 };
	while (ismember(x_val, check))
	{
		x_val = (rand()%(q-1))+1;
	}

	int tau[5] = {alpha, beta, gamma, delta, x_val};

	int Ap_val[6] = { 0 };
	int Bp_val[6] = { 0 };
	int Cp_val[6] = { 0 };
	polyeval_array(Ap, x_val, Ap_val);
	polyeval_array(Bp, x_val, Bp_val);
	polyeval_array(Cp, x_val, Cp_val);
	int Z_val = polyeval(Z, x_val);
	int H_x_val = polyeval(H_x, x_val);
	rmod_1(Ap_val, q);
	rmod_1(Bp_val, q);
	rmod_1(Cp_val, q);
	rmod(&Z_val, q);
	rmod(&H_x_val, q);

	char sigma1_1[3] =
	{
		EC_pmult(alpha, g),
		EC_pmult(beta, g),
		EC_pmult(delta, g)
	};

	char sigma1_3[3] = { 0 };
	int VAL[3] = { 0 };
	for (int i = 0; i < sizeof(Ind_IO) / sizeof(int); i++)
	{
		VAL[i] = beta * Ap_val[i+3] + alpha * Bp_val[i+3] + Cp_val[i+3];
		rmod(&VAL[i], q);
		sigma1_3[i] = EC_pmult(VAL[i], g);
	}

	char sigma1_4[6] = { 0 };
	int val[6] = { 0 };
	for (int i = 0; i < sizeof(Ind_mid) / sizeof(int); i++)
	{
		val[i] = beta * Ap_val[i] + alpha * Bp_val[i] + Cp_val[i];
		rmod(&val[i],q);
		sigma1_4[i] = EC_pmult(val[i], g);
	}

	char sigma2_1[3] =
	{
		EC_pmult(beta, h),
		EC_pmult(gamma, h),
		EC_pmult(delta, h)
	};

/* Prove */
	int r = (rand() % (q - 1)) + 1;
	int s = (rand() % (q - 1)) + 1;
	
	int A_plain_proof = alpha + r * delta;
	for (int i = 0; i < NumWires; i++)
	{
		A_plain_proof += (R[i] * Ap_val[i]);
	}

	int B_plain_proof = beta + s * delta;
	for (int i = 0; i < NumWires; i++)
	{
		B_plain_proof += (R[i] * Bp_val[i]);
	}

	int C_plain_proof = (H_x_val*Z_val) + (s*A_plain_proof*delta) + (r*B_plain_proof*delta);
	int rsdelta = (-r * s*delta*delta);
	rmod(&rsdelta, q);
	C_plain_proof += rsdelta;
	int linear_R_Cpval[3] = { 0 };
	for (int i = 0; i < sizeof(Ind_mid)/sizeof(int); i++)
	{
		C_plain_proof += R[i] * (beta * Ap_val[i] + alpha * Bp_val[i] + Cp_val[i]);
	}

	char proof[3] =
	{
		EC_pmult(A_plain_proof,g),
		EC_pmult(C_plain_proof,g),
		EC_pmult(B_plain_proof,h)
	};
	printf("Proof is ready.\n");


/* Verify */
	two_int LHS = weil(proof[0], proof[2], Group1, Group2);
	cmod(&LHS, p);

	two_int RHS = { 1,0 };
	two_int expon = { 0 };
	RHS = weil(sigma1_1[0], sigma2_1[0], Group1, Group2);
	cmod(&RHS, p);
	for (int i = 0; i < sizeof(Ind_IO) / sizeof(int); i++)
	{
		expon = weil_one(sigma1_3[i], h, Group1, Group2);
		cmod(&expon, p);
		for (int j = 0; j < R[i]; j++)
		{
			RHS.first = RHS.first * expon.first;
			RHS.second = RHS.second * expon.second;
			cmod(&RHS, p);
		}
	}
	two_int last_term = { 0 };
	last_term = weil_one(proof[1], h, Group1, Group2);
	RHS.first = RHS.first * last_term.first;
	RHS.first = RHS.second * last_term.second;
	cmod(&RHS, p);

	int VfyResult = ((LHS.first == RHS.first) && (LHS.second == RHS.second));
	
	if ((RHS.first == 1 && RHS.second == 0) && (LHS.first == 1 && LHS.second == 0))
		printf("Unable to verify\n");
	else if (VfyResult ==1)
		printf("Verification Success\n");
	else
		printf("Verification Failure\n");


}
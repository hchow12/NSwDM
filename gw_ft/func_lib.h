#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_integration.h>
#include "setting.h"

/*
------------------------------------------------------------------------------------------

Function to calculate polytropic constants

------------------------------------------------------------------------------------------
*/

/* The *kap_b function return a pointer of array of size 2 by accepting a polytropic index
 * The pointer contain polytropic constant read from the poly_const_table.csv file 
 */
double *kap_b(double target) // target = polytropic index n = 0.5, 0.6, 0.8, etc.
{
	FILE *fp;
	char row[1000]; 
	const char s[2] = ","; // delimiter ","
	char *token;
	int pos = 1; // The column number to read from poly_const_table.csv
	int i;
	int j = 0;

	double *result = malloc(3 * sizeof(double));

	fp = fopen("poly_const_table.csv", "r"); // Read the whole content of the poly_const_table.csv

	fgets(row, 1000, fp); // Store the content of poly_const_table.csv into row buffer
	token = strtok(row, s); // Get the first token from row string delimited by ","
	
	// Search which column to read from polytropic index
	while (token != NULL && strtod(token, NULL) != target){
		token = strtok(NULL, s); // End the while loop
		pos++;
	}

	// Read polytropic constant kappa_n and b_n (due to DM core)
	while (j < 2){
		fgets(row, 1000, fp);
		token = strtok(row, s);
		for (i = 1; i < pos; i++) // Read until the column for n = target is found
			token = strtok(NULL, s);
		*(result + j) = strtod(token, NULL); // Convert the polytropic constant string (token) to double
		j++;
	}

	fclose(fp); // Close the poly_const_table.csv file

	return result;
}

/*
------------------------------------------------------------------------------------------

PN term Hybrid section

------------------------------------------------------------------------------------------
PN1 + 2, PN25 passed test with python
PN1 + 2 + 3, PN25 passed test with python
*/

/* r = orbital separation (m)
 * rv = first time derivative of r
 * ov = first time derivative of orbital phase 
 * m = binary total mass
 */

double speed(double r, double rv, double ov)
{
	double s1 = pow(rv, 2.0) + pow(r * ov, 2.0);
	return pow(s1, 1.0 / 2.0);
}

double AS(double r, double rv, double ov, double m)
{
	double v = speed(r, rv, ov);
	double s1 = (1.0 - m / r) / pow((1.0 + m / r), 3.0);
	double s2 = -(m / r) * pow(rv, 2.0);
	double s3 = (2.0 - m / r) / (1.0 - pow(m / r, 2.0));
	return s1 + s2 * s3 + pow(v, 2.0);
}

double BS(double r, double m)
{
	double s = -(4.0 - 2.0 * (m / r)) / (1.0 - pow(m / r, 2.0));
	return s;
}

double A1PN(double r, double rv, double ov, double m, double eta)
{
	double v = speed(r, rv, ov);
	double s1 = -3.0 * eta * pow(v, 2.0);
	double s2 = (3.0 / 2.0) * eta * pow(rv, 2.0);
	double s3 = 2.0 * eta * (m / r);
	return s1 + s2 + s3;
}

double B1PN(double eta)
{
	return -2 * eta;
}

double A2PN(double r, double rv, double ov, double m, double eta)
{
	double v = speed(r, rv, ov);
	double s1 = -eta * (3.0 - 4.0 * eta) * pow(v, 4.0);
	double s2 = (1.0 / 2.0) * eta * (13.0 - 4.0 * eta) * (m / r) * pow(v, 2.0);
	double s3 = (3.0 / 2.0) * eta * (3.0 - 4.0 * eta) * pow(rv, 2.0) * pow(v, 2.0);
	double s4 = eta * (25.0 + 2.0 * eta) * (m / r) * pow(rv, 2.0);
	double s5 = -(15.0 / 8.0) * eta * (1.0 - 3 * eta) * pow(rv, 4.0);
	double s6 = -(87.0 / 4.0) * eta * pow((m / r), 2.0);
	return s1 + s2 + s3 + s4 + s5 + s6; 
}

double B2PN(double r, double rv, double ov, double m, double eta)
{
	double v = speed(r, rv, ov);
	double s1 = (1.0 / 2.0) * eta * (15.0 + 4.0 * eta) * pow(v, 2.0);
	double s2 = -(3.0 / 2.0) * eta * (3.0 + 2.0 * eta) * pow(rv, 2.0);
	double s3 = -(1.0 / 2.0) * eta * (41.0 + 8.0 * eta) * (m / r);
	return s1 + s2 + s3;
}

double A3PN(double r, double rv, double ov, double m, double eta)
{
	double v = speed(r, rv, ov);
	double s1 = eta * ((1399.0 / 12.0 - 41.0 * pow(M_PI, 2.0) / 16.0) + (71.0 / 2.0) * eta) * pow((m / r), 3.0);				// certified
	double s2 = eta * ((20827.0 / 840.0 + 123.0 * pow(M_PI, 2.0) / 64.0) - eta * eta) * pow((v * m / r), 2.0);					// certified
	double s3 = -eta * ((22717.0 / 168.0 + 615.0 * pow(M_PI, 2.0) / 64.0) + (11.0 / 8.0) * eta - 7.0 * pow(eta, 2.0)) * pow(rv * m / r, 2.0); // certified
	double s4 = -(1.0 / 4.0) * eta * (11.0 - 49.0 * eta + 52.0 * pow(eta, 2.0)) * pow(v, 6.0);									// certified
	double s5 = (35.0 / 16.0) * eta * (1.0 - 5.0 * eta + 5.0 * pow(eta, 2.0)) * pow(rv, 6.0);									// certified
	double s6 = -(1.0 / 4.0) * eta * (75.0 + 32.0 * eta - 40.0 * pow(eta, 2.0)) * (m / r) * pow(v, 4.0);						// certified
	double s7 = -(1.0 / 2.0) * eta * (158.0 - 69.0 * eta - 60.0 * pow(eta, 2.0)) * (m / r) * pow(rv, 4.0);						// certified
	double s8 = eta * (121.0 - 16.0 * eta - 20.0 * pow(eta, 2.0)) * (m / r) * pow(v * rv, 2.0);									// certified
	double s9 = (3.0 / 8.0) * eta * (20.0 - 79.0 * eta + 60.0 * pow(eta, 2.0)) * pow(v * v * rv, 2.0);							// certified
	double s10 = -(15.0 / 8.0) * eta * (4.0 - 18.0 * eta + 17.0 * pow(eta, 2.0)) * pow(v * rv * rv, 2.0);						// certified

	return s1 + s2 + s3 + s4 + s5 + s6 + s7 + s8 + s9 + s10;
}

double B3PN(double r, double rv, double ov, double m, double eta)
{
	double v = speed(r, rv, ov);
	double s1 = eta * ((5849.0 / 840.0 + 123.0 * pow(M_PI, 2.0) / 32.0) - 25.0 * eta - 8.0 * pow(eta, 2.0)) * pow((m / r), 2.0);// certified
	double s2 = (1.0 / 8.0) * eta * (65.0 - 152.0 * eta - 48.0 * pow(eta, 2.0)) * pow(v, 4.0);									// certified
	double s3 = (15.0 / 8.0) * eta * (3.0 - 8.0 * eta - 2.0 * pow(eta, 2.0)) * pow(rv, 4.0);									// certified
	double s4 = eta * (15.0 + 27.0 * eta + 10.0 * pow(eta, 2.0)) * (m / r) * pow(v, 2.0);										// certified
	double s5 = -(1.0 / 6.0) * eta * (329.0 + 177.0 * eta + 108.0 * pow(eta, 2.0)) * (m / r) * pow(rv, 2.0); 					// certified
	double s6 = -(3.0 / 4.0) * eta * (16.0 - 37.0 * eta - 16.0 * pow(eta, 2.0)) * pow(v * rv, 2.0); 							// certified

	return s1 + s2 + s3 + s4 + s5 + s6;
}


double A25PN(double r, double rv, double ov, double m)
{
	double v = speed(r, rv, ov);
	double s = 18 * pow(v, 2.0) + (2.0 / 3.0) * (m / r) - 25 * pow(rv, 2.0);
	return s;
}

double B25PN(double r, double rv, double ov, double m)
{
	double v = speed(r, rv, ov);
	double s = 6 * pow(v, 2.0) - 2 * (m / r) - 15 * pow(rv, 2.0);
	return s;
}

/*
------------------------------------------------------------------------------------------

The big A index symbol (Ai(a1, a2, a3))

------------------------------------------------------------------------------------------
*/

// The following code make use of the integration method from GSL
struct parameters{double b1; double b2;};

/* The following function f_Ai return the integrand for calculation of the big A index symbol.
 * The formula can be found from section 17 of S.Chandersakar "Ellipsoidal figure of Equilibrium"
 */
double f_A1(double x, void * params)
{	
	struct parameters alpha = * ((struct parameters *) params); // The parameters of integrand is stored in variable alpha
	double y = x / (1 - x); // x is the variable of integration, integrating from x = 0 to x = 1
	double b1, b2;
	b1 = alpha.b1; // Parameter in integrand, which is a function of axes g(a1, a2, a3)
	b2 = alpha.b2; // Parameter in integrand, which is a function of axes g(a1, a2, a3)

	double f2 = (1 + y) * (pow(b1, 2) + y) * (pow(b2, 2) + y);
	double f1 = pow(f2, -1.0 / 2.0);
	double f0 = (1 + y) * pow((1 - x), 2);
	double f = (f1 / f0) * b1 * b2; // The integrand
	return f;
}

double f_A2(double x, void * params)
{
	struct parameters alpha = * ((struct parameters *) params);
	double y = x / (1 - x);
	double b1, b2;
	b1 = alpha.b1;
	b2 = alpha.b2;

	double f2 = (1 + y) * (pow(b1, 2) + y) * (pow(b2, 2) + y);
	double f1 = pow(f2, -1.0 / 2.0);
	double f0 = (pow(b1, 2) + y) * pow((1 - x), 2);
	double f = (f1 / f0) * b1 * b2;
	return f;
}

double f_A3(double x, void * params)
{
	struct parameters alpha = * ((struct parameters *) params);
	double y = x / (1 - x);
	double b1, b2;
	b1 = alpha.b1;
	b2 = alpha.b2;

	double f2 = (1 + y) * (pow(b1, 2) + y) * (pow(b2, 2) + y);
	double f1 = pow(f2, -1.0 / 2.0);
	double f0 = (pow(b2, 2) + y) * pow((1 - x), 2);
	double f = (f1 / f0) * b1 * b2;
	return f;
}

/* The quad function calculate the numerical value of big A index symbol given the length of the ellipsoid axes. 
 * When given the integrand, it evaluate the corresponding value for A1, A2, and A3
 */
double quad(double (* func_ptr)(double, void *), double a1, double a2, double a3)
{
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);

	double result, error;
	double b1 = a2 / a1;
	double b2 = a3 / a1;
	struct parameters alpha = {b1, b2};

	gsl_function F;
	F.function = func_ptr; // Get the integrand and stored it in F
	F.params = &alpha; // Get the parameters of integrand and stored it in F

	gsl_integration_qags (&F, 0, 1, 0, 1e-10, 10000, w, &result, &error); // Integrate integrand stored in F from x = 0 to x = 1 with error 1e-10

	gsl_integration_workspace_free (w);

	return result;
}

// Simiply call the quad function to calculate A1, A2, A3 by inserting suitable integrand
double A1(double a1, double a2, double a3)
{	
	double ans = quad(f_A1, a1, a2, a3);
	return ans;
}

double A2(double a1, double a2, double a3)
{	
	double ans = quad(f_A2, a1, a2, a3);
	return ans;
}

double A3(double a1, double a2, double a3)
{	
	double ans = quad(f_A3, a1, a2, a3);
	return ans;
}

/*
------------------------------------------------------------------------------------------

Function to simplify the evolution equation (the main ODEs)

------------------------------------------------------------------------------------------
*/

/* The variable name used correspond to the Greek letter used in the paper is
 * spin = Omega
 * vort = Lambda
 * kap = kappa_n
 * orbo = theta
 * phi = phi (rotational phase of the ellipsoid)
 * eta = eta (appeared in PN expansion)
 */

// Eq.(10) and (11): a1 * (Omega^2 + Lambda^2) - 2 * a2 * Omega * Lambda
double a12_T(double a1, double a2, double spin, double vort)
{
	double s1 = a1 * (pow(spin, 2.0) + pow(vort, 2.0));
	double s2 = -2.0 * a2 * spin * vort;
	return s1 + s2;
}

// Eq.(10), (11), (12): the self gravitating terms -5/2kappa_n (c_n M + b_n m)
double ns_w(double kap, double c_n, double b_n, double M, double m, double a1, double a2, double a3)
{
	double mass = -(5.0 / (2.0 * kap)) * (c_n * M + b_n * m);
	double axes = 1.0 / (a1 * a2 * a3);
	return mass * axes;
}

// Eq.(10), (11), (12): the internal energy terms (5k1/nkappa_N)(P_c/rho_c)
double ns_u(double kap, double n, double M, double R_0, double a1, double a2, double a3)
{
	double mass = (5.0 * M) / (kap * (5.0 - n));
	double R = pow(a1 * a2 * a3, 1.0 / 3.0);
	double radius = (1.0 / R_0) * pow(R_0 / R, 3.0 / n);
	return mass * radius;
}

// Eq.(13), (14): (a2/a1 - a1/a2)^-1
double spin_pre(double a1, double a2)
{
	double s = a1 * a2 / (a2 * a2 - a1 * a1);
	return s;
}

// Eq.(13), (14): 2(Omega/a2 + Lambda/a1)(first time derivative of a1) - 2(Omega/a1 + Lambda/a2)(first time derivative of a2)
double spinning(double a1, double a2, double a1v, double a2v, double spin, double vort)
{
	double s1 = 2.0 * (spin / a2 + vort / a1) * a1v;
	double s2 = 2.0 * (spin / a1 + vort / a2) * a2v;
	return s1 - s2;
}

// Eq.(13), (14): Tidal interaction terms 3(M' + m')/r^3 * axes * sin (2 alpha) 
double spin_tide(double M, double m, double r, double a1, double a2, double orbo, double phi, char str[])
{
	double mass = (3.0 * (M + m)) / (2.0 * pow(r, 3.0));
	double axes = a1 / a2 + a2 / a1;
	double angle = sin(2.0 * (orbo - phi));

	if (strcmp(str, "spin") == 0)
		return -mass * axes * angle;
	else if (strcmp(str, "vort") == 0)
		return -2.0 * mass * angle;
	else{
		printf("Invalid Input! Check last parameter char");
		return 0;
	}
		
}

// Eq. (15), (16): Tidal interaction terms on orbital evolution (quadrupole order)
double tide_orb(double kap, double M, double m, double a1, double a2, double a3, double orbo, double phi, char str[])
{
	double mass = kap * M / (M + m);
	double angle_r = cos(2.0 * (orbo - phi));
	double axes_r = (3.0 * angle_r + 1.0) * pow(a1, 2.0) - (3.0 * angle_r - 1.0) * pow(a2, 2.0) - 2 * pow(a3, 2.0);
	double angle_o = sin(2.0 * (orbo - phi));
	double axes_o = (pow(a1, 2.0) - pow(a2, 2.0)) * angle_o;

	if (strcmp(str, "r") == 0)
		return mass * axes_r;
	else if (strcmp(str, "o") == 0)
		return mass * axes_o;
	else{
		printf("Invalid Input! Check last parameter char");
		return 0;
	}
}

// Eq.(15), (16): Conservative post-Newtonian terms
double pnc(double r, double rv, double ov, double Mt, double eta, char str[], struct Pa pa)
{
	double ASv = AS(r, rv, ov, Mt);
	double A1 = A1PN(r, rv, ov, Mt, eta);
	double A2 = A2PN(r, rv, ov, Mt, eta);
	double A3 = A3PN(r, rv, ov, Mt, eta);

	double BSv = BS(r, Mt);
	double B1 = B1PN(eta);
	double B2 = B2PN(r, rv, ov, Mt, eta);
	double B3 = B3PN(r, rv, ov, Mt, eta);

	double A = 1.0 - ASv + pa.pn.PN1 * A1 + pa.pn.PN2 * A2 + pa.pn.PN3 * A3;
	double B = -BSv + pa.pn.PN1 * B1 + pa.pn.PN2 * B2 + pa.pn.PN3 * B3;

	if (strcmp(str, "r") == 0) // The Conservative PN terms in equation of acceleration of r
		return A + rv * rv * B;
	else if (strcmp(str, "o") == 0) // The Conservative PN terms in equation of acceleration of theta (orbital phase)
		return rv * ov * B;
	else{
		printf("Invalid Input! Check parameter char");
		return 0;		
	}
}

// Eq.(15), (16): Dissipative post-Newtonian terms
double pnd(double r, double rv, double ov, double Mt, double eta, char str[], struct Pa pa)
{
	double A52 = A25PN(r, rv, ov, Mt);
	double B52 = B25PN(r, rv, ov, Mt);

	double r_eqn = (8.0 / 5.0) * eta * (Mt / pow(r, 2.0)) * (Mt / r) * rv * pa.pn.PN25 * (A52 - B52);
	double o_eqn = -(8.0 / 5.0) * eta * (Mt / pow(r, 2.0)) * (Mt / r) * ov * pa.pn.PN25 * B52;

	if (strcmp(str, "r") == 0) // The Conservative PN terms in equation of acceleration of r
		return r_eqn;
	else if (strcmp(str, "o") == 0) // The Conservative PN terms in equation of acceleration of theta (orbital phase)
		return o_eqn;
	else{
		printf("Invalid Input! Check parameter char");
		return 0;
		}	
}

// Dark repulsion on acceleration of r
double dark_eqn(double r, struct Pa pa)
{
	int on = pa.dark.on;
	double al = pa.dark.alpha; // Dark force coupling constant
	double mV = pa.dark.mV; // Dark force mediator range

	double s1 = on * al * exp(-1.0 * mV * r);
	double s2 = 1.0 + mV * r;

	return 1.0 - s1 * s2;
}
/*
------------------------------------------------------------------------------------------

Function to simplify the equilibrium equation (the 6 equation to solve for axes)

------------------------------------------------------------------------------------------
*/

// Eq.(32), (33): I_11, I_22, I_33
double I_num(double a, double kap, double M)
{
	double s = (1.0 / 5.0) * kap * M * a * a;
	return s;
}

// Eq.(32), (33)
double del(double a1, double a2, double a3, double kap, double M, double m, double r)
{
	double s1 = (3.0 / 2.0) * (2 * I_num(a1, kap, M) - I_num(a2, kap, M) - I_num(a3, kap, M));
	double s2 = (M + m) * r * r;

	return s1 / s2;
}

// Eq.(35)
double g_t(double a1, double a2, double a3, double kap, double M, double r)
{
	double I11 = I_num(a1, kap, M);
	double I22 = I_num(a2, kap, M);
	double I33 = I_num(a3, kap, M);
	double R = pow(a1 * a2 * a3, 1.0 / 3.0);

	return R * (2.0 * I11 - I22 - I33) / pow(r, 3.0);
}

// Eq.(2), (34)
double f(double a1, double a2, double a3)
{
	double I = A1(a1, a2, a3) * a1 * a1 + A2(a1, a2, a3) * a2 * a2 + A3(a1, a2, a3) * a3 * a3;
	double s2 = 2.0 * pow(a1 * a2 * a3, (2.0 / 3.0));

	return I / s2;
}

// Eq.(34)
double T_W_shared(double a1_A, double a2_A, double a3_A,
				  double a1_B, double a2_B, double a3_B,
				  double kap_A, double kap_B, double M_A, 
				  double m_A, double M_B, double m_B, double r)
{
	double mass = (M_A + m_A + M_B + m_B) / pow(r, 3.0);
	double del_A = del(a1_A, a2_A, a3_A, kap_A, M_A, m_A, r);
	double del_B = del(a1_B, a2_B, a3_B, kap_B, M_B, m_B, r);

	return mass * (1.0 + del_A + del_B);
}

// Eq.(34)
double T_W_single(double a1, double a2, double a3,
				  double c_n, double b_n, double kap,
				  double f_R, double M, double m)
{
	double index = kap / (10.0 * (c_n * M + b_n * m));
	double R = pow(a1 * a2 * a3, 1.0 / 3.0);
	double axes1 = (a1 * a1 + a2 * a2);
	double axes2 = (1.0 + pow(f_R * a1 * a2 / axes1, 2.0));
	double axes3 = 4.0 * f_R * a1 * a1 * a2 * a2 / axes1;
	double f_var = f(a1, a2, a3);

	return index * R * (axes1 * axes2 + axes3) / f_var;
}
/*
------------------------------------------------------------------------------------------

The six equilibrium equations for DR(DM) model

------------------------------------------------------------------------------------------
*/
// Eq.(40)
double eqn_1(double a1_A, double a1_B, double a2_A,
			 double a2_B, double a3_A, double a3_B,
			 void * params)
{
	return a1_A + a1_B - 1;
}

// Eq.(40): R/R' - ((a1 a2 a3)/(a1' a2' a3') ^(1/3)) = 0
double eqn_2(double a1_A, double a1_B, double a2_A,
			 double a2_B, double a3_A, double a3_B,
			 void * params)
{
	struct Pa p = * ((struct Pa *) params);

	double R_A_s1 = (1.0 + p.b_n_A * p.m_A / (p.c_n_A * p.M_A));
	double T_W_A_s = T_W_single(a1_A, a2_A, a3_A, p.c_n_A, p.b_n_A, p.kap_A, p.f_R_A, p.M_A, p.m_A);
	double T_W_s = T_W_shared(a1_A, a2_A, a3_A, a1_B, a2_B, a3_B, p.kap_A, p.kap_B, p.M_A, p.m_A, p.M_B, p.m_B, p.r);
	double T_W_A = T_W_A_s * T_W_s;
	double f_var_A = f(a1_A, a2_A, a3_A);
	double R_A_s2 = (p.M_B + p.m_B) / (p.M_A * (p.c_n_A * p.M_A + p.b_n_A * p.m_A));
	double g_t_A = g_t(a1_A, a2_A, a3_A, p.kap_A, p.M_A, p.r);
	double R_A_s3 = R_A_s1 * ((1.0 - 2.0 * T_W_A) * f_var_A - R_A_s2 * g_t_A);

	double R_B_s1 = (1.0 + p.b_n_B * p.m_B / (p.c_n_B * p.M_B));
	double T_W_B_s = T_W_single(a1_B, a2_B, a3_B, p.c_n_B, p.b_n_B, p.kap_B, p.f_R_B, p.M_B, p.m_B);
	double T_W_B = T_W_B_s * T_W_s;
	double f_var_B = f(a1_B, a2_B, a3_B);
	double R_B_s2 = (p.M_A + p.m_A) / (p.M_B * (p.c_n_B * p.M_B + p.b_n_B * p.m_B));
	double g_t_B = g_t(a1_B, a2_B, a3_B, p.kap_B, p.M_B, p.r);
	double R_B_s3 = R_B_s1 * ((1.0 - 2.0 * T_W_B) * f_var_B - R_B_s2 * g_t_B);

	double R_A = pow(R_A_s3, -p.n_A / (3.0 - p.n_A));
	double R_B = pow(R_B_s3, -p.n_B / (3.0 - p.n_B));
	double pre = pow(p.M_A / p.M_B, (1.0 - p.n_A) / (3.0 - p.n_A));

	double RHS_s = a1_A * a2_A * a3_A / (a1_B * a2_B * a3_B);
	double RHS = pow(RHS_s, 1.0 / 3.0);

	return pre * R_A / R_B - RHS;
}

// Eq.(36)
double eqn_3(double a1_A, double a1_B, double a2_A,
			 double a2_B, double a3_A, double a3_B,
			 void * params)
{
	struct Pa p = * ((struct Pa *) params);

	double h_n = 4.0 * p.kap_A / (5.0 * (p.c_n_A * p.M_A + p.b_n_A * p.m_A));
	double mu_R = (p.M_B + p.m_B) / pow(p.r, 3.0);
	double R3 = a1_A * a2_A * a3_A;

	double mass_p = 1.0 + (p.M_A + p.m_A) / (p.M_B + p.m_B);
	double Q12_s1 = p.f_R_A * p.f_R_A * a1_A * a1_A;
	double Q12_s2 = pow(a1_A * a2_A / (a1_A * a1_A + a2_A * a2_A), 2.0);
	double del_A = del(a1_A, a2_A, a3_A, p.kap_A, p.M_A, p.m_A, p.r);
	double del_B = del(a1_B, a2_B, a3_B, p.kap_B, p.M_B, p.m_B, p.r);
	double Q12 = mass_p * Q12_s1 * Q12_s2 * (1.0 + del_A + del_B);

	double Q2_s = p.f_R_A * a2_A * a2_A / (a1_A * a1_A + a2_A * a2_A);
	double Q2 = mass_p * Q2_s * (1.0 + del_A + del_B);

	double eqn_s1 = h_n * mu_R * R3;
	double eqn_s2 = 2.0 + mass_p * (1.0 + del_A + del_B);

	double LHS =  eqn_s1 * (Q12 + (eqn_s2 + 2.0 * Q2) * a1_A * a1_A + a3_A * a3_A);
	double RHS = 2.0 * (A1(a1_A, a2_A, a3_A) * a1_A * a1_A - A3(a1_A, a2_A, a3_A) * a3_A * a3_A);

	return LHS - RHS;
}

// Eq.(37)
double eqn_4(double a1_A, double a1_B, double a2_A,
			 double a2_B, double a3_A, double a3_B,
			 void * params)
{
	struct Pa p = * ((struct Pa *) params);

	double h_n = 4.0 * p.kap_A / (5.0 * (p.c_n_A * p.M_A + p.b_n_A * p.m_A));
	double mu_R = (p.M_B + p.m_B) / pow(p.r, 3.0);
	double R3 = a1_A * a2_A * a3_A;

	double mass_p = (1.0 + (p.M_A + p.m_A) / (p.M_B + p.m_B));
	double Q21_s1 = p.f_R_A * p.f_R_A * a2_A * a2_A;
	double Q21_s2 = pow(a1_A * a2_A / (a1_A * a1_A + a2_A * a2_A), 2.0);
	double del_A = del(a1_A, a2_A, a3_A, p.kap_A, p.M_A, p.m_A, p.r);
	double del_B = del(a1_B, a2_B, a3_B, p.kap_B, p.M_B, p.m_B, p.r);
	double Q21 = mass_p * Q21_s1 * Q21_s2 * (1.0 + del_A + del_B);

	double Q1_s = -p.f_R_A * a1_A * a1_A / (a1_A * a1_A + a2_A * a2_A);
	double Q1 = mass_p * Q1_s * (1.0 + del_A + del_B);

	double eqn_s1 = h_n * mu_R * R3;
	double eqn_s2 = mass_p * (1.0 + del_A + del_B) - 1.0;

	double LHS =  eqn_s1 * (Q21 + (eqn_s2 - 2.0 * Q1) * a2_A * a2_A + a3_A * a3_A);
	double RHS = 2.0 * (A2(a1_A, a2_A, a3_A) * a2_A * a2_A - A3(a1_A, a2_A, a3_A) * a3_A * a3_A);

	return LHS - RHS;
}

// Eq.(36) for companion
double eqn_5(double a1_A, double a1_B, double a2_A,
			 double a2_B, double a3_A, double a3_B,
			 void * params)
{
	struct Pa p = * ((struct Pa *) params);

	double h_n = 4.0 * p.kap_B / (5.0 * (p.c_n_B * p.M_B + p.b_n_B * p.m_B));
	double mu_R = (p.M_A + p.m_A) / pow(p.r, 3.0);
	double R3 = a1_B * a2_B * a3_B;

	double mass_p = (1.0 + (p.M_B + p.m_B) / (p.M_A + p.m_A));
	double Q12_s1 = p.f_R_B * p.f_R_B * a1_B * a1_B;
	double Q12_s2 = pow(a1_B * a2_B / (a1_B * a1_B + a2_B * a2_B), 2.0);
	double del_A = del(a1_A, a2_A, a3_A, p.kap_A, p.M_A, p.m_A, p.r);
	double del_B = del(a1_B, a2_B, a3_B, p.kap_B, p.M_B, p.m_B, p.r);
	double Q12 = mass_p * Q12_s1 * Q12_s2 * (1.0 + del_A + del_B);

	double Q2_s = p.f_R_B * a2_B * a2_B / (a1_B * a1_B + a2_B * a2_B);
	double Q2 = mass_p * Q2_s * (1.0 + del_A + del_B);

	double eqn_s1 = h_n * mu_R * R3;
	double eqn_s2 = 2.0 + mass_p * (1.0 + del_A + del_B);

	double LHS =  eqn_s1 * (Q12 + (eqn_s2 + 2.0 * Q2) * a1_B * a1_B + a3_B * a3_B);
	double RHS = 2.0 * (A1(a1_B, a2_B, a3_B) * a1_B * a1_B - A3(a1_B, a2_B, a3_B) * a3_B * a3_B);

	return LHS - RHS;
}

// Eq.(37) for companion
double eqn_6(double a1_A, double a1_B, double a2_A,
			 double a2_B, double a3_A, double a3_B,
			 void * params)
{
	struct Pa p = * ((struct Pa *) params);

	double h_n = 4.0 * p.kap_B / (5.0 * (p.c_n_B * p.M_B + p.b_n_B * p.m_B));
	double mu_R = (p.M_A + p.m_A) / pow(p.r, 3.0);
	double R3 = a1_B * a2_B * a3_B;

	double mass_p = (1.0 + (p.M_B + p.m_B) / (p.M_A + p.m_A));
	double Q21_s1 = p.f_R_B * p.f_R_B * a2_B * a2_B;
	double Q21_s2 = pow(a1_B * a2_B / (a1_B * a1_B + a2_B * a2_B), 2.0);
	double del_A = del(a1_A, a2_A, a3_A, p.kap_A, p.M_A, p.m_A, p.r);
	double del_B = del(a1_B, a2_B, a3_B, p.kap_B, p.M_B, p.m_B, p.r);
	double Q21 = mass_p * Q21_s1 * Q21_s2 * (1.0 + del_A + del_B);

	double Q1_s = -p.f_R_B * a1_B * a1_B / (a1_B * a1_B + a2_B * a2_B);
	double Q1 = mass_p * Q1_s * (1.0 + del_A + del_B);

	double eqn_s1 = h_n * mu_R * R3;
	double eqn_s2 = mass_p * (1.0 + del_A + del_B) - 1.0;

	double LHS =  eqn_s1 * (Q21 + (eqn_s2 - 2.0 * Q1) * a2_B * a2_B + a3_B * a3_B);
	double RHS = 2.0 * (A2(a1_B, a2_B, a3_B) * a2_B * a2_B - A3(a1_B, a2_B, a3_B) * a3_B * a3_B);

	return LHS - RHS;
}
/*
------------------------------------------------------------------------------------------

Function to evaluate the equilibrium radius and orbital frequency using solved axes length

------------------------------------------------------------------------------------------
*/

// Eq.(34)
double eqn_R(double a1_A, double a2_A, double a3_A,
			 double a1_B, double a2_B, double a3_B,
			 void * params, char a)
{
	struct Pa p = * ((struct Pa *) params);

	double R_A_s1 = (1.0 + p.b_n_A * p.m_A / (p.c_n_A * p.M_A));
	double R_B_s1 = (1.0 + p.b_n_B * p.m_B / (p.c_n_B * p.M_B));
	double T_W_A_s = T_W_single(a1_A, a2_A, a3_A, p.c_n_A, p.b_n_A, p.kap_A, p.f_R_A, p.M_A, p.m_A);
	double T_W_B_s = T_W_single(a1_B, a2_B, a3_B, p.c_n_B, p.b_n_B, p.kap_B, p.f_R_B, p.M_B, p.m_B);
	double T_W_s = T_W_shared(a1_A, a2_A, a3_A, a1_B, a2_B, a3_B, p.kap_A, p.kap_B, p.M_A, p.m_A, p.M_B, p.m_B, p.r);
	double T_W_A = T_W_A_s * T_W_s;
	double T_W_B = T_W_B_s * T_W_s;
	double f_var_A = f(a1_A, a2_A, a3_A);
	double f_var_B = f(a1_B, a2_B, a3_B);
	double R_A_s2 = (p.M_B + p.m_B) / (p.M_A * (p.c_n_A * p.M_A + p.b_n_A * p.m_A));
	double R_B_s2 = (p.M_A + p.m_A) / (p.M_B * (p.c_n_B * p.M_B + p.b_n_B * p.m_B));
	double g_t_A = g_t(a1_A, a2_A, a3_A, p.kap_A, p.M_A, p.r);
	double g_t_B = g_t(a1_B, a2_B, a3_B, p.kap_B, p.M_B, p.r);
	double R_A_s3 = R_A_s1 * ((1.0 - 2.0 * T_W_A) * f_var_A - R_A_s2 * g_t_A);
	double R_B_s3 = R_B_s1 * ((1.0 - 2.0 * T_W_B) * f_var_B - R_A_s2 * g_t_B);

	if (a == 'A')
		return pow(R_A_s3, -p.n_A / (3.0 - p.n_A));
	else if (a == 'B')
		return pow(R_B_s3, -p.n_B / (3.0 - p.n_B));
	else{
		printf("Invalide last char input!");
		return 0;
	}
	
}

// Eq.(31): Modified Kepler's Law
double kepler(double a1_A, double a2_A, double a3_A,
			  double a1_B, double a2_B, double a3_B,
			  void * params)
{
	struct Pa p = * ((struct Pa *) params);
	double r = p.r * (a1_A + a1_B);

	double mass = (p.M_A + p.m_A + p.M_B + p.m_B) / pow(r, 3.0);
	double del_A = del(a1_A, a2_A, a3_A, p.kap_A, p.M_A, p.m_A, r);
	double del_B = del(a1_B, a2_B, a3_B, p.kap_B, p.M_B, p.m_B, r);

	return pow(mass * (1 + del_A + del_B), 1.0 / 2.0);
}
/*
------------------------------------------------------------------------------------------

Solving For Equilibrium configuration section (From the above 6 equlibrium equation)

------------------------------------------------------------------------------------------
*/

// Use the GSL library for the multidimensional newton method
int
rosenbrock_f (const gsl_vector * x, void * params,
              gsl_vector * f)
{

  const double x0 = gsl_vector_get (x, 0);
  const double x1 = gsl_vector_get (x, 1);
  const double x2 = gsl_vector_get (x, 2);
  const double x3 = gsl_vector_get (x, 3);
  const double x4 = gsl_vector_get (x, 4);
  const double x5 = gsl_vector_get (x, 5);

  // Call the 6 equilibrium equations
  const double y0 = eqn_1(x0, x1, x2, x3, x4, x5, params);
  const double y1 = eqn_2(x0, x1, x2, x3, x4, x5, params);
  const double y2 = eqn_3(x0, x1, x2, x3, x4, x5, params);
  const double y3 = eqn_4(x0, x1, x2, x3, x4, x5, params);
  const double y4 = eqn_5(x0, x1, x2, x3, x4, x5, params);
  const double y5 = eqn_6(x0, x1, x2, x3, x4, x5, params);


  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);
  gsl_vector_set (f, 2, y2);
  gsl_vector_set (f, 3, y3);
  gsl_vector_set (f, 4, y4);
  gsl_vector_set (f, 5, y5);

  return GSL_SUCCESS;
}

/* The print_state function is not used */
int 
print_state (size_t iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %3lu x = % .6e % .6e % .6e % .6e % .6e % .6e"
          "f(x) = % .3e % .3e % .3e % .3e % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2),
          gsl_vector_get (s->x, 3),
          gsl_vector_get (s->x, 4),
          gsl_vector_get (s->x, 5),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1),
          gsl_vector_get (s->f, 2),
          gsl_vector_get (s->f, 3),
          gsl_vector_get (s->f, 4),
          gsl_vector_get (s->f, 5));
}

/* This function is called by func_r function below */
struct Eqm newton_6d(struct Pa p)
{   
    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int status;
    size_t i, iter = 0;

    const size_t n = 6;

    gsl_multiroot_function f = {&rosenbrock_f, n, &p}; // Insert the 6 equilibrium equation and parameters

    double x_init[6] = {0.95, 0.94, 0.93, 0.92, 0.91, 0.90}; // Initial Guess of the 6 variable a_i, a'_i
    gsl_vector *x = gsl_vector_alloc (n);

    gsl_vector_set (x, 0, x_init[0]);
    gsl_vector_set (x, 1, x_init[1]);
    gsl_vector_set (x, 2, x_init[2]);
    gsl_vector_set (x, 3, x_init[3]);
    gsl_vector_set (x, 4, x_init[4]);
    gsl_vector_set (x, 5, x_init[5]);

    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc (T, 6);
    gsl_multiroot_fsolver_set (s, &f, x); // This solve the 6 equilibrium equations, stored in variable s

    //print_state (iter, s);

    do
    {
        iter++;
        status = gsl_multiroot_fsolver_iterate (s);

        //print_state (iter, s);

        if (status)   /* check if solver is stuck */
        break;

        status =
        gsl_multiroot_test_residual (s->f, 1e-13);
    }
    while (status == GSL_CONTINUE && iter < 1000);

    //printf ("status = %s\n", gsl_strerror (status));
    //print_state(iter, s);

    struct Eqm retur; // Declare a retur variable to store the equilirbium value of a_i, orbital frequncy and separation
					  // The struct Eqm is declared in the setting.h file

	// Solved ai_A/B value (normalized by a1 + a1')
    double a1_A = gsl_vector_get (s->x, 0);
    double a2_A = gsl_vector_get (s->x, 2);
    double a3_A = gsl_vector_get (s->x, 4);
    double a1_B = gsl_vector_get (s->x, 1);
    double a2_B = gsl_vector_get (s->x, 3);
    double a3_B = gsl_vector_get (s->x, 5);

	// Physical equililbrium radius R = (a1 * a2 * a3)**(1/3) for star A and B
    double R_A = p.R_0_A * eqn_R(a1_A, a2_A, a3_A,
                             a1_B, a2_B, a3_B, &p, 'A');
    double R_B = p.R_0_A * pow(p.M_A / p.M_B, (1.0 - p.n_A) / (3.0 - p.n_A)) * eqn_R(a1_A, a2_A, a3_A,
                                                                                     a1_B, a2_B, a3_B, &p, 'B');
	// Ratio of the axes
    double a21_A = a2_A / a1_A;
    double a31_A = a3_A / a1_A;
    double a21_B = a2_B / a1_B;
    double a31_B = a3_B / a1_B;

	// Physical length of the axes in unit of (m)
    retur.a1_A = R_A * pow(1.0 / (a21_A * a31_A), 1.0 / 3.0);
    retur.a2_A = retur.a1_A * a21_A;
    retur.a3_A = retur.a1_A * a31_A;

    retur.a1_B = R_B * pow(1.0 / (a21_B * a31_B), 1.0 / 3.0);
    retur.a2_B = retur.a1_B * a21_B;
    retur.a3_B = retur.a1_B * a31_B;

	// Physical value of equilibrium orbital frequency without PN terms in unit of (rad/m)
    retur.orb = kepler(retur.a1_A, retur.a2_A, retur.a3_A,
                     retur.a1_B, retur.a2_B, retur.a3_B, &p);

    retur.r_phy = p.r * (retur.a1_A + retur.a1_B); // Physical value of orbital separation r in unit of (m)

    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (x);

    return retur; 
}

/* The func_r provide the equation r/R_0(solved from the above newton_6d) - r/R_0(in binary_conf.c) = 0 
 * which we seek to solve for r/R_0(solved from the above newton_6d) by secant method 
 */
double func_r(double x, void * params)
{
    struct Pa p = * ((struct Pa *) params);
    double aim = p.r; // Set the variable aim = p.r chosen in binary_conf.c
    p.r = x; // Assign a random guess orbital separation x = r/(a1+a1') to p.r
	
	/* Obtain the equilibrium configuration from the random guessed orbital separation x
	 * Compare the solved phyical orbital separation r to aim = p.r chosen in binary_conf.c
	 */
	double r_phy = newton_6d(p).r_phy / p.R_0_A; 
    return r_phy - aim; 
}

/* The func_orb is used to obtain the intial orbital frequency accroding to different PN order */
double func_orb(double x, void * params, void * params2)
{
    struct Pa p = * ((struct Pa *) params); // Binary parameters
    struct Eqm eq = * ((struct Eqm *) params2); // Binary equilibrium value
	
	// Variables to simplify the radial acceleration in orbital motion 
    double Mt = p.M_A + p.m_A + p.M_B + p.m_B;
    double eta = (p.M_A + p.m_A) * (p.M_B + p.m_B) / (Mt * Mt);
    double r = eq.r_phy;
    double r_tide_A = tide_orb(p.kap_A, p.M_A, p.m_A, eq.a1_A, eq.a2_A, eq.a3_A, 0.0, 0.0, "r");
    double r_tide_B = tide_orb(p.kap_B, p.M_B, p.m_B, eq.a1_B, eq.a2_B, eq.a3_B, 0.0, 0.0, "r");
    double ra_c = pnc(r, 0.0, x, Mt, eta, "r", p);
    double dark = dark_eqn(r, p);

	// The equation d^2 r/ dt^2 = 0 (theta = orbital phase)
    double result = -(3.0 / 20.0) * (Mt / (r * r * r * r)) * (r_tide_A + r_tide_B)
                    + r * x * x - (Mt / (r * r)) * dark + p.pn.PNc * (Mt / (r * r)) * ra_c;
    return result;
}

/* Obtain the suitable x = r/(a1 + a1') used in solving the equlibrium axes length of the ellipsoid*/
double secant_r(double func(double, void *), void * params)
{
    struct Pa p = * ((struct Pa *) params);
    double ans, x1, x0, err;
    x1 = 6.0; // Initial guess 1
    x0 = 5.0; // Initial guess 2
    err = fabs(x1 - x0); // Initial error
    int i = 0;
    while (err > 1e-15)
    {
        double f_x1 = func(x1, &p);
        double f_x0 = func(x0, &p);
        double new_x = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0);

        x0 = x1;
        x1 = new_x;
        err = fabs(x1 - x0); // Update error until it is within the desired range (1e-15)
        i++;
    }

    ans = x1;
    //printf("total iteration = %d suitable r scaled is % .5lf\n", i, ans);
    return ans;

}

/* Obtain the initial orbital frequency by secant method */
double secant_orb(double func(double, void *, void *), void * params, void * params2)
{
    struct Pa p = * ((struct Pa *) params); // The binary configuration set in binary_conf.c
    struct Eqm eq = * ((struct Eqm *) params2); // The solved equilibrium axes length
    double ans, x1, x0, err;
    x1 = 2.0e-6;
    x0 = 3.0e-6;
    err = fabs(x1 - x0);
    int i = 0;
    while (err > 1e-15)
    {
        double f_x1 = func_orb(x1, &p, &eq);
        double f_x0 = func_orb(x0, &p, &eq);
        double new_x = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0);

        x0 = x1;
        x1 = new_x;
        err = fabs(x1 - x0);
        i++;
    }

    ans = x1;
    //printf("total iteration = %d orb is % .5e\n", i, ans);
    return ans;

}

/* The initial function return the equilibrium value of a1, a2, a3, orb, and r given the binary configuration set in binary_conf.c */
struct Eqm initial(struct Pa p)
{
    double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
    double sui_r = secant_r(func_r, &p); // Find suitable sui_r = r (physical) / (a1 + a1')
    p.r = sui_r; // Use the suitable sui_r to solve for equlibrium value a1, a2, a3 ,etc.
    struct Eqm a = newton_6d(p);
    a.orb = secant_orb(func_orb, &p, &a); // Use the solved equlibrium value to get initial orbital frequency accroding to the PN terms

    return a;
}

/*
------------------------------------------------------------------------------------------

The main ODEs section

------------------------------------------------------------------------------------------
*/

/* Return the first time derivative as pointer for implementation of rk45 method */
double *eqn(double s[], struct Pa pa)
{
	double * p;
	int dim = 22;
	p = (double *) malloc((dim + 1) * sizeof(double));
	double Mt = pa.M_A + pa.m_A + pa.M_B + pa.m_B;
	double eta = (pa.M_A + pa.m_A) * (pa.M_B + pa.m_B) / (Mt * Mt);

	double ns_w_A = ns_w(pa.kap_A, pa.c_n_A, pa.b_n_A, pa.M_A, pa.m_A, s[1], s[5], s[9]);
	double ns_w_B = ns_w(pa.kap_B, pa.c_n_B, pa.b_n_B, pa.M_B, pa.m_B, s[3], s[7], s[11]);
	double ns_u_A = ns_u(pa.kap_A, pa.n_A, pa.M_A, pa.R_0_A, s[1], s[5], s[9]);
	double ns_u_B = ns_u(pa.kap_B, pa.n_B, pa.M_B, pa.R_0_B, s[3], s[7], s[11]);

	*p =  a12_T(s[1], s[5], s[12], s[16]) + s[1] * A1(s[1], s[5], s[9]) * ns_w_A + (1.0 / s[1]) * ns_u_A
		   + ((pa.M_B + pa.m_B) / (2.0 * pow(s[19], 3.0))) * s[1] * (3 * cos(2.0 * (s[21] - s[13])) + 1);

	*(p + 1) = s[0];

	*(p + 2) = a12_T(s[3], s[7], s[14], s[17]) + s[3] * A1(s[3], s[7], s[11]) * ns_w_B + (1.0 / s[3]) * ns_u_B
		   		+ ((pa.M_A + pa.m_A) / (2.0 * pow(s[19], 3.0))) * s[3] * (3 * cos(2.0 * (s[21] - s[15])) + 1);

	*(p + 3) = s[2];

	*(p + 4) = a12_T(s[5], s[1], s[12], s[16]) + s[5] * A2(s[1], s[5], s[9]) * ns_w_A + (1.0 / s[5]) * ns_u_A
		   		- ((pa.M_B + pa.m_B) / (2.0 * pow(s[19], 3.0))) * s[5] * (3.0 * cos(2.0 * (s[21] - s[13])) - 1.0);

	*(p + 5) = s[4];

	*(p + 6) = a12_T(s[7], s[3], s[14], s[17]) + s[7] * A2(s[3], s[7], s[11]) * ns_w_B + (1.0 / s[7]) * ns_u_B
		   		- ((pa.M_A + pa.m_A) / (2.0 * pow(s[19], 3.0))) * s[7] * (3.0 * cos(2.0 * (s[21] - s[15])) - 1.0);

	*(p + 7) = s[6];

	*(p + 8) = s[9] * A3(s[1], s[5], s[9]) * ns_w_A + (1.0 / s[9]) * ns_u_A - (pa.M_B + pa.m_B) * s[9] / pow(s[19], 3.0);

	*(p + 9) = s[8];

	*(p + 10) = s[11] * A3(s[3], s[7], s[11]) * ns_w_B + (1.0 / s[11]) * ns_u_B - (pa.M_A + pa.m_A) * s[11] / pow(s[19], 3.0);

	*(p + 11) = s[10];

	double s_p_A = spin_pre(s[1], s[5]);
	double s_p_B = spin_pre(s[3], s[7]);
	double s_t_A = spin_tide(pa.M_B, pa.m_B, s[19], s[1], s[5], s[21], s[13], "spin");
	double v_t_A = spin_tide(pa.M_B, pa.m_B, s[19], s[1], s[5], s[21], s[13], "vort");
	double s_t_B = spin_tide(pa.M_A, pa.m_A, s[19], s[3], s[7], s[21], s[15], "spin");
	double v_t_B = spin_tide(pa.M_A, pa.m_A, s[19], s[3], s[7], s[21], s[15], "vort");

	*(p + 12) = s_p_A * (spinning(s[1], s[5], s[0], s[4], s[12], s[16]) + s_t_A);

	*(p + 13) = s[12];

	*(p + 14) = s_p_B * (spinning(s[3], s[7], s[2], s[6], s[14], s[17]) + s_t_B);

	*(p + 15) = s[14];

	*(p + 16) = s_p_A * (spinning(s[5], s[1], s[0], s[4], s[12], s[16]) + v_t_A);

	*(p + 17) = s_p_B * (spinning(s[7], s[3], s[2], s[6], s[14], s[17]) + v_t_B);

	double r_tide_A = tide_orb(pa.kap_A, pa.M_A, pa.m_A, s[1], s[5], s[9], s[21], s[13], "r");
	double orb_tide_A = tide_orb(pa.kap_A, pa.M_A, pa.m_A, s[1], s[5], s[9], s[21], s[13], "o");
	double r_tide_B = tide_orb(pa.kap_B, pa.M_B, pa.m_B, s[3], s[7], s[11], s[21], s[15], "r");
	double orb_tide_B = tide_orb(pa.kap_B, pa.M_B, pa.m_B, s[3], s[7], s[11], s[21], s[15], "o");
	double ra_c = pnc(s[19], s[18], s[20], Mt, eta, "r", pa);
	double orba_c = pnc(s[19], s[18], s[20], Mt, eta, "o", pa);
	double ra_d = pnd(s[19], s[18], s[20], Mt, eta, "r", pa);
	double orba_d = pnd(s[19], s[18], s[20], Mt, eta, "o", pa);

	double dark = dark_eqn(s[19], pa);

	*(p + 18) = -(3.0 / 20.0) * (Mt / pow(s[19], 4.0)) * (r_tide_A + r_tide_B)
				 + s[19] * pow(s[20], 2.0) - (Mt / pow(s[19], 2.0)) * dark + pa.pn.PNc * (Mt / pow(s[19], 2.0)) * ra_c + ra_d;

	*(p + 19) = s[18];

	*(p + 20) = -2.0 * s[18] * s[20] / s[19] - (3.0 * Mt / (10.0 * pow(s[19], 5.0))) * (orb_tide_A + orb_tide_B)
				 + pa.pn.PNc * (Mt / pow(s[19], 2.0)) * orba_c + orba_d;

	*(p + 21) = s[20];

	return p;
}

/*
------------------------------------------------------------------------------------------

The rk45 method (The coefficients and notation are copied from wiki page of Runge–Kutta–Fehlberg method)

------------------------------------------------------------------------------------------
*/

double B[6][5] = {{-1.0, -1.0, -1.0, -1.0, -1.0}
							, {1.0 / 4.0, -1.0, -1.0, -1.0, -1.0}
							, {3.0 / 32.0, 9.0 / 32.0, -1.0, -1.0, -1.0}
							, {1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0, -1.0, -1.0}
							, {439.0 / 216.0, -8.0, 3680.0 / 513.0, -845.0 / 4104.0, -1.0}
							, {-8.0 / 27.0, 2.0, -3544.0 / 2565.0, 1859.0 / 4104.0, -11.0 / 40.0}};

double A[6] = {0.0 , 1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 1.0 / 2.0};

double C[6] = {25.0 / 216.0, 0.0, 1408.0 / 2565.0, 2197.0 / 4104.0, -1.0 / 5.0, -1};

double CH[6] = {16.0 / 135.0, 0.0, 6656.0 / 12825.0, 28561.0 / 56430.0, -9.0 / 50.0, 2.0 / 55.0};

double CT[6] = {1.0 / 360.0, 0.0, -128.0 / 4275.0, -2197.0 / 75240.0, 1.0 / 50.0, 2.0 / 55.0};
							
/* func = the system of ODE ;
 * array = previous point 
 * h = step size
 * dim = number of first order ODE
 * pa = binary configuration setting (set in binary_conf.c)
 * The *rk45 function return a pointer containing numerical value of next point and relative error of each equation
 */
double *rk45(double * func(double *, struct Pa), double array[], double h, int dim, struct Pa pa)
{
	double *p_k1, *p_k2, *p_k3, *p_k4, *p_k5, *p_k6, *p;
	double k2_arg[dim], k3_arg[dim], k4_arg[dim], k5_arg[dim], k6_arg[dim];
	int i;

	p_k1 = func(array, pa);
	for(i = 0; i < dim; i++){
		p_k1[i] = h * p_k1[i];
		k2_arg[i] = array[i] + B[1][0] * p_k1[i];
	}
	

	p_k2 = func(k2_arg, pa);
	for(i = 0; i < dim; i++){
		p_k2[i] = h * p_k2[i];
		k3_arg[i] = array[i] + B[2][0] * p_k1[i] + B[2][1] * p_k2[i];
	}

	p_k3 = func(k3_arg, pa);
	for(i = 0; i < dim; i++){
		p_k3[i] = h * p_k3[i];
		k4_arg[i] = array[i] + B[3][0] * p_k1[i] + B[3][1] * p_k2[i] + B[3][2] * p_k3[i];
	}

	p_k4 = func(k4_arg, pa);
	for(i = 0; i < dim; i++){
		p_k4[i] = h * p_k4[i];
		k5_arg[i] = array[i] + B[4][0] * p_k1[i] + B[4][1] * p_k2[i] + B[4][2] * p_k3[i] + B[4][3] * p_k4[i];
	}

	p_k5 = func(k5_arg, pa);
	for(i = 0; i < dim; i++){
		p_k5[i] = h * p_k5[i];
		k6_arg[i] = array[i] + B[5][0] * p_k1[i] + B[5][1] * p_k2[i] + B[5][2] * p_k3[i] + B[5][3] * p_k4[i] + B[5][4] * p_k5[i];
	}

	p_k6 = func(k6_arg, pa);
	for(i = 0; i < dim; i++)
		p_k6[i] = h * p_k6[i];

	p = (double *) malloc(2 * (dim + 1) * sizeof(double)); // Allocate array to store the numerical value of next point and relative error

	for(i = 0; i < dim; i++)
	{
		p[i] = array[i] + CH[0] * p_k1[i] + CH[1] * p_k2[i] + CH[2] * p_k3[i] + CH[3] * p_k4[i] + CH[4] * p_k5[i] + CH[5] * p_k6[i];
		p[i + dim] = fabs((CT[0] * p_k1[i] + CT[1] * p_k2[i] + CT[2] * p_k3[i] + CT[3] * p_k4[i] + CT[4] * p_k5[i] + CT[5] * p_k6[i]) / p[i]);
	}
	free(p_k1);
	free(p_k2);
	free(p_k3);
	free(p_k4);
	free(p_k5);
	free(p_k6);
	return p;
}

/* return the maxiumum numerical value from an array pointer by p 
 * dim = size of the array
 */
double max(double *p, int dim)
{
	int i;
	double max_num = 0;

	for(i = 0; i < dim; i++)
	{
		if(p[i] > max_num)
			max_num = p[i];
	}

	if(max_num == 0)
	{
		printf("No number greater than 0! returning max_val = 1\n");
		return 1;
	}
	return max_num;
}

/*
------------------------------------------------------------------------------------------

File naming section

------------------------------------------------------------------------------------------
*/

/* return a pointer to the datafile name string accrdoing to the setting in binary_conf.c */
char *string_format(struct Pa p)
{
    char MA[10]; // Star A total mass (M+m) in (solar mass)
    char mA[10]; // Star A m/(M+m) ratio
    char MB[10]; // Star B total mass (M+m) in (solar mass)
    char mB[10]; // Star B m/(M+m) ratio
    char RA[10]; // Star A radius without tide and DM core in unit of (km)
    char nA[10]; // Star A and B polytropic index n
    char spin[10]; // Irrotational Star A and B as Riemann ellipsoid
    char ecc[10]; // Orbital eccentricity (We assume all circular so that e_0 = 0)
    char rp[10]; // Initial orbital separation in unit of (R_0_A)
    char alpha[10]; // Dark force coupling constant
    char mV[10]; // Dark force mediator mass range
    char pn_order[10]; // PN order
    char digit[10] = "%.1f";
    char digit2[10] = "%.2f";
    char digit3[10] = "%.3f";
	char pA[10]; // Star A m/M ratio
	char pB[10]; // Star B m/M ratio

    char *container;
    container = (char *) malloc(200 * sizeof(char)); // buffer to store data file name string

	// PNc mean Conservative PN terms
    if (p.pn.PNc == 0)
        strcpy(pn_order, "C0");
    else if ((p.pn.PNc == 1 && p.pn.PN3 == 0) && (p.pn.PN1 == 1 && p.pn.PN2 == 0))
        strcpy(pn_order, "C1P1");
    else if ((p.pn.PNc == 1 && p.pn.PN3 == 0) && (p.pn.PN1 == 1 && p.pn.PN2 == 1))
        strcpy(pn_order, "C1P12");
    else if ((p.pn.PNc == 1 && p.pn.PN3 == 1) && (p.pn.PN1 == 1 && p.pn.PN2 == 1))
        strcpy(pn_order, "C1P123");

    sprintf(MA, digit2, p.M_A);
    sprintf(mA, digit3, p.m_A);
    sprintf(MB, digit2, p.M_B);
    sprintf(mB, digit3, p.m_B);
    sprintf(RA, digit, p.R_0_A);
    sprintf(nA, digit, p.n_A);
    sprintf(spin, digit, p.f_R_A);
    sprintf(ecc, "e0.0");
    sprintf(rp, digit, p.r);
    sprintf(alpha, digit2, p.dark.alpha);
    sprintf(mV, digit, p.dark.mV);
	sprintf(pA, digit2, p.p_A);
	sprintf(pB, digit2, p.p_B);

	// D2 mean dissipative PN2.5 order 
    if (p.dark.on == 0)
        sprintf(container, "%s_D2_MA%s_mA%s_MB%s_mB%s_RA%s_nA%s_spin%s_%s_rp%s_pA%s_pB%s.txt",
                pn_order, MA, mA, MB, mB, RA, nA, spin, ecc, rp, pA, pB);
    else if (p.dark.on == 1)
        sprintf(container, "%s_D2_MA%s_mA%s_MB%s_mB%s_RA%s_nA%s_spin%s_%s_rp%s_alpha%s_mV%s_pA%s_pB%s.txt",
                pn_order, MA, mA, MB, mB, RA, nA, spin, ecc, rp, alpha, mV, pA, pB);
    else
        sprintf(container, "Invalid dark.on!!");
    return container;

}

/*
------------------------------------------------------------------------------------------

 BH-BH section 

------------------------------------------------------------------------------------------
*/

// Set d^2 r/dt^2 = 0 to solve for initial orbital frequency with PN terms for BH-BH inspiral
double func_orb_BH(double x, void * params)
{
    struct Pa p = * ((struct Pa *) params);
    double Mt = p.M_A + p.M_B;
    double eta = p.M_A * p.M_B / (Mt * Mt);
    double r = p.r * p.R_0_A * 1000;
    double ra_c = pnc(r, 0.0, x, Mt, eta, "r", p);
    double dark = dark_eqn(r, p);

    double result = r * x * x - (Mt / (r * r)) * dark + p.pn.PNc * (Mt / (r * r)) * ra_c;
    return result;
}

// Use secant method to solve for initial orbital frequency
double secant_orb_BH(double func(double, void *), void * params)
{
    struct Pa p = * ((struct Pa *) params);
    double ans, x1, x0, err;
    x1 = 2.0e-6;
    x0 = 3.0e-6;
    err = fabs(x1 - x0);
    int i = 0;
    while (err > 1e-15)
    {
        double f_x1 = func_orb_BH(x1, &p);
        double f_x0 = func_orb_BH(x0, &p);
        double new_x = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0);

        x0 = x1;
        x1 = new_x;
        err = fabs(x1 - x0);
        i++;
    }

    ans = x1;
    //printf("total iteration = %d orb is % .5e\n", i, ans);
    return ans;

}

// Call the above secant method to get initial orbital frequency
struct Eqm initial_BH(struct Pa p)
{
	struct Eqm a;
    a.orb = secant_orb_BH(func_orb_BH, &p);

    return a;
}

// Data file name for BH-BH inspiral
char *string_format_BH(struct Pa p)
{
    char MA[10]; // Star A total mass in (solar mass)
    char MB[10]; // Star B total mass in (solar mass)
    char RA[10]; // Reference radius for initial orbital separation
    char rp[10]; // Initial orbital separation in unit of (R_0_A) 
    char alpha[10]; // Dark force coupling constant
    char mV[10]; // DArk force mediator mass range
    char pn_order[10];
    char digit[10] = "%.1f";
    char digit2[10] = "%.2f";
    char digit3[10] = "%.3f";

    char *container; // buffer to store the data file name string
    container = (char *) malloc(200 * sizeof(char));

    if (p.pn.PNc == 0)
        strcpy(pn_order, "C0");
    else if ((p.pn.PNc == 1 && p.pn.PN3 == 0) && (p.pn.PN1 == 1 && p.pn.PN2 == 0))
        strcpy(pn_order, "C1P1");
    else if ((p.pn.PNc == 1 && p.pn.PN3 == 0) && (p.pn.PN1 == 1 && p.pn.PN2 == 1))
        strcpy(pn_order, "C1P12");
    else if ((p.pn.PNc == 1 && p.pn.PN3 == 1) && (p.pn.PN1 == 1 && p.pn.PN2 == 1))
        strcpy(pn_order, "C1P123");

    sprintf(MA, digit2, p.M_A);
    sprintf(MB, digit2, p.M_B);
    sprintf(RA, digit, p.R_0_A);
    sprintf(rp, digit, p.r);
    sprintf(alpha, digit2, p.dark.alpha);
    sprintf(mV, digit, p.dark.mV);

    if (p.dark.on == 0)
        sprintf(container, "%s_D2_MA%s_MB%s_RA%s_rp%s.txt",
                pn_order, MA, MB, RA, rp);
    else if (p.dark.on == 1)
        sprintf(container, "%s_D2_MA%s_MB%s_RA%s_rp%s_alpha%s_mV%s.txt",
                pn_order, MA, MB, RA, rp, alpha, mV);
    else
        sprintf(container, "Invalid dark.on!!");
    return container;

}

// The system of ODE for BH-BH inspiral
double *eqn_BH(double s[], struct Pa pa)
{
	double *p;
	int dim = 4;
	p = (double *) malloc((dim + 1) * sizeof(double));
	double Mt = pa.M_A + pa.M_B;
	double eta = pa.M_A * pa.M_B / (Mt * Mt);

	double ra_c = pnc(s[1], s[0], s[2], Mt, eta, "r", pa);
	double orba_c = pnc(s[1], s[0], s[2], Mt, eta, "o", pa);
	double ra_d = pnd(s[1], s[0], s[2], Mt, eta, "r", pa);
	double orba_d = pnd(s[1], s[0], s[2], Mt, eta, "o", pa);

	double dark = dark_eqn(s[1], pa);

	*p = s[1] * pow(s[2], 2.0) - (Mt / pow(s[1], 2.0)) * dark + pa.pn.PNc * (Mt / pow(s[1], 2.0)) * ra_c + ra_d;

	*(p + 1) = s[0];

	*(p + 2) = -2.0 * s[0] * s[2] / s[1] + pa.pn.PNc * (Mt / pow(s[1], 2.0)) * orba_c + orba_d;

	*(p + 3) = s[2];

	return p;
}

/*
------------------------------------------------------------------------------------------

 Read data file for GW calculation section

------------------------------------------------------------------------------------------
*/

// Count the number of rows in the data file
int line_counter(char file[], char path[]) {

	FILE *fp;
	int count = 0;
	char c;

	strcat(path, file);
	fp = fopen(path, "r"); // Open the data file

	for (c = getc(fp); c != EOF; c = getc(fp))
		if (c == '\n') // If newline is found
			count = count + 1; // Count how many newline

	fclose(fp); // Close the file

	return count - 1;
}

/* The get_data function return a pointer containing the data stored in the data file 
 * It read the content of the data file 
 */
double *get_data (char file[], char path[], int* lenn, int col) {

	double *p;
	unsigned int i, j;
	int len;
	len = line_counter(file, path); // Get the number of row of the data file
	*lenn = len; // Return the numbe of row as pointer to int pointer *lenn

	// Initialize the 2d array (actually 1D array)
	p = (double *) malloc(len * col * sizeof(double)); // col is the column number to read from the data file

	FILE *fp = fopen(path, "r");
	// Check if file exist, return NULL if error
	if (fp == NULL) {
		fprintf(stderr, "Couldn't open %s: %s\n", path, strerror(errno));
		p = NULL;
		return p;
	}

	// Store the file content in 2d array
	for (i = 0; i < len; i++) {
		for (j = 0; j < col; j++) {
			fscanf(fp, "%lf", p + i * col + j);
		}
	}
	fclose(fp);

	return p;
}

// Cut the string given the start and ending position and return it as *storage
void cut_string(char *s, int start, int end, double *storage){
	
	char dest[10];
	size_t size = strlen(s);
	memcpy(dest, s + start, end);
	sscanf(dest, "%lf", storage);
}

// Cut the string given the ending position only and return it as *dest
char *cut_string2(char *s, int end_reverse) {
	
	char *dest = malloc(200 * sizeof(char));
	size_t size = strlen(s);
	memcpy(dest, s, size - end_reverse);
	return dest;
}

// Return the binary parameters structure variable by given data file name
struct Pa get_const(char file[]){

	double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
	double mass_const = 1.3271244e20 / (c * c);
	struct Pa pa;
	char *token; // Used to spilit the file name accroding to "_"
	char *cut_file; // Store the file name without extension
	cut_file = cut_string2(file, 4); // Ignore the file name extension
	int i;

	// Obtaining PN order from file name
	token = strtok(cut_file, "_");	
	// printf("%s\n", token); // Debug
	if (strcmp(token, "C1P123") == 0) {
		pa.pn.PNc = 1;
		pa.pn.PN1 = 1;
		pa.pn.PN2 = 1;
		pa.pn.PN25 = 1;
		pa.pn.PN3 = 1;
	}
	else if (strcmp(token, "C1P12") == 0) {
		pa.pn.PNc = 1;
		pa.pn.PN1 = 1;
		pa.pn.PN2 = 1;
		pa.pn.PN25 = 1;
		pa.pn.PN3 = 0;
	}
	else if (strcmp(token, "C1P1") == 0) {
		pa.pn.PNc = 1;
		pa.pn.PN1 = 1;
		pa.pn.PN2 = 0;
		pa.pn.PN25 = 1;
		pa.pn.PN3 = 0;
	}
	else if (strcmp(token, "C0") == 0) {
		pa.pn.PNc = 0;
		pa.pn.PN1 = 0;
		pa.pn.PN2 = 0;
		pa.pn.PN25 = 1;
		pa.pn.PN3 = 0;
	}
	else 
		printf("Probably Wrong file name\n");

	// Obtain binary mass and polytropic index
	for (i = 0; i < 10; i++){
		token = strtok(NULL, "_");
		if (i == 1)
			cut_string(token, 2, 5, &pa.M_A);
		if (i == 2)
			cut_string(token, 2, 5, &pa.m_A);
		if (i == 3)
			cut_string(token, 2, 5, &pa.M_B);
		if (i == 4)
			cut_string(token, 2, 5, &pa.m_B);
		if (i == 6)
			cut_string(token, 2, 5, &pa.n_A);
	}

	// Check if the file contain dark force
	char mV[10] = "mV";
	char choice;
	if (strstr(file, mV) != NULL)
		choice = 'm';
	else
		choice = 'd';
	
	// Dark Force Constant
	if (choice == 'm') {
		pa.dark.on = 1;
		for (i = 0; i < 2; i++) {
			token = strtok(NULL, "_");
			if (i == 0) 
				cut_string(token, 5, 9, &pa.dark.alpha);
			if (i == 1)
				cut_string(token, 2, 5, &pa.dark.mV);
		}
	}	
	else if (choice == 'd') {
		pa.dark.on = 0;
		pa.dark.alpha = 0.1;
		pa.dark.mV = 20.0;
	}
	
	double M_A = pa.M_A * (1.0 - pa.m_A);
	double m_A = pa.M_A * pa.m_A;
	double M_B = pa.M_B * (1.0 - pa.m_B);
	double m_B = pa.M_B * pa.m_B;
	// Convert the value obtained from file name to geometric unit
	pa.M_A = M_A * mass_const;
	pa.m_A = m_A * mass_const;
	pa.M_B = M_B * mass_const;
	pa.m_B = m_B * mass_const;

	// Obtain polytropic constant from polytropic index n
	double *kapb;
	kapb = kap_b(pa.n_A);
	pa.n_B = pa.n_A;
	pa.c_n_A = 3.0 / (5.0 - pa.n_A);
	pa.c_n_B = pa.c_n_A;
	pa.b_n_A = *(kapb + 1);
	pa.b_n_B = pa.b_n_A;
	pa.kap_A = *kapb;
	pa.kap_B = pa.kap_A;

	free(kapb);
	free(cut_file);
	return pa;
}

// Orbital equations to obtain acceleration for orbital separation and orbital phase
double *eqn_orb(double s[], struct Pa pa){

	double *p;
	int dim = 2;
	p = malloc(2 * sizeof(double));
	double Mt = pa.M_A + pa.m_A + pa.M_B + pa.m_B;
	double eta = (pa.M_A + pa.m_A) * (pa.M_B + pa.m_B) / (Mt * Mt);
	double r_tide_A = tide_orb(pa.kap_A, pa.M_A, pa.m_A, s[1], s[5], s[9], s[21], s[13], "r");

	double orb_tide_A = tide_orb(pa.kap_A, pa.M_A, pa.m_A, s[1], s[5], s[9], s[21], s[13], "o");
	double r_tide_B = tide_orb(pa.kap_B, pa.M_B, pa.m_B, s[3], s[7], s[11], s[21], s[15], "r");
	double orb_tide_B = tide_orb(pa.kap_B, pa.M_B, pa.m_B, s[3], s[7], s[11], s[21], s[15], "o");
	double ra_c = pnc(s[19], s[18], s[20], Mt, eta, "r", pa);
	double orba_c = pnc(s[19], s[18], s[20], Mt, eta, "o", pa);
	double ra_d = pnd(s[19], s[18], s[20], Mt, eta, "r", pa);
	double orba_d = pnd(s[19], s[18], s[20], Mt, eta, "o", pa);
	
	double dark = dark_eqn(s[19], pa);
	*p = -(3.0 / 20.0) * (Mt / pow(s[19], 4.0)) * (r_tide_A + r_tide_B)
		 + s[19] * pow(s[20], 2.0) - (Mt / pow(s[19], 2.0)) * dark + pa.pn.PNc * (Mt / pow(s[19], 2.0)) * ra_c + ra_d;

 	*(p + 1) = -2.0 * s[18] * s[20] / s[19] - (3.0 * Mt / (10.0 * pow(s[19], 5.0))) * (orb_tide_A + orb_tide_B)
			   + pa.pn.PNc * (Mt / pow(s[19], 2.0)) * orba_c + orba_d;

	return p;	   
}

/* GW h_+ and h_x polarization modes
 * r = orbital separation
 * rv = first time derivative of r
 * ra = second time derivative of r
 * o = orbital phase
 * ov = first time derivative of o
 * oa = second time derivative of o
 */
double *TT_gauge(double r, double rv, double ra, 
				 double o, double ov, double oa)
{
	double *hpc;
	double a;
	double b;
	hpc = malloc(2 * sizeof(double));

	a = rv * rv * cos(2 * o) + r * ra * cos(2 * o)
	  - 4 * r * rv * ov * sin(2 * o)
	  - 2 * r * ov * r * ov * cos(2 * o)
	  - r * r * oa * sin(2 * o);

	b = rv * rv * sin(2 * o) + r * ra * sin(2 * o)
	  + 4 * r * rv * ov * cos(2 * o)
	  - 2 * r * ov * r * ov * sin(2 *o)
	  + r * r * oa * cos(2 * o);

	*hpc = a;
	*(hpc + 1) = b;

	return hpc;
}

//---------------------------------------------------------------------------

#pragma hdrstop

#include "basic_functions.h"

double** Legendre_matrix_GOST(int number_of_harmonics, double phi)
{
	// (number_of_harmonics + 1)*(number_of_harmonics + 1) matrix of Fully Normalized Associated Legendre Functions
	// (FNALF) values.
	// For the increase of speed of the calculations, I store the values of
	// FNALF in matrix, to calculate them only one time for each knot.

	int NOHplus1 = number_of_harmonics + 1;

	double** P = new double *[NOHplus1]; // matrix of FNALF values in the knot concidered

	int n, m;
	for (n = 0; n < NOHplus1; n++)
	{
		P[n] = new double [NOHplus1];
		for (m = 0; m < NOHplus1; m++) //if m <= n, then P[n][m] exists
		{

			double nd = double(n);
			double md = double(m);

			P[n][m] = 0;

			if (n > m && n != 1)
			{
				double z;
				z = nd * nd  ;

				P[n][m] =  P[n - 1][m] * sin(phi) * pow(((4*z - 1)/(pow(nd, 2) - pow(md, 2))), 0.5) -
					P[n - 2][m] *
					pow(((pow(nd - 1.0, 2) - pow(md, 2))*(2*n + 1)/((pow(nd, 2) - pow(md, 2))*(2*n - 3))), 0.5);
			 // it works, but what happens with P[n-2][m], when n = 1 and m = 0?
			 // P[-1][0] ???!?!?
			 // THAT'S THE QUESTION.
			 // One day it might become a realy big shit.
			 // I had to create one more "if" to prevent this. I don't know
			 // if it is right.
			 // Question still exists.
			 // One more question -
			}
			else if (n > m && n == 1)
			{
				P[n][m] =  P[n - 1][m] * sin(phi) * pow(((4*pow(nd, 2) - 1)/(pow(nd, 2) - pow(md, 2))), 0.5);
			}
			else if ((n == m) && (n != 0))
			{
				P[n][m] =  P[n - 1][m - 1] * cos(phi) * pow((2*nd + 1.0)/(2*nd), 0.5);
				//P[n][m] =  P[n - 1][m - 1] * pow(1 - pow(sinphi, 2), 0.5) * pow((2*nd + 1.0)/(2*nd), 0.5);
			}
			else if (n < m)
			{
				P[n][m] = 0.0;  // zero or NAN ?
			}
			else if ((n == m) && (n == 0))
			{
				P[n][m] = 1.0/pow(2, 0.5); //because it is fully normalized
			}
			else cout << "error when tryig to calculate Legendre pol. m or n are NAN " << endl;

			//std::cout << P[n][m] << ", ";  //stupid debugging
		}

		//std::cout << std::endl; //stupid debugging

	}
	//cout << P[number_of_harmonics][number_of_harmonics] << endl;
	return P;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
double** Legendre_derivatives_matrix(int number_of_harmonics, double phi) //ïîòîì ñâåðèòü ñ ìàòëàáîì, ò.ê. ôîðìóëà ìîæåò áûòü íåâåðíîé
{
	int NOHplus1 = number_of_harmonics + 1;
	double** der_P = new double *[NOHplus1]; // matrix of FNALF values in the knot concidered

	double** P = Legendre_matrix_GOST( number_of_harmonics, phi);

	int n, m;
	for (n = 0; n < NOHplus1; n++)
	{
		//cout << "n is " << n << endl;
		der_P[n] = new double [NOHplus1];
		for (m = 0; m < NOHplus1; m++) //if m <= n, then P[n][m] exists
		{
			//cout << "m is " << m << ", " << endl;

			double nd = double(n);
			double md = double(m);

			if (m + 1 < NOHplus1)
			{
				der_P[n][m] = - md * std::tan(phi) * P[n][m] + sqrt(( nd + md + 1 ) * (nd - md )) * P[n][m + 1];
			}
			else
			{
				der_P[n][m] = - md * std::tan(phi) * P[n][m];
			}



			std::cout << "der_P[" << n << "][" << m << "] = " << der_P[n][m] << endl;  //stupid debugging
		}

		std::cout << std::endl; //stupid debugging

	}
	//cout << P[number_of_harmonics][number_of_harmonics] << endl;
	return der_P;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
double Factorial(int n)
{
	double* fact = new double [ n + 1 ];
	fact[0] = 1.0;
	for (int i = 1; i <= n; i++)
	{
		fact[i] = double(i) * fact[ i - 1 ];
	}

	return fact[n];
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
double J(int n)
{
	double alpha = 1.0/298.25784; //table 2, page 7
	double a = 6378136.0;
	double omega = 0.7292115 * 1E-4;
	double fM = 398600.4418 * 1E+9;
	double e = pow(2.0 * alpha - pow(alpha, 2), 0.5);

	double output;


	if (n == 2)
	{
		double e1 = pow(pow(e, 2)/(1-pow(e, 2)), 0.5);
		double m = pow(omega, 2) * pow(a, 3)/fM;
		double q0 = 1.0 / 2.0 * ((1 + 3 / pow(e1, 2)) * atan(e1) - 3 / e1);
		output = pow(e, 2) / 3 - 2.0 / 45.0 * m * pow(e, 3) / q0;
	}
	else if (n == 4 || n == 6 || n == 8)
	{
		output = pow(-1.0, n/2 + 1) * (3 * pow(e, n - 2) /
			(2 * (n + 1) * (n + 3))) * (5 * n * J(2) - (n - 2) * pow(e, 2));
	}
	else
	{
		output = 0;
		cout << "Warning. In J calculation n isn't 2, 4, 6 or 8. Maybe, there is an error." << endl;
	}

	return output;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int EGM_reader(string path, int n_max, double** delta_C, double** F)
{
	//--------------------------------------------------------------------------
	//ñï¸ð êîä îòñþäà
	//http://blog.harrix.org/article/5819#h2_0
	//--------------------------------------------------------------------------
	//	Ïðèìåð èñïîëüçîâàíèÿ ôóíêöèè:
	//  int n_max = 16;
	//	double **C;
	//	double **F;
	//	C = new double* [n_max];
	//	F = new double* [n_max];
	//	for (int n = 2; n <= n_max ; n++)
	//	{
	//		C[n] = new double [n + 1];
	//		F[n] = new double [n + 1];
	//	}
	//  int ier = EGM_reader("..//..//egm96_to360.ascii", n_max, C, F);
	//  cout << C[2][1] << endl;
	//--------------------------------------------------------------------------


	//Ñîçäàåì ôàéëîâûé ïîòîê è ñâÿçûâàåì åãî ñ ôàéëîì
	//÷òî çíà÷èò c_str() ? - ïðåîáðàçîâàíèå str â char
	ifstream file(path.c_str());

	if (file.is_open())//Åñëè îòêðûòèå ôàéëà ïðîøëî óñïåøíî
	{
		double **C;
		C = new double* [n_max + 1];
		for (int n = 2; n <= n_max ; n++)
		{
			C[n] = new double [n + 1];
		}

		cout << "EGM file has been opened." << endl;

		string line;//Ñòðî÷êà òåêñòà

		//Áóäåì ñ÷èòûâàòü èíôîðìàöèþ ïîñòðî÷íî äî òåõ ïîð,
		//ïîêà íå çàêîí÷èòñÿ ôàéë
		while (getline(file, line))
		{
			//Òåïåðü â line õðàíèòñÿ ñîäåðæèìîå ñòðî÷êè èç ôàéëà.
			//Áóäåì å¸ ðàçáèðàòü íà ñîñòàâíûå ÷àñòè.
			//Â íàøåì ôàéëå èäóò äâà int, à ïîòîì 4 ÷èñëà double
			int n, m;
			double Cf, Ff, dCf, dFf;

			//Ñîçäàäèì ïîòîê äëÿ ñ÷èòûâàíèÿ äàííûõ èç ñòðî÷êè
			istringstream iss(line);

			//Òåïåðü ÷åðåç ñòàíäàðòíûé îïåðàòîð >> ñ÷èòàåì äàííûå
			//Ïðîãðàììà ñàìà ïîéìåò ÷òî â êà÷åñòâå ðàçäåëèòåëÿ íàäî
			//èñïîëüçîâàòü çíàê òàáóëÿöèè \t
			//Çà ðàç âñ¸ ñ÷èòàåì
			iss >> n >> m >>  Cf >>  Ff >>  dCf >>  dFf;

			C[n][m] = Cf;
			F[n][m] = Ff;

//			//Âûâåäåì íàøè äàííûå
//			cout << "Data from string:" << endl;
//			cout << "\t n: " << n << endl;
//			cout << "\t m: " << m << endl;
//			cout << "\t C: " << C[n][m] << endl;
//			cout << "\t F: " << F[n][m] << endl;
//			cout << "\t dC: " << dCf << endl;
//			cout << "\t dF: " << dFf << endl;

			if (m >= n_max)
			{
				return 0;
			}

			//calculation of delta_C
			for (int n = 2; n <= n_max; n++)
			{
				for (int m = 0; m <= n; m++)
				{
					if (n == 2 || n == 4 || n == 6 || n == 8)
					{
						delta_C[n][m] = C[n][m] + J(n)/pow(2 * double(n) + 1, 0.5); // is it right?
					}
					else
					{
						delta_C[n][m] = C[n][m];
					}
				}
			}
		}
	}
	else
	{
		cout << "Error: the EGM file hasn't been opened." << endl;
		std::system("pause");
		std::exit (EXIT_FAILURE);
	}

    return 0;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
double*** second_sum(int number_of_harmonics, double* phi, double* lambda, int p_size,
	int l_size)
{
	//C and F reading from the file, creating of empty delta_C array
	double **C;
	double **F;
	double **delta_C;
	C = new double* [number_of_harmonics];
	F = new double* [number_of_harmonics];
	delta_C = new double* [number_of_harmonics];
	for (int n = 2; n <= number_of_harmonics ; n++)
	{
		C[n] = new double [n + 1];
		F[n] = new double [n + 1];
		delta_C[n] = new double [n + 1];
	}
	int ier = EGM_reader("..//..//egm96_to360.ascii", number_of_harmonics, C, F);
	//--------------------------------------------------------------------------
	//calculation of delta_C
	for (int n = 2; n <= number_of_harmonics; n++)
	{
		for (int m = 0; m <= n; m++)
		{
			if (n == 2 || n == 4 || n == 6 || n == 8)
			{
				delta_C[n][m] = C[n][m] + J(n)/pow(2 * double(n) + 1, 0.5); // is it right?
			}
			else
			{
				delta_C[n][m] = C[n][m];
            }
		}
	}
	//--------------------------------------------------------------------------
	//Calculation of FNALF in each phi[p] point
	//And calculationf of sinphi for each phi[p] point
	int p; //number of phi point
	int n; //for FNALF
	int m; //for FNALF

	//3D-matrix of FNALF for each phi
	double*** P;
	P = new double** [p_size];
	//sinphi for speed
	double* sinphi;
	sinphi = new double [p_size];

	for (p = 0; p < p_size ; p++)
	{
		sinphi[p] = sin(phi[p]); //calculation of sin(phi) for each phi value
		double** P_temp = Legendre_matrix_GOST (number_of_harmonics, phi[p]);
		P[p] = new double* [number_of_harmonics + 1];
		for (n = 2; n <= number_of_harmonics; n++)
		{
			P[p][n] = new double [n + 1];
			for (m = 0; m <= n; m++)
			{
				P[p][n][m] = P_temp[n][m];//3D-matrix of FNALF for each phi
				//p - phi number
				//n - FNALF index
				//m - FNALF index. m <= n
			}
		}
	}
	//--------------------------------------------------------------------------
	//calculation of sinlam and coslam
	int l; //number of lambda point
	//coslam[l][m] for speed
	double** coslam;
	//double** sinlam1;
	double** sinlam;
	coslam = new double* [l_size];
	//sinlam[l][m] for speed
	sinlam = new double* [l_size];


	for (l = 0; l < l_size; l++)
	{
		sinlam[l] = new double [number_of_harmonics + 1];
		coslam[l] = new double [number_of_harmonics + 1];
		for (m = 0; m <= number_of_harmonics; m++)
		{
			sinlam[l][m] = sin(lambda[l] * double(m));
			coslam[l][m] = cos(lambda[l] * double(m));
		}
	}
	//--------------------------------------------------------------------------
	//calculation of second sums
	double*** ss = new double** [p_size];
	for (p = 0; p < p_size; p++)
	{
		ss[p] = new double* [l_size];
		for (l = 0; l < l_size; l++)
		{
			ss[p][l] = new double [number_of_harmonics + 1];
			for (n = 2; n <= number_of_harmonics; n++)
			{
				ss[p][l][n] = 0;
//				cout << "n is " << n << endl;
				for (m = 0; m <= n; m++)
				{
					ss[p][l][n] += (delta_C[n][m] * coslam[l][m] + F[n][m] *
						sinlam[l][m]) * P[p][n][m];
//					cout << "\t m : " << m << endl;
//					cout << "\t delta_C[n][m] : " << delta_C[n][m] << endl;
//					cout << "\t F[n][m] : " << F[n][m] << endl;
//					cout << "\t coslam[l][m] : " << coslam[l][m] << endl;
//					cout << "\t sinlam[l][m] : " << sinlam[l][m] << endl;
//					cout << "\t P[p][n][m] : " << P[p][n][m] << endl;
//					cout << "\t ss[p][l][n]: " << ss[p][l][n] << endl;

				}
//				cout << "\t Second sum value:" << endl;
//				cout << "\t p: " << p << endl;
//				cout << "\t l: " << l << endl;
//				cout << "\t n: " << n << endl;
//				cout << "\t phi[p]: " << phi[p]/M_PI*180.0 << endl;
//				cout << "\t lam[l]: " << lambda[l]/M_PI*180.0 << endl;
//				cout << "\t ss[p][l][n]: " << ss[p][l][n] << endl;
			}
		}
	}
	return ss;
}

int kek()
{
	return 0;
}




//---------------------------------------------------------------------------
#pragma package(smart_init)

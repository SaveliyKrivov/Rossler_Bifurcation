/**************************************************************
 * Аттрактор Рёсслера
 * бифуркации для параметра b
 **************************************************************/

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

// параметры модели ////////////////////////////////////
const double t0 = 0.0;      // временной 
const double T = 100.0;     // интервал
const double a = 0.2;
double b = 0.0;
const double c = 5.7;
/////////////////////////////////////////////////////////

// правые части системы ОДУ 
double f_x( double, double, double, double );
double f_y( double, double, double, double );
double f_z( double, double, double, double );

// метод Рунге-Кутты 4-го порядка
std::vector<double>
Runge_Kutta( double(*f_x)( double, double, double, double ), 
	double(*f_y)( double, double, double, double ), 
	double(*f_z)( double, double, double, double ),
	const double&, const double&, const double&,
	const double&, const double&  );


int main( int argc, char** argv ) 
{
	
	// количество точек 
	int n = 5001;
	double dt = (T - t0) / (n - 1);
	
	std::vector<double> x( n, 0 );
	std::vector<double> y( n, 0 );
	std::vector<double> z( n, 0 );
	
	// границы изменения параметра b
	const double b_min = 0.001;
	const double b_max = 2.0;
	int m = 11001;
	double bc = (b_max - b_min) / (m - 1);
	
	std::ofstream  data_file( "Rossler_b.d" ); // файл с числ. данными;
	std::cout.precision(8); 
	
	b = b_min;
	for ( size_t k = 0;  k < m; ++k  )
	{
		// начальные условия
		x[0] = 1.0;
		y[0] = 1.0;
		z[0] = 0.1;
		for ( size_t i = 1; i < n; i++ )
		{
			std::vector<double> 
			p = Runge_Kutta( f_x, f_y, f_z, dt, t0 + dt*(i - 1), 
				x[i - 1], y[i - 1], z[i - 1] );
			x[i] = p[0];
			y[i] = p[1];
			z[i] = p[2];
		}
		//double max_x = -1.0e16; 
		for ( size_t i = 1; i < n - 1; i++ )
		{ 
			// локальные максимумы x
			if ( ( x[i - 1] < x[i] ) && (x[i] > x[i + 1]) ) 
				data_file << std::setw(16) << b << std::setw(16) << x[i] << "\n";	;	
		}	
		
		b += bc;
	}
	std::cout << "Data file: Rossler_b.d\n";		   
	data_file.close();
	return 0;
}

//-------------------------------------------------------------------
// метод Рунге-Кутты 4-го порядка
// решения системы 3-х ОДЕ
// w' = f_w
// z' = f_z
// v' = f_v
//-------------------------------------------------------------------
std::vector<double>
Runge_Kutta( double(*f_x)( double, double, double, double ), 
	double(*f_y)( double, double, double, double ), 
	double(*f_z)( double, double, double, double ),
	const double &dt, 
	const double &t, 
	const double &x,
	const double &y,
	const double &z  )
{
	double	
	k11 = dt * f_x( t, x, y, z ),
	k12 = dt * f_y( t, x, y, z ), 
	k13 = dt * f_z( t, x, y, z ), 
	
	k21 = dt * f_x( t + dt/2, x + k11/2, y + k12/2, z + k13/2 ),
	k22 = dt * f_y( t + dt/2, x + k11/2, y + k12/2, z + k13/2 ),
	k23 = dt * f_z( t + dt/2, x + k11/2, y + k12/2, z + k13/2 ),
	
	k31 = dt * f_x( t + dt/2, x + k21/2, y + k22/2, z + k13/2 ),
	k32 = dt * f_y( t + dt/2, x + k21/2, y + k22/2, z + k13/2 ),
	k33 = dt * f_z( t + dt/2, x + k21/2, y + k22/2, z + k13/2 ),
	
	k41 = dt * f_x( t + dt,   x + k31,   y + k32, z + k33   ),
	k42 = dt * f_y( t + dt,   x + k31,   y + k32, z + k33   ),
	k43 = dt * f_z( t + dt,   x + k31,   y + k32, z + k33   );
	
	return { x + (k11 + 2*k21 + 2*k31 + k41) / 6, 
		y + (k12 + 2*k22 + 2*k32 + k42) / 6,
		z + (k13 + 2*k23 + 2*k33 + k43) / 6  };
}
//-------------------------------------------------------------------
// правая часть уравнения x' = f_x( t, x, y, z ) 
//-------------------------------------------------------------------
double f_x( double t, double x, double y, double z )
{
	return -y - z; 
}
//-------------------------------------------------------------------
// правая часть уравнения y' = f_y( t, x, y, z ) 
//-------------------------------------------------------------------
double f_y( double t, double x, double y, double z )
{
	return  x + a*y; 
}
//-------------------------------------------------------------------
// правая часть уравнения z' = f_z( t, x, y, z ) 
//-------------------------------------------------------------------
double f_z( double t, double x, double y, double z )
{
	return b + z*(x - c); 
}


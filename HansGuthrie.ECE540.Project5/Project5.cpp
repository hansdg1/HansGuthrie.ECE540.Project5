#pragma warning( disable:4996 )
#include <math.h>
#include <stdio.h>
#include <complex>

using namespace std;

#define eps 2.22e-16
int main( );
void part_b( );

// Support Routines for complex class.
double mag( complex<double> a )
{
	double r, i;
	r = a.real( );
	i = a.imag( );
	return( sqrt( r*r + i*i ) );
}

complex<double> conjg( complex<double> a )
{
	double r, i;
	r = a.real( );
	i = a.imag( );
	return( complex<double>( r, -i ) );
}

// Simple function to write out an array of complex numbers.
void WriteToFile( char *name, complex<double> *out, int length, double Dt )
{
	FILE *fout = fopen( name, "w" );

	double T = 0.0;

	for ( int k = 0; k < length; k++ )
		fprintf( fout, "%18.16lg,%18.16lg,%18.16lg\n", out[ k ].real( ), out[ k ].imag( ), T = T + Dt );
	fclose( fout );
}

// Program to simulate the steady state response of
// a circuit to a square wave.
#define SimulationSteps 130000
#define SimulationSteps2 20000

void part_a( )
{
	double w,
		T = 4.0, // Period set at 4.0 seconds.
		Dt = 0.1e-3, // Step length (seconds)
		J = 100e-3,
		Kv = 5.0,
		PI = 4.0*atan( 1.0 );

	complex<double> Hjkw, Xk, Outw, s;
	complex<double> *in = new complex<double>[ SimulationSteps ],
		*out = new complex<double>[ SimulationSteps ]; // Array of complex to hold outputs.

	int k, m, Ki, Kp, i;
	char name[ 32 ];

	FILE *fout = fopen( "FourierCoefficients.csv", "w" );

	w = 2 * PI / T;

	//Dt = 3.0*T / SimulationSteps; // Compute step size to generate 3 periods.

	for ( Kp = 1; Kp <= 4; Kp *= 2 )
	{
		for ( i = 0; i <= 2; i++ )
		{
			if ( i == 0 )
			{
				Ki = 0;
			}
			else
			{
				Ki = Kp*( 1 + i );
			}

			// Initialize output as zero.
			for ( m = 0; m < SimulationSteps; m++ )
			{
				out[ m ] = 0.0;
				in[ m ] = 0.0;
			}

			// Generate output
			k = 1;
			while ( k < 10 )
			{
				if ( k ) // Compute Fourier Coefficients for square wave.
				{
					Xk = 1.0 * sin( k*PI / 2 ) / ( k*PI / 2 );
				}

				/*else // k == 0...
				Xk = 0.0;
				*/

				// Compute Filter response at this frequency
				s = complex<double>( 0.0, k*w );
				Hjkw = ( Kv*Kp*s + Kv*Ki ) / ( J*s*s*s + s*s + Kv*Kp*s + Kv*Ki );

				//fprintf( fout, "%18.16lg, %18.16lg, %18.16lg, %18.16lg\n", Xk.real(), Xk.imag(), Hjkw.real(), Hjkw.imag() );

				// Add in this term into steady state response.
				for ( m = 0; m < SimulationSteps; m++ )
				{
					// Create the value of exp(j*t*w) for frequency (k*w) and time (m*Dt)
					Outw = exp( complex<double>( 0.0, m*Dt*k*w ) );

					// Compute Fourier Series representation of Input.
					in[ m ] = in[ m ] + Xk * Outw;
					// Note k = 0 is a special case
					if ( k ) // which is not a conjugate pair.
						in[ m ] = in[ m ] + conjg( Xk * Outw );

					// Compute Fourier Series representation of Output.
					//out[ m ] = out[ m ] + Xk * Hjkw * Outw;
					out[ m ] += Xk * Hjkw * Outw;
					// Note k = 0 is a special case
					if ( k ) // which is not a conjugate pair.
						out[ m ] += conjg( Xk * Hjkw * Outw );
				}// End of Loop through time steps.

				if ( k >= 9 )
				{
					sprintf( name, "FS_output_Kp_%d_Ki_%d.csv", Kp, Ki );
					WriteToFile( name, out, SimulationSteps, Dt );
				}

				// Move k ahead, in pattern 0,1,3,5,7,...
				//if ( k == 0 )
				//{
				//	k = 1;
				//}
				//else
				k += 2;
			}// End of loop through Fourier Series Components.
		}
	}
	printf( "Turn down for what\n" );
}// End of part_a

void part_b( )
{
	double w,
		T = 4,
		Dt = 0.1e-3,
		J = 100e-3,
		Kv = 5,
		PI = 4.0*atan( 1.0 );

	complex<double> Hjkw, Xk, Outw, w1, w2, s;
	complex<double> *in = new complex<double>[ SimulationSteps2 ],
		*out = new complex<double>[ SimulationSteps2 ]; // Array of complex to hold outputs.

	int k, m, Kp, Ki, n;
	char name[ 32 ];
	FILE *fout = fopen( "FourierCoefficients2.csv", "w" );
	w = 2 * PI / T;

	for ( Kp = 1; Kp <= 4; Kp = 2 * Kp )
	{
		Ki = 2 * Kp;
		for ( m = 0; m < SimulationSteps2; m++ )
		{
			out[ m ] = 0.0;
			in[ m ] = 0.0;
		}
		// Generate input and output
		k = 0;
		while ( k < 10 )
		{
			if ( k ) // Compute Fourier Coefficient k for sawtooth wave.
			{
				w1 = complex<double>( 0.0, 0.2*k*PI );
				w2 = complex<double>( 0.0, 1.8*k*PI );
				Xk = ( 1 / 0.9 )*( -10.0 + 9.0*( 1.0 - w1 )*exp( w1 ) + ( 1.0 + w2 )*exp( -w2 ) ) / ( 4.0*k*k*PI*PI );
			}
			else // k == 0...
				Xk = 0.5;
			// Compute Filter response at this frequency
			s = complex<double>( 0.0, k*w );
			Hjkw = ( Kv*Kp*s + Kv*Ki ) / ( J*s*s*s + s*s + Kv*Kp*s + Kv*Ki );
			fprintf( fout, "%18.16lg, %18.16lg, %18.16lg, %18.16lg\n",
				Xk.real( ), Xk.imag( ), Hjkw.real( ), Hjkw.imag( ) );

			// Loop through time.
			for ( m = 0; m < SimulationSteps2; m++ )
			{
				// Create the value of exp( j*t*w ) for this
				// frequency (k*w) and time (m*Dt);
				Outw = exp( complex<double>( 0.0, m*Dt*k*w ) );

				// Compute Fourier Series representation of Input.
				in[ m ] = in[ m ] + Xk * Outw;
				// Note k = 0 is a special case
				if ( k ) // which is not a conjugate pair.
					in[ m ] = in[ m ] + conjg( Xk * Outw );

				// Compute Fourier Series representation of Output.
				out[ m ] = out[ m ] + Xk * Hjkw * Outw;
				// Note k = 0 is a special case
				if ( k ) // which is not a conjugate pair.
					out[ m ] = out[ m ] + conjg( Xk * Hjkw * Outw );
			} // End of time loop.

			if ( k >= 9 )
			{
				// Save off input and output model at this point.
				sprintf( name, "FS_RST_output_Kp_%i_Ki_%i.csv", Kp, Ki );
				WriteToFile( name, out, SimulationSteps2, Dt );
			}
			k++;
		}
	}
	printf( "More realistic input is done\n" );
}//end of part_b

int main( )
{
	part_a( );
	part_b( );
	printf( "press any key to quit." );
	getchar( );
}
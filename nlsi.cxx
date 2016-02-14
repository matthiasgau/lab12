#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void NLStep( cmplx* const psi,   const double dt, const int Nx);
void init( cmplx* const psi0, const double eta, const double sigma, const double dx,
          const int Nx);
void LinStep( cmplx* const psi, const double dx, const double dt, const int Nx);
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin);
//-----------------------------------
int main(){

	const int Nx = 400;
	const double L = 800;
	const double xmin = 0;
	const double Tend = 50;
	const double dx = L / (Nx - 1);
	const double dt = dx  / 10;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double eta = 0.2;

	stringstream strm;

	cmplx* psi0 = new cmplx[Nx];

	init(psi0, eta, dx, dt,Nx);

	writeToFile(psi0,"psi_0", dx,Nx,xmin);

	LinStep(psi0, dx, dt/2.0, Nx);

	for (int i = 1; i <= Na; i++) {
		LinStep(psi0, dx, dt/2.0, Nx);
	        NLStep(psi0, dt, Nx);
		for (int j = 1; j <= Nk-1; j++) {
			LinStep(psi0, dx, dt, Nx);
			NLStep(psi0, dt, Nx);
		}
		LinStep(psi0, dx, dt/2.0, Nx);
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin);
	}

        delete[] psi0;
	return 0;
}
//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin)
{
	ofstream out(s.c_str());
	for(int i=0; i<Nx; i++){
		double x = xmin + i * dx;
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag() << endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double eta,  const double dx, const double dt,
          const int Nx)
{
	const double x0 = dx*Nx * 0.5;
	const double f = sqrt(2) * eta;
	for(int i=0;i<Nx; i++){
		double x = i*dx - x0;
		psi0[i] = 2*f/cosh(eta * x);
	}
}
//-----------------------------------
void NLStep( cmplx* const psi, const double dt, const int Nx){
	for(int i=0; i<Nx; i++){
		double phi = abs(psi[i]) * abs(psi[i]) * dt ;
	    psi[i] *= cmplx(cos(phi), sin(phi));
	}
}
//-----------------------------------
void LinStep(cmplx* const rhs, const double dx, const double dt,
						 const int Nx)
{
	cmplx* u =new cmplx[Nx];

	cmplx alpha = cmplx(0.0, dt/(dx*dx));

	cmplx d =  (cmplx(1) + cmplx(2.0,0) * alpha) / (-alpha);

	//forward substitution
	u[0] = cmplx(1) / d;
	rhs[0] = rhs[0] / -alpha /d;

	for(int i=1; i<Nx; i++){
		u[i] = cmplx(1) / (d - u[i-1]);
		rhs[i] = (rhs[i]/(-alpha) - rhs[i-1]) / (d - u[i-1]);
	}


	//backward substitution
	for(int i=Nx-2; i>=0; i--)
		rhs[i] = rhs[i] - u[i]*rhs[i+1];

  delete[] u;
}

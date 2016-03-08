#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double eta, const double sigma, const double dx,
          const int Nx);

void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin);
void stepl(cmplx* const f1, cmplx* const f0,const double dt, const double dx, const int Nx);
void stepn(cmplx* const f1,cmplx* const f0,const double dt, const int Nx);
//-----------------------------------
int main(){

	const int Nx = 4000;
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
	cmplx* psi1 = new cmplx[Nx];
	

	init(psi0, eta, dx, dt,Nx);

	writeToFile(psi0,"psi_0", dx,Nx,xmin);
	
	//stepl(psi0, psi0,dt/2.0,dx,Nx);

	for (int i = 1; i <= Na; i++) {
	  //stepl(psi1,psi0,dt/2.0, dx,Nx);
	  //stepn(psi0,psi1,dt/2.0,Nx);

		for (int j = 1; j <= Nk-1; j++) {
		  stepl(psi1,psi0,dt,dx,Nx);
		  stepn(psi0,psi1,dt,Nx);
		}
		//stepl(psi0,psi0,dt/2.0,dx,Nx);
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin);
	}
      delete[] psi0;
      delete[] psi1;
      
	return 0;
}
//-----------------------------------
void stepl(cmplx* const f1, cmplx* const f0,const double dt, const double dx, const int Nx)
{
 const cmplx j=cmplx(0.0,1.0); 
 cmplx a[Nx], b[Nx], d[Nx];
 for(int i=0; i<Nx; i++){
   a[i]=j*dt/(dx*dx);
   b[i]=a[i];
   d[i]=1.0-2.0*j*dt/(dx*dx);
 }
 for(int i=1; i<Nx; i++){
   d[i]-=b[i]/d[i-1]*a[i-1];
   f0[i]-=b[i]/d[i-1]*f0[i-1];
 }
 f1[Nx-1]=f0[Nx-1]/d[Nx-1];
 for(int i=Nx-2; i>=0; i--)
   f1[i]=(f0[i]-a[i]*f1[i+1])/d[i];
}
//-----------------------------------
void stepn(cmplx* const f1,cmplx* const f0,const double dt, const int Nx)
{
const cmplx j=cmplx(0.0,1.0);
  for(int i=0; i<Nx; i++){
    f1[i]=f0[i]*exp(-j*f0[i]*conj(f0[i])*dt);
    
}
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
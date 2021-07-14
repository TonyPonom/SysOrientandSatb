
#include "functionsl3.h"
#include <fstream>
#include <math.h>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;


void main()
{
//исходные данные
double r[9];
r[0]=4226.800251;// rx
r[1]=3085.944251;// ry
r[2]=-4321.376266;// rz
r[3]=0.593397; // Vx
r[4]=5.793711; // Vy
r[5]=4.948645; // Vz
r[6]=1/100000;//0.0174;//pow(10.0,-5); // omx
r[7]=1/1000000;//1.0000000000000001e-006; // omy
r[8]=-1/100000;//-pow(10.0,-5); // omz = 
double A[3][3];
A[0][0]=-0.132203033235;
A[0][1]=0.698974798902;
A[0][2]=0.702820452535;
A[1][0]=-0.966003954060;
A[1][1]=0.068068924482;
A[1][2]=-0.24940525707;
A[2][0]=-0.22216822172;
A[2][1]=-0.71189946763;
A[2][2]=0.66621350124;

//A[0][0]=1;
//A[0][1]=0;
//A[0][2]=0;
//A[1][0]=0;
//A[1][1]=1;
//A[1][2]=0;
//A[2][0]=0;
//A[2][1]=0;
//A[2][2]=1;
double OMmax=628; 
double Lmax=0.01;
double J[3][3];
double Jx,Jy,Jz;
Jx=1;
Jy=2;
Jz=1;
J[0][0]=Jx;J[1][1]=Jy;J[2][2]=Jz;
J[0][1]=0;J[0][2]=0;J[1][2]=0;
J[1][0]=0;J[2][0]=0;J[2][1]=0;
double Idm=1.5/10000;
double Kp,Kd;
Kp=0.01;//1;
Kd=0.215;//2;
double sec,min,hour;
sec=21;
min=40;
hour=10;
double t=0;//с
double dt=0.1;// 100 мс
double tend=500;
int count;
double fi[3];
double sigm[3];
double OM[3] = {0,0,0};
ofstream fout1("res.txt");
fout1 << "\n";
fout1  << setw( 12 ) << setiosflags( ios::right ) << "t ";
fout1  << setw( 12 ) << setiosflags( ios::right ) << "rx ";
fout1  << setw( 12 ) << setiosflags( ios::right ) << "ry ";
fout1  << setw( 12 ) << setiosflags( ios::right ) << "rz ";
fout1  << setw( 12 ) << setiosflags( ios::right ) << "omx ";
fout1  << setw( 12 ) << setiosflags( ios::right ) << "omy ";
fout1  << setw( 12 ) << setiosflags( ios::right ) << "omz ";
fout1  << setw( 12 ) << setiosflags( ios::right ) << "Omx ";
fout1  << setw( 12 ) << setiosflags( ios::right ) << "Omy ";
fout1  << setw( 12 ) << setiosflags( ios::right ) << "Omz ";
fout1  << setw( 12 ) << setiosflags( ios::right ) << "fix ";
fout1  << setw( 12 ) << setiosflags( ios::right ) << "fiy ";
fout1  << setw( 12 ) << setiosflags( ios::right ) << "fiz "<<  endl;
fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << t ;
fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << r[0];
fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << r[1];
fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << r[2];
fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << r[6]*tograd;
fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << r[7]*tograd;
fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << r[8]*tograd;
fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << OM[0]*tograd;
fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << OM[1]*tograd;
fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << OM[2]*tograd;
// перходим кк углам Ёйлера
	fi[0]=atan2(-A[1][2], A[1][1]); //fix
    fi[1]=atan2(-A[2][0], A[0][0]); //fiy
    fi[2]=asin(A[1][0]); //fiz
fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << fi[0]*tograd;
fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << fi[1]*tograd;
fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << fi[2]*tograd << endl;
double h = dt;
while (t<=tend){
	//введение нелинейности типа "насыщени€"
	for ( count = 0; count < 3; count++){
		sigm[count]=-Kp*fi[count]-Kd*r[6+count];
		
		if (abs(sigm[count]) >= Lmax){
			sigm[count]= Lmax*(sign(sigm[count]));
		}
		if (abs(OM[count]) >= OMmax){
			sigm[count]= 0;
		}
		
	}
	//sigm[0]= 0;
	//sigm[1]= 0;
	//sigm[2]= 0;
	//
	double Msum[3];
	Ms(r,A,J,Msum );
	//Msum[0]= 0;
	//Msum[1]= 0;
	//Msum[2]= 0;
	double Mgr[3][3];
	getMgr( sec, min, hour, Mgr);
	double  a_grav[3],a_sun[3];
	getagrav( r, Mgr, a_grav);
	getasun(r,  sec, min, hour, a_sun);
	//RK4OneStepA(RightPartA, t,A,  h,  r);
	//RK4OneStepOm (RightPartOm, t,OM,  h,  Idm,  sigm );
	RK4OneStep (RightPart, t,r,  h, OM,Msum, A, J, a_sun, a_grav, Idm,  sigm );
	double Om[3];
	Om[0]=OM[0]*tograd/10;Om[1]=OM[1]*tograd/10;Om[2]=OM[2]*tograd/10;
	// перходим кк углам Ёйлера
	fi[0]=atan2(-A[1][2], A[1][1]); //fix
    fi[1]=atan2(-A[2][0], A[0][0]); //fiy
    fi[2]=asin(A[1][0]); //fiz
	t+= dt;
	// выводим в файл
	fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << t ;
	fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << r[0];
	fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << r[1];
	fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << r[2];
	fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << r[6]*tograd;
	fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << r[7]*tograd;
	fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << r[8]*tograd;
	fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << Om[0];
	fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << Om[1];
	fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << Om[2];
	fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << fi[0]*tograd;
	fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << fi[1]*tograd;
	fout1  << setw( 12 ) << setiosflags( ios::right )  << setprecision( 3 ) << fixed << fi[2]*tograd << endl;
}
fout1.close();
}
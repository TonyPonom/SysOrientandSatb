#pragma once

#include <math.h>
#include <stdio.h>
#include <string>
using namespace std;


#define PI 3.14159265359
#define PI2 6.28318530718

const double mu= 398600.4415;//гравитационаая постоянная
const double Re = 6371.0;// радиус земли
const double tograd = 180/PI;//перевод в градусы
const double torad = PI/180;//первод в радианы
const double Wer=0.7292115e-4;//угловая скорость вращения земли
const double g0 = 9.80665;
const double OMmax=628;
const double Lmax=0.01;;

double sign(double a){
	if (a > 0) {
		return 1;
	}
	if (a == 0) {
		return 0;
	}
	if (a < 0) {
		return -1;
	}
}
double norm(double vector[3] ) {	//длинна вектора
	return( sqrt (vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]));
};
void transpM(double M[3][3],double MT[3][3]){
	double c;
	MT[0][0]=M[0][0]; MT[0][1]=M[0][1]; MT[0][2]=M[0][2];
	MT[1][0]=M[1][0]; MT[1][1]=M[1][1]; MT[1][2]=M[1][2];
	MT[2][0]=M[2][0]; MT[2][1]=M[2][1]; MT[2][2]=M[2][2];
	for (int count = 0; count < 3; count++){
		for (int count2 = count+1; count2 < 3; count2++){
			c = MT[count][count2];
			MT[count][count2] = MT[count2][count];
			MT[count2][count] = c;
		}
	}
}
void vectpr(double x[3], double y[3], double z[3] ) {	//векторное произведение
	z[0]=x[1]*y[2]-x[2]*y[1];
	z[1]=-x[0]*y[2]+x[2]*y[0];
	z[2]=x[0]*y[1]-x[1]*y[0];
};
void Ms(double r[9], double A[3][3], double J[3][3],double Ms[3] ) {	//внешнего возмущающего гравитационного момента
	double Mu = 3.986004415*pow(10.0,14); // гравитационный параметр Земли [м^3/c^2]
	double AT[3][3];
	transpM(A,AT);
	double f[3];
	f[0]=r[0]*1000;f[1]=r[1]*1000;f[2]=r[2]*1000;// [м]
	double l = norm(f);// [м]
	double eta[3],jeta[3],v[3];
	//for (int count = 0; count < 3; count++){
	//	eta[count]=AT[count][0]*(f[0]/l) + AT[count][1]*(f[1]/l) +AT[count][2]*(f[2]/l);
	//}
	eta[0]=AT[0][0]*(f[0]/l) + AT[0][1]*(f[1]/l) +AT[0][2]*(f[2]/l);
	eta[1]=AT[1][0]*(f[0]/l) + AT[1][1]*(f[1]/l) +AT[1][2]*(f[2]/l);
	eta[2]=AT[2][0]*(f[0]/l) + AT[2][1]*(f[1]/l) +AT[2][2]*(f[2]/l);
	//for (int count = 0; count < 3; count++){
	//	jeta[count]= J[count][0] * eta[0] + J[count][1] * eta[1] + J[count][2] * eta[2];
	//}
	jeta[0]= J[0][0] * eta[0] + J[0][1] * eta[1] + J[0][2] * eta[2];
	jeta[1]= J[1][0] * eta[0] + J[1][1] * eta[1] + J[1][2] * eta[2];
	jeta[2]= J[2][0] * eta[0] + J[2][1] * eta[1] + J[2][2] * eta[2];
	vectpr(eta, jeta, v );
	for (int count = 0; count < 3; count++){
		Ms[count]= (3*Mu/(2*pow(l,3)) * v[count])/1000;
	}
};
void M3x3xM3x3(double M[3][3],double M1[3][3], double Mres[3][3]) {	//произведение двух матриц
	Mres[0][0]=M[0][0]*M1[0][0]+M[0][1]*M1[1][0]+M[0][2]*M1[2][0];
	Mres[0][1]=M[0][0]*M1[0][1]+M[0][1]*M1[1][1]+M[0][2]*M1[2][1];
	Mres[0][2]=M[0][0]*M1[0][2]+M[0][1]*M1[1][2]+M[0][2]*M1[2][2];
	Mres[1][0]=M[1][0]*M1[0][0]+M[1][1]*M1[1][0]+M[1][2]*M1[2][0];
	Mres[1][1]=M[1][0]*M1[0][1]+M[1][1]*M1[1][1]+M[1][2]*M1[2][1];
	Mres[1][2]=M[1][0]*M1[0][2]+M[1][1]*M1[1][2]+M[1][2]*M1[2][2];
	Mres[2][0]=M[2][0]*M1[0][0]+M[2][1]*M1[1][0]+M[2][2]*M1[2][0];
	Mres[2][1]=M[2][0]*M1[0][1]+M[2][1]*M1[1][1]+M[2][2]*M1[2][1];
	Mres[2][2]=M[2][0]*M1[0][2]+M[2][1]*M1[1][2]+M[2][2]*M1[2][2];
};

double JulianDat(double sec,double min,double hour){
	double DayJD,MonthJD,YearJD; 
	DayJD = 22;
    MonthJD = 5;
    YearJD = 2018;
	// вычисление промежуточных коэффициентов
	double aJd,yJd,mJd;
    aJd = floor((14 - MonthJD)/12); 
    yJd = YearJD + 4800 - aJd;
    mJd = MonthJD + 12*aJd - 3;
	// вычисление номера юлианского дня
    double Jdn = DayJD + floor((153*mJd + 2)/5) + 365*yJd + floor(yJd/4) - floor(yJd/100) + floor(yJd/400) - 32045;
    // вычисление юлианской даты
    return(Jdn + (hour - 12)/24 + min/1440 + sec/86400);
}
void getMgr( double sec,double min,double hour,double Mgr[3][3]){//Матрица перехода от J2000 в гринвичскую СК
	double Mu = 398600.4415;
    
    //вычисление гринвичского угла
    double A = -19089.451590;
    double B = 8640184.812866;
    double C = 0.093104;
    double D = -6.2*pow(10.0,-6);
    double JD0 = 2451545; 
    double JDD = 36525; 
    double DS2R = 7.272205216643039903848712*pow(10.0,-5);
    
    double JD = JulianDat(sec, min, hour);
    double t = (JD - JD0)/JDD;
    double f = 86400 *fmod(JD , 1);
    
    double alfa = DS2R * ((A + (B + (C + D*t)*t)*t) + f); 
    alfa = fmod(alfa , PI2);    
    if (alfa < 0){
             alfa = alfa + PI2;
	}
	Mgr[0][0]=cos(alfa); Mgr[0][1]=sin(alfa); Mgr[0][2]=0;
	Mgr[1][0]=-sin(alfa); Mgr[1][1]=cos(alfa); Mgr[1][2]=0;
	Mgr[2][0]=0; Mgr[2][1]=0; Mgr[2][2]=1;

}
void getagrav(double r[9],double Mgr[3][3],double a_grav[3]){
	double re= 6.378137*pow(10.0,6);// средний радиус Земли [м]
    double Mu = 3.986004418*pow(10.0,14); // графитационная постоянная Земли
    double J2 = 1.08262668355*pow(10.0,-3);// коффициент второй зональной гармоники
    double rad[3];
	double a[3];
	double c=0;
	double MgrT[3][3];
	transpM( Mgr,MgrT);
	for (int count = 0; count < 3; count++){
		rad[count] = Mgr[count][0]*r[0]*1000 + Mgr[count][1]*r[1]*1000 + Mgr[count][2]*r[2]*1000 ; //радиус-вектор во вращающейся гринвичской СК в м
	}
	for (int count = 0; count < 3; count++){
		if (count==2){
			c=2;
		}
		// Возмущающее ускорение в гринвичской вращающейся СК - (13)
		a[count] = (-1.5) * J2 * (Mu/pow(norm(rad),2)) * pow((re/norm(rad)),2) * (( 1 + c - 5*pow((rad[2]/norm(rad)),2))*(rad[count]/norm(rad)));
        //в конце делим на 1000 чтобы перейти к [км/с2] 
		//транспонируем и умножаем
		
	}
	for (int count = 0; count < 3; count++){
		a_grav[count] = (MgrT[count][0]*a[0] + MgrT[count][1]*a[1] + MgrT[count][2]*a[2])/1000 ;
	}
   
}
void getasun(double r[9], double sec,double min,double hour,double a_sun[3]){
	double OmandW = 282.940*torad; // долгота восх. узла + арг. перицентра
    double musun = 132712517951; // гравитационный параметр Солнца
    double epsisun = 23.43929111*torad; // наклонение плоскости эклиптики
    
    double JD = JulianDat(sec, min, hour);
    double Tsun = (JD - 2451545.0)/36525; // модифицированная юлианская дата
    double Msun = (357.5226 + 35999.049*Tsun)*torad; // средняя аномалия
    
    // Определяем вектор Земля-Солнце
    // 6892"/206264.8=0.033413359914
    // 72"/206264.8=0.0003490658609(сек в радианной мере)
    double lamdasun = OmandW + Msun + 0.033413359914*sin(Msun) + 0.0003490658609*sin(2*Msun);
    double rsunabs = (149.619 - 2.499*cos(Msun) - 0.021*cos(2*Msun))*pow(10.0,6);
    
    // Вектор rs
	double rs[3];
	rs[0] = rsunabs * cos(lamdasun);
	rs[1] = rsunabs * sin(lamdasun) * cos(epsisun);
	rs[2] = rsunabs * sin(lamdasun) * sin(epsisun);

    double lrs = norm(rs);
	double d[3];
	d[0]= rs[0] - r[0]; d[1]= rs[1] - r[1]; d[2]= rs[2] - r[2];
	double lrazv = norm(d);
    // ускорение
	for ( int count = 0; count < 3; count++){
		a_sun[count] = musun * (d[count]/pow(lrazv,3) - rs[count]/pow(lrs,3));
	}
}							//(double t, double fx[9],double dfx[9],double dOm[3],double a_sun[3],double a_grav[3],double Idm, double sigm[3],double A[3][3],double dA[3][3],double J[3][3],double Om[3],double Ms[3])
void RK4OneStep (void (*f)(double , double [] ,double [],double [] ,double [],double [],double ,double [],double [3][3],double [3][3],double [3][3],double [],double [] ), double t,double xx[9], double h,double Om[3],double Ms[3],double A[3][3],double J[3][3],double a_sun[3],double a_grav[3],double Idm, double sigm[3] ) { // Метод Рунге-Кутта 4-го порядка
    double K1[9],K2[9],K3[9],K4[9],ffx[9],O1[3],O2[3],O3[3],O4[3],A1[3][3],A2[3][3],A3[3][3],A4[3][3];
	double ffOM[3],fA[3][3];
	double h2;
	int count;
	h2=h/2;
	f(t,xx,K1,O1,a_sun, a_grav,Idm,sigm,A, A1, J,Om, Ms);
	for ( count = 0; count < 9; count++){
		ffx[count] = xx[count]+h2*K1[count];
		if (count<3){
			ffOM[count] = Om[count]+h2*O1[count];
		}
	}
	fA[0][0]=A[0][0]+h2*A1[0][0]; fA[0][1]=A[0][1]+h2*A1[0][1]; fA[0][2]=A[0][2]+h2*A1[0][2];
	fA[1][0]=A[1][0]+h2*A1[1][0]; fA[1][1]=A[1][1]+h2*A1[1][1]; fA[1][2]=A[1][2]+h2*A1[1][2];
	fA[2][0]=A[2][0]+h2*A1[2][0]; fA[2][1]=A[2][1]+h2*A1[2][1]; fA[2][2]=A[2][2]+h2*A1[2][2];
	f(t+h2,ffx,K2,O2,a_sun, a_grav,Idm,sigm,fA, A2, J,ffOM, Ms); 
	for (count = 0; count < 9; count++){
		ffx[count] = xx[count]+h2*K2[count];
		if (count<3){
			ffOM[count] = Om[count]+h2*O2[count];
		}
	}
	fA[0][0]=A[0][0]+h2*A2[0][0]; fA[0][1]=A[0][1]+h2*A2[0][1]; fA[0][2]=A[0][2]+h2*A2[0][2];
	fA[1][0]=A[1][0]+h2*A2[1][0]; fA[1][1]=A[1][1]+h2*A2[1][1]; fA[1][2]=A[1][2]+h2*A2[1][2];
	fA[2][0]=A[2][0]+h2*A2[2][0]; fA[2][1]=A[2][1]+h2*A2[2][1]; fA[2][2]=A[2][2]+h2*A2[2][2];
	f(t+h2,ffx,K3,O3,a_sun, a_grav,Idm,sigm,fA, A3, J,ffOM, Ms);
	for (count = 0; count < 9; count++){
		ffx[count] = xx[count]+h2*K3[count];
		if (count<3){
			ffOM[count] = Om[count]+h2*O3[count];
		}
	}
	fA[0][0]=A[0][0]+h2*A3[0][0]; fA[0][1]=A[0][1]+h2*A3[0][1]; fA[0][2]=A[0][2]+h2*A3[0][2];
	fA[1][0]=A[1][0]+h2*A3[1][0]; fA[1][1]=A[1][1]+h2*A3[1][1]; fA[1][2]=A[1][2]+h2*A3[1][2];
	fA[2][0]=A[2][0]+h2*A3[2][0]; fA[2][1]=A[2][1]+h2*A3[2][1]; fA[2][2]=A[2][2]+h2*A3[2][2];
	f(t+h,ffx,K4,O4,a_sun, a_grav,Idm,sigm,fA, A4, J,ffOM, Ms);
	for (count = 0; count < 9; count++){
		xx[count] = xx[count]+h*(K1[count]+K4[count]+2*(K2[count]+K3[count]))/6.;
		if (count<3){
			Om[count] = Om[count]+h*(O1[count]+O4[count]+2*(O2[count]+O3[count]))/6.;
		}
	}
	A[0][0]=A[0][0]+h*(A1[0][0]+A4[0][0]+2*(A2[0][0]+A3[0][0]))/6.;
	A[0][1]=A[0][1]+h*(A1[0][1]+A4[0][1]+2*(A2[0][1]+A3[0][1]))/6.;
	A[0][2]=A[0][2]+h*(A1[0][2]+A4[0][2]+2*(A2[0][2]+A3[0][2]))/6.;
	A[1][0]=A[1][0]+h*(A1[1][0]+A4[1][0]+2*(A2[1][0]+A3[1][0]))/6.;
	A[1][1]=A[1][1]+h*(A1[1][1]+A4[1][1]+2*(A2[1][1]+A3[1][1]))/6.;
	A[1][2]=A[1][2]+h*(A1[1][2]+A4[1][2]+2*(A2[1][2]+A3[1][2]))/6.;
	A[2][0]=A[2][0]+h*(A1[2][0]+A4[2][0]+2*(A2[2][0]+A3[2][0]))/6.;
	A[2][1]=A[2][1]+h*(A1[2][1]+A4[2][1]+2*(A2[2][1]+A3[2][1]))/6.;
	A[2][2]=A[2][2]+h*(A1[2][2]+A4[2][2]+2*(A2[2][2]+A3[2][2]))/6.;

	
	
}
void RightPart(double t, double fx[9],double dfx[9],double dOm[3],double a_sun[3],double a_grav[3],double Idm, double sigm[3],double A[3][3],double dA[3][3],double J[3][3],double Om[3],double Ms[3]){ //Правая часть уравнений
	 double X[3];
	 X[0] = fx[0];
	 X[1] = fx[1];
	 X[2] = fx[2];
	 double  r=norm(X), r3=r*r*r;
	 int count;
	 for (count = 0; count < 3; count++){
		dfx[count] = fx[count+3];
		dOm[count] = -pow(Idm, -1) *sigm[count];
	 }
	 for (count = 3; count < 6; count++){
		dfx[count] = -mu*fx[count-3]/r3 + a_sun[count-3] + a_grav[count-3];
	 }
	 double Jom[3],IOm[3],c[3],v[3],c1[3];
	 Jom[0]=J[0][0]*fx[6];Jom[1]=J[1][1]*fx[7];Jom[2]=J[2][2]*fx[8];
	 IOm[0]=Idm*Om[0];IOm[1]=Idm*Om[1];IOm[2]=Idm*Om[2];
	 c[0]=Jom[0]+IOm[0];c[1]=Jom[1]+IOm[1];c[2]=Jom[2]+IOm[2];
	 c1[0]=fx[6];c1[1]=fx[7];c1[2]=fx[8];
	 vectpr(c1, c, v );
	 for (count = 6; count < 9; count++){
		dfx[count] = (Ms[count-6]-v[count-6]+sigm[count-6]) ;
	 }
	 dfx[6] = dfx[6]/J[0][0] ;
	 dfx[7] = dfx[7]/J[1][1] ; 
	 dfx[8] = dfx[8]/J[2][2] ; 
	 double Aom[3][3];
	 Aom[0][0]=0; Aom[0][1]=-fx[8]; Aom[0][2]=fx[7];
	 Aom[1][0]=fx[8]; Aom[1][1]=0; Aom[1][2]=-fx[6];
	 Aom[2][0]=-fx[7]; Aom[2][1]=fx[6]; Aom[2][2]=0;
	 M3x3xM3x3(A,Aom, dA);
	 //dA[0][0]=0; dA[0][1]=0; dA[0][2]=0;
	 //dA[1][0]=0; dA[1][1]=0; dA[1][2]=0;
	 //dA[2][0]=0; dA[2][1]=0; dA[2][2]=0;

}
void RightPartOm(double t, double fx[3],double dfx[3],double Idm, double sigm[3]){ //Правая часть уравнений
	 for (int count = 0; count < 3; count++){
		dfx[count] = -sigm[count]/Idm;
	 }
}
void RK4OneStepOm (void (*f)(double , double [] ,double [],double  ,double [] ), double t,double xx[3], double h,double Idm, double sigm[3] ) { // Метод Рунге-Кутта 4-го порядка
    double K1[3],K2[3],K3[3],K4[3],ffx[3];
	double h2;
	int count;
	h2=h/2;
	f(t,xx,K1,Idm,sigm);
	for ( count = 0; count < 3; count++){
		ffx[count] = xx[count]+h2*K1[count];
	}
	f(t+h2,ffx,K2,Idm,sigm); 
	for (count = 0; count < 3; count++){
		ffx[count] = xx[count]+h2*K2[count];
	}
	f(t+h2,ffx,K3,Idm,sigm);
	for (count = 0; count < 3; count++){
		ffx[count] = xx[count]+h2*K3[count];
	}
	f(t+h,ffx,K4,Idm,sigm);
	for (count = 0; count < 3; count++){
		xx[count] = xx[count]+h*(K1[count]+K4[count]+2*(K2[count]+K3[count]))/6.;

	}

	
	
}
void RightPartA(double t, double A[3][3],double dA[3][3],double r[9]){ //Правая часть уравнений
	 double Aom[3][3];
	 Aom[0][0]=0; Aom[0][1]=-r[8]; Aom[0][2]=r[7];
	 Aom[1][0]=r[8]; Aom[1][1]=0; Aom[1][2]=-r[6];
	 Aom[2][0]=-r[7]; Aom[2][1]=r[6]; Aom[2][2]=0;
	 M3x3xM3x3(A,Aom, dA);
}
void RK4OneStepA(void (*f)(double , double [3][3] ,double [3][3],double [] ), double t,double A[3][3], double h,double r[9] ) { // Метод Рунге-Кутта 4-го порядка
    double K1[3][3],K2[3][3],K3[3][3],K4[3][3],fA[3][3];
	double h2;
	h2=h/2;
	f(t,A,K1,r);
	fA[0][0]=A[0][0]+h2*K1[0][0]; fA[0][1]=A[0][1]+h2*K1[0][1]; fA[0][2]=A[0][2]+h2*K1[0][2];
	fA[1][0]=A[1][0]+h2*K1[1][0]; fA[1][1]=A[1][1]+h2*K1[1][1]; fA[1][2]=A[1][2]+h2*K1[1][2];
	fA[2][0]=A[2][0]+h2*K1[2][0]; fA[2][1]=A[2][1]+h2*K1[2][1]; fA[2][2]=A[2][2]+h2*K1[2][2];
	f(t+h2,fA,K2,r); 
	fA[0][0]=A[0][0]+h2*K2[0][0]; fA[0][1]=A[0][1]+h2*K2[0][1]; fA[0][2]=A[0][2]+h2*K2[0][2];
	fA[1][0]=A[1][0]+h2*K2[1][0]; fA[1][1]=A[1][1]+h2*K2[1][1]; fA[1][2]=A[1][2]+h2*K2[1][2];
	fA[2][0]=A[2][0]+h2*K2[2][0]; fA[2][1]=A[2][1]+h2*K2[2][1]; fA[2][2]=A[2][2]+h2*K2[2][2];
	f(t+h2,fA,K3,r);
	fA[0][0]=A[0][0]+h2*K3[0][0]; fA[0][1]=A[0][1]+h2*K3[0][1]; fA[0][2]=A[0][2]+h2*K3[0][2];
	fA[1][0]=A[1][0]+h2*K3[1][0]; fA[1][1]=A[1][1]+h2*K3[1][1]; fA[1][2]=A[1][2]+h2*K3[1][2];
	fA[2][0]=A[2][0]+h2*K3[2][0]; fA[2][1]=A[2][1]+h2*K3[2][1]; fA[2][2]=A[2][2]+h2*K3[2][2];
	f(t+h,fA,K4,r);
	A[0][0]=A[0][0]+h*(K1[0][0]+K4[0][0]+2*(K2[0][0]+K3[0][0]))/6.;
	A[0][1]=A[0][1]+h*(K1[0][1]+K4[0][1]+2*(K2[0][1]+K3[0][1]))/6.;
	A[0][2]=A[0][2]+h*(K1[0][2]+K4[0][2]+2*(K2[0][2]+K3[0][2]))/6.;
	A[1][0]=A[1][0]+h*(K1[1][0]+K4[1][0]+2*(K2[1][0]+K3[1][0]))/6.;
	A[1][1]=A[1][1]+h*(K1[1][1]+K4[1][1]+2*(K2[1][1]+K3[1][1]))/6.;
	A[1][2]=A[1][2]+h*(K1[1][2]+K4[1][2]+2*(K2[1][2]+K3[1][2]))/6.;
	A[2][0]=A[2][0]+h*(K1[2][0]+K4[2][0]+2*(K2[2][0]+K3[2][0]))/6.;
	A[2][1]=A[2][1]+h*(K1[2][1]+K4[2][1]+2*(K2[2][1]+K3[2][1]))/6.;
	A[2][2]=A[2][2]+h*(K1[2][2]+K4[2][2]+2*(K2[2][2]+K3[2][2]))/6.;

	
	
}

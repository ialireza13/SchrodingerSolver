#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <stdlib.h>
#include <string.h>
using namespace std;
typedef complex<double> CD;

#define print_potential

double hbar = 1.05457266e-34;
double m_e = 9.10938356e-31;
double _sqrt_2mh = 1.27992151207911e19;
double Pi = acos(-1);
const CD ii(0, 1);

int main(){
	cout<<hbar*hbar*pow(9,4)/2.0/m_e<<endl;
	cout<<25.0*hbar*hbar*pow(9,4)/8.0/m_e<<endl;
	cout<<sqrt(2.0*m_e*4.94446e-37)/hbar*2.0/pow(3,4)<<endl;
	int i,j,l;
	char filename[24];
	strcpy(filename,"potential");
	int N = 0; double data;
	ifstream openn;
	openn.open(filename);
	if(openn.fail()){
		cout<<"!!!!! Error Opening File.\n";exit(1);
	}
	openn.close();
	ifstream inFile;
	inFile.open(filename);
	do{
		// d(i) >> V(i)
		inFile>>data>>data;
		N++;
	}while((!inFile.eof()));
	inFile.close();
	cout<<"Number of barriers= "<<N<<endl;
	double d[N], V[N];
	ifstream inFile1;
	inFile1.open(filename);
	for(i=0;i<N;i++)
		inFile1>>d[i]>>V[i];
	inFile1.close();
	#ifdef print_potential
		double X_m=0,X_M=0;
		for(i=0;i<N;i++){
			X_M += d[i];
			cout<<V[i]<<"\t\t\t"<<X_m<<"<x<="<<X_M<<endl;
			X_m = X_M;
		}
	#endif
	CD D[2][2][N], P[2][2][N];
	CD iD[2][2][N], det;
	int n = int(N/2.0);
	CD M[2][2], MM[2][2][2][N-1], MMM[2][2], k[N-1];
	double V_max=0, E, R, T, dE=1e-37;
	for(i=0;i<N;i++)
		if(fabs(V[i])>V_max)V_max=fabs(V[i]);
	ofstream Tout("ET");
	ofstream Rout("ER");
	double E_min, E_max;
	cout<<"Enter E_min: "; cin>>E_min;
	cout<<"Enter E_max: "; cin>>E_max;
	cout<<"Enter Step: "; cin>>dE;
	#ifdef print_energy
		cout<<"E= "<<V_max+dE;
	#endif
	for(E=E_min+dE; E<=E_max; E+=dE){
		for(i=0;i<N;i++){
			k[i] = sqrt(2.0*m_e*fabs(E-V[i]))/hbar;
		}
		for(i=0;i<N;i++){
			D[0][0][i] = 1.0;
			D[0][1][i] = 1.0;
			D[1][0][i] = k[i];
			D[1][1][i] = -k[i];
			det = D[0][0][i]*D[1][1][i] - D[0][1][i]*D[1][0][i];
			iD[0][0][i] = -k[i]/det;
			iD[0][1][i] = -1.0/det;
			iD[1][0][i] = -k[i]/det;
			iD[1][1][i] = 1.0/det;
		}
		for(i=0;i<N;i++){
			P[0][0][i] = exp(ii*k[i]*d[i]);
			P[0][1][i] = 0;
			P[1][0][i] = 0;
			P[1][1][i] = exp(-ii*k[i]*d[i]);
		}
		for(i=1;i<N-1;i++){
			//Di * Pi
			MM[0][0][0][i]= D[0][0][i]*P[0][0][i] + D[0][1][i]*P[1][0][i];
			MM[0][0][1][i]= D[0][0][i]*P[0][1][i] + D[0][1][i]*P[1][1][i];
			MM[0][1][0][i]= D[1][0][i]*P[0][0][i] + D[1][1][i]*P[1][0][i];
			MM[0][1][1][i]= D[1][0][i]*P[0][1][i] + D[1][1][i]*P[1][1][i];
			// Di*Pi * iDi
			MM[1][0][0][i]= MM[0][0][0][i]*iD[0][0][i] + MM[0][0][1][i]*iD[1][0][i];
			MM[1][0][1][i]= MM[0][0][0][i]*iD[0][1][i] + MM[0][0][1][i]*iD[1][1][i];
			MM[1][1][0][i]= MM[0][1][0][i]*iD[0][0][i] + MM[0][1][1][i]*iD[1][0][i];
			MM[1][1][1][i]= MM[0][1][0][i]*iD[0][1][i] + MM[0][1][1][i]*iD[1][1][i];
		}

		//M1 * M2
		if(N>3){
			MMM[0][0]= MM[1][0][0][1]*MM[1][0][0][2] + MM[1][0][1][1]*MM[1][1][0][2];
			MMM[0][1]= MM[1][0][0][1]*MM[1][0][1][2] + MM[1][0][1][1]*MM[1][1][1][2];
			MMM[1][0]= MM[1][1][0][1]*MM[1][0][0][2] + MM[1][1][1][1]*MM[1][1][0][2];
			MMM[1][1]= MM[1][1][0][1]*MM[1][0][1][2] + MM[1][1][1][1]*MM[1][1][1][2];
		} else{
			MMM[0][0]= MM[1][0][0][1];
			MMM[0][1]= MM[1][0][1][1];
			MMM[1][0]= MM[1][1][0][1];
			MMM[1][1]= MM[1][1][1][1];
		}
		if(N>4){
			for(i=3; i<N-1; i++){
				MM[0][0][0][1]= MMM[0][0]*MM[1][0][0][i] + MMM[0][1]*MM[1][1][0][i];
				MM[0][0][1][1]= MMM[0][0]*MM[1][0][1][i] + MMM[0][1]*MM[1][1][1][i];
				MM[0][1][0][1]= MMM[1][0]*MM[1][0][0][i] + MMM[1][1]*MM[1][1][0][i];
				MM[0][1][1][1]= MMM[1][0]*MM[1][0][1][i] + MMM[1][1]*MM[1][1][1][i];
				MMM[0][0] = MM[0][0][0][1];
				MMM[0][1] = MM[0][0][1][1];
				MMM[1][0] = MM[0][1][0][1];
				MMM[1][1] = MM[0][1][1][1];
			}
		}
		//M2 * DN
			MM[0][0][0][1]= MMM[0][0]*D[0][0][N-1] + MMM[0][1]*D[1][0][N-1];
			MM[0][0][1][1]= MMM[0][0]*D[0][1][N-1] + MMM[0][1]*D[1][1][N-1];
			MM[0][1][0][1]= MMM[1][0]*D[0][0][N-1] + MMM[1][1]*D[1][0][N-1];
			MM[0][1][1][1]= MMM[1][0]*D[0][1][N-1] + MMM[1][1]*D[1][1][N-1];
		//iD0 * M2*DN
			M[0][0]= iD[0][0][0]*MM[0][0][0][1] + iD[0][1][0]*MM[0][1][0][1];
			M[0][1]= iD[0][0][0]*MM[0][0][1][1] + iD[0][1][0]*MM[0][1][1][1];
			M[1][0]= iD[1][0][0]*MM[0][0][0][1] + iD[1][1][0]*MM[0][1][0][1];
			M[1][1]= iD[1][0][0]*MM[0][0][1][1] + iD[1][1][0]*MM[0][1][1][1];
		R = ((M[1][0]*conj(M[1][0]))/(M[0][0]*conj(M[0][0]))).real();
		T = (k[N-1]/(k[0]*M[0][0]*conj(M[0][0]))).real();
		//cout<<"T = "<<T<<endl;
		//cout<<"R = "<<R<<endl;
		//cout<<"\t R+T = "<<R+T<<endl;
		Tout<<sqrt(2.0*m_e*E)/hbar*2.0/pow(3,4)<<'\t'<<T<<endl;
		//Rout<<E/V_max<<'\t'<<R<<endl;
		Rout<<sqrt(2.0*m_e*E)/hbar*2.0/pow(3,4)<<'\t'<<R<<endl;
		#ifdef print_energy
			cout<<"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bE= "<<E;
		#endif
	}
	#ifdef print_energy
		cout<<"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b                                    "<<endl;
	#endif
	Tout.close();
	Rout.close();
	return 0;
}
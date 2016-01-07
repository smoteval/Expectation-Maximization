#include <iostream>
#include <map>
#include <fstream>
#include <vector>
#include <math.h>
#include <time.h>
#include <cmath>
using namespace std;



int Train_Slop[2000];
int Train_Forienn[2000];
int Train_Degar[2000];
int Train_Genes[2000];
int Train_Dunnet[2000];

int Test_Slop[100];
int Test_Forienn[100];
int Test_Degar[100];
int Test_Genes[100];
int Test_Dunnet[100];

double Dunnet[3];
double Gene[2];
double Sloepnea[3][2][2];
double Foriennditis[3][2];
double Degar_Spots[3][2];



double SUMS[2][2][2][2][3];


double mean(double x[], int size)
{
double sum=0;
double mean2;
for ( int i = 0; i <size; i++ )
{
sum += x[i];
}
mean2=sum/ (double) size;
return mean2;
}

double deviation (double x[], int size, double mean)
{
double deviation;
double sum2=0;

for ( int i = 0; i < size; i++ )
{
sum2 += pow((x[i]-mean),2);
}
deviation= sqrt(sum2/(size));
return deviation;
} 


/*
  This function computes P(Dunnet = dunnet , everything else)
*/
double Prob_Dunnet(int Sloep, int Fori, int deg, int gene, int dunnet) {
	double P = Dunnet[dunnet];
	double P_F = Foriennditis[dunnet][Fori];
	double P_D = Degar_Spots[dunnet][deg];
	double P_G = Gene[gene];
	double P_S_Given_G = Sloepnea[dunnet][gene][Sloep];
	double result = P_F*P_D*P_G*P_S_Given_G*P;
	return result;
}



double SumOfStuff() {
	double sum = 0;
	for(int i = 0; i < 2000; i++) {
		if(Train_Dunnet[i] == 1 ) {
			sum = sum + 1;
			SUMS[Train_Slop[i]][Train_Forienn[i]][Train_Degar[i]][Train_Genes[i]][1] = 
			SUMS[Train_Slop[i]][Train_Forienn[i]][Train_Degar[i]][Train_Genes[i]][1] + 1;
		}
		else if( Train_Dunnet[i] == 2) {
			sum = sum + 1;
			SUMS[Train_Slop[i]][Train_Forienn[i]][Train_Degar[i]][Train_Genes[i]][2] = 
			SUMS[Train_Slop[i]][Train_Forienn[i]][Train_Degar[i]][Train_Genes[i]][2] + 1;
		}
		else if( Train_Dunnet[i] == 0 ) {
			sum = sum + 1;
			SUMS[Train_Slop[i]][Train_Forienn[i]][Train_Degar[i]][Train_Genes[i]][0] = 
			SUMS[Train_Slop[i]][Train_Forienn[i]][Train_Degar[i]][Train_Genes[i]][0] + 1;
		}
		else if(Train_Dunnet[i] == -1) {
			double no = Prob_Dunnet(Train_Slop[i], Train_Forienn[i], Train_Degar[i],Train_Genes[i],0);
			double mild = Prob_Dunnet(Train_Slop[i], Train_Forienn[i], Train_Degar[i],Train_Genes[i],1);
			double severe = Prob_Dunnet(Train_Slop[i], Train_Forienn[i], Train_Degar[i],Train_Genes[i],2);
			sum = no+mild+severe;
			double s = no+mild+severe;
			double P_NO = no/s;
			double P_MILD = mild/s;
			double P_SEVERE = severe/s; 
			SUMS[Train_Slop[i]][Train_Forienn[i]][Train_Degar[i]][Train_Genes[i]][0] = 
			SUMS[Train_Slop[i]][Train_Forienn[i]][Train_Degar[i]][Train_Genes[i]][0] + P_NO;
			SUMS[Train_Slop[i]][Train_Forienn[i]][Train_Degar[i]][Train_Genes[i]][1] = 
			SUMS[Train_Slop[i]][Train_Forienn[i]][Train_Degar[i]][Train_Genes[i]][1] + P_MILD;
			SUMS[Train_Slop[i]][Train_Forienn[i]][Train_Degar[i]][Train_Genes[i]][2] = 
			SUMS[Train_Slop[i]][Train_Forienn[i]][Train_Degar[i]][Train_Genes[i]][2] + P_SEVERE;
	    }
	}
	return sum;
}



double SumUp(int b1, int b2, int b3, int b4, int b5) {
	double sum = 0;
	for(int a1 = 0; a1 < 2; a1++) {
		for(int a2 = 0; a2 < 2; a2++) {
			for(int a3 = 0; a3 < 2; a3++) {
				for(int a4 = 0; a4< 2; a4++) {
					for(int a5 = 0; a5 < 3; a5++) {
						if( (a1 == b1 ||  b1 == -1) && (a2 == b2 ||  b2 == -1)  && (a3 == b3 ||  b3 == -1) && 
							(a4 == b4 ||  b4 == -1) && (a5 == b5 ||  b5 == -1)) {
							sum = sum + SUMS[a1][a2][a3][a4][a5];
						}
					}
				}
			}
		}
	}
	return sum;
}



void UpdateDunTable() {
     double mild = SumUp(-1,-1,-1,-1,1);
     double severe = SumUp(-1,-1,-1,-1,2);
     mild = (double) mild/2000;
     severe = (double) severe/2000;
     Dunnet[0] = 1.0 - (mild+severe);
     Dunnet[1] = mild;
     Dunnet[2] = severe;
}


void UpdateGeneTable(){
	double yes = SumUp(-1,-1,-1,1,-1);
	yes = (double) yes/2000;
	Gene[1] = yes;
	Gene[0] = 1.0 - yes;
}


void UpdateSloepneaTable() {
	Sloepnea[0][0][1] = SumUp(1,-1,-1,0,0)/SumUp(-1,-1,-1,0,0);
	Sloepnea[1][0][1] = SumUp(1,-1,-1,0,1)/SumUp(-1,-1,-1,0,1);
	Sloepnea[2][0][1] = SumUp(1,-1,-1,0,2)/SumUp(-1,-1,-1,0,2);
	Sloepnea[0][1][1] = SumUp(1,-1,-1,1,0)/SumUp(-1,-1,-1,1,0);
	Sloepnea[1][1][1] = SumUp(1,-1,-1,1,1)/SumUp(-1,-1,-1,1,1);
	Sloepnea[2][1][1] = SumUp(1,-1,-1,1,2)/SumUp(-1,-1,-1,1,2);
	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 2; j++) {
			Sloepnea[i][j][0] = 1.0 - Sloepnea[i][j][1];
		}
	}
}


void UpdateForidennditisTable() {
	Foriennditis[0][1] = SumUp(-1,1,-1,-1,0)/SumUp(-1,-1,-1,-1,0);
	Foriennditis[1][1] = SumUp(-1,1,-1,-1,1)/SumUp(-1,-1,-1,-1,1);
	Foriennditis[2][1] = SumUp(-1,1,-1,-1,2)/SumUp(-1,-1,-1,-1,2);
	for(int i = 0; i < 3; i++) {
		Foriennditis[i][0] = 1.0 - Foriennditis[i][1];
	}
}


void UpdateDegarTable() {
	Degar_Spots[0][1] = SumUp(-1,-1,1,-1,0)/SumUp(-1,-1,-1,-1,0);
	Degar_Spots[1][1] = SumUp(-1,-1,1,-1,1)/SumUp(-1,-1,-1,-1,1);
	Degar_Spots[2][1] = SumUp(-1,-1,1,-1,2)/SumUp(-1,-1,-1,-1,2);
	for(int i = 0; i < 3; i++) {
		Degar_Spots[i][0] = 1.0 - Degar_Spots[i][1];
	}	
}


void MakeSumZero() {
	for(int a1 = 0; a1 < 2; a1++) {
		for(int a2 = 0; a2 < 2; a2++) {
			for(int a3 = 0; a3 < 2; a3++) {
				for(int a4 = 0; a4< 2; a4++) {
					for(int a5 = 0; a5 < 3; a5++) {
						SUMS[a1][a2][a3][a4][a5] = 0;
					}
				}
			}
		}
	}
}


double TestAccuracy() {
	double correct = 0;
	for(int i = 0; i < 100; i++) {
		double P0 = Prob_Dunnet(Test_Slop[i],Test_Forienn[i],Test_Degar[i],Test_Genes[i],0);
		double P1 = Prob_Dunnet(Test_Slop[i],Test_Forienn[i],Test_Degar[i],Test_Genes[i],1);
		double P2 = Prob_Dunnet(Test_Slop[i],Test_Forienn[i],Test_Degar[i],Test_Genes[i],2);
		int label;
		if(P0 >= P1 && P0 >= P2) {
			label = 0;
		}
		else if(P1 >= P0 && P1 >= P2) {
			label = 1;
		}
		else {
			label = 2;
		}
		if( label == Test_Dunnet[i] ) {
			correct = correct + 1;
		}
	}
	return correct/(double)100;
}


double fRand(double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return f * (fMax );
}


void setinitialvalue(double max) {
	double s;

	Dunnet[0] = 0.5 + fRand(max);
	Dunnet[1] = 0.25 + fRand(max);
	Dunnet[2] = 0.25 + fRand(max);
	s = (Dunnet[0]+Dunnet[1]+Dunnet[2]);
	Dunnet[0] = Dunnet[0]/s;
	Dunnet[1] = Dunnet[1]/s;
	Dunnet[2] = Dunnet[2]/s;



	Gene[0] = 0.9 + fRand(max);
	Gene[1] = 0.1 + fRand(max);
	s = (Gene[0]+Gene[1]);
	Gene[0] = Gene[0] / s;
	Gene[1] = Gene[1] / s;


	Sloepnea[0][0][0] = 0.95 +fRand(max);
	Sloepnea[0][0][1] = 0.05 +fRand(max);
	s = (Sloepnea[0][0][1] + Sloepnea[0][0][0]);
	Sloepnea[0][0][0] = Sloepnea[0][0][0] / s;
	Sloepnea[0][0][1] = Sloepnea[0][0][1] / s;

	Sloepnea[0][1][0] = 1.0 +fRand(max);
	Sloepnea[0][1][1] = 0.0 +fRand(max);
	s = (Sloepnea[0][1][0] + Sloepnea[0][1][1]);
	Sloepnea[0][1][0] = Sloepnea[0][1][0] / s;
	Sloepnea[0][1][1] = Sloepnea[0][1][1] / s;


	Sloepnea[1][0][0] = 0.5 +fRand(max);
	Sloepnea[1][0][1] = 0.5 +fRand(max);
	s = Sloepnea[1][0][0]+ Sloepnea[1][0][1];
	Sloepnea[1][0][0] = Sloepnea[1][0][0]/s;
	Sloepnea[1][0][1] = Sloepnea[1][0][1]/s;

	Sloepnea[2][0][0] = 0.5 +fRand(max);
	Sloepnea[2][0][1] = 0.5 +fRand(max);
	s = Sloepnea[2][0][0]+Sloepnea[2][0][1];
	Sloepnea[2][0][0] = Sloepnea[2][0][0]/s;
	Sloepnea[2][0][1] = Sloepnea[2][0][1]/s;



	Sloepnea[1][1][0] = 0.95 +fRand(max);
	Sloepnea[1][1][1] = 0.05 +fRand(max);
	s = Sloepnea[1][1][0]+ Sloepnea[1][1][1];
	Sloepnea[1][1][0] = Sloepnea[1][1][0]/s;
	Sloepnea[1][1][1] = Sloepnea[1][1][1]/s;

	Sloepnea[2][1][0] = 0.95 +fRand(max);
	Sloepnea[2][1][1] = 0.05 +fRand(max);	
	s = Sloepnea[2][1][0]+Sloepnea[2][1][1];
	Sloepnea[2][1][0] = Sloepnea[2][1][0]/s;
	Sloepnea[2][1][1] = Sloepnea[2][1][1]/s; 



	Foriennditis[0][0] = 0.95 +fRand(max);
	Foriennditis[0][1] = 0.05 +fRand(max);
	s = Foriennditis[0][0]+Foriennditis[0][1];
	Foriennditis[0][0] = Foriennditis[0][0]/s;
	Foriennditis[0][1] = Foriennditis[0][1]/s;


	Foriennditis[1][0] = 0.2 +fRand(max);
	Foriennditis[1][1] = 0.8 +fRand(max);
	s = Foriennditis[1][0]+ Foriennditis[1][1];
	Foriennditis[1][0] = Foriennditis[1][0]/s;
	Foriennditis[1][1] = Foriennditis[1][1]/s;

	Foriennditis[2][0] = 0.6 +fRand(max);
	Foriennditis[2][1] = 0.4 +fRand(max);
	s = Foriennditis[2][0]+Foriennditis[2][1];
	Foriennditis[2][0] = Foriennditis[2][0]/s;
	Foriennditis[2][1] = Foriennditis[2][1]/s;




	Degar_Spots[0][0] = 0.95 +fRand(max);
	Degar_Spots[0][1] = 0.05+ fRand(max);
	s = Degar_Spots[0][0]+Degar_Spots[0][1];
	Degar_Spots[0][0] = Degar_Spots[0][0]/s;
	Degar_Spots[0][1] = Degar_Spots[0][1]/s;


	Degar_Spots[1][0] = 0.6 +fRand(max);
	Degar_Spots[1][1] = 0.4 +fRand(max);
	s = Degar_Spots[1][0]+Degar_Spots[1][1];
	Degar_Spots[1][0] = Degar_Spots[1][0]/s;
	Degar_Spots[1][1] = Degar_Spots[1][1]/s;

	Degar_Spots[2][0] = 0.2 +fRand(max);
	Degar_Spots[2][1] = 0.8 +fRand(max);
	s = Degar_Spots[2][0]+Degar_Spots[2][1];
	Degar_Spots[2][0] = Degar_Spots[2][0]/s;
	Degar_Spots[2][1] = Degar_Spots[2][1]/s;


}




void Start() {
	double before = 1;
	double after = 2;
	while(after-before > 0.1 || after - before < -0.1) {
		MakeSumZero();
		before = SumOfStuff();
		UpdateDunTable();
		UpdateGeneTable();
		UpdateSloepneaTable();
		UpdateForidennditisTable();
		UpdateDegarTable();
		after = SumOfStuff();
	}
}


int main() {
	ifstream my;
	my.open("traindata.txt");
	for(int i = 0; i < 2000; i++) {
		int sloe, fori, degar, gen, dun;
		my >> sloe >> fori >> degar >> gen >> dun;
		Train_Slop[i] = sloe;
		Train_Forienn[i] = fori;
		Train_Degar[i] = degar;
		Train_Genes[i] = gen;
		Train_Dunnet[i] = dun;
	}
	my.close();
	my.open("testdata.txt");
	for(int i = 0; i < 100; i++) {
		int sloe, fori, degar, gen, dun;
		my >> sloe >> fori >> degar >> gen >> dun;
		Test_Slop[i] = sloe;
		Test_Forienn[i] = fori;
		Test_Degar[i] = degar;
		Test_Genes[i] = gen;
		Test_Dunnet[i] = dun;
	}
	my.close();
	/////////////////////
	time_t seconds;
    time(&seconds);
    srand((unsigned int) seconds);
	ofstream deltas,meansbfore,meansafter,sdbefore,sdafter;
	deltas.open("delta");
	meansbfore.open("menbefore");
	meansafter.open("meansafter");
	sdbefore.open("sdbefore");
	sdafter.open("sdafter");


	double add = (double) 4.0/21.0;
	double i = 0;
	while(i < 4) {
		double x[20];
		double y[20];
		for(int j = 0; j < 20; j++) {
			setinitialvalue(i);
			y[j] = TestAccuracy();
			Start();
			x[j] = TestAccuracy();

		}
		double m = mean(x,20);
		double sd = deviation(x,20,m);
		double mbefore = mean(y,20);
		double SDbefore = deviation(y,20,mbefore);
		deltas << i << endl;
		meansbfore  << mbefore << endl;
		meansafter <<  m << endl;
		sdbefore << SDbefore << endl;
		sdafter << sd << endl;
		i = i + add;
	}
	deltas.close();
	deltas.close();
	meansbfore.close();
	meansafter.close();
	sdbefore.close();
	sdafter.close();
}

 

// #include "SingleCell/NeonatalRatAtria.cc"
// #include "SingleCell/RatSAN.cc"
//#include "SingleCell/HumanVentricle_healthy.cpp"
#include "SingleCell/HumanVentricle_HF.cpp"
// #include "SingleCell/TPORd.cc"
// #include "SingleCell/ORd.cc"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include "omp.h"

using namespace std;

// #define SAN
// #define ATRIA
#define VENT
// #define HETEROGENEITY 

int main(int argc, char *argv[])
{
	// --------user configuration list--------
	double dx = 0.15; // mm
	double dt = 0.02; // ms
	int sanCellNum = 0;
	int atrialCellNum = 0;
	int epiCellNum = 0;
	int mCellNum = 0;
	int endoCellNum = 0;
	#if defined(SAN) || defined(HETEROGENEITY)
	sanCellNum = 50; 
	#endif
	#if defined(ATRIA) || defined(HETEROGENEITY) 
	atrialCellNum = 50;
	#endif
	#if defined(VENT)	
	endoCellNum = 25;
	mCellNum = 35;
	epiCellNum = 40;
	#endif
	// following parameters are for 1D or 2D.3D
	double atrialCoeff = 1.0*0.0195; 
	double sanCoeff = 0.1*0.0195; 
//	double ventCoeff = 0.1275; 
	double ventCoeff = 0.1275*0.65; // HF/HF_CO


	// for atria
	#ifdef ATRIA
	int numS1 = 5;
	double BCL = 1000; // ms
	double stopTime = numS1*BCL; //ms
	double stimStrength = -10.0;//8.78; //8.78;//-8.78; // pA
	double stimDuration = 5.0;	// ms
	double stimStart = 0.0; // ms 
	#endif
	
	// for SAN
	#ifdef SAN
	int numS1 = 50;
	double BCL = 300; // ms
	double stopTime = 2000; //ms
	double stimStrength = -0.0; // pA
	double stimDuration = 0.0;	// ms
	double stimStart = 0.0; // ms
	#endif

	#ifdef VENT
	int numS1 = 10;
	double BCL = 1500; // ms 
	double stopTime = numS1*BCL + BCL; //ms
	double stimStrength = -80.0;
	double stimDuration = 1.0;
	double stimStart = 0.0;
	#endif


	#ifdef HETEROGENEITY
	int numS1 = 2;
	double BCL = 1000;
	double stopTime = 1500; //ms
	double stimStrength = -0.0; // pA
	double stimDuration = 0.0;	// ms
	double stimStart = 0.0; // ms 
	#endif

	double cvStartTime = 0;
	double cvEndTime = 0;
	int cvStartFlag = 0; 
	int cvEndFlag = 0;
	double cv = 0;

	int coreNum = 8;
	omp_set_num_threads(2 * coreNum);

	int cellNum = sanCellNum + atrialCellNum + epiCellNum + mCellNum + endoCellNum; 
	typedef CellType* CellPointer;
	ORdHumanVentricle* strand[cellNum];
	double coeff[cellNum]; // 
	double dcoeff_dx[cellNum]; 
	double oldV[cellNum];

	for(int i = 0; i < cellNum; i++)
	{
		#ifdef VENT //
		if(i == 59)
			coeff[i] = 0.2*ventCoeff;
		else
			coeff[i] = ventCoeff;
		#endif
	}

	for(int i = 0; i < cellNum; i++)
	{
		if (i == 0) 
			dcoeff_dx[i] = (coeff[i+1] - coeff[i])/dx;
		else if (i == cellNum-1) 
			dcoeff_dx[i] = (coeff[i] - coeff[i-1])/dx;
		else
			dcoeff_dx[i] = (coeff[i+1] - coeff[i-1])/(2.0*dx);
	}

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < cellNum; i++)
	{
		# ifdef HETEROGENEITY
		if(i < sanCellNum)
		{
			strand[i] = new RatSAN();

		}	
		else
		{
			strand[i] = new NeonatalRatAtria();
			FILE *initfile = fopen("SingleCell/NeonatalRatAtriaInitialValues.dat","r");
			strand[i]->readinAllStates(initfile); 
			fclose(initfile);
		}
		#endif



		#ifdef VENT
		FILE *initfile;
		if(i >= 0 && i < endoCellNum)
		{
			strand[i] = new ORdHumanVentricle(ENDO);
			fclose(initfile);
		}
		else if (i < endoCellNum + mCellNum)
		{
			strand[i] = new ORdHumanVentricle(M);
			fclose(initfile);
		}
		else // i < total cellnum
		{
			strand[i] = new ORdHumanVentricle(EPI);
			fclose(initfile);
		}	
		#endif

		strand[i]->setDt(dt);
	}

	#ifdef ATRIA
	FILE *datafile = fopen("Outputs/AtriaOneDResults.dat","w+");
	#endif

	#ifdef SAN
	FILE *datafile = fopen("Outputs/SANOneDResults.dat","w+");
	#endif
	
	#ifdef VENT
	FILE *datafile = fopen("Outputs/VentOneDResults.dat","w+");
	#endif

	#ifdef HETEROGENEITY
	FILE *datafile = fopen("Outputs/HeterOneDResults.dat","w+");
	#endif

	double time = 0;
	int step = 0;
	for(time = 0.0, step = 0; time <= stopTime; time += dt, step++)
	{
		// 1. progress stats
		if(step%25000 == 0) // 25000 * dt ms = 0.125s 
			cout << "Progress = " << 100.0*time/stopTime << "\%." << endl;

		for(int i = 0; i < cellNum; i++)
		{
			oldV[i] = strand[i]->getV();
		}

		
		for(int i = 0; i < cellNum; i++)
		{

			strand[i]->setIstim(0.0);
			if(time - floor(time/BCL)*BCL >= stimStart && 
		   	   time - floor(time/BCL)*BCL < stimStart + stimDuration)
			{
		    	if(i < 3 && i >= 0)
				{
					strand[i]->setIstim(stimStrength);
				}
			}

			// ---------calculate diffusion, i.e. dVgap---------
		
			double dVgap_dt = 0;
			double first_order;
			double second_order;

			if(i == 0) 
			{
				first_order = (oldV[i+1] - oldV[i])/(1.0*dx);
				second_order = (oldV[i+1] + oldV[i] - 2.0*oldV[i])/(dx*dx);
			}
			else if(i > 0 && i < cellNum - 1) 
			{
				first_order = (oldV[i+1] - oldV[i-1])/(2.0*dx);
				second_order = (oldV[i+1] + oldV[i-1] - 2.0*oldV[i])/(dx*dx);	
			}
			else if(i == cellNum - 1)
			{
				first_order = (oldV[i] - oldV[i-1])/(1.0*dx);
				second_order = (oldV[i] + oldV[i-1] - 2.0*oldV[i])/(dx*dx);	
			}

			dVgap_dt = dcoeff_dx[i]*first_order + coeff[i]*second_order;

			strand[i]->setDVgap_dt(dVgap_dt);
				
			strand[i]->update();
		}


			if(step%10 == 0) // 50*dt = 1 ms once
			{
				for(int j = 0; j < cellNum; j++)
				{
					if(j == 0)
						fprintf(datafile,"%4.10f\t", time);

					fprintf(datafile,"%4.10f\t", strand[j]->getV()); // unit: mV

					if(j == cellNum - 1)
						fprintf(datafile,"\n");
				}
			}


		if (floor(time/BCL) == numS1)
		{ 
			if(strand[10]->getV() >= -30 && cvStartFlag == 0)
			{
				cvStartTime = time;
				cout << "start = " << cvStartTime << endl;
				cvStartFlag = 1;
			}
			if(strand[90]->getV() >= -30 && cvEndFlag == 0)
			{
				cvEndTime = time;
				cout << "end = " << cvEndTime << endl;
				cvEndFlag = 1;
				cv = (dx * 80) / (cvEndTime - cvStartTime);
				cout << "duration = " << cvEndTime - cvStartTime << endl;
			}
		} 

	}
	fclose(datafile);

	if(cvStartFlag == 1 && cvEndFlag == 1)
		cout << "CV = " << cv << " m/s." << endl;
	else
		cout << "Conduction failure!" << endl;
	printf("All done.\n");

	return 0;
}
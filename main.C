#include <iostream>
#include "signal.h"
#include <fstream>

using namespace std;
using namespace std::chrono;

int main(){
	std::cout<<"code6 running..."<<std::endl;
    int N = 1024*1024*128;

    double params[2];//[0] is weight, [1] is Xe mass number 

	Signalcalc signal;

	/*this is the code for binned band hist filling*/
	
	int xBinning = 100;
    double xMin = 0.;
    double xMax = 50;
    double xStep = (xMax - xMin)/(double)xBinning;
    int yBinning = 100;
    double yMin = 0.5;
    double yMax = 4.;
    double yStep = (yMax - yMin)/(double)yBinning;

	double output[20010] = {0};
	
	auto start = high_resolution_clock::now();

	for(int i = 0; i < N ; i++){
		if(i%10000000==0)std::cout<<i<<std::endl;

		params[0] = 1.0;
		double E = signal.GetNR(params);// in keV
		
	    QuantaResult quanta;
		quanta = signal.GetQuanta_NR(E);
		vector<double> S1 = signal.GetS1(quanta);
		vector<double> S2 = signal.GetS2(quanta);

		params[0] *= signal.GetS1Eff(S1[1]);

    	if(params[0]<0.){params[0] = 0.;}

		if(S1[6]>0.&&S2[6]>0.)
		{
    		*(output+0) += 1.;
    		//get values, overflows and underflows in output[1]-output[9]
    		double xvalue = (double)S1[6];
    		double yvalue = (double)log10f(S2[6]/S1[6]);

    		if(xvalue<xMin && yvalue>=yMax)*(output+1) += 1.;
    		else if(xvalue>xMin && xvalue<xMax && yvalue>=yMax)*(output+2) += 1.;
    		else if(xvalue>=xMax && yvalue>=yMax)*(output+3) += 1.;
    		else if(xvalue<xMin && yvalue>yMin && yvalue<yMax)*(output+4) += 1.;
    		else if(xvalue>=xMin && xvalue<xMax && yvalue>=yMin && yvalue<yMax)
    		{
    		    int xbin = (int) ((xvalue - xMin)/xStep) + 1;
    		    int ybin = (int) ((yvalue - yMin)/yStep) + 1;
				int index = 9+(ybin-1)*xBinning+xbin;
				int weightindex = 9+xBinning*yBinning+((ybin-1)*xBinning+xbin);
    		    *(output+index) += params[0];
				*(output+weightindex) += params[0] * params[0];
    		    *(output+5) += 1.;
    		}
    		else if(xvalue>=xMax && yvalue>yMin && yvalue<yMax)*(output+6) += 1.;
    		else if(xvalue<xMin && yvalue<yMin)*(output+7) += 1.;
    		else if(xvalue>xMin && xvalue<xMax && yvalue<yMin)*(output+8) += 1.;
    		else if(xvalue>=xMax && yvalue<yMin)*(output+9) += 1.;
			
		}
		
	}

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	cout << "Time taken by function: "
         << duration.count() << " microseconds" << endl;

	std::ofstream outfile;
    outfile.open("./outputs_cpu.dat");

	for(int i = 0; i<20010; i++){
		outfile <<output[i]<<" ";
	}
	outfile.close();
	

	
	/*
	std::ofstream outfile;
    outfile.open("./outputs_cpu.dat");
	for(int i = 0; i < N; i++){
		params[0] = 1.0;
		double Enr = signal.GetNR(params);// in keV

		std::cout<<Enr<<" "<<params[0]<<std::endl;
		
	    QuantaResult quanta;
		quanta = signal.GetQuanta_NR(Enr);

		vector<double> S1 = signal.GetS1(quanta);
		vector<double> S2 = signal.GetS2(quanta);

		if (outfile.is_open()) 
        {
            outfile << Enr<<" "<<quanta.quanta<<" "
			<<quanta.excitons<<" "<<quanta.ions<<" "
			<<quanta.photons<<" "<<quanta.electrons<<" "
			<<S1[1]<<" "<<S2[2]<<" "<<S1[3]<<" "<<S2[3]<<" "
			<<S1[6]<<" "<<S2[6]<<" "<< params[0] << std::endl;
	    }
	    else{
			std::cout<<"error in file opening"<<std::endl;
		}

		//params[0] *= signal.GetS1Eff(S1[1]);
		//float Nions = (float)(quanta.ions);
		//float Nex = (float)(quanta.excitons);

    	//if(params[0]<0.){params[0] = 0.;}


	}
	outfile.close();
	*/


	return 0;
}
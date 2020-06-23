#ifndef __BETA_DECAY_HPP__
#define __BETA_DECAY_HPP__

#include <cmath>
#include <vector>
#include <string>

class BetaDecay{
	private:
		std::string DecayType;
		double EndPointEnergy;
		int NumBins;
		std::string OutputFile;
		double Charge;

		std::vector<double> Counts;
		std::vector<double> Energy;
		std::vector<double> Kurie;
		std::vector<double> Fermi;
	
		const double LowerBound = 10.0;	
		const double ElectronMass = 510.998950;
		const double FineStructure = 1.0/137.0;

		void CalculateSpectrum();
		void ConvertToW();
		void ConvertToE();
		void WriteData();
		void ApplyShapeFactor();
	public:
		BetaDecay(std::string,std::string,double,double,int);
};

#endif

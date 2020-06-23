/*
BetaDecay Spectrum Calculator
by: Thomas Ruland
*/
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <algorithm>
#include <cctype>

#include <TFile.h>
#include <TH1F.h>

#include "BetaDecay.hpp"

void ShowUsage(const char* );
void ShowHelp();


BetaDecay::BetaDecay(std::string output,std::string decay,double end,double z,int num){
	
	this->DecayType = decay;
	std::transform(this->DecayType.begin(),this->DecayType.end(),this->DecayType.begin(),[](unsigned char c){return std::tolower(c);});

	this->EndPointEnergy = end;
	this->OutputFile = output;
	this->NumBins = num;
	this->Charge = z;

	this->Counts = std::vector<double>(num,0.0);
	this->Energy = std::vector<double>(num,0.0);
	this->Kurie = std::vector<double>(num,0.0);
	this->Fermi = std::vector<double>(num,0.0);

	this->CalculateSpectrum();
	this->WriteData();
}

void BetaDecay::WriteData(){
	if( this->OutputFile.substr(this->OutputFile.find_last_of(".")+1) == "root" ){
		TFile* rootfile = new TFile(this->OutputFile.c_str(),"RECREATE");
		TH1F* energyspectrum = new TH1F("RawBetaSpectrum","RawBetaSpectrum",this->NumBins,this->LowerBound,this->EndPointEnergy);
		TH1F* kuriespectrum = new TH1F("KurieBetaSpectrum","KurieBetaSpectrum",this->NumBins,this->LowerBound,this->EndPointEnergy);
	
		for( size_t ii = 0; ii < this->NumBins; ++ii ){
			energyspectrum->Fill(this->Energy.at(ii),this->Counts.at(ii) * this->Fermi.at(ii) );
			kuriespectrum->Fill(this->Energy.at(ii),this->Kurie.at(ii));
		}
		
		rootfile->Write();
		rootfile->Close();
	}else{
		std::ofstream out;
		out.open(this->OutputFile.c_str());
		for( size_t ii = 0; ii < this->NumBins; ++ii )
			out << this->Energy.at(ii) << '\t' << this->Counts.at(ii)*this->Fermi.at(ii) << '\t' << this->Kurie.at(ii) << '\n';
		out.close();

	}
}

void BetaDecay::CalculateSpectrum(){
	//INITIALIZE ENERGY VECTOR
	double deltaE = (this->EndPointEnergy - this->LowerBound)/static_cast<double>(this->NumBins - 1);
	for( size_t ii = 0; ii < this->NumBins; ++ii )
		this->Energy.at(ii) = ii*deltaE + LowerBound;

	//SHIFT INTO W REGIME
	this->ConvertToW();
	
	//CALCULATE FERMI CORRECTED SPECTRA AND UNCORRECTED ONE AND KURIE PLOT
	double currcount = 1.0;
	double eta = 0.0;
	double beta = 0.0;
	double momentum = 0.0;
	for( size_t ii = 0; ii < this->NumBins; ++ii ){
		currcount = std::sqrt( ( this->Energy.at(ii) *  this->Energy.at(ii) ) - 1.0 );
		momentum = currcount;
		currcount *= this->Energy.at(ii);
		currcount *= ( (this->EndPointEnergy - this->Energy.at(ii)) * (this->EndPointEnergy - this->Energy.at(ii)) );

		beta = std::sqrt( 1.0 - ( 1.0 / ( 1.0 + ( momentum*momentum ) ) ) );
		eta = (M_PI * 2.0 * this->FineStructure * this->Charge) / ( beta );
		this->Counts.at(ii) = currcount;
		this->Fermi.at(ii) = ( eta )/( 1.0 - std::exp(-eta) );
		
		currcount = this->Energy.at(ii);
		currcount *= ( (this->EndPointEnergy - this->Energy.at(ii)) * (this->EndPointEnergy - this->Energy.at(ii)) );
		currcount /= std::sqrt( ( this->Energy.at(ii) *  this->Energy.at(ii) ) - 1.0 );
		this->Kurie.at(ii) = std::sqrt(currcount);
	}
	this->ApplyShapeFactor();

	//NORMALIZE BOTH FERMI CORRECTED AND UNCORRECTED
	double area = std::accumulate(this->Counts.begin(),this->Counts.end(),0.0)/static_cast<double>(this->NumBins);
	for( size_t ii = 0; ii < this->NumBins; ++ii )
		this->Counts.at(ii) /= area;

	//SHIFT BACK INTO REGULAR E REGIME	
	this->ConvertToE();
}

void BetaDecay::ApplyShapeFactor(){
	double psquared = 0.0;
	double qsquared = 0.0;

	if( this->DecayType.find("allowed") != std::string::npos ){
		return;
	}else if( this->DecayType.find("first") != std::string::npos and this->DecayType.find("non") == std::string::npos ){
		
		for( size_t ii = 0; ii < this->NumBins; ++ii ){
			psquared = (this->Energy.at(ii) * this->Energy.at(ii)) - 1.0;
			qsquared = (this->EndPointEnergy - this->Energy.at(ii))*(this->EndPointEnergy - this->Energy.at(ii));
			this->Counts.at(ii) *= (psquared + qsquared);
		}

	}else if( this->DecayType.find("second") != std::string::npos  and this->DecayType.find("non") == std::string::npos ){
		
		for( size_t ii = 0; ii < this->NumBins; ++ii ){
			psquared = (this->Energy.at(ii) * this->Energy.at(ii)) - 1.0;
			qsquared = (this->EndPointEnergy - this->Energy.at(ii))*(this->EndPointEnergy - this->Energy.at(ii));
			this->Counts.at(ii) *= (psquared*psquared + qsquared*qsquared + (10.0/3.0)*psquared*qsquared);
		}

	}else if( this->DecayType.find("third") != std::string::npos  and this->DecayType.find("non") == std::string::npos ){
		
		for( size_t ii = 0; ii < this->NumBins; ++ii ){
			psquared = (this->Energy.at(ii) * this->Energy.at(ii)) - 1.0;
			qsquared = (this->EndPointEnergy - this->Energy.at(ii))*(this->EndPointEnergy - this->Energy.at(ii));
			this->Counts.at(ii) *= (psquared*psquared*psquared + 7.0*psquared*psquared*qsquared + 7.0*qsquared*qsquared*psquared + qsquared*qsquared*qsquared);
		}

	}else if( this->DecayType.find("non") != std::string::npos ){
		
		double a = 0.0;
		double b = 0.0;
		double c = 0.0;
		std::cout << "The shape factor will be calculated as S = 1 + a*W + b/W + c*W^2, where W = 1 + KE[KeV]/m_e[KeV]\nEnter a:";
		std::cin >> a;
		std::cout << "Enter b:";
		std::cin >> b;
		std::cout << "Enter c:";
		std::cin >> c;
		for( size_t ii = 0; ii < this->NumBins; ++ii )
			this->Counts.at(ii) *= ( 1.0 + a*this->Energy.at(ii) + b/this->Energy.at(ii) + c*this->Energy.at(ii)*this->Energy.at(ii) );

	}else{
		std::cout << "Transition type : " << this->DecayType << " is not programmed yet, exiting now\n";
		exit(EXIT_FAILURE);
	}
};

void BetaDecay::ConvertToW(){
	for( size_t ii = 0; ii < this->NumBins; ++ii )
		this->Energy.at(ii) = 1.0 + this->Energy.at(ii)/this->ElectronMass;

	this->EndPointEnergy = 1.0 + this->EndPointEnergy/this->ElectronMass;
}

void BetaDecay::ConvertToE(){
	for( size_t ii = 0; ii < this->NumBins; ++ii )
		this->Energy.at(ii) = (this->Energy.at(ii) - 1.0)*this->ElectronMass;

	this->EndPointEnergy = (this->EndPointEnergy - 1.0)*this->ElectronMass;
}

void ShowUsage(char* programname){
	std::cout << "Usage: " << programname << " -e [endpoint energy in KeV] -o [outputfile] -t [transitiontype] -n [numbins] -z [charge]\n";
	std::cout << "Type: " << programname << " -h for help\n";
}

void ShowHelp(){
		std::cout << " -o outputs data into the root file named [outputfile] \n";
		std::cout << " -z charge, positive for beta- and negative for beta+\n";
		std::cout << " -n numbins in the histograms for root output\n";
		std::cout << " -e endpoint energy of the decay in KeV \n";
		std::cout << " -t transition type, need to put hyphens inbetween words rather than spaces, examples are as follows\n";
		std::cout << "    for allowed transition, pass -t allowed\n";
		std::cout << "    for first forbidden unique, pass -t first-unique, etc.\n";
}

int main( int argc, char* argv[]){
	int opt;
	std::string RootFile;
	double EndPointEnergy = 0.0;
	std::string TransitionType;
	int NumBins = 0;
	double Charge = 0.0;

	while( (opt = getopt(argc,argv,":o:t:e:z:n:h")) != -1 ){
		switch( opt ){
			case 'o':
				RootFile = std::string(optarg);
				break;
			case 'e':
				EndPointEnergy = std::atof(optarg);
				break;
			case 't':
				TransitionType = std::string(optarg);
				break;
			case 'n':
				NumBins = std::atoi(optarg);
				break;
			case 'z':
				Charge = std::atof(optarg);
				break;
			case 'h':
				ShowHelp();
				exit(EXIT_SUCCESS);
				break;
			default:
				std::cout << "Unknown option. See usage" << std::endl;
				ShowUsage(argv[0]);
				break;
		}
	}

	if( RootFile.empty() or TransitionType.empty() or EndPointEnergy <= 0.0 or NumBins <= 0 or Charge == 0.0 ){
		ShowUsage(argv[0]);
		exit(EXIT_FAILURE);
	}

	if( optind > argc ){
		ShowUsage(argv[0]);
		exit(EXIT_FAILURE);
	}

	BetaDecay decay(RootFile,TransitionType,EndPointEnergy,Charge,NumBins);

	exit(EXIT_SUCCESS);
}

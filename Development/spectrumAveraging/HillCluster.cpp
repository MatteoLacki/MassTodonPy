#include <Rcpp.h>
#include <iostream>

using namespace Rcpp;
using namespace std;


bool isNewHill(bool A, bool B)
{
	return A ? false : ( B ? true : false );
}


// [[Rcpp::export]]
IntegerVector HillCluster(
	NumericMatrix 	spectrum,
	double 			tol
){
	int nrow 	= spectrum.nrow();
	int clust  	= 0;
	
	IntegerVector out(nrow);
	out[0] = clust;

	double 	prevMass = spectrum(0,0);
	double 	prevInte = spectrum(0,1);
	double 	massDiff, mass, inte;
	
	bool 	prevClimb 	= true;
	bool  	climb;

	for(int i=1; i<nrow; i++){
		mass = spectrum(i,0);
		inte = spectrum(i,1);

		
		massDiff 	= mass - prevMass;
		climb 		= inte - prevInte >= 0;
		  
		// cout << massDiff << ' ' << tol << ' ' << (massDiff > tol) << ' ';
		if(
			massDiff > tol || // Mass difference too big.
			isNewHill( prevClimb, climb ) 	
		){
		 	clust++;
			prevClimb = true;
		} else {
			prevClimb = climb;
		}

		out[i] = clust;
		// cout << clust << endl;

		prevMass 	= mass;
		prevInte 	= inte;
	}

	return out;
}
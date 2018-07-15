#include <Rcpp.h>
using namespace Rcpp;

// // [[Rcpp::export]]
// NumericVector DirkCluster( 
// 	NumericVector masses, 
// 	double tolerance
// ){
// 	int n = masses.size();
// 	NumericVector clustering(n);

// 	clustering[0] = 1;

// 	for( int i = 1; i < n; ++i){

// 		if( masses[i] - masses[i-1] >  tolerance ){ 
// 			clustering[i] = clustering[i-1] + 1;
// 		} else {
// 			clustering[i] = clustering[i-1];
// 		}
// 	}

// 	return clustering;
// }

// [[Rcpp::export]]
NumericVector DirkCluster(
    NumericVector masses,
	double tolerance,
	double diam
){
	int n = masses.size();
	NumericVector clustering(n);

	clustering[0] = 1;
	double currDiam = 0.0;

	for( int i = 1; i < n; ++i){
		if( 
			(masses[i] - masses[i-1] >  tolerance) ||
			( currDiam > diam )
		){
			clustering[i] = clustering[i-1] + 1;
			currDiam = 0.0;
		} else {
			clustering[i] = clustering[i-1];
			currDiam += masses[i] - masses[i-1];
		}
	}

	return clustering;
}

#include <Rcpp.h>
#include <stdio.h>
#include <algorithm>    // std::min
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector bitonicClusters_cpp(
	NumericVector	mass, 
	NumericVector 	intensity,
	double 			massDiff_threshold
) {
  int 		n = mass.size();
  NumericVector clustering(n);
  int 		clustNo = 0;
  clustering[0] 	= clustNo;
  double 	prevInt = intensity[0];
  double 	prevMas = mass[0];
  double 	curMas, curInt;
  bool 		prevUp 	= true;
  bool		goingUp;

  for(int i=1; i<n; ++i ){
  	curMas = mass[i];
  	curInt = intensity[i];
  	goingUp= curInt > prevInt;

  	if( 
  		curMas-prevMas > massDiff_threshold ||
  		( !prevUp && goingUp )
  	){
  		clustNo++;
  	}

	clustering[i] = clustNo;
	prevMas = curMas;
	prevInt = curInt;	
  }

  return clustering;
}

/*** R
bitonicClusters <- function(
	spectrum, 
	massDiff_threshold
) {
	if( !('data.frame' %in% class(spectrum))) stop('Input is not a data.frame or tbl_df.') 
	if( !all( c('mass', 'intensity') %in% names(spectrum) ) ) stop('There is no mass or intensity among column names.')
	if( nrow(spectrum) == 0 ) return(c())
	
	bitonicClusters_cpp(
		spectrum$mass, 
		spectrum$intensity,
		massDiff_threshold
	)
}		
*/


// [[Rcpp::export]]
NumericVector bitonicClustersEnlighted(
	NumericVector	mass, 
	NumericVector 	intensity,
    NumericVector   coefs
) {
    double intercept = coefs[0]; 
    double linearTerm= coefs[1]; 
    double quadraticTerm = coefs[2];
    int 		n = mass.size();
    NumericVector clustering(n);
    int 		clustNo = 0;
    clustering[0] 	= clustNo;
    double 	prevInt = intensity[0];
    double 	prevMas = mass[0];
    double 	curMas, curInt, masDiff;
    bool 		prevUp 	= true;
    bool		goingUp;



    for(int i=1; i<n; ++i ){
    	curMas = mass[i];
    	curInt = intensity[i];
    	goingUp= curInt > prevInt;
    	masDiff= curMas-prevMas;
    	if( 
    		masDiff > (prevMas*quadraticTerm + linearTerm)*prevMas + intercept   ||
    		( !prevUp && goingUp )
    	){
    		clustNo++;
    	}

        clustering[i] = clustNo;
        prevMas = curMas;
        prevInt = curInt;	
    }

    return clustering;
}


// [[Rcpp::export]]
NumericVector bitonicLocal(
    NumericVector   mass, 
    NumericVector   intensity
) {
    int n       = mass.size();
    int clustNo = 0;    
    
    NumericVector clustering(n);

    double  M_a = mass[0];
    double  M_diff = 1.0;
    double  M_b;

    double  I_a = intensity[0];
    double  I_b;

    bool    up = false; 
    int     i=1;
    while( i<n ){
        M_b = mass[i];
        I_b = intensity[i];
        clustering[i-1] = clustNo;

        if( M_b-M_a > 1.5*M_diff ){ 
            clustNo++;
            up = true;
        }
        else{
            if( !up ){
                if( I_b > I_a ) clustNo++;
            }
            up = I_b > I_a;
        }
        
        M_diff  = M_b - M_a;
        M_a     = M_b;
        I_a     = I_b;
        i++;
    }

    // Classification of last M_b. 
    clustering[i-1] = clustNo;

    return clustering;
}

// [[Rcpp::export]]
IntegerMatrix bitonicLocalMatrix(
    NumericVector   mass, 
    NumericVector   intensity
) {
    int n       = mass.size();
    int clustNo = 0;    
    
    IntegerMatrix clustering(n, 2);

    double  M_a = mass[0];
    double  M_diff = 1.0;
    double  M_b;

    double  I_a = intensity[0];
    double  I_b;

    bool    up = false; 
    int     i=1;
    while( i<n ){
        M_b = mass[i];
        I_b = intensity[i];
        clustering(i-1,0) = clustNo;

        if( M_b-M_a > 1.5*M_diff ){ 
            clustNo++;
            up = true;
            clustering(i-1,1) = 1;
        }
        else{
            clustering(i-1,1) = 0;

            if( !up ){
                if( I_b > I_a ) clustNo++;
            }
            up = I_b > I_a;
        }

        M_diff  = M_b - M_a;
        M_a     = M_b;
        I_a     = I_b;
        i++;
    }

    // Classification of last M_b. 
    clustering(i-1,0) = clustNo;
    clustering(i-1,1) = 0;

    return clustering;
}


// [[Rcpp::export]]
IntegerVector massDiffClustering(
    NumericVector   mass,
    double  thresh = 1.2
) {
    int n       = mass.size();
    int clustNo = 0;    
    IntegerVector clust(n);
    clust[0]    = clustNo;
    double  M_a = mass[0];
    double  M_b;
    double  M_diff;
    for(int i=1; i<n; i++ ){
        M_b = mass[i];
        if( M_b > M_a ){
            if( M_b-M_a > thresh*M_diff ) clustNo++;    
            M_diff  = M_b - M_a;
        }           
        clust[i]= clustNo;
        M_a     = M_b;
    }
    return clust;
}

// [[Rcpp::export]]
IntegerVector massDiffClustAbsolute(
    NumericVector   mass,
    double          thresh 
) {
    int n       = mass.size();
    int clustNo = 0;    
    IntegerVector clust(n);
    clust[0]    = clustNo;
    double  M_a = mass[0];
    double  M_b;
    double  M_diff;
    for(int i=1; i<n; i++ ){
        M_b = mass[i];
        if( M_b > M_a ){
            if( M_b-M_a > thresh ) clustNo++;    
            M_diff  = M_b - M_a;
        }           
        clust[i]= clustNo;
        M_a     = M_b;
    }
    return clust;
}
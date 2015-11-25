#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
int test(int N){
	cout << N << '\n';
	return 0;
}


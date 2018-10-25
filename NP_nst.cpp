// slightly modified version of Model D of Nielsen and Berg (2014), aka SAM
// non-stationary
// cf. ICES WGNSSK REPORT 2015, p. 495
// ini randeff predicted yet unspecified (no stationary density)
// robcode: 0 = ML, 1 = loglog, 2 = ssh
// v0.1
#include <TMB.hpp>

// library needed for the multvariate normal distribution
using namespace density;

// parameter transformation, R -> (-1,+1)
template <class Type>
Type bound11(Type x){
	return (Type(1.0)-exp(-x))/(Type(1.0)+exp(-x));
}

template <class Type> 
Type square(Type x){
	return x*x;
}

// rho function: switch between ML (identity), loglog and ssh
template <class Type>
Type rhofunc(Type x, Type tc, int robcode){
	// tc = tuning constant
	// robcode: 0 = ML, 1 = loglog, 2 = ssh
	switch(robcode){
		case 0: // ML: rho = identity
			return x;
		break;
		case 1: // rho = log-logistic
			return log(Type(1.0)+exp(x+tc))-log(Type(1.0)+exp(tc));
		break;
		case 2: // rho = smothed semi-Huber (ssh)
			return CppAD::CondExpGt(x, -tc, x,
					tc*log((x+tc)/tc+sqrt(1.0+square((x+tc)/tc)))-tc);
		break;
		default:
			// std::cerr<<"Unknown robcode: "<<robcode<<std::endl;
			error("robcode must be 0 (ML), 1 (loglog) or 2 (ssh).");
			return 0;
		break;
	}
}


// objective function
template<class Type>
Type objective_function<Type>::operator() () {

	//--------------------------------------------------------------------------
	// Inputs
	//--------------------------------------------------------------------------

	// Data
	DATA_MATRIX(log_Cat); // log commercial catch (AC x TC)
	DATA_MATRIX(log_Iat); // log survey index (AI x TI)
	DATA_MATRIX(Mat); // natural mortality, same dim as log_Nat (AN x TN)
	DATA_SCALAR(daysprop); // prop into year survey was conducted, days/365
	DATA_VECTOR(tc); // dim=6, one for each contrib/sd
	DATA_INTEGER(robcode); // 0 = ML, 1 = loglog, 2 = ssh

	// Fixed parameters
	PARAMETER(log_sigmaF3); // log sd, proc eq F, age=3
	PARAMETER(log_sigmaF4); // log sd, proc eq F, age>=4
	PARAMETER(t_rho); // transfo corr, proc eq F
	PARAMETER(log_sigmaR); // log sd, proc eq N, age=3
	PARAMETER(log_sigmaN); // log sd, proc eq N, age>3
	PARAMETER(log_sigmaP); // log sd plus group, proc eq N, age=AN
	PARAMETER(log_sigmaC); // log sd, obs eq C, all ages
	PARAMETER(log_q3); // log catchability, obs eq I, age=3
	PARAMETER(log_q4); // log catchability, obs eq I, age=4
	PARAMETER(log_q5); // log catchability, obs eq I, age=5
	PARAMETER(log_q6); // log catchability, obs eq I, age=6
	PARAMETER(log_q7); // log catchability, obs eq I, age=7
	PARAMETER(log_q8); // log catchability, obs eq I, age=8
	PARAMETER(log_sigmaI); // log sd, obs eq I, all ages
	
	// Random effects
	PARAMETER_MATRIX(log_Fat); // log fishing mortality rate (AF x TF)
	PARAMETER_MATRIX(log_Nat); // log abundance thousands numbers (AN x TN)
	

	//--------------------------------------------------------------------------
	// Setup, procedures and init
	//--------------------------------------------------------------------------

	// dim of F proc
	int AF = log_Fat.rows(); // age = 3, ..., 9+
	int TF = log_Fat.cols(); // t = 1967, ..., 2016

	// dim of N proc
	int AN = log_Nat.rows(); // age = 3, ..., 10+
	int TN = log_Nat.cols(); // t = 1967, ..., 2016

	// dim of C obs
	// int AC = log_Cat.rows(); // age = 3, ..., 10+ // not used
	int TC = log_Cat.cols(); // t = 1967, ..., 2015

	// dim of I obs
	// int AI = log_Iat.rows(); // age = 3, ..., 8+ // not used
	int TI = log_Iat.cols(); // t = 1992, ..., 2016
	int t1992 = TN-TI; // time offset for variables ranging 1967-2016

	// Transform back all param and randeff
	Type sigmaF3 = exp(log_sigmaF3);
	Type sigmaF4 = exp(log_sigmaF4);
	Type rho = bound11(t_rho);
	Type sigmaR = exp(log_sigmaR);
	Type sigmaN = exp(log_sigmaN);
	Type sigmaP = exp(log_sigmaP);
	Type sigmaC = exp(log_sigmaC);
	Type q3 = exp(log_q3);
	Type q4 = exp(log_q4);
	Type q5 = exp(log_q5);
	Type q6 = exp(log_q6);
	Type q7 = exp(log_q7);
	Type q8 = exp(log_q8);
	Type sigmaI = exp(log_sigmaI);

	matrix<Type> Fat(AF,TF); // same dim as log_Fat (AF x TF)
	for (int t = 0; t < TF; t++){
		for (int a = 0; a < AF; a++){
			Fat(a,t) = exp(log_Fat(a,t));
		}
	}

	matrix<Type> Zat(AN,TN); // total mortality, same dim as Mat (AN x TN)
	for (int t = 0; t < TN; t++){
		// 0<=a<=AF-1
		for (int a = 0; a < AF; a++){
			Zat(a,t) = Mat(a,t) + Fat(a,t);
		}
		// a=AF=AN-1
		Zat(AF,t) = Mat(AF,t) + Fat(AF-1,t); // constant F for a>AF-1
	}

	matrix<Type> Nat(AN,TN); // same dim as log_Nat (AN x TN)
	for (int t = 0; t < TN; t++){
		for (int a = 0; a < AN; a++){
			Nat(a,t) = exp(log_Nat(a,t));
		}
	}

	Type nll = 0.0; // ini (robustified) neg loglik

	//--------------------------------------------------------------------------
	// Proc eq F: RW on log scale for all ages, AR(1) varcov across ages
	//--------------------------------------------------------------------------

	matrix<Type> Sigmaxi(AF,AF); // varcov matrix of noise
	Sigmaxi(0,0) = square(sigmaF3); // variance for recruits
	for (int i = 1; i < AF; i++){
		Sigmaxi(i,0) = pow(rho,Type(i))*sigmaF4*sigmaF3; // AR(1) across ages
		Sigmaxi(0,i) = Sigmaxi(i,0); // symmetry
		for (int j = 1; j < AF; j++){
			Sigmaxi(i,j) = pow(rho,Type(abs(i-j)))*square(sigmaF4); // AR(1)
		}
	}
	MVNORM_t<Type> xidist(Sigmaxi); // multinorm with cov mat Sigmaxi

	// dynamics for log F
	for (int t = 1; t < TF; t++){
		vector<Type> diff_logFat(AF);
		for (int a = 0; a < AF; a++){
			diff_logFat(a) = log_Fat(a,t) - log_Fat(a,t-1); // RW
		}
		nll -= rhofunc(-xidist(diff_logFat),tc(0),robcode);
	}

	//--------------------------------------------------------------------------
	// Proc eq N: RW for age=3, survival for other ages
	//--------------------------------------------------------------------------

	// dynamics for log N
	for (int t = 1; t < TN; t++){
		// a=0, RW for recruits
		nll -= rhofunc(dnorm(log_Nat(0,t), log_Nat(0,t-1), sigmaR, true),
					   tc(1),robcode);

		// 1<=a<(AN-1), survival
		for (int a = 1; a < (AN-1); a++){
			Type mu_logNat = log_Nat(a-1,t-1)-Fat(a-1,t-1)-Mat(a-1,t-1);
			nll -= rhofunc(dnorm(log_Nat(a,t), mu_logNat, sigmaN, true),
						   tc(2),robcode);

		}

		// a=AN-1=AF, survival for plus group, Fat fixed at (AF-1)
		Type mu_logNAt = log(Nat(AN-2,t-1)*exp(-Fat(AF-1,t-1)-Mat(AN-2,t-1))
						 + Nat(AN-1,t-1)*exp(-Fat(AF-1,t-1)-Mat(AN-1,t-1)));
		nll -= rhofunc(dnorm(log_Nat(AN-1,t), mu_logNAt, sigmaP, true),
					   tc(3),robcode); // tc(2) in v0.4
	}


	//--------------------------------------------------------------------------
	// Obs eq C: Baranov catch eq, same as in NP_st
	//--------------------------------------------------------------------------

	for (int t = 0; t < TC; t++){
		// 0<=a<=AF-1
		for (int a = 0; a < AF; a++){
			Type mu_logCat = log_Fat(a,t) - log(Zat(a,t))
				+ log(Type(1.0) - exp(-Zat(a,t))) + log_Nat(a,t);
			nll -= rhofunc(dnorm(log_Cat(a,t), mu_logCat, sigmaC, true),
						   tc(4),robcode); // tc(3) in v0.4
		}
		// a=AF=AN-1, Fat fixed at a=AF-1
		Type mu_logCAt = log_Fat(AF-1,t) - log(Zat(AF,t))
				+ log(Type(1.0) - exp(-Zat(AF,t))) + log_Nat(AF,t);
		nll -= rhofunc(dnorm(log_Cat(AF,t), mu_logCAt, sigmaC, true),
						tc(4),robcode); // tc(3) in v0.4
	}
	
	//--------------------------------------------------------------------------
	// Obs eq I: abundance proportional to survey index, same as in NP_st
	//--------------------------------------------------------------------------

	for (int t = 0; t < TI; t++){
		Type mu_logIat = 0.0; // ini mean logIat
		// a=0
		mu_logIat = log_q3 - Zat(0,t+t1992)*daysprop + log_Nat(0,t+t1992);
		nll -= rhofunc(dnorm(log_Iat(0,t),mu_logIat,sigmaI,true),tc(5),robcode);
		// a=1
		mu_logIat = log_q4 - Zat(1,t+t1992)*daysprop + log_Nat(1,t+t1992);
		nll -= rhofunc(dnorm(log_Iat(1,t),mu_logIat,sigmaI,true),tc(5),robcode);
		// a=2
		mu_logIat = log_q5 - Zat(2,t+t1992)*daysprop + log_Nat(2,t+t1992);
		nll -= rhofunc(dnorm(log_Iat(2,t),mu_logIat,sigmaI,true),tc(5),robcode);
		// a=3
		mu_logIat = log_q6 - Zat(3,t+t1992)*daysprop + log_Nat(3,t+t1992);
		nll -= rhofunc(dnorm(log_Iat(3,t),mu_logIat,sigmaI,true),tc(5),robcode);
		// a=4
		mu_logIat = log_q7 - Zat(4,t+t1992)*daysprop + log_Nat(4,t+t1992);
		nll -= rhofunc(dnorm(log_Iat(4,t),mu_logIat,sigmaI,true),tc(5),robcode);
		// a=5
		mu_logIat = log_q8 - Zat(5,t+t1992)*daysprop + log_Nat(5,t+t1992);
		nll -= rhofunc(dnorm(log_Iat(5,t),mu_logIat,sigmaI,true),tc(5),robcode);
	}
	// ^ tc(4) in v0.4

	//--------------------------------------------------------------------------
	// Outputs
	//--------------------------------------------------------------------------

	// Reports misc quantities
	// REPORT(Sigmaxi);

	// Reports simulated quantities
	// SIMULATE{
	// 	REPORT(log_Nat);
	// 	REPORT(log_Fat);
	// 	REPORT(log_Cat);
	// 	REPORT(log_Iat);
	// }

	// Reports on transformed parameters
	ADREPORT(sigmaF3);
	ADREPORT(sigmaF4);
	ADREPORT(rho);
	ADREPORT(sigmaR);
	ADREPORT(sigmaN);
	ADREPORT(sigmaP);
	ADREPORT(sigmaC);
	ADREPORT(q3);
	ADREPORT(q4);
	ADREPORT(q5);
	ADREPORT(q6);
	ADREPORT(q7);
	ADREPORT(q8);
	ADREPORT(sigmaI);

	// Reports on randeff
	ADREPORT(Fat);
	ADREPORT(Nat);

	// Reports on derived quantities of interest

	return nll;
}

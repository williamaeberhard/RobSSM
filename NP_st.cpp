// modified version of Model D of Nielsen and Berg (2014), aka SAM
// stationary
// cf. ICES WGNSSK REPORT 2015, p. 495
// ini randeff set at stationary density, mean baseline level for F and N
// robcode: 0 = ML, 1 = loglog, 2 = ssh
// v0.1
#include <TMB.hpp>

// library needed for the multivariate normal distribution
using namespace density;

// parameter transformation, R -> (-scale,+scale)
template <class Type>
Type bound11(Type x, Type scale){
	return scale*(Type(1.0)-exp(-x))/(Type(1.0)+exp(-x));
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
	DATA_SCALAR(boundARcoef); // bounds AR coef phiF, phiR and phiP
	DATA_VECTOR(tc); // dim=6, one for each contrib/sd
	DATA_INTEGER(robcode); // 0 = ML, 1 = loglog, 2 = ssh
	DATA_INTEGER(ntrunc); // truncation threshold infinite sums Var logN a=AN-1
	
	// Fixed parameters
	PARAMETER(meanlogF3); // st mean log scale, proc eq F, age=3
	PARAMETER(meanlogF4); // st mean log scale, proc eq F, age=4
	PARAMETER(meanlogF5); // st mean log scale, proc eq F, age=5
	PARAMETER(meanlogF6); // st mean log scale, proc eq F, age=6
	PARAMETER(meanlogF7); // st mean log scale, proc eq F, age=7
	PARAMETER(meanlogF8); // st mean log scale, proc eq F, age=8
	PARAMETER(meanlogF9); // st mean log scale, proc eq F, age=9+
	PARAMETER(t_phiF); // transfo AR(1) coef, proc eq F, all ages
	PARAMETER(log_sigmaF3); // log sd, proc eq F, age=3
	PARAMETER(log_sigmaF4); // log sd, proc eq F, age>=4
	PARAMETER(t_rho); // transfo corr across ages, proc eq F
	PARAMETER(meanlogN3); // st mean log recruits, proc eq N, age=3
	PARAMETER(t_phiR); // transfo AR(1) coef, proc eq N, age=3
	PARAMETER(log_sigmaR); // log sd, proc eq N, age=3
	PARAMETER(t_phiN); // transfo AR(1) coef survival, proc eq N, 3<age<AN
	PARAMETER(log_sigmaN); // log sd survival, proc eq N, 3<age<AN
	PARAMETER(t_phiP); // transfo AR(1) coef plus group, proc eq N, age=AN
	PARAMETER(log_sigmaP); // log sd plus group, proc eq N, age=AN
	PARAMETER(log_sigmaC); // log sd, obs eq C, all ages
	PARAMETER(log_q3); // log catchability, obs eq I, age=3
	PARAMETER(log_q4); // log catchability, obs eq I, age=4
	PARAMETER(log_q5); // log catchability, obs eq I, age=5
	PARAMETER(log_q6); // log catchability, obs eq I, age=6
	PARAMETER(log_q7); // log catchability, obs eq I, age=7
	PARAMETER(log_q8); // log catchability, obs eq I, age=8+
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
	Type phiF = bound11(t_phiF,boundARcoef);
	Type sigmaF3 = exp(log_sigmaF3);
	Type sigmaF4 = exp(log_sigmaF4);
	Type rho = bound11(t_rho,Type(1.0));
	Type phiR = bound11(t_phiR,boundARcoef);
	Type sigmaR = exp(log_sigmaR);
	Type phiN = bound11(t_phiN,Type(1.0)); // stationary even at phiN=1
	Type sigmaN = exp(log_sigmaN);
	Type phiP = bound11(t_phiP,boundARcoef);
	Type sigmaP = exp(log_sigmaP);
	Type sigmaC = exp(log_sigmaC);
	Type q3 = exp(log_q3);
	Type q4 = exp(log_q4);
	Type q5 = exp(log_q5);
	Type q6 = exp(log_q6);
	Type q7 = exp(log_q7);
	Type q8 = exp(log_q8);
	Type sigmaI = exp(log_sigmaI);

	vector<Type> meanlogF(AF); // simplifies proc eq below
	meanlogF(0) = meanlogF3;
	meanlogF(1) = meanlogF4;
	meanlogF(2) = meanlogF5;
	meanlogF(3) = meanlogF6;
	meanlogF(4) = meanlogF7;
	meanlogF(5) = meanlogF8;
	meanlogF(6) = meanlogF9; // AF-1=6

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
	// Proc eq F: AR(1) on log scale for all ages, AR(1) varcov across ages
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

	matrix<Type> Sigmaxist = Sigmaxi/(Type(1.0)-square(phiF)); // varcov st dist
	MVNORM_t<Type> xidistst(Sigmaxist); // multinorm st dist

	// ini st dist + mean baseline level
	vector<Type> diff_logFat(AF);
	for (int a = 0; a < AF; a++){
		diff_logFat(a) = log_Fat(a,0) - meanlogF(a);
	}
	nll -= rhofunc(-xidistst(diff_logFat),tc(0),robcode); // st dist
	// vector<Type> xi(AF); // ini vector of AR(1) proc error for simul
	// SIMULATE{
	// 	xi = MVNORM(Sigmaxi).simulate(); // sim multinorm
	// 	for (int a = 0; a < AF; a++){
	// 		log_Fat(a,0) = log_meanF(a) + xi(a);
	// 		Fat(a,0) = exp(log_Fat(a,0)); // update sim
	// 	}
	// }

	// dynamics for log F
	for (int t = 1; t < TF; t++){
		for (int a = 0; a < AF; a++){
			diff_logFat(a) = log_Fat(a,t) - (Type(1.0)-phiF)*meanlogF(a)
							 - phiF*log_Fat(a,t-1);
		}
		nll -= rhofunc(-xidist(diff_logFat),tc(0),robcode);
		// SIMULATE{
		// 	xi = MVNORM(Sigmaxi).simulate(); // sim multi Gaussian
		// 	for (int a = 0; a < AF; a++){
		// 		log_Fat(a,t) = (Type(1.0)-phiF)*log_meanF(a)
		// 					   + phiF*log_Fat(a,t-1) + xi(a);
		// 		Fat(a,t) = exp(log_Fat(a,t)); // update sim
		// 	}
		// }
	}

	// SIMULATE{ // update Zat in simulated data
	// 	for (int t = 0; t < TN; t++){
	// 		// 0<=a<=AF-1
	// 		for (int a = 0; a < AF; a++){
	// 			Zat(a,t) = Mat(a,t) + Fat(a,t);
	// 		}
	// 		// a=AF=AN-1
	// 		Zat(AF,t) = Mat(AF,t) + Fat(AF-1,t); // constant F beyond AF
	// 	}
	// }

	// vector<Type> meanF = exp(log_meanF); // wrong, multiv log-normal mean
	vector<Type> meanF(AF); // st mean over time, log-normal
	for (int a = 0; a < AF; a++){
		meanF(a) = exp(meanlogF(a)+Sigmaxist(a,a)/Type(2.0));
	}


	//--------------------------------------------------------------------------
	// Proc eq N: AR(1) for age=3, modified survival for other ages
	//--------------------------------------------------------------------------

	vector<Type> meanlogN(AN); // st mean log N for 0<=a<AN
	vector<Type> varlogN(AN); // st var
	vector<Type> sdlogN(AN); // st sd

	// ini st dist: a=0, recruits with meanlogN3
	meanlogN(0) = meanlogN3; // st mean AR(1) recruits
	varlogN(0) = square(sigmaR)/(Type(1.0)-square(phiR)); // st var
	sdlogN(0) = sqrt(varlogN(0)); // st sd
	nll -= rhofunc(dnorm(log_Nat(0,0),meanlogN(0),sdlogN(0),true),tc(1),robcode);
	// SIMULATE{
	// 	log_Nat(0,0) = rnorm(log_meanN3, sigmaR); // ini st dist
	// 	Nat(0,0) = exp(log_Nat(0,0)); // update sim
	// }

	// ini st dist: 1<=a<(AN-1), modified survival
	for (int a = 1; a < (AN-1); a++){ // plus-group treated separately below
		Type sumphiFM = 0.0; // ini
		for (int j = 1; j <= a; j++){
			sumphiFM += pow(phiN,j)*(meanF(a-j)+Mat(a-j,0)); // Mat invariant
		}
		meanlogN(a) = pow(phiN,a)*meanlogN3 - sumphiFM; // st mean logN

		Type sumphivarF = 0.0; // ini
		for (int j = 1; j <= a; j++){
			for (int i = 1; i <= a; i++){
				sumphivarF += pow(phiN,j)*pow(phiN,i)
				  *exp(meanlogF(a-j) + meanlogF(a-i)
					  +(Sigmaxist(a-j,a-j)+Sigmaxist(a-i,a-i))/Type(2.0))
				  *(exp(pow(phiF,Type(abs(i-j)))*Sigmaxist(a-j,a-i))-Type(1.0));
			}
		}
		varlogN(a) = pow(phiN,Type(2*a))*varlogN(0)
					 + square(sigmaN)*(Type(1.0)-pow(phiN,Type(2*a)))
						 /(Type(1.0)-square(phiN)) + sumphivarF; // st var logN
		
		sdlogN(a) = sqrt(varlogN(a)); // st sd

		nll -= rhofunc(dnorm(log_Nat(a,0), meanlogN(a), sdlogN(a), true),
					   tc(2),robcode);
		// SIMULATE{
		// 	log_Nat(a,0) = rnorm(log_meanN(a), sigmaN); // ini stationary dist
		// 	Nat(a,0) = exp(log_Nat(a,0)); // update sim
		// }
	}

	// ini st dist: a=AN-1=AF, modified survival for plus group
	Type sumphiFMA = 0.0; // ini
	for (int j = 1; j <= AF; j++){ // j<AN
		sumphiFMA += pow(phiN,j)*(meanF(AF-j)+Mat(AF-j,0)); // Mat invariant
		// ^ same meanF for a=AN-2=AF-1 and a=AN-1=AF
	}
	meanlogN(AN-1) = (pow(phiN,AF)*meanlogN3 - phiP*(meanF(AF-1)+Mat(AF,0))
						- sumphiFMA)/(Type(1.0)-phiP); // st mean logN

	Type sumcov1 = 0.0; // ini, equivalent to uV1 in notes
	for (int j = 1; j <= ntrunc; j++){
		for (int i = 1; i <= ntrunc; i++){
			sumcov1 += pow(phiP,Type(i+j))*(exp(pow(phiF,Type(abs(i-j)))*
											Sigmaxist(AF-1,AF-1))-Type(1.0));
		}
	}
	sumcov1 = sumcov1*exp(Type(2.0)*meanlogF(AF-1) + Sigmaxist(AF-1,AF-1));

	Type sumcov2 = 0.0; // ini, equivalent to uV2 in notes
	for (int j = 0; j <= ntrunc; j++){
		for (int i = 0; i <= ntrunc; i++){
			for (int k = 1; k < AN; k++){
				for (int l = 1; l < AN; l++){
					sumcov2 += pow(phiP,Type(i+j))*pow(phiN,Type(k+l))
					*exp(meanlogF(AN-1-k)+meanlogF(AN-1-l)
						+(Sigmaxist(AN-1-k,AN-1-k)
							+Sigmaxist(AN-1-l,AN-1-l))/Type(2.0))
					*(exp(pow(phiF,Type(abs(j+k-i-l)))
						*Sigmaxist(AN-1-k,AN-1-l))-Type(1.0));
				}
			}
		}
	}

	Type sumcov3 = 0.0; // ini, equivalent to uV3.alt in notes
	for (int j = 1; j <= ntrunc; j++){
		for (int i = 0; i <= ntrunc; i++){
			for (int k = 1; k < AN; k++){
				sumcov3 += pow(phiP,Type(i+j))*pow(phiN,k)
				*exp(meanlogF(AF-1)+meanlogF(AN-1-k)
						+(Sigmaxist(AF-1,AF-1)
							+Sigmaxist(AN-1-k,AN-1-k))/Type(2.0))
					*(exp(pow(phiF,Type(abs(j-i-k)))
						*Sigmaxist(AN-1-k,AF-1))-Type(1.0));
					// ^ no issue with plus-group different between F and N?
			}
		}
	}

	varlogN(AN-1) = square(sigmaP)/(Type(1.0)-square(phiP))
		+ square(sigmaN)*(square(phiN)*(Type(1.0)-pow(phiN,Type(2*(AN-2)))))
			/((Type(1.0)-square(phiP))*(Type(1.0)-square(phiN)))
		+ square(sigmaR)*pow(phiN,Type(2*(AN-1)))/square(phiR-phiP)
			*(square(phiR)/(Type(1.0)-square(phiR))
				+square(phiP)/(Type(1.0)-square(phiP))
				-Type(2.0)*phiR*phiP/(Type(1.0)-phiR*phiP))
		+ sumcov1 + sumcov2 + Type(2.0)*sumcov3; // st var

	sdlogN(AN-1) = sqrt(varlogN(AN-1)); // st sd
	nll -= rhofunc(dnorm(log_Nat(AN-1,0), meanlogN(AN-1), sdlogN(AN-1), true),
					   tc(3),robcode); // tc(2) in v0.4
	// SIMULATE{
	// 	log_Nat(AN-1,0) = rnorm(log_meanN(AN-1), sigmaP); // ini stationary dist
	// 	Nat(AN-1,0) = exp(log_Nat(AN-1,0)); // update sim
	// }

	vector<Type> meanN(AN); // st mean over time, log-normal
	for (int a = 0; a < AN; a++){
		meanN(a) = exp(meanlogN(a)+varlogN(a)/Type(2.0));
	}

	// dynamics for log N
	for (int t = 1; t < TN; t++){
		// a=0, AR(1) for recruits
		Type mu_logNat = (Type(1.0)-phiR)*meanlogN3 + phiR*log_Nat(0,t-1);
		nll -= rhofunc(dnorm(log_Nat(0,t), mu_logNat, sigmaR, true),
					   tc(1),robcode);
		// SIMULATE{
		// 	log_Nat(0,t) = rnorm((Type(1.0)-phiR)*meanlogN3
		// 						 + phiR*log_Nat(0,t-1), sigmaR);
		// 	Nat(0,t) = exp(log_Nat(0,t)); // update sim
		// }

		// 1<=a<(AN-1), modified survival
		for (int a = 1; a < (AN-1); a++){
			Type mu_logNat = phiN*(log_Nat(a-1,t-1)-Fat(a-1,t-1)-Mat(a-1,t-1));
			nll -= rhofunc(dnorm(log_Nat(a,t), mu_logNat, sigmaN, true),
						   tc(2),robcode);
			// SIMULATE{
			// 	mu_logNat = phiN*(log_Nat(a-1,t-1)-Fat(a-1,t-1)-Mat(a-1,t-1));
			// 	log_Nat(a,t) = rnorm(mu_logNat, sigmaN);
			// 	Nat(a,t) = exp(log_Nat(a,t)); // update sim
			// }
		}

		// a=AN-1=AF, modified survival for plus group, Fat fixed at (AF-1)
		Type mu_logNAt = phiN*(log_Nat(AN-2,t-1)-Fat(AF-1,t-1)-Mat(AN-2,t-1))
						 + phiP*(log_Nat(AN-1,t-1)-Fat(AF-1,t-1)-Mat(AN-1,t-1));
		nll -= rhofunc(dnorm(log_Nat(AN-1,t), mu_logNAt, sigmaP, true),
					   tc(3),robcode); // tc(2) in v0.4
		// SIMULATE{
		// 	mu_logNAt = phiN*(log_Nat(AN-2,t-1)-Fat(AF-1,t-1)-Mat(AN-2,t-1))
		// 				+ phiP*(log_Nat(AN-1,t-1)-Fat(AF-1,t-1)-Mat(AN-1,t-1));
		// 	log_Nat(AN-1,t) = rnorm(mu_logNAt, sigmaP);
		// 	Nat(AN-1,t) = exp(log_Nat(AN-1,t)); // update sim
		// }
	}

	//--------------------------------------------------------------------------
	// Obs eq C: Baranov catch eq, same as in NP_nst
	//--------------------------------------------------------------------------
 
	for (int t = 0; t < TC; t++){
		// 0<=a<=AF-1
		for (int a = 0; a < AF; a++){
			Type mu_logCat = log_Fat(a,t) - log(Zat(a,t))
				+ log(Type(1.0) - exp(-Zat(a,t))) + log_Nat(a,t);
			nll -= rhofunc(dnorm(log_Cat(a,t), mu_logCat, sigmaC, true),
						   tc(4),robcode); // tc(3) in v0.4
			// SIMULATE{
			// 	log_Cat(a,t) = rnorm(log_Fat(a,t)-log(Zat(a,t))
			// 		+log(Type(1.0)-exp(-Zat(a,t)))+log_Nat(a,t),sigmaC);
			// }
		}
		// a=AF=AN-1, Fat fixed at a=AF-1
		Type mu_logCAt = log_Fat(AF-1,t) - log(Zat(AF,t))
				+ log(Type(1.0) - exp(-Zat(AF,t))) + log_Nat(AF,t);
		nll -= rhofunc(dnorm(log_Cat(AF,t), mu_logCAt, sigmaC, true),
					   tc(4),robcode); // tc(3) in v0.4
		// SIMULATE{
		// 	log_Cat(AF,t) = rnorm(log_Fat(AF-1,t)-log(Zat(AF,t))
		// 			+log(Type(1.0)-exp(-Zat(AF,t)))+log_Nat(AF,t),sigmaC);
		// }
	}
	
	//--------------------------------------------------------------------------
	// Obs eq I: abundance proportional to survey index, same as in NP_nst
	//--------------------------------------------------------------------------

	for (int t = 0; t < TI; t++){
		Type mu_logIat = 0.0; // ini mean logIat
		// a=0
		mu_logIat = log_q3 - Zat(0,t+t1992)*daysprop + log_Nat(0,t+t1992);
		nll -= rhofunc(dnorm(log_Iat(0,t),mu_logIat,sigmaI,true),tc(5),robcode);
		// SIMULATE{
		// 	log_Iat(0,t) = rnorm(log_q3-Zat(0,t+t1992)*daysprop
		// 			+log_Nat(0,t+t1992),sigmaI);
		// }
		// a=1
		mu_logIat = log_q4 - Zat(1,t+t1992)*daysprop + log_Nat(1,t+t1992);
		nll -= rhofunc(dnorm(log_Iat(1,t),mu_logIat,sigmaI,true),tc(5),robcode);
		// SIMULATE{
		// 	log_Iat(1,t) = rnorm(log_q4-Zat(1,t+t1992)*daysprop
		// 			+log_Nat(1,t+t1992),sigmaI);
		// }
		// a=2
		mu_logIat = log_q5 - Zat(2,t+t1992)*daysprop + log_Nat(2,t+t1992);
		nll -= rhofunc(dnorm(log_Iat(2,t),mu_logIat,sigmaI,true),tc(5),robcode);
		// SIMULATE{
		// 	log_Iat(2,t) = rnorm(log_q5-Zat(2,t+t1992)*daysprop
		// 			+log_Nat(2,t+t1992),sigmaI);
		// }
		// a=3
		mu_logIat = log_q6 - Zat(3,t+t1992)*daysprop + log_Nat(3,t+t1992);
		nll -= rhofunc(dnorm(log_Iat(3,t),mu_logIat,sigmaI,true),tc(5),robcode);
		// SIMULATE{
		// 	log_Iat(3,t) = rnorm(log_q6-Zat(3,t+t1992)*daysprop
		// 			+log_Nat(3,t+t1992),sigmaI);
		// }
		// a=4
		mu_logIat = log_q7 - Zat(4,t+t1992)*daysprop + log_Nat(4,t+t1992);
		nll -= rhofunc(dnorm(log_Iat(4,t),mu_logIat,sigmaI,true),tc(5),robcode);
		// SIMULATE{
		// 	log_Iat(4,t) = rnorm(log_q7-Zat(4,t+t1992)*daysprop
		// 			+log_Nat(4,t+t1992),sigmaI);
		// }
		// a=5
		mu_logIat = log_q8 - Zat(5,t+t1992)*daysprop + log_Nat(5,t+t1992);
		nll -= rhofunc(dnorm(log_Iat(5,t),mu_logIat,sigmaI,true),tc(5),robcode);
		// SIMULATE{
		// 	log_Iat(5,t) = rnorm(log_q8-Zat(5,t+t1992)*daysprop
		// 			+log_Nat(5,t+t1992),sigmaI);
		// }
	}
	// ^ tc(4) in v0.4

	//--------------------------------------------------------------------------
	// Outputs
	//--------------------------------------------------------------------------

	// Reports misc quantities for testing
	// REPORT(Sigmaxi);
	REPORT(meanlogF); // st mean vector log F
	REPORT(Sigmaxist); // st varcov matrix log F
	// REPORT(meanF); // st mean vector F
	REPORT(meanlogN); // st mean vector log N
	REPORT(varlogN); // st var vector log N
	// REPORT(meanN); // st mean vector N

	// // Reports simulated quantities
	// SIMULATE{
	// 	REPORT(log_Nat);
	// 	REPORT(log_Fat);
	// 	REPORT(log_Cat);
	// 	REPORT(log_Iat);
	// }

	// Reports on transformed parameters
	ADREPORT(meanlogF3);
	ADREPORT(meanlogF4);
	ADREPORT(meanlogF5);
	ADREPORT(meanlogF6);
	ADREPORT(meanlogF7);
	ADREPORT(meanlogF8);
	ADREPORT(meanlogF9);
	ADREPORT(phiF);
	ADREPORT(sigmaF3);
	ADREPORT(sigmaF4);
	ADREPORT(rho);
	ADREPORT(meanlogN3);
	ADREPORT(phiR);
	ADREPORT(sigmaR);
	ADREPORT(phiN);
	ADREPORT(sigmaN);
	ADREPORT(phiP);
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
	ADREPORT(meanF); // vector of st log-normal mean F, original scale
	ADREPORT(meanN); // vector of approx st log-normal mean N, original scale

	return nll;
}

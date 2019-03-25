/***
	NAME
		GM_DEBPCoD.c
	PURPOSE
		Escalator Boxcar Train implementation of the Pilot Whale DEB PCoD model as described in:
		
		"Bio-energetic modeling of medium-sized cetaceans shows high sensitivity 
		to disturbance in seasons of low resource supply" 
		
		by

		Vincent Hin, John Harwood & André M. de Roos
		Ecological Applications 2019

	HISTORY
		Original implementation: AMdR
		Revised last: VH - Mar 25, 2019
***/

#include "escbox.h"

/*
*====================================================================================================================================
*
*  LABELLING ENVIRONMENT AND I-STATE VARIABLES
*
*====================================================================================================================================
*/

#define resource                  env[1]

#define age                       (i_state(0))
#define reserves                  (i_state(1))
#define totalmort                 (i_state(2))

#define IDpregstart               (i_const( 0))
#define IDpregclock               (i_const( 1))
#define IDcalved				  (i_const( 2))
#define IDweaned				  (i_const( 3))
#define IDpregnant                (i_const( 4))
#define IDlactating               (i_const( 5))

#define IDlength                  (i_const( 6))
#define IDbones                   (i_const( 7))
#define IDweight                  (i_const( 8))
#define IDweightM                 (i_const( 9))
#define IDfatratio                (i_const(10))

#define IDingest                  (i_const(11))
#define IDringest                 (i_const(12))
#define IDmingest                 (i_const(13))

#define IDmaint                   (i_const(14))
#define IDgrowth                  (i_const(15))
#define IDpregcosts               (i_const(16))
#define IDlactcosts               (i_const(17))
#define IDnet_energy              (i_const(18))

#define IDmortality               (i_const(19))
#define IDbackground              (i_const(20))
#define IDstarvation              (i_const(21))

#define IDagefirstrepro	          (i_const(22))
#define IDagefirstweaned          (i_const(23))
#define IDreprofemale		      (i_const(24))
#define IDR0estimated			  (i_const(25))
#define IDearlydeaths			  (i_const(26))

#define IDminloglifespan		  (i_const(27))

#define ispregnant(n)             (popIDcard[MAMMALS][(n)][IDpregnant] > TINY)
#define islactating(n)            (popIDcard[MAMMALS][(n)][IDlactating] > TINY)
#define isdying(n)                ((pop[MAMMALS][(n)][reserves] < MINFATRATIO * popIDcard[MAMMALS][n][IDweight]) || (pop[MAMMALS][(n)][totalmort] > popIDcard[MAMMALS][(n)][IDminloglifespan]))

/*
*====================================================================================================================================
*
*  DEFINING AND LABELLING CONSTANTS AND PARAMETERS
*
*====================================================================================================================================
*/

#define RANDOM_DEATH			  0

#define MAMMALS                   0
#define MOTHER                    0
#define CHILD                     1

#define YEAR                      365.0                                             // Total year length
#define MAXAGE                    1.0E9
#define MAXMORTALITY              5.0
#define MINFATRATIO               0.005                                             // Minimum fat level below individuals die with certainty
#define TINY                      1.0E-3
#define MINSURV					  1.0E-9

#define Tg                        parameter[ 0]                                     // Gestation period
#define Tl                        parameter[ 1]                                     // Lactation period
#define Tf                        parameter[ 2]                                     // Age at independent foraging
#define Tr                        parameter[ 3]                                     // Age at 100% foraging capacity

#define Lb                        parameter[ 4]                                     // Length at birth
#define Linf                      parameter[ 5]                                     // Asymptotic length
#define K                         parameter[ 6]                                     // Laird/Gompertz growth rate

#define OmegaS                    parameter[ 7]                                     // Weight-length scalar
#define OmegaE                    parameter[ 8]                                     // Weight-length exponent
#define OmegaM                    parameter[ 9]                                     // Relative maintenance reserves
#define OmegaC                    parameter[10]                                     // Relative maintenance foetus

#define Rho                       parameter[11]                                     // Target reserve-weight ratio
#define RhoS                      parameter[12]                                     // Starvation reserve-weight ratio

#define EtaF                      parameter[13]                                     // Steepness reserves-assimilation response
#define Gamma                     parameter[14]                                     // Shape parameter of resource assimilation-age response

#define PhiM                      parameter[15]                                     // Energy provisioning ratio through lactation

#define XiM                       parameter[16]                                     // Non-linearity in mother reserves-lactation rate
#define XiF                       parameter[17]                                     // Non-linearity in resource assimilation - calve age rate
#define XiC                       parameter[18]                                     // Non-linearity in milk assimilation-calve age rate

#define SigmaM                    parameter[19]                                     // Field metabolic rate scalar
#define SigmaG                    parameter[20]                                     // Structural mass growth costs
#define SigmaL                    parameter[21]                                     // Reserve-milk conversion ratio

#define Alpha1                    parameter[22]                                     // Background mortality parameter
#define Beta1					  parameter[23]                                     // Background mortality parameter
#define Alpha2                    parameter[24]                                     // Background mortality parameter
#define Beta2                     parameter[25]                                     // Background mortality parameter

#define MuS                       parameter[26]                                     // Starvation mortality scalar

#define EpsPlus                   parameter[27]                                     // Anabolic conversion efficiency
#define EpsMin                    parameter[28]                                     // Catabolic conversion efficiency

#define RESOURCE                  parameter[29]                                     // Resource density
#define RESOURCEAMPLITUDE         parameter[30]                                     // Amplitude of annual resource fluctuations
#define DISTURBANCESTART          parameter[31]                                     // Start date of annual disturbance
#define DISTURBANCEPERIOD         parameter[32]                                     // Duration of annual disturbance

#define PREGNANCYDELAY			  parameter[33]                                     // Average waiting time for female to get pregnant

static double                     SizeAtBirth, TotalNeonateCosts;
static void                       UpdateIDcards(double *env, population *pop);
static void                       UpdateStats(double value, double *mean, double *sum_sq, double *minval, double *maxval, long n);
static double                     bexp(double x);
static double					  MotherEinit = 78.81;
static long unsigned			  obsIBI = 0L;
static long unsigned			  obsIWI = 0L;
static double                     meanIBI, ssIBI, minIBI, maxIBI;
static double					  meanIWI, ssIWI, minIWI, maxIWI;

#if (RANDOM_DEATH == 1)
#define UNIDEVSEED				  -1
static double                     UniDev(void);
static long unsigned              UniDevSeed	= 0L;
#endif

#define StdDev(ss, n)             sqrt((ss)/((double)(n - 1)))

FILE *R0file;

/*
*====================================================================================================================================
*
*  USER INITIALIZATION ROUTINE ALLOWS OPERATIONS ON INITIAL POPULATIONS
*
*====================================================================================================================================
*/

void UserInit(int argc, char **argv, double *env, population *pop)

{
	double  NeonateWeightCosts, NeonateMaintCosts, NeonateLactCosts;
		
	switch (argc)
	{
		case 4:
			RESOURCE = atof(argv[3]);
		case 3:
			DISTURBANCEPERIOD = atof(argv[2]);
		default:
		break;
	}
	
#if (RANDOM_DEATH == 1)
	// Initialize the random number generator
    if (UNIDEVSEED > 0.0) UniDevSeed = (long)(UNIDEVSEED + TINY);
    UniDev();
#endif
	
	resource = RESOURCE*(1.0 - RESOURCEAMPLITUDE*cos(2*M_PI*91.25 / YEAR));
	
	SizeAtBirth         = OmegaS*pow(Lb, OmegaE);
	NeonateWeightCosts  = SigmaG*SizeAtBirth;                                         // Structural mass costs in MJ
	NeonateWeightCosts /= EpsMin;                                                     // Convert into kg stored reserves
	NeonateWeightCosts += RhoS*SizeAtBirth/(1-RhoS);                                  // Add the reserves transferred at birth to the calve

	NeonateMaintCosts   = SigmaM*(4/(3*OmegaE+4))*pow(OmegaC*OmegaS*pow(Lb, OmegaE), 0.75)*Tg;      // Maintenance costs in MJ
	NeonateMaintCosts  /= EpsMin;                                                     // Convert into kg stored reserves

	NeonateLactCosts    = 0.0;                                                        // DO NOT KNOW YET HOW MUCH THIS IS
	NeonateLactCosts   /= EpsMin;                                                     // Convert into kg stored reserves
  
	TotalNeonateCosts   = NeonateWeightCosts;			  							  // TotalNeonateCosts is used in pregnancy threshold
	
	meanIBI = ssIBI = minIBI = maxIBI = 0.0;
	meanIWI = ssIWI = minIWI = maxIWI = 0.0;
	
	// Delete all existing cohorts
	if (cohort_no[0])
	{
		int             i;
		extern void     SievePop(void);

		for (i = 0; i < cohort_no[0]; i++)
			pop[MAMMALS][i][number] = -1.0;

		pop[MAMMALS][MOTHER][number] = 1.0;
		SievePop();
	}
	else                                                                              // Add a newborn individual
		AddCohorts(pop, MAMMALS, 1);


	if (MotherEinit < 0)
	{
		pop[MAMMALS][MOTHER][number]    = 1.0;
		pop[MAMMALS][MOTHER][age]       = 0.0;
		pop[MAMMALS][MOTHER][reserves]  = RhoS * SizeAtBirth/(1-RhoS);
		pop[MAMMALS][MOTHER][totalmort] = 0.0;
	}
	else
	{
		pop[MAMMALS][MOTHER][number]    = 1.0;
		pop[MAMMALS][MOTHER][age]       = Tl;
		pop[MAMMALS][MOTHER][reserves]  = (MotherEinit > 0) ? MotherEinit : (RhoS * SizeAtBirth/(1-RhoS));
		pop[MAMMALS][MOTHER][totalmort] = 0.0;
	}
	
	popIDcard[MAMMALS][MOTHER][IDpregstart] = MAXAGE;
	popIDcard[MAMMALS][MOTHER][IDpregclock] = MAXAGE;
	popIDcard[MAMMALS][MOTHER][IDcalved] = MAXAGE;
	popIDcard[MAMMALS][MOTHER][IDweaned] = MAXAGE;
	popIDcard[MAMMALS][MOTHER][IDpregnant]  = 0;
	popIDcard[MAMMALS][MOTHER][IDlactating] = 0;

	popIDcard[MAMMALS][MOTHER][IDagefirstrepro] = 0;
	popIDcard[MAMMALS][MOTHER][IDagefirstweaned] = 0;
	popIDcard[MAMMALS][MOTHER][IDreprofemale] = 0;
	popIDcard[MAMMALS][MOTHER][IDR0estimated] = 0;
	popIDcard[MAMMALS][MOTHER][IDearlydeaths] = 0;

#if (RANDOM_DEATH == 1)
	double tmp;
	
    tmp = UniDev();
    tmp = max(tmp, DBL_MIN);
    tmp = min(tmp, 1.0 - DBL_MIN);
    tmp = -log(tmp);
	
	popIDcard[MAMMALS][MOTHER][IDminloglifespan] = tmp;
#else
	popIDcard[MAMMALS][MOTHER][IDminloglifespan] = 15.3;
#endif
  
	ReportNote("%-67s:  S=%.1f kg (W=%.1f kg)", "Calve size at birth", SizeAtBirth, SizeAtBirth/(1.0 - RhoS));
	ReportNote("%-67s:  TotalNeonateCosts=%.2f kg", "Total costs of neonate determines pregnancy threshold", TotalNeonateCosts);
  
	char tmpstr[2048];
	sprintf(tmpstr, "%s%s", argv[1], "_R0.dat");

	R0file = fopen(tmpstr, "a");
	  
	UpdateIDcards(env, pop);
	
	return;
}


/*
*====================================================================================================================================
*
*  SPECIFICATION OF THE NUMBER AND VALUES OF BOUNDARY POINTS
*
*====================================================================================================================================
*/

void SetBpointNo(double *env, population *pop, int *bpoint_no)

{
	bpoint_no[0] = 0;

	return;
}


/*==================================================================================================================================*/

void SetBpoints(double *env, population *pop, population *bpoints)

{
	return;
}


/*
*====================================================================================================================================
*
*  SPECIFICATION OF DERIVATIVES
*
*====================================================================================================================================
*/

void Gradient(double *env, population *pop, population *ofs, double *envgrad, population *popgrad, population *ofsgrad, population *bpoints)

{
	int           i;


	UpdateIDcards(env, pop);

	for (i = 0; i < cohort_no[MAMMALS]; i++)
	{
		popgrad[MAMMALS][i][number]		= 0.0;
		popgrad[MAMMALS][i][age]        = 1.0;
		popgrad[MAMMALS][i][totalmort]	= popIDcard[MAMMALS][i][IDmortality];
		popgrad[MAMMALS][i][reserves]   = popIDcard[MAMMALS][i][IDnet_energy];
		popgrad[MAMMALS][i][reserves]  /= (popIDcard[MAMMALS][i][IDnet_energy] > 0) ? EpsPlus : EpsMin;
	}

	envgrad[0] = 1;
	envgrad[1] = 0;

	return;
}


/*
*====================================================================================================================================
*
*  SPECIFICATION OF EVENT LOCATION AND DYNAMIC COHORT CLOSURE
*
*====================================================================================================================================
*/

void EventLocation(double *env, population *pop, population *ofs, population *bpoints, double *events)

{
	double eval;
	
	UpdateIDcards(env, pop);

	events[0] = -1.0;										// Event signalling both giving birth and conception
	events[1] = -1.0;										// Event signalling end of lactation (weaning)

	if (ispregnant(MOTHER) && islactating(MOTHER))  		// A pregnant and lactating mother must first wean
	  {
          events[1] = max(events[1], (pop[MAMMALS][MOTHER][age] - popIDcard[MAMMALS][MOTHER][IDcalved]) - Tl);
	  }
	else if (ispregnant(MOTHER)) 							// A pregnant, but non-lactating mother gives birth
	  {
          events[0] = max(events[0], (pop[MAMMALS][MOTHER][age] - popIDcard[MAMMALS][MOTHER][IDpregstart]) - Tg);
	  }
    else						 							// A non-pregnant mother can wean (if lactating) or become pregnant (two events!)
	  {
          if (islactating(MOTHER)) events[1] = max(events[1], (pop[MAMMALS][MOTHER][age] - popIDcard[MAMMALS][MOTHER][IDcalved]) - Tl);		// A lactating mother weans

		  if (popIDcard[MAMMALS][MOTHER][IDpregclock] == MAXAGE)																			// If pregnancy clock is not set
		    {
          	  eval      = min((pop[MAMMALS][MOTHER][reserves] - RhoS*popIDcard[MAMMALS][MOTHER][IDweight]) - TotalNeonateCosts, pop[MAMMALS][MOTHER][age] - Tl - 365); // The female must by older than Tl + 365 to become receptive
			  if (islactating(MOTHER)) eval = min(eval, (pop[MAMMALS][MOTHER][age] - popIDcard[MAMMALS][MOTHER][IDcalved]) - Tl + Tg - 1.0);// New pregnancy can be initiated 1 day after entering the last year of lactation
			  events[0] = max(events[0], eval);
	  		}
		  else																																// If pregnancy clock is already set
		    {
		      events[0] = max(events[0], (pop[MAMMALS][MOTHER][age] - popIDcard[MAMMALS][MOTHER][IDpregclock]) - PREGNANCYDELAY);			
		    }
	  }

	return;
}


/*==================================================================================================================================*/

int ForceCohortEnd(double *env, population *pop, population *ofs, population *bpoints)

{
	UpdateIDcards(env, pop);

	if (ispregnant(MOTHER) && iszero(pop[MAMMALS][MOTHER][age] - popIDcard[MAMMALS][MOTHER][IDpregstart] - Tg))	return COHORT_END;	// Add calve
	  
    if (islactating(MOTHER) && iszero(pop[MAMMALS][MOTHER][age] - popIDcard[MAMMALS][MOTHER][IDcalved] - Tl))						// Deal with lactating mother (pregnant or non-pregnant) that ends lactation
	  {
		  popIDcard[MAMMALS][MOTHER][IDlactating] = 0;
		  popIDcard[MAMMALS][MOTHER][IDreprofemale] += 0.5;
		  popIDcard[MAMMALS][MOTHER][IDR0estimated] += 0.5 * exp(-pop[MAMMALS][MOTHER][totalmort]) * exp(-pop[MAMMALS][CHILD][totalmort]);
		  if (popIDcard[MAMMALS][MOTHER][IDagefirstweaned] == 0) popIDcard[MAMMALS][MOTHER][IDagefirstweaned] = pop[MAMMALS][MOTHER][age];
		  
  		  if (popIDcard[MAMMALS][MOTHER][IDweaned] < MAXAGE)
  		  	{
				obsIWI++;
  				UpdateStats(pop[MAMMALS][MOTHER][age] - popIDcard[MAMMALS][MOTHER][IDweaned], &meanIWI, &ssIWI, &minIWI, &maxIWI, obsIWI);
  		  	}

		  popIDcard[MAMMALS][MOTHER][IDweaned] = pop[MAMMALS][MOTHER][age];
	   	  pop[MAMMALS][CHILD][number] = -1.0;																						// Throw away child after weaning.

	  
	  	  if ((!ispregnant(MOTHER)) && (popIDcard[MAMMALS][MOTHER][IDpregclock] == MAXAGE) && (pop[MAMMALS][MOTHER][age] > (Tl + 365)) && pop[MAMMALS][MOTHER][reserves] > (RhoS*popIDcard[MAMMALS][MOTHER][IDweight] + TotalNeonateCosts)) // Check whether new pregnancy can start if not done so already
		    {
				popIDcard[MAMMALS][MOTHER][IDpregclock] = pop[MAMMALS][MOTHER][age];
				popIDcard[MAMMALS][MOTHER][IDpregnant]  = 0;
		    }
		return NO_COHORT_END;
	  }
	else																															// All other cases are individuals (lactating or non-lactating) that become pregnant
	  {
		  if (popIDcard[MAMMALS][MOTHER][IDpregclock] == MAXAGE)																	// Set pregnancy clock if not done so
		    {
		  		popIDcard[MAMMALS][MOTHER][IDpregclock] = pop[MAMMALS][MOTHER][age];
				popIDcard[MAMMALS][MOTHER][IDpregnant]  = 0;
		    }
		  else																														// Set pregnancy when clock is already set
			{
				popIDcard[MAMMALS][MOTHER][IDpregnant]  = 1;
				popIDcard[MAMMALS][MOTHER][IDpregstart] = pop[MAMMALS][MOTHER][age];
				popIDcard[MAMMALS][MOTHER][IDpregclock] = MAXAGE;
			}
		  return NO_COHORT_END;
	  }
}


/*
*====================================================================================================================================
*
*  SPECIFICATION OF BETWEEN COHORT CYCLE DYNAMICS
*
*====================================================================================================================================
*/

void InstantDynamics(double *env, population *pop, population *ofs)

{
	int i;

	UpdateIDcards(env, pop);
	
	// A pregnant mother giving birth
	if (ispregnant(MOTHER) && iszero(pop[MAMMALS][MOTHER][age] - popIDcard[MAMMALS][MOTHER][IDpregstart] - Tg))
	  {
		if (cohort_no[0] < 2) AddCohorts(pop, MAMMALS, 1);

		pop[MAMMALS][CHILD][number]    = 1.0;
		pop[MAMMALS][CHILD][age]       = 0.0;
		pop[MAMMALS][CHILD][reserves]  = RhoS*SizeAtBirth/(1-RhoS);
		pop[MAMMALS][CHILD][totalmort] = 0.0;

		popIDcard[MAMMALS][CHILD][IDpregstart] = MAXAGE;
		popIDcard[MAMMALS][CHILD][IDpregnant]  = 0;
		popIDcard[MAMMALS][CHILD][IDcalved] = MAXAGE;
		popIDcard[MAMMALS][CHILD][IDweaned] = MAXAGE;
		popIDcard[MAMMALS][CHILD][IDlactating] = 0;
		popIDcard[MAMMALS][CHILD][IDagefirstrepro] = 0;
		popIDcard[MAMMALS][CHILD][IDagefirstweaned] = 0;
		popIDcard[MAMMALS][CHILD][IDreprofemale] = 0;
		popIDcard[MAMMALS][CHILD][IDR0estimated] = 0;
		popIDcard[MAMMALS][CHILD][IDearlydeaths] = 0;
		
#if (RANDOM_DEATH == 1)
		double tmp;
	
		tmp = UniDev();
		tmp = max(tmp, DBL_MIN);
		tmp = min(tmp, 1.0 - DBL_MIN);
		tmp = -log(tmp);
	
		popIDcard[MAMMALS][CHILD][IDminloglifespan] = tmp;
#else
		popIDcard[MAMMALS][CHILD][IDminloglifespan] = 15.3;
#endif
		
		if (popIDcard[MAMMALS][MOTHER][IDcalved] < MAXAGE)
		{
			obsIBI++;
			UpdateStats(pop[MAMMALS][MOTHER][age] - popIDcard[MAMMALS][MOTHER][IDcalved], &meanIBI, &ssIBI, &minIBI, &maxIBI, obsIBI);
		}
		
		pop[MAMMALS][MOTHER][reserves]         	-= pop[MAMMALS][CHILD][reserves];
		popIDcard[MAMMALS][MOTHER][IDcalved] 	= pop[MAMMALS][MOTHER][age];
		popIDcard[MAMMALS][MOTHER][IDpregnant]  = 0;
		popIDcard[MAMMALS][MOTHER][IDpregstart] = MAXAGE;
		popIDcard[MAMMALS][MOTHER][IDlactating] = 1;
		if (popIDcard[MAMMALS][MOTHER][IDagefirstrepro] == 0) popIDcard[MAMMALS][MOTHER][IDagefirstrepro] = pop[MAMMALS][MOTHER][age];
	  }

	// Delete all superfluous cohorts
	for (i = 2; i < cohort_no[0]; i++) pop[MAMMALS][i][number] = -1.0;

	// Kill the mother if she is dying and not lactating!
	if (isdying(MOTHER) && !islactating(MOTHER))
	{
		fprintf(R0file, "%16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\t %16.8G\n", 
				DISTURBANCEPERIOD, RESOURCE, 
				pop[MAMMALS][MOTHER][age], pop[MAMMALS][MOTHER][reserves], pop[MAMMALS][MOTHER][totalmort], 
				popIDcard[MAMMALS][MOTHER][IDminloglifespan], popIDcard[MAMMALS][MOTHER][IDagefirstrepro], popIDcard[MAMMALS][MOTHER][IDagefirstweaned],
				popIDcard[MAMMALS][MOTHER][IDreprofemale], popIDcard[MAMMALS][MOTHER][IDR0estimated], popIDcard[MAMMALS][MOTHER][IDearlydeaths], 
				meanIBI, ((obsIBI > 1) ? StdDev(ssIBI, obsIBI) : 0), minIBI, maxIBI, (double)obsIBI, meanIWI, ((obsIWI > 1) ? StdDev(ssIWI, obsIWI) : 0), minIWI, maxIWI, (double)obsIWI);

		pop[MAMMALS][MOTHER][number]  = -1.0;
	}
	
	// If the calf has lost all fat or when calf survival drops below MINSURV, throw away the calf
	if ((cohort_no[0] > 1) && isdying(CHILD))
	{
		pop[MAMMALS][CHILD][number]   = -1.0;
		popIDcard[MAMMALS][MOTHER][IDlactating]  = 0;
		// If calf death is early (before end of lactation) add 1 to earlydeaths
		if(pop[MAMMALS][MOTHER][age] < popIDcard[MAMMALS][MOTHER][IDcalved] + Tl) popIDcard[MAMMALS][MOTHER][IDearlydeaths] += 1.0;

		// Is there sufficient reserves to start pregnancy of the mother immediately again, if not done so already?
  	  	if ((!isdying(MOTHER)) && (!ispregnant(MOTHER)) && (popIDcard[MAMMALS][MOTHER][IDpregclock] == MAXAGE) && (pop[MAMMALS][MOTHER][age] > (Tl + 365)) && pop[MAMMALS][MOTHER][reserves] > (RhoS*popIDcard[MAMMALS][MOTHER][IDweight] + TotalNeonateCosts)) // Check whether new pregnancy can start
	      {
			  popIDcard[MAMMALS][MOTHER][IDpregclock] = pop[MAMMALS][MOTHER][age];
			  popIDcard[MAMMALS][MOTHER][IDpregnant]  = 0;
	      }
	}

	UpdateIDcards(env, pop);

	return;
}


/*
*====================================================================================================================================
*
*  SPECIFICATION OF OUTPUT VARIABLES
*
*====================================================================================================================================
*/

void DefineOutput(double *env, population *pop, double *output)

{
	int nr = 0;

	UpdateIDcards(env, pop);
	
	output[nr++] = env[1];

	output[nr++] = pop[MAMMALS][MOTHER][number];                                      // Column 3
	output[nr++] = pop[MAMMALS][MOTHER][age]/365.0;
	output[nr++] = pop[MAMMALS][MOTHER][reserves];
	output[nr++] = pop[MAMMALS][MOTHER][totalmort];
	
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDlength];                              // Column 7
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDbones];
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDweight];
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDweightM];
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDfatratio];

	output[nr++] = popIDcard[MAMMALS][MOTHER][IDingest];                              // Column 12
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDringest];
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDmingest];

	output[nr++] = popIDcard[MAMMALS][MOTHER][IDmaint];                               // Column 15
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDgrowth];
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDpregcosts];
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDlactcosts];
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDnet_energy];

	output[nr++] = popIDcard[MAMMALS][MOTHER][IDmortality];                           // Column 20
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDbackground];
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDstarvation];
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDminloglifespan];

	output[nr++] = popIDcard[MAMMALS][MOTHER][IDpregnant];                            // Column 24
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDpregstart];
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDpregclock];
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDlactating];
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDcalved];
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDweaned];

	output[nr++] = popIDcard[MAMMALS][MOTHER][IDagefirstrepro];							// Column 30
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDagefirstweaned];
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDreprofemale];
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDR0estimated];
	output[nr++] = popIDcard[MAMMALS][MOTHER][IDearlydeaths];

	output[nr++] = (pop[MAMMALS][MOTHER][reserves] - TotalNeonateCosts - RhoS*popIDcard[MAMMALS][MOTHER][IDweight]); // Column 35
  										
	if (cohort_no[0] > 1)
	{
		output[nr++] = pop[MAMMALS][CHILD][number];                                   // Column 36
		output[nr++] = pop[MAMMALS][CHILD][age]/365.0;
		output[nr++] = pop[MAMMALS][CHILD][reserves];
		output[nr++] = pop[MAMMALS][CHILD][totalmort];

		output[nr++] = popIDcard[MAMMALS][CHILD][IDlength];                           // Column 40
		output[nr++] = popIDcard[MAMMALS][CHILD][IDbones];
		output[nr++] = popIDcard[MAMMALS][CHILD][IDweight];
		output[nr++] = popIDcard[MAMMALS][CHILD][IDweightM];
		output[nr++] = popIDcard[MAMMALS][CHILD][IDfatratio];

		output[nr++] = popIDcard[MAMMALS][CHILD][IDingest];                           // Column 45
		output[nr++] = popIDcard[MAMMALS][CHILD][IDringest];
		output[nr++] = popIDcard[MAMMALS][CHILD][IDmingest];

		output[nr++] = popIDcard[MAMMALS][CHILD][IDmaint];                            // Column 48
		output[nr++] = popIDcard[MAMMALS][CHILD][IDgrowth];
		output[nr++] = popIDcard[MAMMALS][CHILD][IDpregcosts];
		output[nr++] = popIDcard[MAMMALS][CHILD][IDlactcosts];
		output[nr++] = popIDcard[MAMMALS][CHILD][IDnet_energy];

		output[nr++] = popIDcard[MAMMALS][CHILD][IDmortality];                        // Column 53
		output[nr++] = popIDcard[MAMMALS][CHILD][IDbackground];
		output[nr++] = popIDcard[MAMMALS][CHILD][IDstarvation];                       
		output[nr++] = popIDcard[MAMMALS][CHILD][IDminloglifespan];					  // Column 56
	}


	return;
}


/*
*====================================================================================================================================
*
*  ADDITIONAL, PROBLEM-SPECIFIC FUNCTIONS
*
*====================================================================================================================================
*/

static void UpdateIDcards(double *env, population *pop)

	/*
	* UpdateIDcards - Updates the values of all IDconst variables for the current time and environmental conditions
	*/

{
	int     i, notdisturbed;
	double  len, bones, fat, weight, weightM, fatratio;
	double  tau, foetuslen, foetusbones;
	double  milkage;
	double  ringest, mingest, tmp;
	double  maint, growth, pregcosts;
	double  background, starvation;
	double	timeinyear; 

	// Assume that the focal female is born in spring!
	resource = RESOURCE*(1.0 - RESOURCEAMPLITUDE*cos(2*M_PI*((env[0] + 91.25) / YEAR)));
    
	timeinyear   = env[0] - floor(env[0]/YEAR)*YEAR;
	notdisturbed = ((timeinyear >= DISTURBANCESTART) && (timeinyear <= (DISTURBANCESTART + DISTURBANCEPERIOD))) ? 0 : 1;

	for (i = 0; i < cohort_no[MAMMALS]; i++)
	{
		if (pop[MAMMALS][i][number] < 0.0) continue;

		len       = Linf - (Linf - Lb)*bexp(-K*pop[MAMMALS][i][age]);
		bones     = OmegaS*pow(len, OmegaE);
		fat       = pop[MAMMALS][i][reserves];
		weight    = bones + fat;
		weightM   = bones + OmegaM*fat;

		// Add foetus' weight to mother's weight
		pregcosts = 0.0;
		if (ispregnant(i))
		{
			tau           = pop[MAMMALS][MOTHER][age] - popIDcard[MAMMALS][MOTHER][IDpregstart];
			foetuslen     = Lb*tau/Tg;
			foetusbones   = OmegaS*pow(foetuslen, OmegaE);
			weight       += foetusbones;
			weightM      += OmegaC*foetusbones;
			pregcosts     = SigmaG*OmegaS*OmegaE*pow(Lb/Tg, OmegaE)*pow(tau, OmegaE-1);
		}

		fatratio  = fat/weight;

		// Calculate resource ingestion
		ringest   = notdisturbed*resource*pow(bones, 2.0/3.0);
		ringest  /= 1 + bexp(-EtaF*(Rho/fatratio - 1.0));
		ringest  *= pow(pop[MAMMALS][i][age], Gamma) / (pow(Tr, Gamma) + pow(pop[MAMMALS][i][age], Gamma));

		// Calculate milk ingestion
		mingest = 0.0;
		
		if ((!isdying(MOTHER)) && (pop[MAMMALS][i][age] < Tl))
		{
			mingest  = PhiM*pow(bones, 2.0/3.0);

			if (i == CHILD)                                                           // For the initial individual maximal milk supply
			{
				tmp  = max((1 - XiM)*(pop[MAMMALS][MOTHER][reserves] - RhoS*popIDcard[MAMMALS][MOTHER][IDweight]), 0.0);
				tmp /= (Rho - RhoS)*popIDcard[MAMMALS][MOTHER][IDweight] - XiM*(pop[MAMMALS][MOTHER][reserves] - RhoS*popIDcard[MAMMALS][MOTHER][IDweight]);
				mingest *= max(tmp, 0.0);
			}
			
			mingest  /= 1 + bexp(-EtaF*(Rho/fatratio - 1.0));
			milkage   = max(0, (1 - (pop[MAMMALS][i][age] - Tf)/(Tl-Tf))/(1 - XiC*(pop[MAMMALS][i][age]-Tf)/(Tl-Tf)));
			mingest  *= min(1, milkage);
		}

		// Compute the maintenance costs
		maint = SigmaM*pow(weightM, 0.75);

		// Compute the growth costs
		growth = OmegaE*SigmaG*K*(Linf - len)*OmegaS*pow(len,OmegaE-1);

		// Determine the survival probablity
		background = Alpha1 * bexp(-Beta1*pop[MAMMALS][i][age]) + Alpha2 * bexp(Beta2*pop[MAMMALS][i][age]);
		starvation = 0.0;

		if (fatratio < MuS*RhoS/(MuS + MAXMORTALITY)) 
			starvation = MAXMORTALITY;
		else
			starvation = max(MuS*(RhoS/fatratio - 1.0), 0);

		// Initialize the IDconst to be reset with default empty values
		popIDcard[MAMMALS][i][IDlength]     = len;
		popIDcard[MAMMALS][i][IDbones]      = bones;
		popIDcard[MAMMALS][i][IDweight]     = weight;
		popIDcard[MAMMALS][i][IDweightM]    = weightM;
		popIDcard[MAMMALS][i][IDfatratio]   = fatratio;

		popIDcard[MAMMALS][i][IDingest]     = ringest + mingest;
		popIDcard[MAMMALS][i][IDringest]    = ringest;
		popIDcard[MAMMALS][i][IDmingest]    = mingest;

		popIDcard[MAMMALS][i][IDmaint]      = maint;
		popIDcard[MAMMALS][i][IDgrowth]     = growth;
		popIDcard[MAMMALS][i][IDpregcosts]  = pregcosts;
		popIDcard[MAMMALS][i][IDlactcosts]  = 0.0;
		popIDcard[MAMMALS][i][IDnet_energy] = ringest + mingest - maint - growth - pregcosts;

		popIDcard[MAMMALS][i][IDmortality]  = background + starvation;
		popIDcard[MAMMALS][i][IDbackground] = background;
		popIDcard[MAMMALS][i][IDstarvation] = starvation;
	}


	// If mother is lactating deal with lactation costs
	if (!isdying(MOTHER) && islactating(MOTHER) && (cohort_no[0] > 1) && (pop[MAMMALS][CHILD][number] > 0.0))
	{
		popIDcard[MAMMALS][MOTHER][IDlactcosts]   = max(popIDcard[MAMMALS][CHILD][IDmingest]/SigmaL, 0.0);
		popIDcard[MAMMALS][MOTHER][IDnet_energy] -= popIDcard[MAMMALS][MOTHER][IDlactcosts];
	}

	return;
}


#define MAX_EXP                   50.0

static double bexp(double x)

{
	double  pw = x;

	pw = max(pw, -MAX_EXP);
	pw = min(pw, MAX_EXP);

	return exp(pw);
}

static void UpdateStats(double value, double *mean, double *sum_sq, double *minval, double *maxval, long n)

/*
 *  UpdateStats - Routine that updates the value of the mean, the sum of the
 *                squared deviations from the mean, the maximum and minimum
 *                for a series of number of which "value" is the next element.
 *
 *  Arguments   - value   : The current element in the series of numbers.
 *                mean    : The mean value of all preceding number in the series.
 *                sum_sq  : The sum of the squared deviations of all preceding numbers.
 *                minval  : The minimum value of all preceding number in the series.
 *                maxval  : The maximum value of all preceding number in the series.
 *                n       : The order number of the current value in the series.
 *                          It hence equals the number of observations on which
 *                          the resulting (!!) statistics are based.
 */

{
  double  deviation;
  
  deviation  = value -*mean;
  *mean     += (deviation/n);
  *sum_sq   += (deviation*(value -*mean));

  if ((n == 1) || (value >*maxval)) *maxval = value;
  if ((n == 1) || (value <*minval)) *minval = value;

  return;
}

// From here onwards we have to undefine time as an alias for env[0]
#if (RANDOM_DEATH == 1)

#undef time
#include <sys/time.h>
#include <time.h>

#define MAX_STORED                100
#ifndef RAND_MAX
#define RAND_MAX                  32767
#endif

#define M1                        259200
#define IA1                       7141
#define IC1                       54773
#define RM1                       (1.0/M1)
#define M2                        134456
#define IA2                       8121
#define IC2                       28411
#define RM2                       (1.0/M2)
#define M3                        243000
#define IA3                       4561
#define IC3                       51349

double UniDev(void)

/*
 * UniDev - Routine is adapted from the routine ran1() in "Numerical
 *          Recipes in C". It returns a uniform random deviate between
 *          0.0 and 1.0 and is entirely portable.
 *
 * Return - randev : the random number.
 */

{
  static long       ix1, ix2, ix3;
  static double     prevrans[MAX_STORED];
  static int        iff = 0;
  int               j;
  double            randev;
  struct timeval    tv;
  struct timezone   tz;

  if (iff == 0)
    {
      // time(NULL) returns seconds since some day
      if (!UniDevSeed)
        {
          gettimeofday(&tv, &tz);
          UniDevSeed = (unsigned)(tv.tv_usec % RAND_MAX);
          srand(UniDevSeed);
          UniDevSeed = (unsigned)(rand() % RAND_MAX);
        }
      ix1 = (IC1 - UniDevSeed) % M1;
      ix1 = (IA1*ix1 + IC1) % M1;
      ix2 = ix1 % M2;
      ix1 = (IA1*ix1 + IC1) % M1;
      ix3 = ix1 % M3;
      for (j = 0; j < MAX_STORED; j++)
        {
          ix1         = (IA1*ix1 + IC1) % M1;
          ix2         = (IA2*ix2 + IC2) % M2;
          prevrans[j] = (ix1 + ix2*RM2)*RM1;
        }
      iff = 1;
    }

  ix1 = (IA1*ix1 + IC1) % M1;
  ix2 = (IA2*ix2 + IC2) % M2;
  ix3 = (IA3*ix3 + IC3) % M3;
  j   = (int)((MAX_STORED*ix3)/(double)M3);

  if ((j < 0) || (j >= MAX_STORED)) ErrorExit(1, "UniDev: Invalid index.");

  randev      = prevrans[j];
  prevrans[j] = (ix1 + ix2*RM2)*RM1;

  return randev;
}

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3

#endif
/*==================================================================================================================================*/

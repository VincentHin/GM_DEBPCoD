/***
   NAME
     GM_DEBPCoD.h
   PURPOSE

	Escalator Boxcar Train header file for the life history simulation of
	the Pilot Whale DEB PCoD model as described in:
	
	"Bio-energetic modeling of medium-sized cetaceans shows high sensitivity 
	to disturbance in seasons of low resource supply" 
		
	by

	Vincent Hin, John Harwood & André M. de Roos
	Ecological Applications 2019
     
   HISTORY
	Original implementation: AMdR
	Revised last: VH - Mar 25, 2019
***/

#define POPULATION_NR             1
#define I_STATE_DIM               3
#define I_CONST_DIM               28
#define ENVIRON_DIM               2
#define OUTPUT_VAR_NR             55
#define PARAMETER_NR              34
#define TIME_METHOD               RKCK
#define EVENT_NR                  2
#define DYNAMIC_COHORTS           0
#define BIFURCATION               0

/*================================================================================*/

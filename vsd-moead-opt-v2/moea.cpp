/*==========================================================================
//  C++ Implementation of MOEA/D Based on Differential Evolution (DE) for Contest Multiobjective
//  Problems in CEC2009
//
//  Author: Hui Li
//
//  See the details of MOEA/D-DE and test problems in the following papers
//
//  1) H. Li and Q. Zhang, Comparison Between NSGA-II and MOEA/D on a Set of Multiobjective Optimization
//  Problems with Complicated Pareto Sets, Technical Report CES-476, Department of Computer Science,
//  University of Essex, 2007
//
//  2) H. Li and Q. Zhang, Multiobjective Optimization Problems with Complicated Pareto Sets, MOEA/D and NSGA-II,
//  IEEE Transaction on Evolutionary Computation, 2008, to appear.
//
//  If you have any questions about the codes, please contact
//  Dr. Hui Li       at   hzl@cs.nott.ac.uk   or
//  Dr. Qingfu Zhang at   qzhang@essex.ac.uk
//
//  Date: 14/NOV/2008
//
// ===========================================================================*/



#include "algorithm.h"
#include <omp.h>

void InitializeBounds(int nvar, char * Instance)
{
	if( !strcmp("UF1", Instance) || !strcmp("UF2", Instance) || !strcmp("UF3", Instance) || !strcmp("UF4", Instance) || !strcmp("UF5", Instance) || !strcmp("UF6", Instance) || !strcmp("UF7", Instance) || !strcmp("UF8", Instance) || !strcmp("UF9", Instance) || !strcmp("UF10", Instance))
	{
		for(int i = 0 ;  i < nvar; i++)
		{
		   vlowBound[i]=0.0;
		   vuppBound[i]=1.0;//2.0*(i+1.0);
		}
	}
	
	if( !strcmp("WFG1", Instance) || !strcmp("WFG2", Instance) || !strcmp("WFG3", Instance) || !strcmp("WFG4", Instance) || !strcmp("WFG5", Instance) || !strcmp("WFG6", Instance) || !strcmp("WFG7", Instance) || !strcmp("WFG8", Instance) || !strcmp("WFG9", Instance))
	{
		for(int i = 0 ;  i < nvar; i++)
		{
		   vlowBound[i]=0.0;
		   vuppBound[i]=2.0*(i+1.0);
		}
	}
	if( !strcmp("DTLZ1", Instance) || !strcmp("DTLZ2", Instance) || !strcmp("DTLZ3", Instance) || !strcmp("DTLZ4", Instance) || !strcmp("DTLZ5", Instance) || !strcmp("DTLZ6", Instance) || !strcmp("DTLZ7", Instance) )
	{
		for(int i = 0 ;  i < nvar; i++)
		{
		   vlowBound[i]=0.0;
		   vuppBound[i]=1.0;
		}
	}
	if( !strcmp("RWP1", Instance))
        {
                for(int i = 0 ;  i < nvar; i++)
                {
                   vlowBound[i]=0.0;
                   vuppBound[i]=1.0;
                }
        }
        if( !strcmp("RWP2", Instance))
        {
                for(int i = 0 ;  i < nvar; i++)
                {
                   vlowBound[i]=1.0;
                   vuppBound[i]=3.0;
                }
        }

}
int main(int argc, char *argv[])
{

	int pop;

	//std::ifstream readf("TestInstance.txt");
	//std::ifstream readf(argv[1]);

	//int numOfInstance;
	//readf>>numOfInstance;

	//printf("\n -- %d Instances are being tested ---\n\n", numOfInstance);

//	for(int inst=1; inst<=numOfInstance; inst++)
	{
		// the parameter setting of test instance
		//readf>>strTestInstance;
		//readf>>nvar;
		//readf>>nobj;
		int index = 1;
		int run = 1;
		strcpy(strpath, argv[index++]);
		strcpy(strTestInstance, argv[index++]);
		run= atoi(argv[index++]);
		nobj = atoi(argv[index++]);
		pops = atoi(argv[index++]);
		max_nfes= atoi(argv[index++]);
		niche =  atoi(argv[index++]);
		prob = atof(argv[index++]);
		if(argc <= index+2)//Two nvar and Di
		{
		   nvar = atoi(argv[index++]);
		 
		}
		else  //WFG instances..
		{
		   param_l = atoi(argv[index++]);
		   param_k = atoi(argv[index++]);
		   nvar = param_l + param_k;
		}
		Di = sqrt(nvar)*atof(argv[index++]);
		InitializeBounds(nvar, strTestInstance);

		//printf("\n -- Instance: %s, %d variables %d objectives di: %f\n\n", strTestInstance, nvar, nobj, Di);


		clock_t start, temp, finish;
		double  duration, last = 0;
		start = clock();

		std::fstream fout;
		char logFilename[1024];
		sprintf(logFilename, "LOG/LOG_MOEAD_%s.dat", strTestInstance);
		fout.open(logFilename,std::ios::out);
		fout<<"Inst: "<<strTestInstance<<endl;
	    fout<<"Time: \n\n";
		//#pragma omp parallel for	
		//for(int run=1; run<=35; run++)
		{
			//printf("\n -- %d-th run  -- \n", run);
			CMOEAD MOEAD;
		//	MOEAD.load_parameter();
		//	pop = MOEAD.pops;
			MOEAD.exec_emo(run);
			temp = clock();
			duration = (double)(temp - start) / CLOCKS_PER_SEC;
			fout<<run<<":"<<duration - last<<" ";
			last = duration;
			if(run%10==0) fout<<"\n";
			//break;
		}

		fout<<"\n\n";

		finish = clock();
		duration = (double)(finish - start) / CLOCKS_PER_SEC;


		fout<<"Mean  CPU Time  "<<duration/30<<" seconds"<<endl;
		fout<<"Total CPU Time  "<<duration<<" seconds"<<endl;
		fout.close();

	}
}

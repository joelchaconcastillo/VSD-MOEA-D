#ifndef __EVOLUTION_H_
#define __EVOLUTION_H_

#include <queue>
#include <iomanip>
#include "global.h"
#include "recomb.h"
#include "common.h"
#include "individual.h"

class CMOEAD
{

public:
	CMOEAD();
	virtual ~CMOEAD();

	void init_neighbourhood();               // calculate the neighbourhood of each subproblem
	void init_population();                  // initialize the population

	void load_parameter();

	void update_reference(CIndividual &ind);                 // update ideal point which is used in Tchebycheff or NBI method
	void update_problem(CIndividual &ind, int &id, int &type); // update current solutions in the neighbourhood
	void update_memory();
	int max_near_point(vector<int> &reference, vector<CIndividual> &candidates, vector<int> &penalized, vector<bool> &active);
//	int max_near_point(vector< CIndividual> &reference, vector<CIndividual> &candidates, vector<int> &penalized, vector<bool> &active);
//	double get_distance_near_point( vector<CIndividual> &SetA, CIndividual &ind);
	double get_distance_near_point( vector<int> &SetA, int index,  vector<CIndividual> &candidates);
	void replacement_phase();

	void evol_population();                                      // DE-based recombination
	void mate_selection(vector<int> &list, int cid, int size, int type);  // select mating parents
	

	// execute MOEAD
	void exec_emo(int run);

	void save_front(char savefilename[4024]);       // save the pareto front into files
	void save_pos(char savefilename[4024]);

	void tour_selection(int depth, vector<int> &selected);
	void comp_utility();

	void update_parameterD();

	double DCN(CIndividual &ind, int &index);
	double distance( vector<double> &a, vector<double> &b);
	vector <CSubproblem> population;
	vector<CIndividual> child_pop, best;	// memory solutions
	vector<int> indexSeeds;

	vector <double>      utility;
	void operator=(const CMOEAD &moea);

public:

	// algorithm parameters
	int		max_gen, curren_gen;       //  the maximal number of generations and current gen
//	int     pops;          //  the population size
    //int	    niche;         //  the neighborhood size
	int     limit;         //  the maximal number of solutions replaced
	//double  prob;          //  the neighboring selection probability
	double  rate;          //  the differential rate between solutions

	int     nfes;          //  the number of function evluations
	double	D;	//Current minimum distance

};

CMOEAD::CMOEAD()
{

}

CMOEAD::~CMOEAD()
{

}
void CMOEAD::update_parameterD()
{
	double TElapsed = nfes;
        double TEnd = max_nfes;

        D = Di - Di * ( 2.0*TElapsed / TEnd  );
  	//this->D = DI - (DF - DI )* ( 2*TElapsed / TEnd  );
	D = ( -1 > D)? -1: D ;
}
double CMOEAD::distance( vector<double> &a, vector<double> &b)
{
	double dist = 0 ;
   for(int i = 0; i < a.size(); i++)
	{
	   double factor = (a[i]-b[i])/(vuppBound[i]-vlowBound[i]);
	   dist += factor*factor;
	}
   return sqrt(dist);
}
double CMOEAD::DCN(CIndividual &ind, int &index)
{
	double min_dist = 10000000;
   for(int i = 0; i < pops; i++)
	if(i == index) continue;
	else
	min_dist  = min(min_dist, distance(population[i].indiv.x_var, ind.x_var));
   return min_dist;
}

double CMOEAD::get_distance_near_point( vector<int> &SetA, int index,  vector<CIndividual> &candidates)
{

   double min_distance = INFINITY;

   if(SetA.empty()) return min_distance;
   for(int i = 0 ; i < SetA.size(); i++)
	min_distance = min(min_distance, distance( candidates[SetA[i]].x_var, candidates[index].x_var) );

   return min_distance;
}

int CMOEAD::max_near_point(vector<int> &reference, vector<CIndividual> &candidates, vector<int> &penalized, vector<bool> &active)
{
   int index=-1, kpenalized=-1;
   double max_dcn = -INFINITY;
   for(int i = 0 ; i  < penalized.size(); i++)
   {
					int k = penalized[i];
					if( ! active[k/3] ) continue;
					double distance = get_distance_near_point( reference, k , candidates);
					if( distance > max_dcn )	
					{
						 index = k;
						 kpenalized = i;
						 max_dcn = distance;
					}
   }
   //Eliminar al k-ésimo elemento, intercambiarlo por el último
   iter_swap(penalized.begin() + kpenalized, penalized.end()-1);
   penalized.pop_back();
//   penalized.erase(penalized.begin()+kpenalized);
   return index;
}

//Se implementa una fase de reemplazo con el mismo procedimiento que la especiacion
//en el caso mono-objetivo para funciones multimodales
void CMOEAD::replacement_phase()
{
   if(D<=0) 
	 {

   for(int i = 0, j = 0 ; i < pops; i++) population[i].indiv = best[i];
	 return;
	 }

   vector<CIndividual > Candidates;
   vector<int> selected_pop;
   
   priority_queue< pair<double, int>> pq;
   double f1, f2, f3;
   for(int i = 0, j = 0 ; i < pops; i++)
	{
		Candidates.push_back(population[i].indiv);
		f1 = fitnessfunction( population[i].indiv.y_obj , population[i].namda);
		pq.push(make_pair(-f1, j++));

		Candidates.push_back(best[i]);
		f2 = fitnessfunction( best[i].y_obj , population[i].namda);
		pq.push(make_pair(-f2, j++));

		Candidates.push_back(child_pop[i]);
		f3 = fitnessfunction( child_pop[i].y_obj , population[i].namda);
		pq.push( make_pair(-f3, j++));
	}

  vector<int> penalized;//( Candidates.size(), 1);
	vector<int> weigth_indexes;
	vector<bool> active(pops, true); //Los subgrupos que están activos
	//Agregar al mejor individuo y penalizar a los que tengan una distancia mínima....
	//Max numbee of seeds....

///  //Extremes are seeds.....
///  for(int m = 0; m < nobj; m++)
///	{
///	   double minM=INFINITY;
///		 int index=-1;
///	   for(int i= 0; i < Candidates.size(); i++)
///		 {
///		 	double fm = Candidates[i].y_obj[m];
///			  if(fm < minM)
///				{
///			     minM = fm;
///					 index = i;
///				}
///		 }
///		  if(active[index/3])
///			{
///	     selected_pop.push_back(index);		 
///		   weigth_indexes.push_back(index/3);
///		   active[index/3]=false;
///			 }
///	}


	while(!pq.empty())
	{
	   pair<double, int> data = pq.top();
	   int index = data.second;
	   pq.pop(); 

	   if( !active[index/3]  ) continue;

	   if(get_distance_near_point(selected_pop, index, Candidates)  < D)
	   {
			penalized.push_back(index);
	

	   }
	   else
		{
		   selected_pop.push_back(index);
		   weigth_indexes.push_back(index/3);
		   indexSeeds.push_back(index/3);	
		   active[index/3]=false;
		}
	}	


//
//	int max_seeds =  (  2*curren_gen* pops )  / ( max_gen) ;
//	 max_seeds = max(1, max_seeds);


//  while(!pq.empty() )
//	{
//	   pair<double, int> data = pq.top();
//	   int index = data.second;
//	   pq.pop(); 
//	   if( !active[index/3]  ) continue;
//
//	   //if(get_distance_near_point(selected_pop, Candidates[index])  < D)
//		if((index+2)%3 == 0 &&  countbestseeds <= max_seeds)
//			{
//			   countbestseeds++;
//				 selected_pop.push_back(Candidates[index]);
//		     weigth_indexes.push_back(index/3);
//		     active[index/3]=false;
//			}
//	   else
//		{
//			   penalized.push_back(index);
//		}
//	}




        for(int i = penalized.size()-1; i >= 0 ;  i--)
	{
	   if( (penalized[i]+2)%3 ==0 ||  active[penalized[i]/3] == false ) 
	   {
		iter_swap(penalized.begin() + i, penalized.end()-1);
		penalized.pop_back();
	   }
	}


     vector<double> v_distances(penalized.size(), INFINITY);
     for(int i = 0 ;  i < penalized.size(); i++)
        {
           for(int j = 0; j< selected_pop.size(); j++)
           {
	      
              v_distances[i] = min(v_distances[i], distance( Candidates[penalized[i]].x_var, Candidates[selected_pop[j]].x_var));
           }
        }
//	if(selected_pop.size() < pops)
//	Advanced_Crowding_Procedure(selected_pop, Candidates, penalized, active, weigth_indexes);
	while(selected_pop.size() < pops)
	{
//	  int index = max_near_point(selected_pop, Candidates, penalized, active); 
	double maxd = -INFINITY;
        int index = -1;

        //find the max near point from survivor to children....
        for(int i = 0 ; i < penalized.size(); i++)
        {

		if( ! active[penalized[i]/3] )
		{
		   continue;
		}
                if( v_distances[i] > maxd)
                {
                        maxd = v_distances[i];
                        index = i;//penalized[i];
                }
        }

	  selected_pop.push_back(penalized[index]);
	  weigth_indexes.push_back(penalized[index]/3);
	  active[penalized[index]/3] = false;
	  for(int i = 0 ; i < penalized.size(); i++)
          {
             if( active[penalized[i]/3]==false) continue;
             v_distances[i] = min(v_distances[i] , distance(Candidates[penalized[index]].x_var, Candidates[penalized[i]].x_var ));
          }
	   iter_swap(penalized.begin() + index, penalized.end()-1);
	   penalized.pop_back();
	   iter_swap(v_distances.begin() + index, v_distances.end()-1);
	   v_distances.pop_back();
	}
	//Asignar a los individuos padres..
	for(int i = 0; i < pops; i++)
	{
	   population[weigth_indexes[i]].indiv = Candidates[selected_pop[i]];
	}

}
void CMOEAD::init_population()
{

	idealpoint = vector<double>(nobj, 1.0e+30);
	utility    = vector<double>(pops, 1.0);

	char filename[1024];
	// Read weight vectors from a data file
//	sprintf(filename,"/home/joel.chacon/Current/MyResearchTopics/MOEA-D-Diversity/MOEAD-DE/vsd-moead-opt/ParameterSetting/Weight/W%dD_%d.dat", nobj, pops);
	sprintf(filename,"%s/ParameterSetting/Weight/W%dD_%d.dat", strpath, nobj, pops);
	std::ifstream readf(filename);

    for(int i=0; i<pops; i++)
	{

		CSubproblem sub;

		// Randomize and evaluate solution
		sub.indiv.rnd_init();
		sub.indiv.obj_eval();

		sub.saved = sub.indiv;

		// Initialize the reference point
		update_reference(sub.indiv);

		// Load weight vectors
		for(int j=0; j<nobj; j++)
		{
		    readf>>sub.namda[j];
		}
	//	printf("\n"); getchar();

		// Save in the population
		population.push_back(sub);
		child_pop.push_back(sub.indiv);
		best.push_back(sub.indiv);
		nfes++;
	}
	readf.close( );
}

void CMOEAD::operator=(const CMOEAD &alg)
{
	//population = alg.population;
}

void CMOEAD::init_neighbourhood()
{
    vector<double> dist   = vector<double>(pops, 0);
	vector<int>    indx   = vector<int>(pops, 0);

	for(int i=0; i<pops; i++)
	{
		// calculate the distances based on weight vectors
		for(int j=0; j<pops; j++)
		{
		    dist[j]    = dist_vector(population[i].namda,population[j].namda);
			indx[j]  = j;
		}

		// find 'niche' nearest neighboring subproblems
		minfastsort(dist,indx,population.size(),niche);

		// save the indexes of the nearest 'niche' neighboring weight vectors
		for(int k=0; k<niche; k++)
		{
			population[i].table.push_back(indx[k]);
		}

	}
	dist.clear();
	indx.clear();
}

void CMOEAD::tour_selection(int depth, vector<int> &selected)
{
	// selection based on utility
	vector<int> candidate;
	for(int k=0;    k<nobj; k++)    selected.push_back(k);   // select first m weights
	for(int n=nobj; n<pops; n++)    candidate.push_back(n);  // set of unselected weights

	while(selected.size()<int(pops/5.0))
	{
	    int best_idd = int(rnd_uni(&rnd_uni_init)*candidate.size()), i2;
		int best_sub = candidate[best_idd], s2;
		for(int i=1; i<depth; i++)
		{
		    i2  = int(rnd_uni(&rnd_uni_init)*candidate.size());
			s2  = candidate[i2];
			if(utility[s2]>utility[best_sub])
			{
				best_idd = i2;
			    best_sub = s2;
			}
		}
		selected.push_back(best_sub);
		candidate.erase(candidate.begin()+best_idd);
	}
}

void CMOEAD::comp_utility()
{
	double f1, f2, uti, delta;
    for(int n=0; n<pops; n++)
	{
		f1 = fitnessfunction(population[n].indiv.y_obj, population[n].namda);
		f2 = fitnessfunction(population[n].saved.y_obj, population[n].namda);
		delta = f2 - f1;
        if(delta>0.001)  utility[n] = 1.0;
		else{
            uti        = 0.95*(1.0+delta/0.001)*utility[n];
			utility[n] = uti<1.0?uti:1.0;
		}
        population[n].saved = population[n].indiv;
	}
}

void CMOEAD::update_problem(CIndividual &indiv, int &id, int &type)
{
	// indiv: child solution
	// id:   the id of current subproblem
	// type: update solutions in - neighborhood (1) or whole population (otherwise)
	int size, time = 0;
	if(type==1)	size = population[id].table.size();  // from neighborhood
	else        size = population.size();            // from whole population

	// a random order to update
	std::vector<int> perm(std::vector<int>(size, 0));
	for(int k=0; k<size; k++) perm[k] = k;
	random_shuffle(perm.begin(), perm.end());

		child_pop[id] = indiv;
		//Sólo actualizar el problema actual... luego se revisa el efecto de actualizar información entre los vecindarios....
//		memory[id].push_back(indiv);

	//Check the best problem option......
///	double MinF = INFINITY;
///	int index = -1;
    for(int i=0; i<size; i++)
	{
		// Pick a subproblem to update
		int k;
		if(type==1) k = population[id].table[perm[i]];
		else        k = perm[i];

		////double fitness_K = fitnessfunction(indiv.y_obj, population[k].namda);
		////if( MinF >  fitness_K )
		////{
		////   MinF =  fitness_K;
		////   index = k;
		////}
		//memory[k].push_back(indiv);
		// calculate the values of objective function regarding the current subproblem
		double f1, f2;
		//if( DCN_best(indiv, id) < D) continue;


		f1 = fitnessfunction(best[k].y_obj, population[k].namda);
		f2 = fitnessfunction(indiv.y_obj, population[k].namda);
		if(f2<f1)
		{
			best[k] = indiv;
			time++;
		}
	//	return;	

		// //the maximal number of solutions updated is not allowed to exceed 'limit'
	//	if(time>=limit)
	//	{
	//		return;
	//	}
	}
	////////Update the selected problem...

//		double f1, f2;
//
//		f1 = fitnessfunction(best[index].y_obj, population[index].namda);
//		f2 = fitnessfunction(indiv.y_obj, population[index].namda);
//		if(f2<f1)
//			best[index]= indiv;


	perm.clear();
}

void CMOEAD::update_reference(CIndividual &ind)
{
	//ind: child solution
	for(int n=0; n<nobj; n++)
	{
		if(ind.y_obj[n]<idealpoint[n])
		{
			idealpoint[n] = ind.y_obj[n];
		}
	}
}

void CMOEAD::mate_selection(vector<int> &list, int cid, int size, int type){
	// list : the set of the indexes of selected mating parents
	// cid  : the id of current subproblem
	// size : the number of selected mating parents
	// type : 1 - neighborhood; otherwise - whole population
	int ss   = population[cid].table.size(), id, parent;
    while(list.size()<size)
	{
		if(type==1){
		    id      = int(ss*rnd_uni(&rnd_uni_init));
			parent  = population[cid].table[id];
		}
		else
			parent  = int(population.size()*rnd_uni(&rnd_uni_init));

		// avoid the repeated selection
		bool flag = true;
		for(int i=0; i<list.size(); i++)
		{
			if(list[i]==parent) // parent is in the list
			{
				flag = false;
				break;
			}
		}

		if(flag) list.push_back(parent);
	}
}

void CMOEAD::evol_population()
{

	// random order of subproblems at each generation
	//vector<int> order(vector<int>(pops,0));
	//for(int i=0; i<pops; i++)  order[i] = i;
	//random_shuffle(order.begin(), order.end());

	///vector<int> order;	this->tour_selection(10, order);
	vector<int> order;
	for(int i = 0; i < pops; i++) order.push_back(i);
	

    for(int sub=0; sub<order.size(); sub++)
	{

		int c_sub = order[sub];    // random order
		//  printf("%d ", c_sub);

		int type;
		double rnd = rnd_uni(&rnd_uni_init);

		// mating selection based on probability
		if(rnd<prob)     type = 1;   // from neighborhood
		else             type = 2;   // from population

		// select the indexes of mating parents
		vector<int> plist;
		mate_selection(plist, c_sub, 2, type);  // neighborhood selection


		// produce a child solution
		CIndividual child1, child2;
		double rate2 = 0.5; //rate + 0.25*(rnd_uni(&rnd_uni_init) - 0.5);
		//double rate2 = rate + 0.25*(rnd_uni(&rnd_uni_init) - 0.5);
		diff_evo_xoverB(population[c_sub].indiv,population[plist[0]].indiv,population[plist[1]].indiv, child1, rate2);

		//real_sbx_xoverA(population[plist[0]].indiv, population[plist[1]].indiv, child1, child2);
		plist.clear();

		// apply polynomial mutation
		realmutation(child1, 1.0/nvar);
	//	realmutation(child2, 1.0/nvar);
		// repair method
		//repait
		// evaluate the child solution
		child1.obj_eval();
	//	child2.obj_eval();

		// update the reference points and other solutions in the neighborhood or the whole population
		update_reference(child1);
	//	update_reference(child2);

		//child_pop[c_sub] = child;
		update_problem(child1, c_sub, type);
	//	update_problem(child2, c_sub, type);
	//	nfes++;
	}
		replacement_phase();

}


void CMOEAD::exec_emo(int run)
{
    char filename1[5024];
    char filename2[5024];
		seed = run;
	seed = (seed + 23)%1377;
	rnd_uni_init = -(long)seed;

	// initialization
	nfes      = 0;
	init_population();
        init_neighbourhood();
//	sprintf(filename1,"/home/joel.chacon/Current/MyResearchTopics/MOEA-D-Diversity/MOEAD-DE/vsd-moead-opt/POS/POS_MOEAD_%s_RUN%d_seed_%d_nobj_%d.dat_bounded",strTestInstance,run, seed, nobj);
	sprintf(filename1,"%s/POS/POS_MOEAD_%s_RUN%d_seed_%d_nobj_%d_niche_%d.dat_bounded",strpath, strTestInstance,run, seed, nobj, niche);
	//sprintf(filename2,"/home/joel.chacon/Current/MyResearchTopics/MOEA-D-Diversity/MOEAD-DE/vsd-moead-opt/POF/POF_MOEAD_%s_RUN%d_seed_%d_nobj_%d.dat_bounded",strTestInstance,run, seed, nobj);
	sprintf(filename2,"%s/POF/POF_MOEAD_%s_RUN%d_seed_%d_nobj_%d_niche_%d.dat_bounded",strpath, strTestInstance,run, seed, nobj, niche);
	//for(int gen=1; gen<=max_gen; gen++)
//	for(nfes=1; nfes<=max_nfes
	while(nfes<max_nfes)
	{
		//curren_gen = gen;	
		update_parameterD();
		evol_population();
	        nfes += pops;
	}
		save_pos(filename1);
		save_front(filename2);
	population.clear();
	idealpoint.clear();

}

void CMOEAD::load_parameter()
{
	char filename[4024];

	//sprintf(filename,"/home/joel.chacon/Current/MyResearchTopics/MOEA-D-Diversity/MOEAD-DE/vsd-moead-opt/ParameterSetting/%s.txt", strTestInstance);
	sprintf(filename,"%s/ParameterSetting/%s.txt", strpath,strTestInstance);

	char temp[4024];
	std::ifstream readf(filename);
	readf.getline(temp, 4024);
	//puts(temp);
int trash;	
	readf>>pops;
	readf>>max_gen;
	readf>>trash;
	//readf>>niche;
	readf>>limit;
	readf>>prob;
	readf>>rate;
	niche=10;
//	printf("\n Parameter Setting in MOEA/D \n");
//	printf("\n pop %d, gen %d, niche %d, limit %d, prob %f, rate %f\n\n", pops, max_gen, niche, limit, prob, rate);

	readf.close();
}


void CMOEAD::save_front(char saveFilename[4024])
{

    std::fstream fout;
	//fout.open(saveFilename,std::ios::out);
	fout.open(saveFilename,fstream::app|fstream::out );
	for(int n=0; n<pops; n++)
	{
		for(int k=0;k<nobj;k++)
			fout<<best[n].y_obj[k]<<"  ";
		for(int k=0;k<nobj;k++)
			fout<<population[n].indiv.y_obj[k]<<"  ";
	for(int k=0;k<nobj;k++)
			fout<<child_pop[n].y_obj[k]<<"  ";
		fout<<"\n";
	}
	fout.close();
}

void CMOEAD::save_pos(char saveFilename[4024])
{
    std::fstream fout;
	//fout.open(saveFilename,std::ios::out);
	fout.open(saveFilename, fstream::app|fstream::out);
	for(int n=0; n<pops; n++)
	{
		for(int k=0;k<nvar;k++)
			fout<<population[n].indiv.x_var[k] << "  ";
			//fout<<population[n].indiv.x_var[k]<< fixed << setprecision(30) << "  ";
//	  for(int k=0;k<nvar;k++)
//			fout<<best[n].x_var[k]<<"  ";
//	  for(int k=0;k<nvar;k++)
//			fout<<child_pop[n].x_var[k]<<"  ";
		fout<<"\n";
	}
	fout.close();
}



#endif

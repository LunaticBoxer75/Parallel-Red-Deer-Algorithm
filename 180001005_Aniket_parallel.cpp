/*
Author: Aniket Sangwan (180001005)
*/

#include <bits/stdc++.h>
#include<omp.h>
using namespace std;

#define double long double
#define pb push_back
#define F first
#define S second

const int inf=1e6;

mt19937 rnd(chrono::system_clock::now().time_since_epoch().count());
const double RANDOM_MAX = 4294967295;

// Returns random value
double get_rand(){
   return (double)rnd() / RANDOM_MAX;
}


vector<double>fitness_scores; // Stores the fitness values of every Red Deer

double LB=-1, UB=1; // Lower and Upper Bound used in the algorithm

// Structure for City in TSP problem
struct City{

    int id; // Stores the id of the city (Unique for every city)
    int x,y; // Coordinates of the city

    // Calculates distance between two cities
    double distance(City c){

        double xdis=abs(x-c.x);
        double ydis=abs(y-c.y);
        double distance=sqrt(xdis*xdis+ydis*ydis);

        return distance;

    }

    // Prints the coordinates of the city
    void print_coord(){
        cout<<"("<<x<<", "<<y<<")\n";
    }

};

// Calculates the total distance of the route taken
double routeDistance(vector<double>route, City index[]){

    double path_dist=0,fitness=0.0;

    vector<double>sorted_RD=route;
    sort(sorted_RD.begin(), sorted_RD.end());
    map<double,int>RD_to_city;

    for(int i=0;i<sorted_RD.size();i++){
        RD_to_city[sorted_RD[i]]=i;
    }

    // TSP starts and ends at the same place. So the initial city is inserted again
    route.pb(route[0]);

    // Calculates pairwise distance
    for(int i=0;i<route.size()-1;i++){

        City fromCity= index[RD_to_city[route[i]]];
        City toCity= index[RD_to_city[route[i+1]]];

        path_dist+=fromCity.distance(toCity);

    }

    return path_dist;

}

/* Returns fitness value of given route.
We aim to minimize the distance. Here, I am trying to minimize the fitness too
*/
double routeFitness(vector<double>route, City index[]){

    double fitness=routeDistance(route, index);
    return fitness;

}


// Randomly creates initial population of Red deer using randomizer
vector<vector<double>> initialPopulation(int pop_size,int no_of_cities){

    vector<vector<double>>RD_list;

    RD_list.resize(pop_size);

    #pragma omp parallel for
    for(int i=0;i<pop_size;i++){
        RD_list[i].resize(no_of_cities);
        for(int j=0;j<no_of_cities;j++){
            RD_list[i][j]=get_rand();
        } 
    }

    return RD_list;

}

// Sorts the population in increasing fitness score and returns in a new vector
vector<vector<double>> rankRDs(vector<vector<double>>population, City index[]){

    vector<pair<double,int>>fitnessResults; // Stores fitness value and inital route id
    fitnessResults.resize(population.size());
    fitness_scores.resize(population.size());

    #pragma omp parallel for
    for(int i=0;i<population.size();i++){
        fitnessResults[i].F=routeFitness(population[i],index);
        fitnessResults[i].S=i;
    }

    // Sorts in ascending order wrt fitness value
    sort(fitnessResults.begin(),fitnessResults.end());

    vector<vector<double>>sorted_population; // Will store population sorted in increasing fitness value
    sorted_population.resize(population.size());

    #pragma omp parallel for
    for(int i=0;i<fitnessResults.size();i++){

        sorted_population[i]=population[fitnessResults[i].S];
        fitness_scores[i]=(fitnessResults[i].F);

    }

    return sorted_population;

}

// Separates males and hinds by separating 'no_of_males' Red deers with minimum fitness values
vector<vector<vector<double>>> Separate_male_and_hinds(vector<vector<double>>population, int no_of_males, int no_of_hinds, City index[]){

    vector<vector<vector<double>>>Separate;
    Separate.resize(2);

    Separate[0].resize(no_of_males);
    Separate[1].resize(no_of_hinds);

    for(int i=0;i<no_of_males;i++)
        Separate[0][i]=population[i]; // Store males
    

    for(int i=no_of_males;i<population.size();i++)
        Separate[1][i-no_of_males]=population[i]; // Stores Hinds

    return Separate;
}

/* Simulating the roaring process.
Finds the neighbours of male RD and if their OF is better, we replace them with the male_new
*/
vector<vector<double>> roar(vector<vector<double>>males, City index[]){

    double a1=get_rand(),a2=get_rand(),a3=get_rand();

    for(int i=0;i<males.size();i++){

        vector<double>male_new(males[i].size());

        // Updating the position of males
        for(int j=0;j<males[i].size();j++){
            if(a3>=0.5)
                male_new[j]=males[i][j]+a1*((UB-LB)*a2)+LB;
            else
                male_new[j]=males[i][j]-a1*((UB-LB)*a2)+LB;

        }
        if(routeFitness(male_new,index)<routeFitness(males[i],index))
            males[i]=male_new;
    }

    return males;

}

// Separates commanders and stags by separating gamma*(no_of_males) with minimum fitness into commanders
vector<vector<vector<double>>> Separate_commanders_and_stags(vector<vector<double>>males, double gamma){

    int no_of_commanders=gamma*males.size();
    vector<vector<vector<double>>>Separate;
    Separate.resize(2);

    Separate[0].resize(no_of_commanders);
    Separate[1].resize(males.size()-no_of_commanders);

    for(int i=0;i<no_of_commanders;i++)
        Separate[0][i]=males[i];
    

    for(int i=no_of_commanders;i<males.size();i++)
        Separate[1][i-no_of_commanders]=males[i];

    return Separate;
}

/*
Each commander fights with every stag and two new solutions are generated. The best among the 
commander and generated solutions is assigned the new commander
*/
void Fight(vector<vector<double>>&commanders, vector<vector<double>>stags, City index[]){

    int no_of_stags=stags.size();
    int sz=commanders[0].size();

    for(int i=0;i<commanders.size();i++){
       
        for(int j=0;j<stags.size();j++){
        vector<double> New1(sz),New2(sz);   // Two new solutions
        double b1=get_rand(),b2=get_rand();

            for(int k=0;k<sz;k++){
                New1[k]=(commanders[i][k]+stags[j][k])/2.0+b1*((UB-LB)*b2)+LB;
                New2[k]=(commanders[i][k]+stags[j][k])/2.0-b1*((UB-LB)*b2)+LB;
            }

            // Calculating minimum fitness and 
            double fitness1=routeFitness(commanders[i], index);
            double fitness2=routeFitness(New1, index);
            double fitness3=routeFitness(New2, index);
            double minfitness=min(fitness1,min(fitness2,fitness3));

            if(fitness2==minfitness)commanders[i]=New1;
            else if(fitness3==minfitness)commanders[i]=New2;
        }
    }

}

// Forming harems. Hinds are divided proportionally among male commanders according to their power
vector<vector<vector<double>>> form_harems(vector<vector<double>>&commanders, vector<vector<double>>hinds, City index[]){

    double tot_fitness=0;   // Stores total fitness
    shuffle(hinds.begin(),hinds.end(),std::default_random_engine(rand()));

    // Stores fitness values
    vector<double>fitness;
    fitness.resize(commanders.size());

    for(int i=0;i<commanders.size();i++){
        fitness[i]=routeFitness(commanders[i], index);
        tot_fitness+=fitness[i];
    }

    // 
    vector<vector<vector<double>>>harems; // Stores the hinds in each harem
    harems.resize(commanders.size());

    int sz=commanders[0].size();
    int hinds_taken=0;  // No of hinds assigned to harems
    int no_of_hinds=hinds.size();   // Total no of hinds
    int harem_size=(fitness[0]/tot_fitness)*no_of_hinds;    // No of hinds in next harem

    // Proportionally assigning harems
    for(int i=0;i<commanders.size();i++){
        harems[i].resize(harem_size);
        for(int j=hinds_taken;j<hinds_taken+harem_size;j++){

            harems[i][j-hinds_taken]=(hinds[j]);

        }
        hinds_taken+=harem_size;
        if(i<commanders.size()-2)
            harem_size=(fitness[i+1]/tot_fitness)*no_of_hinds;
        else
            harem_size=no_of_hinds-hinds_taken;
    }

return harems;

}

/*
Mating commander of a harem with alpha percent of hinds in his harem
Mating commander of a harem with beta percent of hinds in another random harem
*/
vector<vector<double>> Mate_alpha_beta(vector<vector<double>>commanders, vector<vector<vector<double>>>harems, double alpha, double beta){

    int sz=commanders[0].size();
    int c=get_rand();
    vector<vector<double>>offs;

    #pragma omp parallel
    {
        vector<vector<double>>offs_temp;
        for(int i=0;i<commanders.size();i++){
            for(int j=0;j<alpha*harems[i].size();j++){

                c=get_rand();
                vector<double>child(sz);
                for(int k=0;k<sz;k++){
                    child[k]=(commanders[i][k]+harems[i][j][k])/2+(UB-LB)*c;
                }
                offs_temp.pb(child);
            }
        }

        #pragma omp critical
        offs.insert(offs.end(),offs_temp.begin(),offs_temp.end());
    }

    #pragma omp parallel
    {
        vector<vector<double>>offs_temp;
        for(int i=0;i<commanders.size();i++){

            int w=rand()%harems.size(); // Randomly selected harem
            while(w==i) w=rand()%harems.size();

            for(int j=0;j<beta*harems[w].size();j++){
                c=get_rand();
                vector<double>child(sz);
                for(int k=0;k<sz;k++){
                    child[k]=(commanders[i][k]+harems[w][j][k])/2+(UB-LB)*c;
                }
                offs_temp.pb(child);
            }
        }

        #pragma omp critical
        offs.insert(offs.end(),offs_temp.begin(),offs_temp.end());
    }

return offs;
}

// Mating stag with the nearest hind
vector<vector<double>> Mate_stag_hind(vector<vector<double>>stags, vector<vector<double>>hinds, City index[]){

    vector<vector<double>>offs;

    // Finding hind with minimum distance from 'stags[i]' stag and mating with it
    #pragma omp parallel
    {
        vector<vector<double>>offs_temp;

        #pragma omp for nowait
        for(int i=0;i<stags.size();i++){
            int cindx=-1,cdist=inf;
            for(int j=0;j<hinds.size();j++){

                vector<double>diff(stags[i].size());

                for(int k=0;k<stags[i].size();k++)
                    diff[k]=abs(stags[i][k]-hinds[j][k]);

                int diff_fitness=routeFitness(diff, index);
                if(diff_fitness < cdist)
                    cdist=diff_fitness,cindx=j;

            }

            // Generating offspring with the selected hind
            int j=cindx;
            double c=get_rand();
            vector<double>child(stags[i].size());
            for(int k=0;k<stags[i].size();k++){
                child[k]=(stags[i][k]+hinds[j][k])/2+(UB-LB)*c;
            }
            offs_temp.pb(child);

        }
        #pragma omp critical
        offs.insert(offs.end(), offs_temp.begin(), offs_temp.end());
    }

return offs;
}

// Selecting the next generation using elitism and roulette wheel mechanism
vector<vector<double>> NextGen(vector<vector<double>>offs, vector<vector<double>>males, int pop_size, City index[]){

    vector<vector<double>>nextGen;

    // All the males are selected in the next generation
    for(int i=0;i<males.size();i++){
        nextGen.pb(males[i]);
    }

    vector<double>fitnessscores;
    fitnessscores.resize(offs.size());

    double tot_fitness=0; // Stores total fitness

    #pragma omp parallel for reduction(+:tot_fitness)
    for(int i=0;i<offs.size();i++){
        fitnessscores[i]=routeFitness(offs[i], index);
        tot_fitness+=fitnessscores[i];
    }
        
    // Assigning fitness weighed probability 
    for(int i=0;i<fitnessscores.size();i++){

        fitnessscores[i]=100*((fitnessscores[i]/tot_fitness));
        if(i!=0) fitnessscores[i]+=fitnessscores[i-1];

    }

    // Using fitness proportionate selection (Roulette Wheel mechanism)
    #pragma omp parallel
    {   
        vector<vector<double>>nextGen_temp;

        #pragma omp for nowait
        for(int i=0;i<pop_size-males.size();i++){

            double pick=get_rand();
            pick*=100;
            for(int j=0;j<fitnessscores.size();j++){

                if(j==fitnessscores.size()-1){
                    nextGen_temp.pb(offs[j]);
                    break;
                }
                if(pick>100-fitnessscores[j]){
                    // cout<<j<<endl;
                    nextGen_temp.pb(offs[j]);
                    break;
                } 

            }
        }

        #pragma omp critical
        nextGen.insert(nextGen.end(),nextGen_temp.begin(),nextGen_temp.end());

    }

return nextGen;
}

// Using all the above functions to get next generation
vector<vector<double>> nextGeneration(vector<vector<double>>population, int popSize, int no_of_males, double alpha, double beta, double gamma, City index[]){

    vector<vector<double>>popRanked=rankRDs(population, index);

    // cout<<routeFitness(popRanked[0], index)<<", ";

    vector<vector<vector<double>>>separate = Separate_male_and_hinds(popRanked, no_of_males, popSize-no_of_males, index);
    vector<vector<double>>Males=separate[0];
    vector<vector<double>>Hinds=separate[1];

    Males=roar(Males, index);
    Males=rankRDs(Males, index);

    separate=Separate_commanders_and_stags(Males, gamma);
    vector<vector<double>>Commanders=separate[0];
    vector<vector<double>>Stags=separate[1];

    Fight(Commanders, Stags, index);

    vector<vector<vector<double>>>harems=form_harems(Commanders, Hinds, index);

    vector<vector<double>>offs=Mate_alpha_beta(Commanders, harems, alpha, beta);

    vector<vector<double>>offs1=Mate_stag_hind(Stags, Hinds, index);

    offs.insert(offs.end(),offs1.begin(),offs1.end());

    offs.insert(offs.end(),Commanders.begin(),Commanders.end());
    offs.insert(offs.end(),Stags.begin(),Stags.end());
    offs.insert(offs.end(),Hinds.begin(),Hinds.end());

    vector<vector<double>>nextGen=NextGen(offs, Males, popSize, index);

return nextGen;
}

void RD_Algo( vector<City>cityList, int popSize, int no_of_males, int generations,double alpha, double beta, double gamma, City index[]){

    // Creating initial population from city list
    vector<vector<double>>population;

    population=initialPopulation(popSize, cityList.size());

    cout<<"Initial Distance was: "<<routeDistance(population[0], index)<<endl;

    for(int i=0;i<generations;i++){

        fitness_scores.clear();
        population=nextGeneration(population, popSize, no_of_males, alpha, beta, gamma, index);

    }
    cout<<"Final Distance is: "<<routeDistance(population[0], index)<<endl;

}

int main(int argc, char **argv){
    ios_base::sync_with_stdio(false);cin.tie(NULL);cout.tie(NULL);
    
    omp_set_num_threads(atoi(argv[1]));

    int no_of_cities=40; // No of cities
    int popSize=100;    // Population Size
    int no_of_males=30; // No of males
    int generations=100;   // No of generations

    double alpha=0.9;
    double beta=0.4;
    double gamma=0.7;

    int no_of_hinds = popSize-no_of_males;

    

    vector<City>cityList; // Stores the initial list of cities
    City index[no_of_cities];

    // Assigning random coordinates to cities
    for(int i=0;i<no_of_cities;i++){
        double x=200*((double)rand()/RAND_MAX);
        double y=200*((double)rand()/RAND_MAX);
        City c;
        c.x=x,c.y=y,c.id=i;
        cityList.pb(c);
        index[i]=c; // Assigning index to cities
    }

    auto begin = chrono::high_resolution_clock::now();

    RD_Algo(cityList, popSize, no_of_males, generations, alpha, beta, gamma, index);

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end-begin);
    cout<<fixed<<setprecision(5) <<"Time Taken by serial program: "<< duration.count()/1000000.0 <<" s "<<endl;
    
return 0;
}

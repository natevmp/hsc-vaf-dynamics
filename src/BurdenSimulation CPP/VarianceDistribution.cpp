
#include "MersenneTwister.h"


# include <iostream>
# include <cmath>
# include <ctime>
# include <cstdlib>
# include <cstdio>
# include <vector>
# include <fstream>
#include <math.h>
#include <random>

using std::vector;
using std::cout;
using std::cin;



int cells = 50000;

int divisions = 350;

double asymprob = 1;

double mutationrate = 3;

int runs = 1000;



time_t seconds ;
MTRand mtrand1(time(NULL));

std::random_device rd;
std::mt19937 gen(rd());


std::default_random_engine generator;
std::poisson_distribution<int> distribution(mutationrate);



vector <double> State (cells,0);

vector <double> Mean (runs,0);

vector <double> Variance (runs,0);





int main()
{
    int tick = 0;
    
  while(tick < runs)
  {
    
    int count = 0;
    
    

    double s0 = 0;
    double s1 = 0;
    double s2 = 0;
      
      
      for ( int i=0; i<cells; i++)
      {
          State.at(i)=0;
      }
      
    
      
    
    while(count < cells * divisions)
    {
        
     
    
    s0 = mtrand1.randExc();
        
        if(s0 < asymprob)
        {
            s1 = mtrand1.randInt(cells-1);
            int number = distribution(generator);
            State.at(s1) = State.at(s1)+number;
        }
        else
        {
            s1 = mtrand1.randInt(cells-1);
            s2 = mtrand1.randInt(cells-1);
            int number = distribution(generator);
            State.at(s1) = State.at(s2)+number;
            
        }
    
        

        count++;
        
    }
      
    
                double MeanData = 0;
                double SumOfSquares = 0;
                
                
                for (int i=0; i<cells; i++)
                {
                    MeanData = MeanData + State.at(i);
                }
                 
                MeanData = MeanData/(cells-1) ;
                
                cout << tick << "\n";
                
                Mean.at(tick) = MeanData;
                
                for (int i=0; i< cells; i++)
                {
                    SumOfSquares = SumOfSquares + pow((State.at(i) - MeanData),2);
                }
                
                Variance.at(tick) = SumOfSquares/(cells-1);
        
      
    
   
      tick++;
  }
    
    std::fstream datei ;
    datei.open ("MeanDistribution.txt" ,std::ios::out);
    for ( int i=0 ; i<1 ; i++)
    { if (i!=0)
    { datei << "\n" ;}
        for ( int j=0 ; j< runs ; j++)
            datei << Mean.at(j) << " " ;}
    datei.close() ;
    
    datei.open ("VarianceDistribution.txt" ,std::ios::out);
    for ( int i=0 ; i<1 ; i++)
    { if (i!=0)
    { datei << "\n" ;}
        for ( int j=0 ; j< runs ; j++)
            datei << Variance.at(j) << " " ;}
    datei.close() ;
    
    datei.open ("VarianceMeanDistribution_N50000_p1_mu3.txt" ,std::ios::out);
    for ( int i=0 ; i<1 ; i++)
    { if (i!=0)
    { datei << "\n" ;}
        for ( int j=0 ; j< runs ; j++)
            datei << Variance.at(j)/Mean.at(j) << " " ;}
    datei.close() ;
    
}












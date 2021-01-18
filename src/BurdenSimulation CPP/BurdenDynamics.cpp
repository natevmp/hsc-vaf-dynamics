

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



int cells = 100;

int divisions = 3000;

double asymprob = 0.9;

int runs = 150;


double exprand( double lambda );

time_t seconds ;
MTRand mtrand1(time(NULL));

std::random_device rd;
std::mt19937 gen(rd());


std::default_random_engine generator;
std::poisson_distribution<int> distribution(1.3);



vector <double> State (cells,0);

vector <double> Mean (divisions,0);

vector <double> Variance (divisions,0);

vector < vector <double> > SummaryMean (runs ,vector <double> (divisions,0));

vector < vector <double> > SummaryVariance (runs ,vector <double> (divisions,0));



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
      
      for ( int i=0; i<divisions; i++)
      {
          Mean.at(i)=0;
          Variance.at(i)=0;
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
    
        
        if ((count % cells) == 0)
        {
            double MeanData = 0;
            double SumOfSquares = 0;
            
            
            for (int i=0; i<cells; i++)
            {
                MeanData = MeanData + State.at(i);
            }
             
            MeanData = MeanData/(cells-1) ;
            
            cout << MeanData << "\n";
            
            Mean.at((double)count/cells) = MeanData;
            
            for (int i=0; i< cells; i++)
            {
                SumOfSquares = SumOfSquares + pow((State.at(i) - MeanData),2);
            }
            
            Variance.at((double)count/cells) = SumOfSquares/(cells-1);
        }
        
        
        
        
    

        count++;
        
    }
      
      for (int i=0; i<divisions;i++)
      {
          SummaryMean.at(tick).at(i)=Mean.at(i);
          SummaryVariance.at(tick).at(i)=Variance.at(i);
      }
   
      tick++;
  }
    
    std::fstream datei ;
    datei.open ("Burden.txt" ,std::ios::out);
    for ( int i=0 ; i<1 ; i++)
    { if (i!=0)
    { datei << "\n" ;}
        for ( int j=0 ; j< cells ; j++)
            datei << State.at(j) << " " ;}
    datei.close() ;
    
    datei.open ("Mean.txt" ,std::ios::out);
    for ( int i=0 ; i<1 ; i++)
    { if (i!=0)
    { datei << "\n" ;}
        for ( int j=0 ; j< divisions ; j++)
            datei << Mean.at(j) << " " ;}
    datei.close() ;
    
    datei.open ("Variance.txt" ,std::ios::out);
    for ( int i=0 ; i<1 ; i++)
    { if (i!=0)
    { datei << "\n" ;}
        for ( int j=0 ; j< divisions ; j++)
            datei << Variance.at(j) << " " ;}
    datei.close() ;
    
    datei.open ("MeanSummary.txt" ,std::ios::out);
    for ( int i=0 ; i<runs ; i++)
    { if (i!=0)
    { datei << "\n" ;}
        for ( int j=0 ; j< divisions ; j++)
            datei << SummaryMean.at(i).at(j) << " " ;}
    datei.close() ;
    
    datei.open ("VarianceSummary.txt" ,std::ios::out);
    for ( int i=0 ; i<runs ; i++)
    { if (i!=0)
    { datei << "\n" ;}
        for ( int j=0 ; j< divisions ; j++)
            datei << SummaryVariance.at(i).at(j) << " " ;}
    datei.close() ;
    
}


double exprand( double lambda)
{
    double y = mtrand1.randExc();
    return (-log(1.-y)/lambda);
}










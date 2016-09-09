#ifndef __JACK_STATS_H_INCLUDED__
#define __JACK_STATS_H_INCLUDED__

#include <cmath>
#include <vector>

// Some (not very clever) utility functions for jackknifing numeric data.
// Assumes an Nj x Nd matrix where the rows index the jackknife 
// sample and the columns index the data.


double jack_avg(std::vector<double>& data)
{
  int Njks = data.size();
  double avg(0.0);
  
  for(int i=0; i<Njks; ++i){ avg += data[i]; }
  
  return avg / static_cast<double>(Njks);
}

std::vector<double> jack_avg(std::vector<std::vector<double>>& data)
{
  int Njks = data.size();
  int Ndat = data[0].size();
  std::vector<double> avg(Ndat, 0.0);
  
  for(int i=0; i<Njks; ++i){ 
    for(int j=0; j<Ndat; ++j){ 
      avg[j] += data[i][j] / static_cast<double>(Njks);
    }
  }
  
  return avg;
}

double jack_std(std::vector<double>& data, double& avg, bool resample)
{
  int Njks = data.size();
  double std(0.0);
  std::vector<double> jack_samples(Njks);
  
  // Compute jackknife samples: 
  // jack_samples[i] is the mean of data without data[i]
  if(resample){
    for(int i=0; i<Njks; ++i){
      std::vector<double> this_jack(Njks-1);
      for(int j=0; j<Njks; ++j){
        if(j < i){ this_jack[j] = data[j]; }
        else if(j > i){ this_jack[j-1] = data[j]; }
      }
      jack_samples[i] = jack_avg(this_jack);
    }
  } else { 
    jack_samples = data; 
  }
  
  // Compute jackknife variance estimate from jack_samples and avg
  double Z = static_cast<double>(Njks) - 1.0; 
  Z /= static_cast<double>(Njks);
  for(int i=0; i<Njks; ++i){
    std += Z*pow(jack_samples[i]-avg, 2.0);
  }
  
  // We want the std
  return sqrt(std);
}

std::vector<double> jack_std(std::vector<std::vector<double>>& data, std::vector<double>& avg, bool resample)
{
  int Njks = data.size();
  int Ndat = data[0].size();
  std::vector<double> std(Ndat, 0.0);
  std::vector<std::vector<double>> jack_samples(Njks, std::vector<double>(Ndat,0.0));
  
  // Compute jackknife samples: 
  // jack_samples[i] is the row-wise mean of data without row data[i]
  if(resample){
    for(int i=0; i<Njks; ++i){
      std::vector<std::vector<double>> this_jack(Njks-1,std::vector<double>(Ndat,0.0));
      for(int j=0; j<Njks; ++j){
        if(j < i){ this_jack[j] = data[j]; }
        else if(j > i){ this_jack[j-1] = data[j]; }
      }
      jack_samples[i] = jack_avg(this_jack);
    }
  } else { 
    jack_samples = data; 
  }
  
  // Compute jackknife variance estimate from jack_samples and avg
  double Z = static_cast<double>(Njks) - 1.0; 
  Z /= static_cast<double>(Njks);
  for(int i=0; i<Njks; ++i){
    for(int j=0; j<Ndat; ++j){
      std[j] += Z*pow(jack_samples[i][j]-avg[j], 2.0);
    }
  }
  
  // We want the std
  for(int i=0; i<Ndat; ++i){ std[i] = sqrt(std[i]); }
  return std;
}

#endif
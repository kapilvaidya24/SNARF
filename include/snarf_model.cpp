#include<iostream>
#include<algorithm>
#include<cmath>
#include <random>
#include <fstream>
#include <chrono>
#include <set>
#include <ctime>
#include <cassert>
#include <cstring>
using namespace std;
using namespace std::chrono; 

//Implementation of the model used in SNARF
template <class T>
struct snarf_model
{
  double level_0_slope,level_0_bias;
  int num_models=1000000;
  vector<T> first_level;
  vector<double> level_1_slope,level_1_bias;

  // Generates Slopes and Biases of linear models in level 1
  void generate_slope_bias_level_1(vector<T> &keys,vector<double> &ecdf)
  {

    uint64_t N=keys.size();
    
    level_1_bias.resize(num_models,0.0);
    level_1_slope.resize(num_models,0.0);

    vector<T> max_val_vec(num_models,0),min_val_vec(num_models,0);
    vector<int> count_items_model(num_models,0);
    vector<double> max_cdf_vec(num_models,0.0),min_cdf_vec(num_models,0.0);

    //find min/max elements per linear model 
    for(int i=0;i<N;i++)
    {
      double est_cdf;
      T consider;
      int bin_index=binary_search(keys[i]);

      if(keys[i]>=first_level[bin_index])
      {
        consider=0;
        est_cdf=num_models-1; 
      }
      else
      {
        consider=first_level[bin_index]-keys[i];
        est_cdf=bin_index;
      }

      int est_pos=floor(est_cdf);
      est_pos=max(0,est_pos);
      est_pos=min(num_models-1,est_pos);

      if(count_items_model[est_pos]==0)
      {
        max_val_vec[est_pos]=consider;
        min_val_vec[est_pos]=consider;
        max_cdf_vec[est_pos]=ecdf[i];
        min_cdf_vec[est_pos]=ecdf[i];
      }
      else
      {
        min_val_vec[est_pos]=min(min_val_vec[est_pos],consider);
        max_val_vec[est_pos]=max(max_val_vec[est_pos],consider);
        max_cdf_vec[est_pos]=max(max_cdf_vec[est_pos],ecdf[i]);
        min_cdf_vec[est_pos]=min(min_cdf_vec[est_pos],ecdf[i]);
      }

      count_items_model[est_pos]++;
    }
    

    vector<T> new_max_val_vec(num_models,0),new_min_val_vec(num_models,0);
    //needed to handle empty models 
    for(int i=0;i<num_models;i++)
    {
      if(i==0)
      {
        new_max_val_vec[i]=first_level[i]-keys[0];
      }
      else
      {
        new_max_val_vec[i]=first_level[i]-first_level[i-1];
      }
    }

    //get final slope bias values
    double temp_slope,temp_bias;
    for(int i=0;i<num_models;i++)
    {
      uint64_t diff=(new_max_val_vec[i]-new_min_val_vec[i]);
      temp_slope=(max_cdf_vec[i]-min_cdf_vec[i])*1.00/((diff)*1.00);

      temp_bias=max_cdf_vec[i]-(new_min_val_vec[i]*temp_slope);
      level_1_slope[i]=temp_slope;
      level_1_bias[i]=temp_bias;
     
    }

    return ;
  }

  //Builds snarf model given array of keys
  void snarf_model_builder(vector<T> &keys)
  {
    uint64_t N=keys.size();

    bool testbool = (N>10000);
    assert(("The number of keys are smaller than 10000, difficult to build a model for this!", testbool));

    //Number of models used is num_keys/10000.0. VARY THIS PARAMETER TO GET BETTER PRECISION
    num_models=ceil(N/10000.0);

    sort(keys.begin(),keys.end());

    vector<double> ecdf(N,0.0);

    //Build level 0 and the cdf
    for(int i=0;i<keys.size();i++)
    {
      ecdf[i]=i*1.00/N;
    }

    first_level.resize(num_models,0.0);

    for(int i=0;i<num_models;i++)
    {
      first_level[i]=keys[int(((i+1)*N*1.00)/num_models)-1];
    }

    first_level[num_models-1]=keys[N-1];

    //build level 1
    generate_slope_bias_level_1(keys,ecdf);

    return ;
  }

  //binary searching the first level to obtain the index of the linear model in level 1.
  int binary_search(T key)
  {
    int start=0,end=first_level.size()-1;
    int mid=(start+end)/2;

    while((end-start)>10)
    {
      if (first_level[mid]<key)
      {
        start=mid;
      }
      else
      {
        end=mid;
      }
      mid=(start+end)/2;

    }

    for(int i=start;i<=end;i++)
    {
      if(first_level[i]>=key)
      {
        return i;
      }
    }

    return first_level.size()-1;


  }

  //get the estimated cdf for a key
  double infer(T key)
  {

    double est_cdf;
    T consider;

    // Get the index of the level 1 model
    int index=binary_search(key);

    if(key>=first_level[index])
    {
      consider=0;
      est_cdf=num_models-1; 
    }
    else
    {
      consider=first_level[index]-key;
      est_cdf=index;
    }


    int est_pos=floor(est_cdf);
    est_pos=max(0,est_pos);
    est_pos=min(num_models-1,est_pos);

    // Use the level 1 model
    double ans=level_1_bias[est_pos]-level_1_slope[est_pos]*consider;
   
    ans=max(0.0,ans);
    ans=min(1.0,ans);

    return ans;
  }

  // Returns the size used by the snarf_model in bytes
  int return_size()
  {
    double a=10.09,b;

    int size_of_template=sizeof(first_level[0]);
    int total_size=0;
    total_size+=sizeof(num_models);
    total_size+=sizeof(size_of_template);
    total_size+=num_models*sizeof(first_level[0]);
    total_size+=2*num_models*sizeof(a);
    
    return total_size; 
  }

  // writes the model contents into a char array
  void serialize(unsigned char* arr)
  {

    double a=0.0;
    int size_of_template=sizeof(first_level[0]);

    int total_size=0;
    total_size+=sizeof(num_models);
    total_size+=sizeof(size_of_template);
    total_size+=num_models*sizeof(first_level[0]);
    total_size+=2*num_models*sizeof(a);

    int offset=0;
    memcpy(arr+offset,&num_models,sizeof(num_models));
    offset+=sizeof(num_models);

    memcpy(arr+offset,&size_of_template,sizeof(num_models));
    offset+=sizeof(size_of_template);

    for(int i=0;i<num_models;i++)
    {
      memcpy(arr+offset,&first_level[i],sizeof(first_level[0]));
      offset+=sizeof(first_level[0]);
    }

    for(int i=0;i<num_models;i++)
    {
      memcpy(arr+offset,&level_1_slope[i],sizeof(level_1_slope[0]));
      offset+=sizeof(level_1_slope[0]);
    }

    for(int i=0;i<num_models;i++)
    {
      memcpy(arr+offset,&level_1_bias[i],sizeof(level_1_bias[0]));
      offset+=sizeof(level_1_bias[0]);
    }


    return ;
  }

  // reads the model contents from a char array
  void deserialize(unsigned char* arr)
  {
    double a=10.09,b;
    int temp_int=0;
    int offset=0;
    int size_of_template=0;

    memcpy(&num_models,arr+offset,sizeof(temp_int));
    offset+=sizeof(temp_int);

    memcpy(&size_of_template,arr+offset,sizeof(temp_int));
    offset+=sizeof(temp_int);

    first_level.resize(num_models,0);
    level_1_slope.resize(num_models,0);
    level_1_bias.resize(num_models,0);



    for(int i=0;i<num_models;i++)
    {
      memcpy(&first_level[i],arr+offset,sizeof(first_level[0]));
      offset+=sizeof(first_level[0]);
    }

    first_level.resize(1,0);

    for(int i=0;i<num_models;i++)
    {
      memcpy(&level_1_slope[i],arr+offset,sizeof(level_1_slope[0]));
      offset+=sizeof(level_1_slope[0]);
    }

    for(int i=0;i<num_models;i++)
    {
      memcpy(&level_1_bias[i],arr+offset,sizeof(level_1_bias[0]));
      offset+=sizeof(level_1_bias[0]);
    }

    return ;
  }


};


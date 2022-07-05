
#include<iostream>
#include<algorithm>
#include<cmath>
#include <random>
#include <fstream>
#include <chrono>
#include <set>
#include <ctime>
#include <cstring>
using namespace std;
using namespace std::chrono; 
#include "include/snarf.cpp"

int main()
{
  //----------------------------------------
  //GENERATING DATA
  //----------------------------------------

  uint64_t N=10'000'000;  
  vector<uint64_t> v_keys(N,0);

  //Key are multiples of 10000. S= {10000, 20000, 30000,...}
  for(uint64_t i=0;i<N;i++)
  {
      v_keys[i]=i*10000;
  }

  //----------------------------------------
  //SNARF CONSTRUCTION
  //----------------------------------------

  //Initializing variables needed for SNARF construction
  //Bits per Key indicates how many bits per key should the snarf instance use
  double bits_per_key=10.00;
  //batch_size indicates what should be the split of the compressed bit array
  //Smaller batch sizes increases the size of the snrf instance while making range queries faster.
  // We recommend using batch_sz ~ 100-200
  int batch_size=100;

  //declare and initialize a snarf instance
  snarf_updatable_gcs<uint64_t> snarf_instance;
  snarf_instance.snarf_init(v_keys,bits_per_key,batch_size);
    
  //get the size of the snarf instance
  int snarf_sz=snarf_instance.return_size();
  cout<<"Bits per key used by SNARF: "<<snarf_sz*8.00/v_keys.size()<<endl;
  



  //----------------------------------------
  //QUERYING SNARF
  //----------------------------------------

  // SNARF suppports 3 main operations
  // Range Query: checks the existence of keys in a range
  // Insert: inserts the key
  // Delete: deletes the key

  uint64_t left_end=15000;
  uint64_t right_end=16000;

  //range query example
  if(snarf_instance.range_query(left_end,right_end))
  {
    cout<<"False Positive for: [ "<<left_end<<", "<<right_end<<"]"<<endl;
  }
  else
  {
     cout<<"True Negative for: [ "<<left_end<<", "<<right_end<<"]"<<endl;
  }

  //Insert example
  snarf_instance.insert_key(15000);
  if(snarf_instance.range_query(left_end,right_end))
  {
    cout<<"True Positive for: [ "<<left_end<<", "<<right_end<<"]"<<endl;
  }
  else
  {
    //We won't reach here as false negatives do not occur in filters.
     cout<<"False Negative for: [ "<<left_end<<", "<<right_end<<"]"<<endl;
  }

  //Delete example
  snarf_instance.delete_key(15000);
  if(snarf_instance.range_query(left_end,right_end))
  {
    cout<<"False Positive for: [ "<<left_end<<", "<<right_end<<"]"<<endl;
  }
  else
  {
     cout<<"True Negative for: [ "<<left_end<<", "<<right_end<<"]"<<endl;
  }

  uint64_t false_positive_left_end=10001;

  //false positive example
  if(snarf_instance.range_query(false_positive_left_end,right_end))
  {
    cout<<"False Positive for: [ "<<false_positive_left_end<<", "<<right_end<<"]"<<endl;
  }
  else
  {
     cout<<"True Negative for: [ "<<false_positive_left_end<<", "<<right_end<<"]"<<endl;
  }



    return 0;
}
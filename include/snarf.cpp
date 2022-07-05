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

#include "snarf_model.cpp"
#include "snarf_bitset.cpp"

//SNARF implementation which is updatable(handles deletes and inserts) and uses Golomb Coding(GCS)
template <class T>
struct snarf_updatable_gcs
{
  //Snarf model
  snarf_model<T> rmi;

  //snarf bit array
  vector< snarf_bitset > bb_bitset_vec;

  //Parameters used in snarf
  uint64_t N,P,block_size,bit_size,total_blocks;
  uint64_t gcs_size;

  //Stores the number of keys in each bit array block
  vector<int> vec_num_keys;
  

  //Used for Inference
  uint64_t delta_query_index=0;
  uint64_t delta_query_remainder=0;
  uint64_t delta_query_val_until_now=0;
  uint64_t delta_query_data_size=0;
  uint64_t delta_query_offset_1=0;
  uint64_t delta_query_offset_2=0;
  uint64_t delta_query_zero_count=0;
  uint64_t delta_query_one_count=0;
  double query_cdf1,query_cdf2=0.0;

  //Name of snarf instance
  char name_curr;

  snarf_updatable_gcs<T>(){

    N=0;
    P=0;
    block_size=0;
    bit_size=0;
    total_blocks=0;
    gcs_size=0;
    name_curr=' ';
  }

 


  //Get the bit locations of the bits that need to be set to 1
  void get_locations(vector<T> &keys,vector<uint64_t> &temp_locations)
  {
    double cdf;
    uint64_t temp_loc;

    sort(keys.begin(),keys.end());

    //Iterate over the key array to find the bit location corresponding to each key
    uint64_t past_loc=0;
    for(int i=0;i<keys.size();i++)
    { 
      cdf=rmi.infer(keys[i]);
      
      temp_loc=floor(cdf*N*P);
      temp_loc=max((uint64_t)0,temp_loc);
      temp_loc=min(N*P-1,temp_loc);
      temp_locations.push_back(temp_loc);
      past_loc=temp_loc;
    }


    sort(temp_locations.begin(),temp_locations.end());

    return ;
  }

  //Create a new bit block(var bb_temp) for a batch of values(curr_batch).
  void create_new_gcs_block(vector<uint64_t> &curr_batch, snarf_bitset &bb_temp)
  {
    // resize the bit array for this new batch
    bb_temp.init((bit_size+1)*curr_batch.size()+block_size);

    int offset_bits=0;

    //Write the binary code in the bit array
    for(int i=0;i<curr_batch.size();i++)
    {
      bb_temp.bitset_write_bits(offset_bits,curr_batch[i]%P,bit_size);
      offset_bits+=bit_size;
    }

    //Write the unary code in the bit array
    int delta_zero_count=0,delta_one_count=0;
    for(int i=0;i<curr_batch.size();i++)
    {

      uint64_t temp=curr_batch[i];

      temp=temp/P;


      if(temp>delta_zero_count)
      {
        bb_temp.bitset_write_bits(offset_bits,0,1);
        offset_bits++;
        delta_zero_count++;
        i--;
      }
      else
      {
        bb_temp.bitset_write_bits(offset_bits,1,1);
        offset_bits++;
      }

    }

    return ;
  }


  //Creates the array of bit vectors
  //It batches values corresponding to each bit block and then calls "create_new_gcs_block" to create the actual blocks
  uint64_t build_bb(vector<uint64_t> &temp_locations)
  {
    vector<uint64_t> curr_batch;
    uint64_t num_batches=ceil(temp_locations.size()*1.00/block_size);

    uint64_t curr_index=0;
    uint64_t total_bits_used=0;
    uint64_t curr_size=0;


    for(int i=0;i<num_batches;i++)
    {
      curr_batch.resize(0);

      int j=curr_index;
      uint64_t lower=(i*block_size*P);
      uint64_t upper=((i+1)*block_size*P);

      while(j<temp_locations.size() && lower<=temp_locations[j] && temp_locations[j]<upper )
      {
        curr_batch.push_back(temp_locations[j]-lower);
        j++;
      }
      curr_index=j;
      create_new_gcs_block(curr_batch,bb_bitset_vec[i]);

      vec_num_keys.push_back(curr_batch.size());

    }



    return total_bits_used;

  }

  
  //initialize snarf 
  void snarf_init(vector<T> &keys,double bits_per_key,int num_ele_per_block)
  {

    bool testbool = (bits_per_key>3);
    assert(("Bits per Key are too low!", testbool));

    double target_fpr=pow(0.5,bits_per_key-3.0);
    //Set parameter values
    N=keys.size();
    P=pow(2,ceil(log2(1.00/target_fpr)));
    bit_size=ceil(log2(1.00/target_fpr));
    block_size=num_ele_per_block;
    total_blocks=ceil(N*1.00/block_size);
    bb_bitset_vec.resize(total_blocks);
    
    //build snarf model
    rmi=snarf_model<T>();
    rmi.snarf_model_builder(keys);

    N=keys.size();
    P=pow(2,ceil(log2(1.00/target_fpr)));
    bit_size=ceil(log2(1.00/target_fpr));
    block_size=num_ele_per_block;
    total_blocks=ceil(N*1.00/block_size);
    bb_bitset_vec.resize(total_blocks);

    //Get bit locations of set bits
    vector<uint64_t> temp_locations;
    get_locations(keys,temp_locations);

    //Build bit blocks using the set bit location values
    gcs_size=build_bb(temp_locations);


    return ;
  }


  //Inserts a value(var val) into bit block at certain index(var bb_index)
  //Current implementation is not performant.
  // It simply read the block to get a list of values and adds the new value to this list. Then created a new value for this list
  void  insert_in_block(uint64_t val,int bb_index)
  {
    
    int num_keys=vec_num_keys[bb_index];
    int delta_zero_count=0,delta_one_count=0;
    int offset_dense_itr=num_keys*bit_size,offset_sparse_itr=0;

    int bit_val;
    uint64_t temp;
    vector<uint64_t> val_list;
    for(int i=0;i<num_keys;i++)
    {
      bit_val=bb_bitset_vec[bb_index].bitset_read_bit(offset_dense_itr,1);
      
      if(bit_val==0)
      {
        delta_zero_count++;
        offset_dense_itr++;
        i--;
      }
      else
      {
        temp=delta_zero_count*P+bb_bitset_vec[bb_index].bitset_read_bits(offset_sparse_itr,bit_size);
        val_list.push_back(temp);
        offset_dense_itr++;
        offset_sparse_itr+=bit_size;
      }
    }

    val_list.push_back(val);
    sort(val_list.begin(),val_list.end());

    vec_num_keys[bb_index]++;
    create_new_gcs_block(val_list,bb_bitset_vec[bb_index]);
    

    return ;

  }

  //Deletes a value(var val) from a bit block at certain index(var bb_index)
  //Current implementation is not performant.
  // It simply read the block to get a list of values and removes the value from this list. Then created a new value for this list
  void  delete_from_block(uint64_t val,int bb_index)
  {
    int num_keys=vec_num_keys[bb_index];
    int delta_zero_count=0,delta_one_count=0;
    int offset_dense_itr=num_keys*bit_size,offset_sparse_itr=0;

    int bit_val;
    uint64_t temp;
    int done=0;
    vector<uint64_t> val_list;
    for(int i=0;i<num_keys;i++)
    {
      bit_val=bb_bitset_vec[bb_index].bitset_read_bit(offset_dense_itr,1);
      
      if(bit_val==0)
      {
        delta_zero_count++;
        offset_dense_itr++;
        i--;
      }
      else
      {
        
        temp=delta_zero_count*P+bb_bitset_vec[bb_index].bitset_read_bits(offset_sparse_itr,bit_size);
        
        offset_dense_itr++;
        if(done==0 && temp==val)
        {
          offset_sparse_itr+=bit_size;  
          done=1;
        }
        else
        {
          val_list.push_back(temp);
          offset_sparse_itr+=bit_size;  
        }
        
        
      }
    }
   
    
    sort(val_list.begin(),val_list.end());

    bool testbool = (done==1);
    assert(("The key to delete was not present!", testbool));


    vec_num_keys[bb_index]-=done;
    create_new_gcs_block(val_list,bb_bitset_vec[bb_index]);
    

    return ;

  }


  //checks if there is a value in a certain block(var bb_temp) that is between low_val and upper_val
  bool range_query_in_block(uint64_t low_val,uint64_t upper_val,snarf_bitset &bb_temp,int num_keys_read)
  {
    int num_keys=num_keys_read;
    int delta_zero_count=0,delta_one_count=0;
    int offset_dense_itr=num_keys*bit_size,offset_sparse_itr=0;

    int bit_val;
    uint64_t temp;
    for(int i=0;i<num_keys;i++)
    {
      bit_val=bb_temp.bitset_read_bit(offset_dense_itr,1);
      
      if(bit_val==1 && ((delta_zero_count+1)*P>=low_val) && ((upper_val)>=delta_zero_count*P))
      {
        //calculate the value
        temp=delta_zero_count*P+bb_temp.bitset_read_bits(offset_sparse_itr,bit_size);
        
        //value is between the range
        if(temp>=low_val && temp<=upper_val)
        {
          return true;
        }
        
      }
      delta_zero_count+=(1-bit_val);
      offset_dense_itr++;
      i-=(1-bit_val);
      offset_sparse_itr+=(bit_val*bit_size);

    }


    return false;


  }


  //finds the bit location corresponding to the key and inserts it in the corresponding block
  void insert_key(T key)
  {
    query_cdf1=rmi.infer(key);  

    uint64_t temp_loc_upper=floor(query_cdf1*N*P);
    temp_loc_upper=min(N*P-1,max((uint64_t)0,temp_loc_upper));

    delta_query_index=floor(temp_loc_upper*1.00/(block_size*P));
    delta_query_remainder=temp_loc_upper-delta_query_index*block_size*P;
   
    insert_in_block(delta_query_remainder,delta_query_index);
   
    return;

  }

  //finds the bit location corresponding to the key and deletes it from the corresponding block
  void delete_key(T key)
  {
    
    query_cdf1=rmi.infer(key); 
    
    uint64_t temp_loc_upper=floor(query_cdf1*N*P);
    temp_loc_upper=min(N*P-1,max((uint64_t)0,temp_loc_upper));

    delta_query_index=floor(temp_loc_upper*1.00/(block_size*P));
    delta_query_remainder=temp_loc_upper-delta_query_index*block_size*P;

    delete_from_block(delta_query_remainder,delta_query_index);

    return;

  }
 
  //finds the bit location corresponding to the query endpoints and checks the corresponding block or blocks for a value
  bool range_query(T lower_val,T upper_val)
  {
    uint64_t temp_loc_lower,temp_loc_upper,large_delta_query_index;

    query_cdf1=rmi.infer(upper_val);
    query_cdf2=rmi.infer(lower_val);

    temp_loc_upper=floor(query_cdf1*N*P);
    temp_loc_upper=min(N*P-1,max((uint64_t)0,temp_loc_upper));

     
    temp_loc_lower=floor(query_cdf2*N*P);
    temp_loc_lower=min(N*P-1,max((uint64_t)0,temp_loc_lower));

    delta_query_index=floor(temp_loc_lower*1.00/(block_size*P));
    large_delta_query_index=floor(temp_loc_upper*1.00/(block_size*P));

    //in case the query endpoints are  in two different blocks we need to query multiple times
    if(delta_query_index==large_delta_query_index)
    {
      return range_query_in_block(temp_loc_lower-delta_query_index*block_size*P,temp_loc_upper-delta_query_index*block_size*P,bb_bitset_vec[delta_query_index],vec_num_keys[delta_query_index]);
    }
    else
    {
      if(range_query_in_block(temp_loc_lower-delta_query_index*block_size*P-1,block_size*P+1,bb_bitset_vec[delta_query_index],vec_num_keys[delta_query_index]))
      {
        return true;
      }

      if(range_query_in_block(0,temp_loc_upper-large_delta_query_index*block_size*P,bb_bitset_vec[large_delta_query_index],vec_num_keys[large_delta_query_index]))
      {
        return true;
      } 

      for(int i=delta_query_index+1;i<large_delta_query_index;i++)
      {
        if(range_query_in_block(0,block_size*P-1,bb_bitset_vec[i],vec_num_keys[i]))
        {
          return true;
        }
      }

      return false;



    }


    return false;

  }

  //returns the space used by snarf overall
  int return_size()
  {
    int total_size=0;
    total_size+=sizeof(name_curr);
    total_size+=rmi.return_size();
    total_size+=7*sizeof(N);

    for(int i=0;i<vec_num_keys.size();i++)
    {
      total_size+=sizeof(vec_num_keys[i]);
    }

    for(int i=0;i<bb_bitset_vec.size();i++)
    {
      total_size+=bb_bitset_vec[i].return_size();
    }

    return total_size;
  }


};



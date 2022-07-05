#include<iostream>
#include<algorithm>
#include<cmath>
#include <random>
#include <fstream>
#include <chrono>
#include <cassert>
#include <set>
#include <ctime> // time_t
#include <boost/dynamic_bitset.hpp>
#include <cstring>
using namespace std;
using namespace std::chrono; 



// Bitset implementation to store the Golomb Coded values
// It relies on dynamic_bitset library in boost
struct snarf_bitset
{
  boost::dynamic_bitset<> bb_bitset;

  //initialize a dynamic bitset of particular size
  void init(uint64_t size)
  {
    bb_bitset.resize(size);
    bb_bitset.reset();
    return;
  }

  //writes certain amount of bits(var num_bits) from a value (var val) at an offset (var offset)
  void bitset_write_bits(uint64_t offset,uint64_t val, uint64_t num_bits)
  {
    uint32_t temp;
    while(num_bits>0)
    {
      temp=val%2;
      val=val/2;
      if(temp)
      {
        bb_bitset.set(offset);
      }
      else
      {
        bb_bitset.reset(offset);
      }
      num_bits--;
      offset++;
    }

    return ;
  }

  //returns certain amount of bits(var num_bits) at an offset (var offset)
  uint64_t bitset_read_bits(uint64_t offset,uint64_t num_bits)
  {
    uint64_t ans=0;
    uint32_t factor=1;
    while(num_bits>0)
    {
      ans+=bb_bitset.test(offset)*factor;
      factor*=2;
      num_bits--;
      offset++;
    }


    return ans;
  }

  // reads a single bit at an offset(var offset)
  uint64_t bitset_read_bit(uint64_t offset,uint64_t num_bits)
  {
    return bb_bitset.test(offset);
  }


  //returns space used by the structure in bytes
  int return_size()
  {
    double a=bb_bitset.size();
    int total_size=sizeof(a);

    total_size=ceil(a/8.00);

    return total_size; 
  }

  // writes the model contents into a char array
  void serialize(unsigned char* arr)
  {

    double a=bb_bitset.size();

    int offset=0;
    memcpy(arr+offset,&a,sizeof(a));
    offset+=sizeof(a);


    for(int i=0;i<ceil(a/8.00);i++)
    {
      unsigned char temp=0;
      for(int j=0;j<8 && (i*8+j<a);j++)
      {
        if(bb_bitset.test(i*8+j))
        {
          temp+=pow(2,j);
        }
      }

      memcpy(arr+offset,&temp,sizeof(temp));
      offset+=sizeof(temp);
    }

    return ;
  }

  // writes the model contents into a char array
  void deserialize(unsigned char* arr)
  {

    double bitset_size=10.09;
    int offset=0;

    memcpy(&bitset_size,arr+offset,sizeof(bitset_size));
    offset+=sizeof(bitset_size);

    bb_bitset.resize(bitset_size);
    bb_bitset.reset();

    for(int i=0;i<ceil(bitset_size/8.00);i++)
    {
      unsigned char temp;
      memcpy(&temp,arr+offset,sizeof(temp));
      offset+=sizeof(temp);

      for(int j=0;j<8 && i*8+j<bitset_size;j++)
      {
        if(temp%2==1)
        {
          bb_bitset.set(i*8+j);
        }
        else
        {
          bb_bitset.reset(i*8+j);
        }
        temp=temp/2;
      }
    }

    return ;
  }


};


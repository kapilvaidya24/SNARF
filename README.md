# SNARF(Sparse Numerical Array Based Range Filter)
Code repository for SNARF


**SNARF** is a fast and compact updatable range filter . This is the source code for our
[SIGMOD paper](https://dl.acm.org/doi/10.14778/3529337.3529347).


## Simple Example
A simple example can be found [here](https://github.com/kapilvaidya24/SNARF/blob/main/example.cpp). To run the example:
```
g++ -std=c++17 -O3 -w -fpermissive example.cpp -o example.out
./example.out
```

Please be carful with using the code for numbers close to integer represenation limit (>2^60). Integer overflow might occur which might result in inaccurate results. 

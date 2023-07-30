#include<iostream>
#include <cstring>
#include <string>
#include <set>
#include <math.h>
#include <fstream>
#include <stack>
#include <algorithm>
#include <string>
#include <chrono>
#include <ctime>
#include <unordered_map>
#include <vector>
#include <map>
#include <cmath>

#include <divsufsort64.h>    
#include <sdsl/bit_vectors.hpp>           
#include <sdsl/rmq_support.hpp>	
#include <sdsl/io.hpp>

using namespace std;
using namespace sdsl;

typedef int64_t INT;

int64_t LCParray ( unsigned char * text, INT n, INT * SA, INT * ISA, INT * LCP )
{
        INT i=0, j=0;

        LCP[0] = 0;
        for ( i = 0; i < n; i++ ) // compute LCP[ISA[i]]
                if ( ISA[i] != 0 )
                {
                        if ( i == 0) j = 0;
                        else j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]]-1 : 0;
                        while ( text[i+j] == text[SA[ISA[i]-1]+j] )
                                j++;
                        LCP[ISA[i]] = j;
                }
        return ( 1 );
}


int main(int argc, char **argv)
{
	
 	/* Read in list of sparse suffixes */
  	
	
	
	ifstream seq(argv[1], ios::in | ios::binary);
 	seq.seekg(0, ios::end);
   	uint64_t file_size = seq.tellg();
   	
   	uint64_t text_size = 0;
   	unsigned char * sequence2 =  ( unsigned char * ) calloc( file_size , sizeof( unsigned char ) );
   	
	char c = 0;
	seq.seekg (0, ios::beg);
	int64_t pos2 = 1;
		
	for (int64_t i = 0; i < file_size; i++)
	{
		seq.read(reinterpret_cast<char*>(&c), 1);

		/*if( ( char) c == '\n' || ( char) c == ' ' )
			continue;
		else
		{*/	
			sequence2[text_size] = c; 
			text_size++;
		//}
		
	}
	seq.close();
	
	vector<int64_t> * suffix_list = new vector<int64_t>();
	
	 ifstream suff_list(argv[2], ios::in | ios::binary);
  	suff_list.seekg(0, ios::end);
  	file_size = suff_list.tellg();
  	
  	string line="";
	suff_list.seekg (0, ios::beg);
	uint64_t b = 0;
	
	for (uint64_t i = 0; i < file_size; i++)
	{
		suff_list.read(reinterpret_cast<char*>(&c), 1);
		
		if( ( char) c == '\n' || ( char) c == ' ' )
		{
		
			if( (uint64_t)stoi(line) >= text_size )
			{
				cout<<"Suffix list contains a suffix larger than sequence length"<<endl;
				return 1;
			}	
			
			suffix_list->push_back(stoi(line));
			
			b++;
			line = "";
			
		}	
		else line += c;
	}
	suff_list.close();	
	
	int64_t n = text_size;
	
	int64_t * SA = ( int64_t * ) malloc( ( n ) * sizeof( int64_t ) );
	if( divsufsort64( sequence2, SA,  n ) != 0 )
  	{
  		fprintf(stderr, " Error: SA computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}
  	
  	int64_t * invSA = ( int64_t * ) calloc( n , sizeof( int64_t ) );
 	if( ( invSA == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
        	return ( 0 );
  	}
  	
  	for ( int64_t i = 0; i < n; i ++ )
  	{
  		invSA [SA[i]] = i;
  	}


  	int64_t * LCP = ( int64_t * ) calloc  ( n, sizeof( int64_t ) );
  	if( ( LCP == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
        	return ( 0 );
  	}
  	if( LCParray( sequence2, n, SA, invSA, LCP ) != 1 )
  	{
		fprintf(stderr, " Error: LCP computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}
  	
  	unordered_map<uint64_t,uint64_t> suffixes;
  	
  	vector<uint64_t> * real_SA = new vector<uint64_t>();
  	vector<uint64_t> * real_LCP = new vector<uint64_t>();
  	
  	
  	for(int i = 0; i<suffix_list->size(); i++)
  	{
 		suffixes[ suffix_list->at(i)] = 1;
 	}
 	
 		 	
  	for(int i = 0; i<text_size; i++)
  	{
  		if (suffixes.find( SA[i] ) != suffixes.end()) 
  		{
  		
  			real_SA->push_back( SA[i] );
  		}
  	
  	}

  	string output = argv[3];
  	string suff_lcp = output + "_accurate.lcp";
	string suff_ssa = output + "_accurate.ssa";
	
	ofstream output_ssa(suff_ssa);
	ofstream output_lcp(suff_lcp);


  	for(int i = 0; i<real_SA->size(); i++ )
  		output_ssa<<real_SA->at(i)<<endl;

  	
  	for(int i = 1; i<real_SA->size(); i++ )
  	{
  		int j = 0;
  		while( sequence2[real_SA->at(i)+j] == sequence2[real_SA->at(i-1)+j] )
  			j++;
  		
  		real_LCP->push_back(j);
  	}

	output_lcp<<0<<endl;
	
	for(int i = 0; i<real_LCP->size(); i++ )
  		output_lcp<<real_LCP->at(i)<<endl;

  	
	
	return 0;
}


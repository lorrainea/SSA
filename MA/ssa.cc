#include<iostream>
#include <cstring>
#include <string>
#include <set>
#include <math.h>
#include <fstream>
#include <stack>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <unordered_map>
#include <vector>
#include <map>
#include <cmath>
#include "krfp.h"
#include "unordered_dense.h"
#include <assert.h>
#include <sys/time.h>
#include <numeric>
#include <memory>
#include <sstream>
#define IPS4 false

#if IPS4 == true
#include "../include/ips4o.hpp"
#endif

#define DEBUG false
#define THRESHOLD 1500000

using namespace std;


struct SSA
{
	uint64_t lcp;
	vector<uint64_t> L;
	SSA() : lcp(0), L() {}
};

double gettime( void )
{
	struct timeval ttime;
	gettimeofday( &ttime , 0 );
	return ttime.tv_sec + ttime.tv_usec * 0.000001;
};
double prep_total;
double hash_total;
double kr_total;
double gr_total;
double order_total;
double sort_total;

auto compare(unsigned char * sequence, vector<uint64_t> * A, uint64_t lcp )
{
	return [sequence, A, lcp](uint64_t a, uint64_t b) 
	{
		return sequence[ A->at(a)+lcp ] < sequence[ A->at(b)+lcp ];
	};
}

/* Computing the fingerprint of an SSA */
uint64_t fingerprint( uint64_t ssa, uint64_t * FP, uint64_t fp_len, uint64_t l, unsigned char * sequence, uint64_t text_size )
{
	uint64_t fp = 0;
	uint64_t ssa_end = (text_size >= ssa+l) ? ssa+l : text_size; //this is the end of the substring we are interested in PLUS 1

	if( l > fp_len )
	{
		uint64_t fp_short = 0;
		if( ssa == 0 ) // then there is not short fragment before the substring
		{
			; // do nothing
		} 
		else if( ssa % fp_len  !=  0 ) // we are in-between two stored prefixes
		{
			uint64_t prefix = ssa / fp_len; // this is the prefix we are in but we have something stored on the left
			uint64_t start = 0;
			
			if( prefix != 0 )
			{
				fp_short = FP[prefix - 1];
				start = prefix * fp_len;
			}
			
			for(uint64_t i = start; i<ssa; i++)	fp_short = karp_rabin_hashing::concat( fp_short, sequence[i], 1);

		}
		else 	// we have the fp_short stored and we read it from FP
		{	
			uint64_t prefix = ssa / fp_len;
			fp_short = FP[prefix - 1];
		}	

		uint64_t fp_long = 0;

		if( ssa_end % fp_len  != 0 )
		{
			
                        uint64_t prefix = ssa_end / fp_len; // this is the prefix we are in but we have something stored on the left
                        uint64_t start = 0;

                        if( prefix != 0 )
                        {
                        	fp_long = FP[prefix - 1];
                        	start = prefix * fp_len;
                        }

                        for(uint64_t i = start; i< ssa_end; i++)	fp_long = karp_rabin_hashing::concat( fp_long, sequence[i] , 1 );
                    }
                else
                { 
                	uint64_t prefix = ssa_end / fp_len;
                	fp_long = FP[prefix - 1];
                }

                fp = karp_rabin_hashing::subtract(fp_long, fp_short, ssa_end - ssa);

            }
            else 
            {
            	for(uint64_t i=ssa; i< ssa_end; ++i)	fp =  karp_rabin_hashing::concat( fp, sequence[i], 1 );
            }

        return fp;
    }

uint64_t group( vector<SSA> * B, vector<uint64_t> * A, uint64_t * FP, uint64_t fp_len, uint64_t l, unsigned char * sequence, uint64_t text_size, uint64_t b, uint64_t &m, uint64_t &z )
{
    	vector<SSA> * B_prime = new vector<SSA>();
	(*B_prime).reserve(b); 

	const auto Bsz = B->size();
	
	for(uint64_t i = 0; i<Bsz; ++i )
	{
		const uint64_t s = (*B)[i].L.size();
		uint64_t k = 0;	
		vector<vector<uint64_t>> vec;
		vector<uint64_t> tmp;

		if( s <= (const uint64_t) z )
		{
			double start = gettime();
			vector<pair<uint64_t,uint64_t> > vec_to_sort;
			for(auto it=(*B)[i].L.begin();it!=(*B)[i].L.end(); ++it)
			{
				if ( (*A)[(*it)]+(*B)[i].lcp + l > text_size ) { tmp.push_back((*it)); continue; }
				uint64_t fp = fingerprint( (*A)[(*it)]+(*B)[i].lcp, FP, fp_len, l, sequence, text_size );
				vec_to_sort.push_back( make_pair(fp,*it) );
			}

			#if IPS4 == true
				ips4o::sort(vec_to_sort.begin(),vec_to_sort.end());
			#else
				sort(vec_to_sort.begin(),vec_to_sort.end());
			#endif

			const auto vsz=vec_to_sort.size();
			for(uint64_t i=0;i<vsz;++i)
			{
				if (i == 0 || (i>0 && vec_to_sort[i].first != vec_to_sort[i - 1].first)) //if we have new group 
				{
					vector<uint64_t> new_vec;		
					new_vec.push_back(vec_to_sort[i].second);
					vec.push_back(new_vec);	 //adds into vec a new vector
					k++;
				}
				else	vec.back().push_back(vec_to_sort[i].second); //adds the id into the last added vector in vec
			}
			double end = gettime();
			sort_total += end - start;
		}
		else
		{	
			double start = gettime();
    			auto groups = ankerl::unordered_dense::map<uint64_t, uint64_t >();
			for(auto it=(*B)[i].L.begin();it!=(*B)[i].L.end(); ++it)
			{
				if ( (*A)[(*it)]+(*B)[i].lcp + l > text_size ) { tmp.push_back((*it)); continue; }
				uint64_t fp = fingerprint( (*A)[(*it)]+(*B)[i].lcp, FP, fp_len, l, sequence, text_size );
				auto itx = groups.find(fp);
				if(itx == groups.end())
				{
					groups[fp] = k;
					vector<uint64_t> v;
					v.push_back(*it);
					vec.push_back(v);
					k++;
				}
				else	vec[itx->second].push_back(*it);
			}
			groups.clear();
			double end = gettime();
			hash_total += end - start;
		}

		double start = gettime();
		vector<uint64_t>().swap((*B)[i].L);

	        
		for( uint64_t j = 0; j < k; j++ )
		{	
			const auto itsz=vec[j].size();

			if( itsz == s )
			{
				(*B)[i].lcp += l;
				for(const auto& value: vec[j])	(*B)[i].L.push_back(value);				

			}
			else if( itsz >= 2 )
			{
				m++; 
				SSA new_ssa;

				for(const auto& value: vec[j])	new_ssa.L.push_back(value);				
				
				new_ssa.lcp = (*B)[i].lcp + l; 
				B_prime->push_back( new_ssa );
				(*B)[i].L.push_back(m);

				A->push_back( (*A)[ vec[j][0] ] );											
			}
			else if ( itsz == 1 )
			{
				(*B)[i].L.push_back( vec[j][0] );
			}
		}
		for (auto& it : tmp) (*B)[i].L.push_back( it );
                vector<uint64_t>().swap(tmp);
		vector<vector<uint64_t>>().swap(vec);
		double end = gettime();
		gr_total += end - start;
	}

	double start = gettime();
	B->insert(std::end(*B), std::begin(*B_prime), std::end(*B_prime));
	delete( B_prime);
	double end = gettime();
	gr_total += end - start;

	return 0;
}

uint64_t order( vector<uint64_t> * final_ssa, vector<uint64_t> * final_lcp, vector<SSA> * B, vector<uint64_t> * A, unsigned char * sequence, uint64_t text_size, uint64_t b )
{

	const uint64_t Bsz=B->size();
	for(uint64_t i = 0; i<Bsz; i++)
	{	
		#if IPS4 == true
			ips4o::sort((*B)[i].L.begin(), (*B)[i].L.end(), compare(sequence,A,(*B)[i].lcp));
		#else
			sort((*B)[i].L.begin(), (*B)[i].L.end(), compare(sequence,A,(*B)[i].lcp));
		#endif
	}
	stack<pair<uint64_t,uint64_t>> S; 
	
	S.push( make_pair<uint64_t, uint64_t>((uint64_t)b, 0) );  //b is the correct first index, not b+1
	
	uint64_t l = 0;
	uint64_t mymax=numeric_limits<uint64_t>::max();
	while( !S.empty() )
	{
		uint64_t i = S.top().first;
		uint64_t l_prime = S.top().second;
		S.pop();
		
		if( l_prime < l )
			l = l_prime;

		if(i>=b) //it is not one of the initial groups
		{
			uint64_t lcp = (*B)[i-b].lcp;
			auto myl=(*B)[i-b].L;
			for(vector<uint64_t>::reverse_iterator it=myl.rbegin();it!=myl.rend();++it)
			{				  					
				S.push(make_pair<uint64_t, uint64_t>( (uint64_t) *it, (uint64_t) lcp ));
			}
		}
		else
		{
			final_ssa->push_back( (*A)[i] );
			final_lcp->push_back( l );
			l = mymax;
		}	
	}

	return 0;
}


int main(int argc, char **argv)
{

	if( argc < 4 )
        {
                cout<<"Check arguments!\n";
                cout<<"./ssa <sequence_file> <suffix_list> <output_filename> [sort/hash_hreshold_z]\n";
                exit(-1);
        }

        uint64_t z;
        if ( argc > 4 )
        {
                std::string str4(argv[4]);
                std::stringstream(str4)>>z;
        }
        else z = THRESHOLD;
        cout<<"Threshold z="<<z<<endl;
	
	/* Read in sequence file */
	ifstream seq(argv[1], ios::in | ios::binary);

	//vector<char> input_seq_vec;
	unsigned char *input_seq_char=NULL;  
	uint64_t cnt=0;
	char c;
  	while (seq.get(c))     
  	{
  		if(cnt==0 || cnt%1048576)
  			input_seq_char = ( unsigned char * ) realloc ( input_seq_char,   ( cnt + 1048576 ) * sizeof ( unsigned char ) );		
		input_seq_char[ cnt ] = c;
		cnt++;
  	}	
  	seq.close();

  	unsigned char * sequence=input_seq_char;
	uint64_t text_size=cnt;	

	cout<<"Text length n = " << text_size << endl;

	/* Read in list of sparse suffixes */
	ifstream suff_list(argv[2], ios::in | ios::binary);
	suff_list.seekg(0, ios::end);
	uint64_t file_size = suff_list.tellg();

	vector<uint64_t> * ssa_list = new vector<uint64_t>();

	c = 0;
	string line="";
	suff_list.seekg (0, ios::beg);
	
	if( file_size == 0 )
	{
		cout<<"Empty suffixes list"<<endl;
		return 1;
	}
	
	for (uint64_t i = 0; i < file_size; i++)
	{
		suff_list.read(reinterpret_cast<char*>(&c), 1);
		
		if( ( char) c == '\n' || ( char) c == ' ' )
		{
			ssa_list->push_back(stoull(line));
			
			if( stoull(line) >= text_size )
			{
				cout<<"Suffix list contains a suffix larger than sequence length"<<endl;
				return 1;
			}

			line = "";
			
		}	
		else line += c;
	}
	suff_list.close();	
	
	karp_rabin_hashing::init();
	uint64_t b = ssa_list->size();
	cout<<"Number of suffixes b = " << b << endl;
	uint64_t l = 1ULL << static_cast<uint64_t>(log2(text_size));
	uint64_t s = 2*b;
	if ( s > text_size ) s = text_size;
        uint64_t fp_len = text_size / s;
        if ( (text_size - fp_len * s) > text_size/s )
                s = text_size/fp_len;
        cout<<"Block length = "<<fp_len<<endl;
        cout<<"Size s of FP table = " << s <<endl<<endl;
	
	// computing fingerprints
	uint64_t * FP =  ( uint64_t * ) calloc( s , sizeof( uint64_t ) );
	uint64_t fp = 0;
	
	std::chrono::steady_clock::time_point start_total = std::chrono::steady_clock::now();
	
	prep_total = 0;
	double start = gettime();
	cout<<"Preprocessing starts"<<endl;
	uint64_t i = 0;
	for(uint64_t j = 0; j < s; ++j )
	{
		for (uint64_t k = 0; k < fp_len; k ++ )	fp = karp_rabin_hashing::concat(fp, sequence[i+k], 1);
			FP[j]=fp;
		i += fp_len;
		if ( i + fp_len > text_size ) break;
	}
	cout<<"Preprocessing ends"<<endl<<endl;
	double end = gettime();
	prep_total = end - start;

	vector<SSA> * B = new vector<SSA>();
	(*B).reserve(b); 

	vector<uint64_t> * A = new vector<uint64_t>();
	vector<uint64_t> L;

	const auto ssa_list_sz=ssa_list->size();
	for(uint64_t i = 0; i<ssa_list_sz; ++i )
	{	
		A->push_back( (*ssa_list)[i] );
		L.push_back(i);
	}
	
	SSA initial;
	initial.lcp = 0;
	initial.L = L;
	
	B->push_back( initial );
	A->push_back( (*ssa_list)[0] );
	uint64_t m = ssa_list->size();
	
	delete( ssa_list );
	
	hash_total = 0;
	kr_total = 0;
	gr_total = 0;
	cout<<"Grouping starts"<<endl;
	while( l > 0 )
	{
		cout<< "l: " << l <<", nodes: "<< m <<endl;
		group( B, A, FP, fp_len, l, sequence, text_size, b, m, z );
		l=l>>1;
	}
	cout<<"Grouping ends"<<endl;
	
	free( FP );
	
	vector<uint64_t> * final_ssa = new vector<uint64_t>();
	vector<uint64_t> * final_lcp = new vector<uint64_t>();
	
	cout<<"Sorting starts"<<endl;

	order_total = 0;
	start = gettime();
	order( final_ssa, final_lcp, B, A, sequence, text_size, b);	
	end = gettime();
	order_total = end - start;
	
	cout<<"Sorting ends"<<endl;
	
	std::chrono::steady_clock::time_point end_total = std::chrono::steady_clock::now();
	std::cout<<"Time taken "<<std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total).count()/(double)1000 << "[s]" << std::endl;
	std::cout<<" Time taken by preprocessing "<<prep_total << "[s]" << std::endl;
	std::cout<<" Time taken by sorting fingerpints "<<sort_total << "[s]" << std::endl;
	std::cout<<" Time taken by hashing fingeprints "<<hash_total << "[s]" << std::endl;
	std::cout<<" Time taken by the rest of grouping "<<gr_total << "[s]" << std::endl;
	std::cout<<" Time taken by ordering "<<order_total << "[s]" << std::endl;

	std::cout<<"Writing the output"<< std::endl;
	string output = argv[3];
	
	string suff_lcp = output + ".lcp";
	string suff_ssa = output + ".ssa";
	
	ofstream output_ssa(suff_ssa);
	
	for(uint64_t i = 0; i<final_ssa->size(); i++)
	{
		output_ssa<<final_ssa->at(i)<<endl;
		if(DEBUG)	cout<<(*final_ssa)[i]<<" ";
	}
	ofstream output_lcp(suff_lcp);
	for(uint64_t i = 0; i<final_lcp->size(); i++)
	{
		output_lcp<<final_lcp->at(i)<<endl;
		if(DEBUG)	cout<<(*final_lcp)[i]<<" ";
	}
	cout<<endl;
	delete( final_lcp );
	delete( final_ssa );
	delete( B );
	delete( A );
	free ( input_seq_char);
	input_seq_char = NULL;

	return 0;
}

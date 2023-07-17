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

#define DEBUG false

using namespace std;

struct SSA
{
	
	uint64_t m;
	uint64_t lcp;
	vector<uint64_t> L;
};


auto compare(unsigned char * sequence, vector<uint64_t> * A, uint64_t lcp )
{
	return [sequence, A, lcp](uint64_t a, uint64_t b) 
	{
		
		
		return sequence[ A->at(a)+lcp ] <= sequence[ A->at(b)+lcp ];
	};
}

	
/* Computing the fingerprint of an SSA */
uint64_t fingerprint( uint64_t ssa, uint64_t * FP, uint64_t fp_len, uint64_t l, unsigned char  * sequence, uint64_t text_size, uint64_t b )
{

	uint64_t fp = 0;

	if( l <= text_size / b )
	{ 
	    const auto msz=(text_size >= ssa+l) ? ssa+l : text_size; 

		//for(uint64_t i = ssa; i<min( text_size, ssa+l ); i++)
		for(uint64_t i=ssa; i< msz; ++i)
		{
			fp =  karp_rabin_hashing::concat( fp, sequence[i], 1 );
		}
	
	}
	else
	{
		//ssa = min( ssa, text_size );
		
		uint64_t closest_start = ssa/fp_len;
		uint64_t end =  text_size >= ssa+ l ? ssa+l : text_size;
		uint64_t closest_end = end/fp_len;	
		uint64_t fp_short = 0;
		uint64_t fp_long = 0;
		
		//cout<<" closest start "<<closest_start<<endl;
		
		if( ssa == 0 )
		{
			fp_short = 0;
		}
		else if( ssa % fp_len  !=  0 )
		{
			const uint64_t diff = ssa - ( fp_len * closest_start ) ;
			uint64_t c_start = 0;
			
			if( closest_start != 0 )
			{
				fp_short = FP[closest_start-1];
				c_start = closest_start;
			}
			else
			{ 
				fp_short = 0;
				c_start = 0;
			}
			
			
			const auto x=c_start*fp_len; 
			for(uint64_t i = 0; i<diff; i++)
			{
				fp_short = karp_rabin_hashing::concat( fp_short, sequence[x + i], 1);
			}
			
		}
		else fp_short = FP[closest_start-1]; 

		if( end %  fp_len  != 0 )
		{
			uint64_t diff = end - ( (fp_len * closest_end  ) );
			uint64_t c_end = 0;
			
			
			if( closest_end == 0 )
			{
				fp_long = 0;
				c_end = 0;
			}
			else 
			{
				fp_long = FP[closest_end-1];
				c_end = closest_end ;
			}
			const auto x=c_end*fp_len; 

			for(uint64_t i = 0; i<diff; i++)
			{
				fp_long = karp_rabin_hashing::concat( fp_long, sequence[x + i] , 1 );
			}
		}
		else fp_long = FP[closest_end-1];
		
		
		fp = karp_rabin_hashing::subtract(fp_long, fp_short, end-ssa);

	}
	
	if(DEBUG)
	{	for(uint64_t x=ssa; x<ssa+l;++x)
			cout<<sequence[x];
		cout<<"| "<<fp<<endl;
	}
return fp;
}

/* Grouping SSAs according to shared fingerprints */
uint64_t group( vector<SSA> * B, vector<uint64_t> * A, uint64_t * FP, uint64_t fp_len, uint64_t l, unsigned char * sequence, uint64_t text_size, uint64_t b, uint64_t &m )
{
	//unordered_map<uint64_t, vector<uint64_t>> groups; //= new unordered_map<uint64_t, vector<uint64_t>>();
	auto groups = ankerl::unordered_dense::segmented_map<uint64_t, vector<uint64_t> >();

	vector<SSA> * B_prime = new vector<SSA>();
	
	const auto Bsz=B->size();
	
	for(uint64_t i = 0; i<Bsz; ++i)
	{
	
		const uint64_t s = (*B)[i].L.size();
		
	    for(auto it=(*B)[i].L.begin();it!=(*B)[i].L.end();)
	    {
			
			uint64_t fp = fingerprint( (*A)[(*it)]+(*B)[i].lcp, FP, fp_len, l, sequence, text_size, b );
			auto itx=groups.find(fp);
			if(itx==groups.end())
			{
				vector<uint64_t> vec;
				vec.push_back(*it);
					groups[fp]=vec;
			}
			else
			{
				groups[fp].push_back(*it);
			}
			it=((*B)[i].L).erase(it);
		}
		
		for( auto it = groups.begin(); it!= groups.end(); ++it)
		{	
			const auto itsz=it->second.size();

			if( itsz >= 2 )
			{
				m++; 
			 	
				SSA new_ssa;
				
				new_ssa.m = m;
			
		
		   		for(auto &it2 : it->second)
					new_ssa.L.push_back(it2);				
				
				new_ssa.lcp =  (*B)[i].lcp + l; //greg
				
				B_prime->push_back( new_ssa );
				
				(*B)[i].L.push_back(m);
				
				A->push_back( (*A)[ it->second[0] ] );											
			}
			else if( itsz == s )
			{
				(*B)[i].lcp += l; //greg
				
				for(auto &it2 : (it->second))
					(*B)[i].L.push_back(it2);
			 				
			}
			else if ( itsz == 1 )
			{
				(*B)[i].L.push_back( (it->second)[0] );
			}
		}
	
		
		groups.clear();
	}
		const uint64_t Bpsz=B_prime->size();
		for(uint64_t i = 0; i<Bpsz; i++)	
			B->push_back( (*B_prime)[i] );
		
	
	delete( B_prime);	
	
return 0;
}


uint64_t order( vector<uint64_t> * final_ssa, vector<uint64_t> * final_lcp, vector<SSA> * B, vector<uint64_t> * A, unsigned char * sequence, uint64_t text_size, uint64_t b )
{

	const uint64_t Bsz=B->size();
	for(uint64_t i = 0; i<Bsz; i++)
	{
			
		sort((*B)[i].L.begin(), (*B)[i].L.end(), compare(sequence,A,(*B)[i].lcp));
	}
	if(DEBUG)
	{
		cout<<"\n print B:\n";
		for(auto it : *B)
		{
			 SSA x=(SSA)(it);
		 	cout<<"i:";
		 	for(auto &it2 : x.L)
		 	{
				 cout<<it2<<" ";
		 	}
		 	cout<<endl;
		}
	}
	
	stack<pair<uint64_t,uint64_t>> S; //= new stack<pair<uint64_t,uint64_t>>();
	
	S.push( make_pair<uint64_t, uint64_t>((uint64_t)b, 0) );  //b is the correct first index, not b+1
	if(DEBUG)
		cout<<"first_added: "<<b<<" "<<0<<endl;
	
	uint64_t l = 0;
	uint64_t mymax=numeric_limits<uint64_t>::max();
	while( !S.empty() )
	{
		uint64_t i = S.top().first;
		uint64_t l_prime = S.top().second;
		S.pop();
		
		if( l_prime < l )
			l = l_prime;
				
		if(DEBUG)
			cout<<"("<<i<<" "<<b<<")"<<endl;      //this is what I got from top
		
		if(i>=b) //it is not one of the initial groups
		{
			uint64_t lcp = (*B)[i-b].lcp;
			auto myl=(*B)[i-b].L;
			for(vector<uint64_t>::reverse_iterator it=myl.rbegin();it!=myl.rend();++it)
			{				  					
						S.push(make_pair<uint64_t, uint64_t>( (uint64_t) *it, (uint64_t) lcp ));
						if(DEBUG) cout<<"Stack add (i,lcp): "<<*it<<" "<<lcp<<endl;		
			}
		}
		else
		{
			final_ssa->push_back( (*A)[i] );
			final_lcp->push_back( l );
			l = mymax;
			if(DEBUG)
			{
				cout<<"SSA add:"<<A->at(i)<<endl;
				cout<<"SLCP add:"<<l<<endl;
			}
		}	
	}
		
return 0;
}


int main(int argc, char **argv)
{

	

	if( argc < 4 )
 	{
        	cout<<"Check arguments!\n";
 		cout<<"./ssa <sequence_file> <suffix_list> <output_filename>\n";
 		exit(-1);
 	}

	/* Read in sequence file */
 	ifstream seq(argv[1], ios::in | ios::binary);
 	
	
 	seq.seekg(0, ios::end);
   	uint64_t file_size = seq.tellg();
   	
   	uint64_t text_size = 0;
   	unsigned char * sequence =  ( unsigned char * ) calloc( file_size , sizeof( unsigned char ) );
   	
	char c = 0;
	seq.seekg (0, ios::beg);
	uint64_t pos = 1;
	
	
	if( file_size == 0 )
	{
		cout<<"Empty sequence file"<<endl;
		return 1;
	}
		
	for (uint64_t i = 0; i < file_size; i++)
	{
		seq.read(reinterpret_cast<char*>(&c), 1);

		if( ( char) c == '\n' || ( char) c == ' ' )
			continue;
		else
		{	
			sequence[text_size] = static_cast<int>(c); 
			text_size++;
		}
		
	}
	seq.close();
	//return 0;
	//cout<<text_size<<endl;
 	/* Read in list of sparse suffixes */
  	ifstream suff_list(argv[2], ios::in | ios::binary);
  	suff_list.seekg(0, ios::end);
  	file_size = suff_list.tellg();
   
  	vector<uint64_t> * ssa_list = new vector<uint64_t>();
 	vector<uint64_t> * slcp_list = new vector<uint64_t>();
 	
 	
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
			ssa_list->push_back(stoi(line));
			
			if( stoi(line) >= text_size )
			{
				cout<<"Suffix list contains a suffix larger than sequence length"<<endl;
				return 1;
			}	
			slcp_list->push_back(0);
			
			line = "";
			
		}	
		else line += c;
		
		
	}
	suff_list.close();	
	
	karp_rabin_hashing::init();
	uint64_t b = ssa_list->size();
	//uint64_t l = pow(2, (uint64_t) log2(text_size));
	uint64_t l = 1ULL << static_cast<uint64_t>(log2(text_size));

	//uint64_t fp_blocks = (text_size/(text_size/b)) ;
	uint64_t fp_len = text_size / b;
	
	// computing fingerprints
	uint64_t * FP =  ( uint64_t * ) calloc( b , sizeof( uint64_t ) );
	uint64_t fp = 0;
	pos = 0;
	
	std::chrono::steady_clock::time_point start_total = std::chrono::steady_clock::now();
	
	uint64_t cnt=0;
    for(uint64_t i=0;i<text_size;++i)
    {
		fp = karp_rabin_hashing::concat(fp, sequence[i], 1);
		cnt++;
		if(cnt%fp_len==0) //start new block
		{
			FP[pos]=fp;
			pos++;
		}
		
	}
	if(DEBUG)
	{		for(int i=0;i<pos;++i)
		{
			cout<<i<<"->"<<FP[i]<<endl;
		}
	}
	vector<SSA> * B = new vector<SSA>();
	vector<uint64_t> * A = new vector<uint64_t>();
	vector<uint64_t> L;

	const auto ssa_list_sz=ssa_list->size();
	for(uint64_t i = 0; i<ssa_list_sz; ++i )
	{	
		A->push_back( (*ssa_list)[i] );
		L.push_back(i);
	}
	
	SSA initial;
	initial.m = ssa_list->size();
	initial.lcp = 0;
	initial.L = L;
	
	B->push_back( initial );
	A->push_back( (*ssa_list)[0] );
	uint64_t m = ssa_list->size();
	
	delete( ssa_list );
	delete( slcp_list );
	

	while( l > 0 )
	{
		group( B, A, FP, fp_len, l, sequence, text_size, b, m );
		l=l>>1;
	}
	
	free( FP );
	
	vector<uint64_t> * final_ssa = new vector<uint64_t>();
	vector<uint64_t> * final_lcp = new vector<uint64_t>();
	
	//cout<<" ORDERING "<<endl;
	
	//for(uint64_t i = 0; i < B->size(); i++ )
	//{
		//cout<<"i "<<i<<endl;
		//cout<<B->at(i).lcp<<endl;
		
	//	for(uint64_t j = 0; j<B->at(i).L.size(); j++)
	//		cout<<B->at(i).L.at(j)<<endl;
		
//	}
	order( final_ssa, final_lcp, B, A, sequence, text_size, b);	
	
	std::chrono::steady_clock::time_point end_total = std::chrono::steady_clock::now();
	std::cout<<"Time taken "<<std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total).count() << "[ms]" << std::endl;
	

	string output = argv[3];
	
	string suff_lcp = output + ".lcp";
	string suff_ssa = output + ".ssa";
	
	ofstream output_ssa(suff_ssa);
	
	cout<<"SSA:\n";

	for(uint64_t i = 0; i<final_ssa->size(); i++)
	{
		output_ssa<<final_ssa->at(i)<<endl;
		cout<<(*final_ssa)[i]<<" ";
	}
	cout<<"\nSLCP:\n";
	ofstream output_lcp(suff_lcp);
	for(uint64_t i = 0; i<final_lcp->size(); i++)
	{
		output_lcp<<final_lcp->at(i)<<endl;
		cout<<(*final_lcp)[i]<<" ";
	}
	cout<<endl;
	delete( final_lcp );
	delete( final_ssa );
	delete( B );
	delete( A );
	
	
	
	free( sequence );
	
	return 0;
}


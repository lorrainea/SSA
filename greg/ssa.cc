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
	//cout<<"l:"<<l<<" text_size/b:"<<text_size/b<<endl;
	if( l <= text_size / b )
	{
	
		for(uint64_t i = ssa; i<min( text_size, ssa+l ); i++)
		{
			fp =  karp_rabin_hashing::concat( fp, sequence[i], 1 );
		}
	
	}
	else
	{
		//ssa = min( ssa, text_size );
		
		uint64_t closest_start = ssa/fp_len;
		uint64_t end = min( text_size , ssa+ l );
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
			uint64_t diff = ssa - ( fp_len * closest_start ) ;
			uint64_t c_start = 0;
			
			if( closest_start == 0 )
			{ 
				fp_short = 0;
				c_start = 0;
			}
			else
			{
				fp_short = FP[closest_start-1];
				c_start = closest_start;
			}
			
			
			for(uint64_t i = 0; i<diff; i++)
			{
				fp_short = karp_rabin_hashing::concat( fp_short, sequence[c_start * fp_len + i], 1);
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
					
			for(uint64_t i = 0; i<diff; i++)
			{
				fp_long = karp_rabin_hashing::concat( fp_long, sequence[c_end * fp_len + i] , 1 );
			}
		}
		else fp_long = FP[closest_end-1];
		
		
		fp = karp_rabin_hashing::subtract(fp_long, fp_short, end-ssa);

	}
	
	for(uint64_t x=ssa; x<ssa+l;++x)
		cout<<sequence[x];
	cout<<"| "<<fp<<endl;
return fp;
}

/* Grouping SSAs according to shared fingerprints */
uint64_t group( vector<SSA> * B, vector<uint64_t> * A, uint64_t * FP, uint64_t fp_len, uint64_t l, unsigned char * sequence, uint64_t text_size, uint64_t b, uint64_t &m )
{
	unordered_map<uint64_t, vector<uint64_t>> groups; //= new unordered_map<uint64_t, vector<uint64_t>>();
	
	vector<SSA> * B_prime = new vector<SSA>();
	
	for(uint64_t i = 0; i<B->size(); i++)
	{
	
		uint64_t s = B->at(i).L.size();
		
		/*for(uint64_t j = 0; j<B->at(i).L.size(); j++)		
		{
		    
			uint64_t fp = fingerprint( A->at(B->at(i).L.at(j))+B->at(i).lcp, FP, fp_len, l, sequence, text_size, b );
						
			auto itx=groups.find(fp);
			if(itx==groups.end())
			{
				vector<uint64_t> vec;
				vec.push_back(B->at(i).L.at(j));
					groups[fp]=vec;
			}
			else
			{
				groups[fp].push_back(B->at(i).L.at(j));
			}			
		}
		*/
		
		//B->at(i).L.clear();  //greg 
	// ************************************** Delete one by one above instead **************************************** 
		
	    for(auto it=(*B)[i].L.begin();it!=(*B)[i].L.end();)
	    {
			
			uint64_t fp = fingerprint( A->at(*it)+(*B)[i].lcp, FP, fp_len, l, sequence, text_size, b );
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
		
		for( auto it = groups.begin(); it!= groups.end(); it++)
		{	
			
			if( it->second.size() >= 2 )
			{
				m = m + 1;
			 	
				SSA new_ssa;
				
				new_ssa.m = m;
			
		
				for(uint64_t j = 0; j<it->second.size(); j++)
				{				
					//cout<<" new node "<<m<<" "<<it->second.at(j)<<" l "<<l<<endl;
					new_ssa.L.push_back( it->second.at(j) );
				}				
				
				new_ssa.lcp =  B->at(i).lcp + l; //greg
				
				B_prime->push_back( new_ssa );
				
				B->at(i).L.push_back(m);
				
				A->push_back( A->at( it->second.at(0) ) );											
			}
			else if( it->second.size() == s )
			{
				B->at(i).lcp += l; //greg
				
				for(uint64_t j = 0; j<it->second.size(); j++)
			 		B->at(i).L.push_back( it->second.at(j) );			
			}
			else if ( it->second.size() == 1 )
			{
				B->at(i).L.push_back( it->second.at(0) );
			}
		}
	
		
		groups.clear();
	}
	
	//cout<<"B_prime:"<<endl;
	for(uint64_t i = 0; i<B_prime->size(); i++)
	{	
		B->push_back( B_prime->at(i) );
		//cout<<"B'[i] ";
		/*for(auto &it2 : B_prime->at(i).L)
			cout<<it2<<" ";
		cout<<endl;*/
	}
	
	
	delete( B_prime);
	//delete( groups );
	

	
return 0;
}


uint64_t order( vector<uint64_t> * final_ssa, vector<uint64_t> * final_lcp, vector<SSA> * B, vector<uint64_t> * A, unsigned char * sequence, uint64_t text_size, uint64_t b )
{

	
	for(uint64_t i = 0; i<B->size(); i++)
	{
			
		sort( B->at(i).L.begin(), B->at(i).L.end(), compare(sequence,A,B->at(i).lcp));	
	
	}
	
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
	
	
	stack<pair<uint64_t,uint64_t>> S; //= new stack<pair<uint64_t,uint64_t>>();
	
	S.push( make_pair<uint64_t, uint64_t>((uint64_t)b, 0) );  //b is the correct first index, not b+1
	cout<<"first_added: "<<b<<" "<<0<<endl;
	
	uint64_t l = 0;
	
	while( !S.empty() )
	{
		uint64_t i = S.top().first;
		uint64_t l_prime = S.top().second;
		S.pop();
		
		if( l_prime< l )
			l = l_prime;
				
		cout<<"("<<i<<" "<<b<<")"<<endl;      //this is what I got from top
		if(i>=b) //it is not one of the initial groups
		{
			uint64_t lcp = B->at(i-b).lcp;
			auto myl=B->at(i-b).L;
			for(vector<uint64_t>::reverse_iterator it=myl.rbegin();it!=myl.rend();++it)
			{				  					
						S.push(make_pair<uint64_t, uint64_t>( (uint64_t) *it, (uint64_t) lcp ));
						cout<<"Stack add (i,lcp): "<<*it<<" "<<lcp<<endl;		
			}
		}
		else
		{
			final_ssa->push_back( A->at(i) );
			final_lcp->push_back( l );
			cout<<"SSA add:"<<A->at(i)<<endl;
			cout<<"SLCP add:"<<l<<endl;
			l = numeric_limits<uint64_t>::max();
		}	
		/*cout<<"Sz: "<<S.size()<<" i:"<<i<<" b:"<<b<<" i-b "<<i-b<<endl;	
		
		uint64_t lcp = B->at(i-b-1).lcp;
		
		if(B->at(i-b).L.empty() )
		{
			final_ssa->push_back( A->at(i) );
			final_lcp->push_back( l );
			cout<<"ifSSA add:"<<A->at(i)<<endl;
			cout<<"ifSLCP add:"<<l<<endl;
			l = numeric_limits<uint64_t>::max();
		}
		else
		{
			cout<<"B->at(i-b).L.size():"<<B->at(i-b).L.size();
			/*for(uint64_t j = B->at(i-b).L.size()-1; j>=0; j--)
			{
				;//S.push( make_pair<uint64_t, uint64_t>( (uint64_t) j, (uint64_t) lcp ) );
			}*/
		/*	auto myl=B->at(i-b-1).L;
			for(vector<uint64_t>::reverse_iterator it=myl.rbegin();it!=myl.rend();++it)
			{				  					
						S.push(make_pair<uint64_t, uint64_t>( (uint64_t) *it, (uint64_t) lcp ));
						cout<<"added: "<<*it<<" "<<lcp<<endl;									
			}
		
		}*/
		//cout<<"Sz: "<<S.size()<<endl;
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
	uint64_t l = pow(2, (uint64_t) log2(text_size));
	uint64_t fp_blocks = (text_size/(text_size/b)) ;
	uint64_t fp_len = text_size / b;
	
	// computing fingerprints
	uint64_t * FP =  ( uint64_t * ) calloc( fp_blocks , sizeof( uint64_t ) );
	uint64_t fp = 0;
	pos = 0;
	
	std::chrono::steady_clock::time_point start_total = std::chrono::steady_clock::now();
	
	/*for(uint64_t i = 0; i<fp_blocks *(text_size/b) ; i++)
	{
		fp = karp_rabin_hashing::concat(fp, sequence[i], 1);
		if( i > 0 &&  ( i + 1 ) % fp_len == 0 )
		{
			FP[pos] = fp ;
			pos++;
		}
		
		
	}
	*/
	cout<<"!"<<fp_len<<endl;
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
	for(int i=0;i<pos;++i)
	{
		cout<<i<<"->"<<FP[i]<<endl;
	}
	
	vector<SSA> * B = new vector<SSA>();
	vector<uint64_t> * A = new vector<uint64_t>();
	vector<uint64_t> L;

	
	for(uint64_t i = 0; i<ssa_list->size(); i++ )
	{	
		A->push_back( ssa_list->at(i) );
		L.push_back(i);
		
	}
	
	SSA initial;
	initial.m = ssa_list->size();
	initial.lcp = 0;
	initial.L = L;
	
	B->push_back( initial );
	A->push_back( ssa_list->at(0) );
	uint64_t m = ssa_list->size();
	
	delete( ssa_list );
	delete( slcp_list );
	

	while( l > 0 )
	{
		//cout<<" l "<<l<<endl;
		group( B, A, FP, fp_len, l, sequence, text_size, b, m );
		l = l/2;
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
	
	string output = argv[3];
	
	string suff_lcp = output + ".lcp";
	string suff_ssa = output + ".ssa";
	
	ofstream output_ssa(suff_ssa);
	
	for(uint64_t i = 0; i<final_ssa->size(); i++)
	{
		output_ssa<<final_ssa->at(i)<<endl;
	}
	
	ofstream output_lcp(suff_lcp);
	for(uint64_t i = 0; i<final_lcp->size(); i++)
	{
		output_lcp<<final_lcp->at(i)<<endl;
	}
	
	delete( final_lcp );
	delete( final_ssa );
	delete( B );
	delete( A );
	
	std::chrono::steady_clock::time_point end_total = std::chrono::steady_clock::now();
	std::cout<<"Time taken "<<std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total).count() << "[ms]" << std::endl;
	
	
	free( sequence );
	
	return 0;
}


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
#include "unordered_dense.h"
#include "krfp.h"

using namespace std;

struct SSA
{
	
	uint64_t ssa;
	uint64_t lcp;
	uint64_t leaf;
	uint64_t leaf_value;
};

auto compare(unsigned char * sequence )
{
	return [sequence](SSA a, SSA b) 
	{
		return sequence[ a.leaf_value+a.lcp ] <= sequence[ b.leaf_value+b.lcp ];
	};
}

/* Computing the fingerprint of an SSA */
uint64_t fingerprint( uint64_t ssa, uint64_t * FP, uint64_t fp_len, uint64_t l, unsigned char  * sequence, uint64_t text_size, uint64_t b )
{

	uint64_t fp = 0;
	
	if( l <= text_size / b )
	{
	
		for(uint64_t i = ssa; i<min( text_size, ssa+l ); i++)
		{
			fp =  karp_rabin_hashing::concat( fp, sequence[i], 1 );
		}
	
	}
	else
	{
		ssa = min( ssa, text_size );
		
		uint64_t closest_start = ssa/fp_len;
		uint64_t end = min( text_size , ssa+ l );
		uint64_t closest_end = end/fp_len;	
		uint64_t fp_short = 0;
		uint64_t fp_long = 0;
		
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
	
return fp;
}

/* Grouping SSAs according to shared fingerprints */
uint64_t group( vector<vector<SSA>> * ssa_struct, uint64_t * FP, uint64_t fp_len, uint64_t l, unsigned char * sequence, uint64_t text_size, uint64_t b )
{
	auto groups = ankerl::unordered_dense::map<uint64_t, vector<SSA>>();
	
	uint64_t size = ssa_struct->size();
	
	for(uint64_t i = 0; i<size; i++)
	{
		uint64_t original_size = ssa_struct->at(i).size();
	
		for(uint64_t j = 0; j<original_size; j++)
		{
			uint64_t ssa = 0;
			uint64_t lcp = 0;
			
			// Check if the node is a leaf or not
			if( ssa_struct->at(i).at(j).leaf == 1 )
			{
				ssa = ssa_struct->at(i).at(j).ssa;
				lcp = ssa_struct->at(i).at(j).lcp;
			}
			else
			{
				SSA ssa_pos =  ssa_struct->at( ssa_struct->at(i).at(j).ssa ).at(0);
				
				while( ssa_pos.leaf == 0 )
				{
					ssa_pos = ssa_struct->at( ssa_pos.ssa ).at(0);	
				}
				
				ssa = ssa_pos.ssa;
				lcp = ssa_struct->at(i).at(0).lcp;
			}
			
		
			// Find the fingerprint of the suffix and add to hash map
			uint64_t fp = fingerprint( ssa+lcp, FP, fp_len, l, sequence, text_size, b );
			
			
			if( groups.count( fp ) == 0 )
				groups[ fp ] =  {ssa_struct->at(i).at(j)};
			else groups[ fp ].push_back( ssa_struct->at(i).at(j) );
			
		}
		
		// Only create new node if some nodes were grouped together
		if( groups.size() != original_size )
		{
			ssa_struct->at(i).clear();
			
			for (auto const& [key, val] : groups)
			{
				if( val.size() > 1 )
				{
					// these nodes have been grouped together, place them in a new node and remove them from current node and replace
					// with new node ID.
						
					uint64_t original_lcp = val.at(0).lcp;
					
					vector<SSA> new_inp;
					for( int64_t k=0; k<val.size(); k++)
					{
						SSA new_ssa = {};
						new_ssa.leaf = val.at(k).leaf;
						new_ssa.ssa = val.at(k).ssa;
						new_ssa.lcp = val.at(k).lcp + l;
						
						new_inp.push_back( new_ssa );
							
					}
					ssa_struct->push_back( new_inp );
					
					
					//create new internal node and add to original ssa vector
					SSA new_internal = {};
					new_internal.leaf = 0;
					//new_internal.leaf_value = it->second.at(0).ssa;
					new_internal.ssa = ssa_struct->size()-1;
					new_internal.lcp = original_lcp;
						
					ssa_struct->at(i).push_back( new_internal );
					// now put new node value back into original.
				}
				else
				{
					SSA original = {};
					original.leaf =  val.at(0).leaf;
					original.ssa =  val.at(0).ssa;
					original.lcp =  val.at(0).lcp;
						
					ssa_struct->at(i).push_back( original );
				
			
				}
			}
			
		}
		groups.clear();
	}
	
	

return 0;
}


uint64_t order( vector<vector<SSA>> * ssa_struct, unsigned char * sequence, uint64_t text_size, string output )
{
	vector<SSA> * to_order = new vector<SSA>();
	
	string suff_lcp = output + ".lcp";
	string suff_ssa = output + ".ssa";
	
	ofstream output_ssa(suff_ssa);
	ofstream output_lcp(suff_lcp);
	
	for(uint64_t i = 0; i<ssa_struct->size(); i++)
	{
		// Depth first search to identify leaf/ssa values for each internal node
		for(uint64_t j = 0; j<ssa_struct->at(i).size(); j++)
		{
			if( ssa_struct->at(i).at(j).leaf == 1 )
			{
				ssa_struct->at(i).at(j).leaf_value = ssa_struct->at(i).at(j).ssa;	
			}
			else
			{
				SSA ssa_pos =  ssa_struct->at( ssa_struct->at(i).at(j).ssa ).at(0);
				
				while( ssa_pos.leaf == 0 )
				{
					ssa_pos = ssa_struct->at( ssa_pos.ssa ).at(0);	
				}
				
				ssa_struct->at(i).at(j).leaf_value = ssa_pos.ssa;
			} 	
		
		}
		
		// Sort all of the nodes	
		sort( ssa_struct->at(i).begin(), ssa_struct->at(i).end(), compare(sequence));	

	}
	
	delete( to_order );
	

	stack<pair<uint64_t,uint64_t>> * ssa_stack = new stack<pair<uint64_t,uint64_t>>();
	
	vector<uint64_t> * final_lcp = new vector<uint64_t>();
	
	uint64_t i = 0;
	uint64_t j = 0;
	bool found = false;
	
	// Depth-first search to find first occurrence of internal node
	for(i = 0; i<ssa_struct->size(); i++)
	{
		for(j = 0; j<ssa_struct->at(i).size(); j++)
		{
			if( ssa_struct->at(i).at(j).leaf == 1 )
			{
				output_ssa<<ssa_struct->at(i).at(j).ssa<<endl;
				output_lcp<<ssa_struct->at(i).at(j).lcp<<endl;
			}
			else
			{
				ssa_stack->push( make_pair<uint64_t, uint64_t>( (uint64_t) i,  (uint64_t)j ) );
				found = true;
				break;
			}
		}
		
		if( found == true )
			break;
	}
	
	uint64_t lcp = 0;
	bool leaf = false;
	
	// Add nodes to final array using stack
	while( ssa_stack->size() > 0 )
	{
		if( ssa_struct->at(i).at(j).leaf == 1 )
		{
		
			if( leaf == true )
				lcp = ssa_struct->at(i).at(j).lcp;
			
			output_lcp<<lcp<<endl;
			output_ssa<<ssa_struct->at(i).at(j).ssa<<endl;
			
			leaf = true;
			
			j++;
			
			while( ssa_stack->size() > 0 && j > ssa_struct->at(i).size()-1 )
			{
				i = ssa_stack->top().first; 
				j = ssa_stack->top().second+1;
					
				ssa_stack->pop();
				
			}	
		}
		else
		{
			
			if( leaf == true )
				lcp = ssa_struct->at(i).at(j).lcp;
	
			leaf = false;
			
			ssa_stack->push( make_pair<uint64_t, uint64_t>( (uint64_t) i, (uint64_t) j) );
		
			i = ssa_struct->at(i).at(j).ssa;
			j = 0;
				
		}
	}
	
	delete( ssa_stack );

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
	
 	/* Read in list of sparse suffixes */
  	ifstream suff_list(argv[2], ios::in | ios::binary);
  	suff_list.seekg(0, ios::end);
  	file_size = suff_list.tellg();
   
 	vector<vector<SSA>> * ssa_struct = new vector<vector<SSA>>();
	
 	
 	c = 0;
	string line="";
	suff_list.seekg (0, ios::beg);
	uint64_t b = 0;
	
	if( file_size == 0 )
	{
		cout<<"Empty suffixes list"<<endl;
		return 1;
	}
	
	vector<SSA> initial;
	
	for (uint64_t i = 0; i < file_size; i++)
	{
		suff_list.read(reinterpret_cast<char*>(&c), 1);
		
		if( ( char) c == '\n' || ( char) c == ' ' )
		{
		
			SSA in = {};
			in.ssa = stoi(line);
			
			in.leaf = 1;

			
			if( stoi(line) >= text_size )
			{
				cout<<"Suffix list contains a suffix larger than sequence length"<<endl;
				return 1;
			}	
			
			in.lcp = 0;
			initial.push_back(in);
			
			b++;
			line = "";
			
		}	
		else line += c;
	}
	suff_list.close();	
	
	ssa_struct->push_back( initial );
	
	karp_rabin_hashing::init();
	uint64_t l = pow(2, (uint64_t) log2(text_size));
	uint64_t fp_blocks = (text_size/(text_size/b)) ;
	uint64_t fp_len = text_size / b;
	
	// computing fingerprints
	uint64_t * FP =  ( uint64_t * ) calloc( fp_blocks , sizeof( uint64_t ) );
	uint64_t fp = 0;
	pos = 0;
	
	std::chrono::steady_clock::time_point start_total = std::chrono::steady_clock::now();
	
	for(uint64_t i = 0; i<fp_blocks *(text_size/b) ; i++)
	{
		fp = karp_rabin_hashing::concat(fp, sequence[i], 1);
		if( ( i + 1 ) % fp_len == 0 )
		{
			FP[pos] = fp ;
			pos++;
		}
		
		
	}
	
	while( l > 0 )
	{
		group( ssa_struct, FP, fp_len, l, sequence, text_size, b );
		l = l/2;
	}
	
	order( ssa_struct, sequence, text_size, argv[3] );

	delete( ssa_struct );
	free( sequence );
	free( FP );
	
	std::chrono::steady_clock::time_point end_total = std::chrono::steady_clock::now();
	std::cout<<"Time taken "<<std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total).count() << "[ms]" << std::endl;


	return 0;
}


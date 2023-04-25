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

using namespace std;

struct SSA
{
	
	uint64_t ssa;
	uint64_t lcp;
	uint64_t leaf;
	uint64_t leaf_value;
};


void mergeSort(vector<SSA> * array, int l, int r, uint64_t * sequence );

// Merge two subarrays L and M into array
void merge(vector<SSA> * array, int p, int q, int r, uint64_t * sequence) 
{
  	int n1 = q - p + 1;	
	int n2 = r - q;
	
	vector<SSA> * L = new vector<SSA>();
	vector<SSA> * M = new vector<SSA>();

	for (int i = 0; i < n1; i++)
		L->push_back(array->at(p + i));
	for (int j = 0; j < n2; j++)
		M->push_back(array->at(q + 1 + j));

	int i, j, k;
	i = 0;
	j = 0;
	k = p;

	while (i < n1 && j < n2) 
	{
		if ( sequence[ L->at(i).leaf_value+L->at(i).lcp ]  <= sequence[ M->at(j).leaf_value+M->at(j).lcp ]  ) 
		{
			array->at(k) = L->at(i);
			i++;
		} 
		else 
		{
	      		array->at(k) = M->at(j);
	      		j++;
	      		
	    	}
	   	k++;
	}

	while (i < n1) 
	{
	 	array->at(k) = L->at(i);
	 	i++;
  		k++;
	}


	while (j < n2) 
	{
	  	array->at(k) = M->at(j);
	    	j++;
		k++;
	}
	
	delete( L );
	delete( M );
}


void mergeSort(vector<SSA> * array, int l, int r, uint64_t * sequence  ) 
{
	if (l < r) 
	{
    		int m = l + (r - l) / 2 ;
    		
		mergeSort(array, l, m, sequence );
	    	mergeSort(array, m + 1, r, sequence );
	    	merge(array, l, m, r, sequence );
  	}
}
 
uint64_t modular(uint64_t base,  uint64_t exp,  uint64_t mod)
{
    uint64_t x = 1;
    uint64_t i;
    uint64_t power = base % mod;

    for (i = 0; i < sizeof(int) * 8; i++) {
        uint64_t least_sig_bit = 0x00000001 & (exp >> i);
        if (least_sig_bit)
            x = (x * power) % mod;
        power = (power * power) % mod;
    }

    return x;
}

/* Computing the fingerprint of an SSA */
uint64_t fingerprint( uint64_t ssa, uint64_t * FP, uint64_t fp_len, uint64_t l, uint64_t  * sequence, uint64_t text_size, uint64_t r, uint64_t q )
{

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
			fp_short =  ( q + ( ( fp_short * r +  sequence[c_start * fp_len + i] ) % q ) )  % q  ;
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
			fp_long = ( q  + ( ( fp_long * r +  sequence[c_end * fp_len + i] )  % q ) ) % q ;
		}
	}
	else fp_long = FP[closest_end-1];

	uint64_t power = modular( r, (end - 1) - (ssa - 1), q) ;
	uint64_t mult = fp_short * power  % q ;
	uint64_t fp =  ( q + ( fp_long - mult) ) % q  ;
	
	
return fp;
}

/* Grouping SSAs according to shared fingerprints */
uint64_t group( vector<vector<SSA>> * ssa_struct, uint64_t * FP, uint64_t fp_len, uint64_t l, uint64_t * sequence, uint64_t text_size, uint64_t r, uint64_t q )
{
	unordered_map<uint64_t, vector<SSA>> * groups = new unordered_map<uint64_t, vector<SSA>>();
	
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
			uint64_t fp = fingerprint( ssa+lcp, FP, fp_len, l, sequence, text_size, r, q );
				
			if( groups->count( fp ) == 0 )
				groups->insert(pair<uint64_t, vector<SSA>>(fp, {ssa_struct->at(i).at(j)}));
			else groups->at( fp ).push_back( ssa_struct->at(i).at(j) );
			
		}
		
		// Only create new node if some nodes were grouped together
		if( groups->size() != original_size )
		{
			ssa_struct->at(i).clear();
			
			for( auto it = groups->begin(); it!= groups->end(); it++)
			{
				if( it->second.size() > 1 )
				{
					// these nodes have been grouped together, place them in a new node and remove them from
					// current node and replace with new node ID.
				
					uint64_t original_lcp = it->second.at(0).lcp;
					
					vector<SSA> new_inp;
					for( int64_t k=0; k<it->second.size(); k++)
					{
						SSA new_ssa = {};
						new_ssa.leaf = it->second.at(k).leaf;
						new_ssa.ssa = it->second.at(k).ssa;
						new_ssa.lcp = it->second.at(k).lcp + l;
						
						new_inp.push_back( new_ssa );
							
					}
					ssa_struct->push_back( new_inp );
					
					
					// create new internal node and add to original ssa vector
					SSA new_internal = {};
					new_internal.leaf = 0;
					//new_internal.leaf_value = it->second.at(0).ssa;
					new_internal.ssa = ssa_struct->size()-1;
					new_internal.lcp = original_lcp;
						
					ssa_struct->at(i).push_back( new_internal );
				}
				else
				{
					// Nodes that were not grouped together can just be re-inserted into the array
					SSA original = {};
					original.leaf =  it->second.at(0).leaf;
					original.ssa =  it->second.at(0).ssa;
					original.lcp =  it->second.at(0).lcp;
						
					ssa_struct->at(i).push_back( original );
				}
			}
			
		}
		groups->clear();
	}
	
	delete( groups );
	

	
return 0;
}


uint64_t order( vector<uint64_t> * final_ssa, vector<vector<SSA>> * ssa_struct, uint64_t * sequence, uint64_t text_size )
{
	vector<SSA> * to_order = new vector<SSA>();
	
	for(uint64_t i = 0; i<ssa_struct->size(); i++)
	{
		// Depth first search to identify leaf/ssa values for each internal node
		for(uint64_t j = 0; j<ssa_struct->at(i).size(); j++)
		{
			SSA new_ssa = {};
			
			if( ssa_struct->at(i).at(j).leaf == 1 )
			{
				new_ssa.ssa = ssa_struct->at(i).at(j).ssa;
				new_ssa.lcp = ssa_struct->at(i).at(j).lcp;
				new_ssa.leaf = ssa_struct->at(i).at(j).leaf;
				new_ssa.leaf_value = ssa_struct->at(i).at(j).ssa;
				
			}
			else
			{
				new_ssa.ssa = ssa_struct->at(i).at(j).ssa;
				new_ssa.lcp = ssa_struct->at(i).at(j).lcp;
				new_ssa.leaf = ssa_struct->at(i).at(j).leaf;
				
				SSA ssa_pos =  ssa_struct->at( ssa_struct->at(i).at(j).ssa ).at(0);
				
				while( ssa_pos.leaf == 0 )
				{
					ssa_pos = ssa_struct->at( ssa_pos.ssa ).at(0);	
				}
				
				new_ssa.leaf_value = ssa_pos.ssa;
			} 	
			
			to_order->push_back( new_ssa );
		
		
		}
		
		// Merge-sort all of the nodes		
		mergeSort( to_order, 0, to_order->size()-1, sequence );
		
		ssa_struct->at(i).clear();
		
		for(uint64_t j = 0; j<to_order->size(); j++)
		{
			
			ssa_struct->at(i).push_back( to_order->at(j) );
		
		}
		
		to_order->clear();
	
	}
	
	delete( to_order );
	

	stack<pair<uint64_t,uint64_t>> * ssa_stack = new stack<pair<uint64_t,uint64_t>>();
	//vector<uint64_t> * final_ssa = new vector<uint64_t>();
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
				
				final_ssa->push_back( ssa_struct->at(i).at(j).ssa );
				final_lcp->push_back( ssa_struct->at(i).at(j).lcp );
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
			
			final_lcp->push_back( lcp);
				
			leaf = true;
			final_ssa->push_back( ssa_struct->at(i).at(j).ssa );
		
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
	
	for(uint64_t i = 0; i<final_ssa->size(); i++)
		cout<<final_ssa->at(i)<<" ";
		
	cout<<endl;
	
	for(uint64_t i = 0; i<final_lcp->size(); i++)
		cout<<final_lcp->at(i)<<" ";
	
	cout<<ssa_struct->size()<<endl;
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
   	uint64_t * sequence =  ( uint64_t * ) calloc( file_size , sizeof( uint64_t ) );
   	
	char c = 0;
	seq.seekg (0, ios::beg);
	uint64_t pos = 1;
	
	
	if( file_size == 0 )
		cout<<"Empty sequence file"<<endl;
		
	for (uint64_t i = 0; i < file_size; i++)
	{
		seq.read(reinterpret_cast<char*>(&c), 1);

		if( ( char) c == '\n' || ( char) c == ' ' )
			continue;
		else
		{	
			sequence[text_size] =  (int) c; 
			text_size++;
		}
		
	}
	seq.close();
	

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
		cout<<"Empty suffixes list"<<endl;
	
	for (uint64_t i = 0; i < file_size; i++)
	{
		suff_list.read(reinterpret_cast<char*>(&c), 1);
		
		if( ( char) c == '\n' || ( char) c == ' ' )
		{
			ssa_list->push_back(stoi(line));
			
			slcp_list->push_back(0);
			
			line = "";
			
		}	
		else line += c;
		
		
	}
	suff_list.close();	
	

	uint64_t b = ssa_list->size();
	uint64_t l = pow(2, (uint64_t) log2(text_size));
	uint64_t fp_blocks = (text_size/(text_size/b)) ;
	uint64_t fp_len = text_size / b;
	
	// computing fingerprints
	uint64_t * FP =  ( uint64_t * ) calloc( fp_blocks , sizeof( uint64_t ) );
	uint64_t q = 5000298533; // q chosen at random, constant for now
	srand (time (0));
	uint64_t r = 2 + (rand() % q-1);

	cout<<" q and r "<<q<<" "<<r<<endl;
	uint64_t fp = 0;
	pos = 0;
	
	std::chrono::steady_clock::time_point start_total = std::chrono::steady_clock::now();
	
	// FP of first block
	for(uint64_t i = 0; i<fp_len; i++)
	{
		fp =  ( q + ( fp * r + sequence[i]  % q ) )  % q ;
	}
	FP[0] = fp ;
	pos++;
	
	uint64_t fp2 = 0;
	for(uint64_t i = fp_len; i<fp_blocks *(text_size/b) ; i++)
	{
		
		fp =  ( q + ( fp * r + sequence[i]  % q ) )  % q ; 
		

		if( i > 0 &&  ( i + 1 ) % fp_len == 0 )
		{
			FP[pos] = fp ;
			pos++;
		}
		
		
	}

	vector<vector<SSA>> * ssa_struct = new vector<vector<SSA>>();
	
	vector<SSA> initial;
	
	for(uint64_t i = 0; i<ssa_list->size(); i++ )
	{	
		SSA in = {};
		in.ssa = ssa_list->at(i);
		in.lcp = slcp_list->at(i);
		in.leaf = 1;
		initial.push_back(in);
		
	}
	
	ssa_struct->push_back( initial );
	
	delete( ssa_list );
	delete( slcp_list );
	
	while( l > 0 )
	{
		group( ssa_struct, FP, fp_len, l, sequence, text_size, r, q ); 
		l = l/2;
		
	}
	
	vector<uint64_t> * final_ssa = new vector<uint64_t>();
	
	order( final_ssa, ssa_struct, sequence, text_size );
	
	std::chrono::steady_clock::time_point end_total = std::chrono::steady_clock::now();
	std::cout << endl<<"Time taken "<<std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total).count() << "[ms]" << std::endl;
	

	free( sequence );
	return 0;
}


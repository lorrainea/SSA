/*
 *  This file is part of rk-lce.
 *  Copyright (c) by
 *  Nicola Prezza <nicolapr@gmail.com>
 *
 *   rk-lce is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   rk-lce is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details (<http://www.gnu.org/licenses/>).
 */

/*
 * Modification: This code is modified from the original to read the input 
 * from text files and to compute and output the sparse LCP array from the 
 * LCE data structure. These modifications were done by the authors of SSA.
 */


#include <includes.hpp>
#include <rk_lce.hpp>

using namespace std;
using namespace rklce;

string output_file = string();
string input_text = string();
string input_pos = string();

void help(){
	cout << "sa-rk: build subset of suffix array." << endl << endl;
	cout << "Usage: sa-rk <input text> <output> [input positions]" << endl;
	cout << "   <input text>        Input text file name. Required."<<endl;
	cout << "   <output>            Output file name (suffix-sorted positions. 64-bit unsigned integers). Required." << endl;
	cout << "   [input positions]   List of text positions to sort, stored as 64-bit unsigned integers. If not specified, " << endl;
	cout <<	"                       compute the full suffix array. NOTE: this file must contain s+1 64-bit unsigned " << endl;
	cout <<	"                       integers; the first integer is s (size), the others are the subset of s text positions " << endl;
	cout <<	"                       to suffix-sort." << endl << endl;
	exit(0);
}

int main(int argc, char** argv){

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;

	//parse options

	if(argc<3 or argc > 4) help();

	input_text = string(argv[1]);
	output_file = string(argv[2]);

	if(argc == 4)
		input_pos = string(argv[3]);

	if(input_pos.compare("")==0){

		cout << "Computing full suffix array of the file " << input_text << endl;
		cout << "Output will be stored in " << output_file << endl;

	}else{

		cout << "Computing partial suffix array of the file " << input_text << endl;
		cout << "Reading positions to sort from file " << input_pos << endl;
		cout << "Output will be stored in " << output_file << endl;

	}

	cout << endl;

    auto t1 = high_resolution_clock::now();

    cout << "Building LCP data structure and allocating memory for SA ... " << endl;

    auto lce = rk_lce(input_text);

    //test_LCP(T,1000000);exit(0);

    cout << " File size is " << lce.size() << endl;
    cout << " Alphabet size is (rounded to a power of 2) " << lce.alphabet_size() << endl;
    cout << " Size of LCE structure (Bytes) " << lce.bit_size()/8 << endl;

    //Suffix Array Subset
    vector<uint64_t> SA;

    if(input_pos.compare("")==0){

    	//compute full suffix array: begin with vector <0, 1, 2, ..., n-1>

    	SA = vector<uint64_t>(lce.length());
    	for(uint64_t i=0;i<SA.size();++i) SA[i] = i;

    }else{

	cout<<input_pos.c_str()<<endl;
	std::ifstream input(input_pos.c_str());

		SA = vector<uint64_t>();
		
		for( std::string line; getline( input, line ); )
		{
		   	uint64_t s;
		   	std::istringstream iss(line);
  			  iss >> s;
			SA.push_back( s);
			
		}

		
	

    }

    auto t2 = high_resolution_clock::now();


    cout << "Computing suffix array ... " << flush;

    //SA-RK algorithm
    std::sort(SA.begin(), SA.end(), lce.lex_less_than());



    auto t3 = high_resolution_clock::now();

    cout << "Saving suffix array to "  << output_file << " ... " << flush;

	vector<uint64_t> slcp;
    
    string output_file_ssa = output_file + "_rk_lce.ssa";
    string output_file_lcp = output_file + "_rk_lce.lcp";
    	
   ofstream ofs(output_file_ssa);
    	for(int i = 0; i<SA.size(); i++)
    		ofs << SA[i]<<"\n";

    	//ofs.write((char*)SA.data(),SA.size()*sizeof(uint64_t));
    
    
    slcp.push_back(0);
    for(int i = 1; i<SA.size(); i++)
    
    {
    	uint64_t lcp = lce.LCE( SA[i-1], SA[i] );
    	slcp.push_back(lcp);
    }

	ofstream ofs2(output_file_lcp);
    	
    		
	for(int i = 0; i<slcp.size(); i++)
    		ofs2 << slcp[i]<<"\n";
    		

    auto t4 = high_resolution_clock::now();

	uint64_t constr = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	uint64_t comp = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();
	uint64_t save = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
	uint64_t tot = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t1).count();

	cout << "Time to build LCE structure: " << (double)constr << " ms" << endl;
	cout << "Time to suffix-sort text positions: " << (double)comp << " ms" << endl;
	cout << "Save time: " << (double)save << " ms" << endl;
	cout << "Total time: " << (double)tot << " ms" << endl;

}

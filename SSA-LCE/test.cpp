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

#include <includes.hpp>
#include "internal/rk_lce.hpp"
#include "internal/rk_lce_bin.hpp"
#include "internal/bitv.hpp"

using namespace std;
using namespace rklce;

string output_file = string();
string input_text = string();
string input_pos = string();

void test_bin_lce(){
	srand(time(NULL));
	uint64_t n = 10000*127;

	vector<bool> B(n);

	//probability of a 1
	double p = 0.8;

	for(int i=0;i<n;++i){

		if(i<10)
			B[i] = 0;
		else
			B[i] = double(rand())/RAND_MAX < p;

	}

	auto lce = rk_lce_bin(B);

	cout << "Size of the input: " << n/8 << " bytes" << endl;
	cout << "Size of the structure: " << lce.bit_size()/8 << " bytes" << endl;

	cout << "testing access ... " << flush;

	for(int i=0;i<B.size();++i){

		assert(B[i]==lce[i]);

	}

	cout << "done.\ntesting LCE ... "<<flush;

	for(int t =0;t<500000;++t){

		auto i = rand()%n;
		auto j = rand()%n;

		auto lce_n = lce.LCE_naive(i,j);
		auto lce_f = lce.LCE(i,j);

		if(lce_n!=lce_f){

			cout << "i = " << i << ", j = " << j << endl;
			cout << "naive = " << lce_n << endl;
			cout << "fast = " << lce_f << endl;

			for(int ki = i;ki<i+lce_f+1;++ki) cout << B[ki];cout << endl;
			for(int kj = j;kj<j+lce_f+1;++kj) cout << B[kj];cout << endl;

			exit(0);

		}

	}

	cout << "\n\ndone" << endl;
}

int main(int argc, char** argv){

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;

	srand(time(NULL));

	string text_path = string(argv[1]);
	int rep = atoi(argv[2]);

    auto t0 = high_resolution_clock::now();
	time_t x=0;

    /*
     * we measure the time taken to access the string stored as an array of chars to
     * make a comparison with the structure
     */
	{

		ifstream ifs(text_path);

		stringstream strStream;
		strStream << ifs.rdbuf();//read the file
		string str = strStream.str();//str holds the content of the file

	    t0 = high_resolution_clock::now();

	    for(int k=0;k<rep;++k){

			auto c = str[rand()%str.size()];
			x |= c;

		}

		/* since we use x, the compiler is forced to execute the above operations
		* (otherwise some compilers may optimize and not execute the accesses)		
		*/
		srand(time(&x));

	}

	cout << "Building structure ... " << endl;

    auto t1 = high_resolution_clock::now();

    rk_lce lce(text_path);

	cout << "Size of the input: " << lce.size() << " Bytes" << endl;
	cout << "Size of the structure: " << lce.bit_size()/8 << " Bytes" << endl;
	cout << "Difference (structure-input): " << int64_t(lce.bit_size()/8) -  int64_t(lce.size()) << " Bytes" << endl;

	cout << "Testing LCE ... " << endl;

	double sum_lce=0;

    	auto t2 = high_resolution_clock::now();

	for(int k=0;k<rep;++k){

		uint64_t i = rand()%lce.size();
		uint64_t j = rand()%lce.size();

		auto lce_fast = lce.LCE(i,j);

		sum_lce += lce_fast;

		/*if(lce_fast!= lce.LCE_naive(i,j)){
			cout << "error."<<endl;
			exit(0);
		}*/

	}

	cout << "Average LCE: " << sum_lce/rep << endl;

    auto t3 = high_resolution_clock::now();

	cout << "Testing access ... " << endl;

	for(int k=0;k<rep;++k){

		x |= lce[rand()%lce.size()];

	}

    srand(time(&x));

    auto t4 = high_resolution_clock::now();

	uint64_t access_str = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
	uint64_t build = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	uint64_t compute = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();
	uint64_t access = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();

	cout << "Build time: " << (double)build << " ms" << endl << endl;
	cout << "LCE time: " << (double)compute << " ms" << endl;
	cout << "Time per LCE query: " << (double)compute/rep << " ms" << endl<<endl;
	cout << "Access time: " << (double)access << " ms" << endl;
	cout << "Time per access query: " << (double)access/rep << " ms" << endl<<endl;
	cout << "Access time for a string object: " << (double)access_str << " ms" << endl;
	cout << "Time per access query on a string object: " << (double)access_str/rep << " ms" << endl<<endl;

}

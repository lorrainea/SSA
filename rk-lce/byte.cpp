#include<fstream>
#include <iostream>
#include <inttypes.h>
#include <stdio.h>
using namespace std;


int main(void)
{
  FILE *fout;
   fout = fopen("output", "wt");
	uint64_t first = 3;
	uint64_t second = 1;
	uint64_t third = 2;
	uint64_t fourth = 3;
	fputc((first & 0xFF), fout);
	 first >>= 8;
	fputc((second & 0xFF), fout);
	 second >>= 8;
	fputc((third & 0xFF), fout);
	 third >>= 8;
	fputc((fourth & 0xFF), fout);
    fourth >>= 8;
    
}


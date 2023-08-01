#include <iostream>
#include <unordered_set>
#include <algorithm>
#include <set>
using namespace std;
  
int main(int argc, char **argv)
{
	if( argc < 3 )
        {
                cout<<"Check arguments!\n";
                cout<<"./random <n> <#suffixes>\n";
                exit(-1);
    	}
    	uint64_t N = (uint64_t) atoi(argv[1]);
    	uint64_t suf = (uint64_t) atoi(argv[2]);
	
	if (suf > N )
        {
                cout<<"Check arguments!\n";
                cout<<"n must be greater than #suffixes!!!\n";
                exit(-1);
    	}

	unordered_set<uint64_t> set;
    	uint64_t * S;
    	S = (uint64_t *)calloc(suf, sizeof(uint64_t));
    	
	for (uint64_t i = 0; i < suf; i++)	
	{
		while(1)
		{
			uint64_t key = rand() % N;
			if ( set.find(key) == set.end() )	
			{
				S[i] = key;
				set.insert(key);
				break;
			}
		}
	}

	std::sort(S, S + suf);
    	
	for (uint64_t i = 0; i < suf; i++)	cout<<S[i]<<endl;
  
    	free(S);

    	return 0;
}

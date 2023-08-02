#include <iostream>
#include <fstream>
#include <vector>
#include <string>

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <file1> <file2>" << std::endl;
        return 1;
    }

    std::string file1_name = argv[1];
    std::string file2_name = argv[2];
    std::ifstream file1(file1_name);
    std::ifstream file2(file2_name);

    if (!file1 || !file2) {
        std::cerr << "Error opening files." << std::endl;
        return 1;
    }

    std::vector<int> data1;
    std::vector<int> data2;
    int value;

    // Read integers from file1 and store them in data1 vector
    int line_number = 1;
    while (file1 >> value) {
        data1.push_back(value);
        line_number++;
    }

    // Read integers from file2 and store them in data2 vector
    line_number = 1; // Reset line_number for the second file
    while (file2 >> value) {
        data2.push_back(value);
        line_number++;
    }

    // Close the files after reading
    file1.close();
    file2.close();

    int diff_lines=0;
    // Compare the two vectors and display the lines at which they are different
    bool are_equal = true;
    size_t min_size = std::min(data1.size(), data2.size());
    for (size_t i = 0; i < min_size; ++i) {
        if (data1[i] != data2[i]) {
            std::cout << "Files are different at line " << i + 1 << std::endl;
            are_equal = false;
	    diff_lines++;	
        }
    }

    if (!are_equal) {
        std::cout << "Files have "<<diff_lines<<" different number of lines." << std::endl;
        
    }

    if (are_equal) {
        std::cout << "The two files are identical." << std::endl;
    }

    return 0;
}

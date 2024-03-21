#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <endian.h>

using namespace std;

void binRead(float outputArray[], const string& fileName, const unsigned int length) {
	ifstream inputData;
  inputData.open(fileName, ios::binary);
	if (inputData) {
		inputData.read(reinterpret_cast<char*>(outputArray), 4);
		inputData.close();
	} 
}


int main() {

	int ignoreInt;
	float ignorefloat;
	double value;
  
	string file = "../h2.r4dscf.r4dorbene";
//	string file = "../h2o.631g.r4dscf.r4dorbene";
	ifstream eneFile;
  eneFile.open(file.c_str(), ios::binary);

  if (!eneFile) {
    cerr << "The file can't be open." << endl;
    exit(EXIT_FAILURE);
  }
	
	eneFile.read((char*)&ignoreInt, 4);
	cout << " fucicj: " << ignoreInt << endl;
	
	int i;
	i = 0;
	while (eneFile.read((char*)&value, 8)) {
		cout << i << " " << value << endl;
		i++;
	}
  eneFile.close();

	return 0;
}


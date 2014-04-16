#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include "pmsig.h"


int main(int argc, char** argv) {

	if (argc < 4) {
		std::cerr << "Usage: " << argv[0] << " mutation_infomation_file output_directory number_of_clusters" << std::endl;
		return -1;
	}

	std::string mut_file = argv[1];
	std::string out_dir = argv[2];
	int cluster_num = atoi(argv[3]);

	long unsigned int MCMCNum = 1000;
	long unsigned int burnInNum = 100;	

	srand(static_cast<std::size_t>(time(NULL)));

	PMSig pmsig(mut_file);
	pmsig.preparation(cluster_num, static_cast<std::size_t>(MCMCNum));

	// ********************
	// burn-in step
	for (long unsigned int i = 0; i < burnInNum; i++) {
		pmsig.gibbsUpdate();
		if ((i % 10 == 0) & (i > 0)) {
			std::cout << i << " times finished." << "\n";
		}
	}
	// *********************
	std::cerr << "The burn-in cycles are finished.\n";

	for (long unsigned int i = 0; i < MCMCNum; i++) {
		pmsig.gibbsUpdate();
		pmsig.incrementParam();
		pmsig.updateBayesianDeviance(i);
		if (i % 100 == 0) {
			std::cout << i << " times finished." << "\n";
		}
	}

     
	pmsig.printMean_Phi(out_dir);
	pmsig.printMean_Psi(out_dir);
	pmsig.printMean_Theta(out_dir);
	pmsig.printPenalizedMeanBayesianDeviance(out_dir);
	
	return(0);

}

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <unistd.h>
#include "pmsig.h"


int main(int argc, char** argv) {

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " mutation_infomation_file output_directory number_of_clusters" << std::endl;
        return -1;
    }

    std::string mut_file = argv[1];
    std::string out_dir = argv[2];

    int cluster_num = 5;
    long unsigned int MCMCNum = 10000;
    long unsigned int burnInNum = 1000;

    int c = 0;
    while ((c = getopt (argc, argv, "k:m:n:")) != -1)
        switch (c) {
            case 'k':
                cluster_num = atoi(optarg);
                break;
            case 'm': 
                burnInNum = atoi(optarg);
                break;
            case 'n':
                MCMCNum = atoi(optarg); 
                break;
           case '?':
             if (optopt == 'k')
               fprintf (stderr, "Option -%c requires an argument.\n", optopt);
             else if (isprint (optopt))
               fprintf (stderr, "Unknown option `-%c'.\n", optopt);
             else
               fprintf (stderr,
                        "Unknown option character `\\x%x'.\n",
                        optopt);
             return 1;
           default:
             abort ();
           }
     

    srand(static_cast<std::size_t>(time(NULL)));

    std::string log_file = out_dir + ".log.txt";
    std::ofstream lofs(log_file.c_str());
    if (lofs.fail()) {
        std::cerr << "Error: Could not open the output file " << log_file << "\n";
        exit(8);
    }

    std::cerr << "Reading input data from: " << mut_file << "\n";
    lofs << "Reading input data from: " << mut_file << "\n";
    lofs.flush();

    PMSig pmsig(mut_file);
    std::cerr << "Input data information" << "\n";
    lofs << "Input data information" << "\n";
    std::cerr << "# of mutations : " << pmsig.getMutationNum() << "\n";
    lofs << "# of mutations : " << pmsig.getMutationNum() << "\n";
    std::cerr << "# of samples : " << pmsig.getSampleNum() << "\n";
    lofs << "# of samples : " << pmsig.getSampleNum() << "\n";
    std::cerr << "size of adjacent sequence : " << pmsig.getSeqSize() << "\n\n";
    lofs << "size of adjacent sequence : " << pmsig.getSeqSize() << "\n\n";
    lofs.flush();

    std::cerr << "Model parameters" << "\n";
    lofs << "Model parameters" << "\n";
    std::cerr << "# of mutation signatures : " <<  cluster_num << "\n\n";
    lofs << "# of mutation signatures : " <<  cluster_num << "\n\n";
    std::cerr << "Learning parameters" << "\n";
    lofs << "Learning parameters" << "\n";
    std::cerr << "# of burn-in cycles : " << burnInNum << "\n";
    lofs << "# of burn-in cycles : " << burnInNum << "\n";
    std::cerr << "# of MCMC cycles : " << MCMCNum << "\n\n";
    lofs << "# of MCMC cycles : " << MCMCNum << "\n\n";
    lofs.flush();

    pmsig.preparation(cluster_num, static_cast<std::size_t>(MCMCNum));

    // ********************
    // burn-in step
    for (long unsigned int i = 0; i < burnInNum; i++) {
        pmsig.gibbsUpdate();
        if ((i % 1000 == 0) & (i > 0)) {
            std::cerr << i << " times finished." << "\n";
            lofs << i << " times finished." << "\n";
            lofs.flush();
        }
    }
    // *********************
    std::cerr << "The burn-in cycles are finished.\n";

    for (long unsigned int i = 0; i < MCMCNum; i++) {
        pmsig.gibbsUpdate();
        pmsig.incrementParam();
        pmsig.updateBayesianDeviance(i);
        if (i % 1000 == 0) {
            std::cout << i << " times finished." << "\n";
            lofs << i << " times finished." << "\n";
            lofs.flush();
        }
    }

     
    pmsig.printMean_Phi(out_dir);
    pmsig.printMean_Psi(out_dir);
    pmsig.printMean_Theta(out_dir);
    pmsig.printPenalizedMeanBayesianDeviance(out_dir);
    
    lofs.close();
    return(0);

}

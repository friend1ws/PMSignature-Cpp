#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <unistd.h>
#include <iomanip>
#include "pmsig.h"

int main(int argc, char** argv) {

    std::clock_t start;
    double duration;
    start = std::clock();

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " method mutation_infomation_file output_directory" << std::endl;
        return -1;
    }

    // int method = 0;
    // if (static_cast<std::string>(argv[1]) == "full") {
    //     method = 1;
    // }

    std::string mut_file = argv[1];
    std::string out_dir = argv[2];

    int cluster_num = atoi(argv[3]);
    int maxRepeatNum = 1000000;

    int c = 0;
    while ((c = getopt(argc, argv, "k:r:")) != -1) {
        switch (c) {
            case 'k':
                cluster_num = atoi(optarg);
                break;
            case 'r': 
                maxRepeatNum = atoi(optarg);
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
    std::cerr << "# of max repeat cycles : " << maxRepeatNum << "\n";
    lofs << "# of max repeat cycles : " << maxRepeatNum << "\n\n";
    lofs.flush();

    pmsig.preparation(cluster_num);

    // pmsig.printF(out_dir);
    // pmsig.printQ(out_dir);

    // double diffParam = 0;
    double prevL = 0;
    double currL = 0;
    int cycle = 0;
    // ********************
    for (cycle = 0; cycle < maxRepeatNum; cycle++) {
        pmsig.EstepUpdate();
        pmsig.MstepUpdate();
       
        // diffParam = pmsig.getDiffF() + pmsig.getDiffQ();
        currL = pmsig.getLikelihood();
        if (cycle > 0) {
            // std::cout << cycle << "\t" << currL << "\t" << currL - prevL << "\t" << pmsig.getDiffF() << "\t" << pmsig.getDiffQ() << "\n";
            if (currL - prevL < 1e-4) {
                break;
            }
        }
        prevL = pmsig.getLikelihood();
        // std::cout << "printTheta" << "\n";
        // pmsig.printTheta(out_dir);
        // std::cout << "printF" << "\n";
        // pmsig.printF(out_dir);
        // std::cout << "printQ" << "\n";
        // pmsig.printQ(out_dir);
    }
    // *********************

    pmsig.printF(out_dir);
    pmsig.printQ(out_dir);   

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    std::cerr << "Learning Result" << "\n";
    lofs << "Learning Result" << "\n";
    std::cerr << std::fixed << std::setprecision(4) << "log-likelihood : " << currL << "\n";
    lofs << std::fixed << std::setprecision(4) << "log-likelihood : " << currL << "\n";
    std::cerr << "# of cycles for estimation : " << cycle << "\n";
    lofs << "# of cycles for estimation : " << cycle << "\n";
    std::cerr << "execution time : " << duration << "\n";
    lofs << "execution time : " << duration << "\n";

    lofs.flush();
 
    lofs.close();
    return(0);

}

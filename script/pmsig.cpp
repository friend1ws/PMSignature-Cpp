#include "pmsig.h"
#include <cmath>


PMSig::PMSig(const std::string& file) {

    std::ifstream ifs(file.c_str());
    if (ifs.fail()) {
        std::cerr << "Error: Could not open the input file.\n";
        exit(8);
    }

    // ********************
    // first row, the number of samples and mutation features
    std::string line;
    std::getline(ifs, line);
    std::istringstream is1(line);
  
    is1 >> S_; 
    is1 >> L_;
    mutFeatureDim_ = new int[L_];
    // ********************

    // ********************
    // the second row, the possible number of each mutation feature
    std::getline(ifs, line);
    std::istringstream is2(line);

    int tL = 0;
    int possibleMutationPatterns = 1;
    while(is2 >>mutFeatureDim_[tL]) {
        possibleMutationPatterns = possibleMutationPatterns * mutFeatureDim_[tL];
        tL = tL + 1;
    }
    P_ = possibleMutationPatterns;

    if (L_ != tL) {
        std::cerr << "At the first raw, the number of mutation feature is inconsistent.\n";
    }

    mutData_ = new int*[S_];
    for (int s = 0; s < S_; ++s) {
        mutData_[s] = new int[P_];
        for (int p = 0; p < P_; p++) {
            mutData_[s][p] = 0;
        }
    }
    // ********************


    // ********************
    // for each row, sample ID and mutation features for each mutation
    int curNum = 0;
    int sampleID;
    int* mutationFeature;
    int mutPatternID;
    int tempDigit;
    mutationFeature = new int[L_];

    while(std::getline(ifs, line)) {

        std::istringstream is(line);
        curNum = curNum + 1;

        is >> sampleID;
        if (sampleID <= 0 or sampleID > S_) {
            std::cerr << "At the " << curNum << " th key, the input data for sample ID is inconsistent.\n";
            exit(8);
        }   

        mutPatternID = 0;
        tempDigit = 1;
        for (int i = 0; i < L_; i++) {
            is >> mutationFeature[i];
            if (mutationFeature[i] <= 0 or mutationFeature[i] > mutFeatureDim_[i]) {
                std::cerr << "At the " << curNum << " th key, the input data for mutation features is inconsistent.\n";
                exit(8);
            }
            mutPatternID = mutPatternID + tempDigit * (mutationFeature[i] - 1);
            tempDigit = tempDigit * mutFeatureDim_[i]; 
        }

        mutData_[sampleID - 1][mutPatternID] = mutData_[sampleID - 1][mutPatternID] + 1;

    }
    // ********************

    // the number of mutations
    N_ = curNum;

    delete [] mutationFeature;
}


PMSig::~PMSig() {

    delete [] mutFeatureDim_;

    for (int s = 0; s < S_; ++s) {
        delete [] mutData_[s];
    }
    delete [] mutData_;


    for (int s = 0; s < S_; ++s) {
        for (int k = 0; k < K_; ++k) {
            delete [] Theta_[s][k];
        }
        delete [] Theta_[s];
    }
    delete [] Theta_;


    for (int k = 0; k < K_; ++k) {
        for (int l = 0; l < L_; ++l) {
            delete [] F_[k][l];
        }
        delete [] F_[k];
    }
    delete [] F_;

    for (int k = 0; k < K_; ++k) {
        for (int l = 0; l < L_; ++l) {
            delete [] newF_[k][l];
        }
        delete [] newF_[k];
    }
    delete [] newF_;


    for (int s = 0; s < S_; ++s) {
        delete [] Q_[s];
    }
    delete [] Q_;

    for (int s = 0; s < S_; ++s) {
        delete [] newQ_[s];
    }   
    delete [] newQ_;

    delete [] currentDigits_;


}






void PMSig::preparation(const int K) {

    K_ = K; 

    // ********************
    // Allocate memories for parameters
    try {

        Theta_ = new double**[S_];
        for (int s = 0; s < S_; ++s) {
            Theta_[s] = new double*[K_];
            for (int k = 0; k < K_; ++k) {
                Theta_[s][k] = new double[P_];
                for (int p = 0; p < P_; p++) {
                    Theta_[s][k][p] = 0;
                }
            }
        }

        F_ = new double**[K_];
        for (int k = 0; k < K_; ++k) {
            F_[k] = new double*[L_];
            for (int l = 0; l < L_; ++l) {
                F_[k][l] = new double[mutFeatureDim_[l]];
                for (int m = 0; m < mutFeatureDim_[l]; m++) {
                // F_[k][l] = new double[6];
                // for (int m = 0; m < 6; m ++) {
                    F_[k][l][m] = 0;
                }
            }
        }

        Q_ = new double*[S_];
        for (int s = 0; s < S_; ++s) {
            Q_[s] = new double[K_];
            for (int k = 0; k < K_; ++k) {
                Q_[s][k] = 0;
            }   
        }   
    
        newF_ = new double**[K_];
        for (int k = 0; k < K_; ++k) {
            newF_[k] = new double*[L_];
            for (int l = 0; l < L_; ++l) {
                newF_[k][l] = new double[mutFeatureDim_[l]];
                for (int m = 0; m < mutFeatureDim_[l]; m++) {
                // newF_[k][l] = new double [6];
                // for (int m = 0; m < 6; m++) {
                    newF_[k][l][m] = 0;
                }
            }
        }

        newQ_ = new double*[S_];
        for (int s = 0; s < S_; ++s) {
            newQ_[s] = new double[K_];
            for (int k = 0; k < K_; ++k) {
                newQ_[s][k] = 0;
            }
        }

        currentDigits_ = new int[L_];
        for (int l = 0; l < L_; l++) {
            currentDigits_[l] = 0;
        }

    
    }



    catch (...) {
        std::cerr << "init(): Out of memmory" << std::endl;
        exit(EXIT_FAILURE);
    }

    double tempRand;
    double tempSum;
    for (int k = 0; k < K_; ++k) {
        for (int l = 0; l < L_; ++l) {
            
            tempSum = 0;
            for (int m = 0; m < mutFeatureDim_[l]; m++) {
                tempRand = gen_rand(1.0);
                F_[k][l][m] = tempRand;
                tempSum = tempSum + tempRand;
            }

            for (int m = 0; m < mutFeatureDim_[l]; m++) {
                F_[k][l][m] = F_[k][l][m] / tempSum;
            }

        }
    }


    for (int s = 0; s < S_; ++s) {
        
        tempSum = 0;
        for (int k = 0; k < K_; ++k) {
            tempRand = gen_rand(1.0);
            Q_[s][k] = tempRand;
            tempSum = tempSum + tempRand;
        }   

        for (int k = 0; k < K_; ++k) {
            Q_[s][k] = Q_[s][k] / tempSum;
        }

    }   


}



void PMSig::EstepUpdate() {

    for (int s = 0; s < S_; ++s) {
        for (int k = 0; k < K_; ++k) {
            for (int p = 0; p < P_; p++) {
                Theta_[s][k][p] = Q_[s][k];
            }
        }
    }


    for (int l = 0; l < L_; l++) {
        currentDigits_[l] = 0;
    }

    for (int k = 0; k < K_; k++) {
        for (int p = 0; p < P_; p++) {
        
            for (int l = 0; l < L_; l++) {
                for (int s = 0; s < S_; s++) {
                    Theta_[s][k][p] = Theta_[s][k][p] * F_[k][l][currentDigits_[l]];
                }
            }

            int tl = 0;
            while(currentDigits_[tl] + 1 >= mutFeatureDim_[tl]) {
                currentDigits_[tl] = 0;
                tl = tl + 1;
            }
            currentDigits_[tl] = currentDigits_[tl] + 1;
       
            // for debug
            // std::cerr <<  currentDigits[0];
            // for (int l = 1; l < L_; l++) {
            //     std::cerr << "\t" << currentDigits[l]; 
            // }
            // std::cerr << "\n";

        }     
    }

    double tempSum;
    likelihood_ = 0;
    for (int s = 0; s < S_; s++) {
        for (int p = 0; p < P_; p++) {

            tempSum = 0;
            for (int k = 0; k < K_; k++) {
                tempSum = tempSum + Theta_[s][k][p];
            }

            if (mutData_[s][p] > 0) {
                likelihood_ = likelihood_ + mutData_[s][p] * log(tempSum);
            }

            if (tempSum > 0) {
                for (int k = 0; k < K_; k++) {
                    Theta_[s][k][p] = Theta_[s][k][p] / tempSum;
                }
            } else {
                for (int k = 0; k < K_; k++) {
                    Theta_[s][k][p] = 1.0 /static_cast<double>(K_);
                }
            }

        }
    }

}


void PMSig::MstepUpdate() {


    for (int k = 0; k < K_; k++) {

        for (int l = 0; l < L_; l++) {
            currentDigits_[l] = 0;
        }

        for (int p = 0; p < P_; p++) {

            for (int s = 0; s < S_; s++) {

                newQ_[s][k] = newQ_[s][k] + mutData_[s][p] * Theta_[s][k][p];
                for (int l = 0; l < L_; l++) {
                    newF_[k][l][currentDigits_[l]] = newF_[k][l][currentDigits_[l]] + mutData_[s][p] * Theta_[s][k][p];
                }

            }

            int tl = 0;
            while(currentDigits_[tl] + 1 >= mutFeatureDim_[tl]) {
                currentDigits_[tl] = 0;
                tl = tl + 1;
            }
            currentDigits_[tl] = currentDigits_[tl] + 1;

            // for debug
            // std::cout <<  currentDigits[0];
            // for (int l = 1; l < L_; l++) {
            //     std::cout << "\t" << currentDigits[l]; 
            // }
            // std::cout << "\n";

        }

    }

    // std::cout << "print newF" << "\n";
    // for (int k = 0; k < K_; k++) {
    //     // ofs << Q_[s][0];
    //     for (int l = 0; l < L_; l++) {
    //         std::cout << newF_[k][l][0];
    //         for (int m = 1; m < mutFeatureDim_[l]; m++) {
    //             // ofs << "\t" << F_[k][l][m];
    //             std::cout << "\t" << newF_[k][l][m];
    //         }
    //         // ofs << "\n";
    //         std::cout << "\n";
    //     }
    // }


    double tempSum;
    diffF_ = 0;
    for (int k = 0; k < K_; k++) {
        for (int l = 0; l < L_; l++) {

            tempSum = 0;
            for (int m = 0; m < mutFeatureDim_[l]; m++) {
                tempSum = tempSum + newF_[k][l][m];
            }
            
            if (tempSum > 0) {
                for (int m = 0; m < mutFeatureDim_[l]; m++) {
                    newF_[k][l][m] = newF_[k][l][m] / tempSum;
                    diffF_ = diffF_ + (F_[k][l][m] - newF_[k][l][m]) * (F_[k][l][m] - newF_[k][l][m]);
                    F_[k][l][m] = newF_[k][l][m];
                }
            } else {
                for (int m = 0; m < mutFeatureDim_[l]; m++) {
                    newF_[k][l][m] = 1.0 /static_cast<double>(mutFeatureDim_[l]);
                    diffF_ = diffF_ + (F_[k][l][m] - newF_[k][l][m]) * (F_[k][l][m] - newF_[k][l][m]);
                    F_[k][l][m] = newF_[k][l][m];
                }
            }
        }
    }


    diffQ_ = 0;
    for (int s = 0; s < S_; s++) {

        tempSum = 0;
        for (int k = 0; k < K_; k++) {
            tempSum = tempSum + newQ_[s][k];
        }

        if (tempSum > 0) {
            for (int k = 0; k < K_; k++) {
                newQ_[s][k] = newQ_[s][k] / tempSum;
                diffQ_ = diffQ_ + (Q_[s][k] - newQ_[s][k]) * (Q_[s][k] - newQ_[s][k]);
                Q_[s][k] = newQ_[s][k];
            }
        } else {
            for (int k = 0; k < K_; k++) {
                newQ_[s][k] = 1.0 /static_cast<double>(K_); 
                diffQ_ = diffQ_ + (Q_[s][k] - newQ_[s][k]) * (Q_[s][k] - newQ_[s][k]);
                Q_[s][k] = newQ_[s][k];

            }
        }

    }    

}



std::size_t PMSig::getMutationNum () {
    return(N_);
}

double PMSig::getLikelihood() {
    return(likelihood_);
}

double PMSig::getDiffF() {
    return(diffF_);
}

double PMSig::getDiffQ() {
    return(diffQ_);
}

int PMSig::getSampleNum() {
    return(S_);
}

int PMSig::getSeqSize() {
    return(L_);
}


void PMSig::printF(const std::string& out_dir) {

    std::string out_file = out_dir + ".F.txt";
    std::ofstream ofs(out_file.c_str());
    if (ofs.fail()) {
        std::cerr << "Error: Could not open the output file.\n";
        exit(8);
    }

    for (int k = 0; k < K_; k++) {
        for (int l = 0; l < L_; l++) {
            ofs << F_[k][l][0];
            // std::cout << F_[k][l][0];
            for (int m = 1; m < mutFeatureDim_[l]; m++) {
                ofs << "\t" << F_[k][l][m];
                // std::cout << "\t" << F_[k][l][m];
            }
            ofs << "\n";
            // std::cout << "\n";
        }
    }
    ofs.close();
}

void PMSig::printQ(const std::string& out_dir) {

    std::string out_file = out_dir + ".Q.txt";
    std::ofstream ofs(out_file.c_str());
    if (ofs.fail()) {
        std::cerr << "Error: Could not open the output file.\n";
        exit(8);
    }

    for (int s = 0; s < S_; s++) {
        ofs << Q_[s][0];
        // std::cout << Q_[s][0];
        for (int k = 1; k < K_; k++) {
            ofs << "\t" << Q_[s][k];
            // std::cout << "\t" << Q_[s][k];
        }
        ofs << "\n";
        // std::cout << "\n";
    }
    ofs.close();
}


void PMSig::printTheta(const std::string& out_dir) {

    std::string out_file = out_dir + ".Theta.txt";
    std::ofstream ofs(out_file.c_str());
    if (ofs.fail()) {
        std::cerr << "Error: Could not open the output file.\n";
        exit(8);
    }

    for (int s = 0; s < S_; ++s) {
        for (int p = 0; p < P_; ++p) {
            std::cout << Theta_[s][0][p];
            for (int k = 0; k < K_; k++) {
                std::cout << "\t" << Theta_[s][k][p];
            }
            std::cout << "\n";
        }
    }

}

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <map>
#include "random.h"

#ifndef __PMSIG_H__
#define __PMSIG_H__

class PMSig {

    private:
        std::size_t N_; // # of mutation
        int S_; // # of samples
        int L_; // # of mutation features 
        int K_; // # of cluster
        int P_; // # of possible mutation patterns

        int** mutData_; // data matrix representing the number of mutation for each mutation pattern for each sample
        int* mutFeatureDim_; // dimension representing the number of possibilities for each mutation feature 

        std::size_t maxRepNum_; // # of repeats after burn-in

        double*** Theta_;
        double*** F_;
        double** Q_;
        double*** newF_;
        double** newQ_;

        int* currentDigits_; 
        double diffF_;
        double diffQ_;

        double likelihood_;

    public:
        PMSig(const std::string& file);
        virtual ~PMSig();
        void preparation(const int K);
        void EstepUpdate();
        void MstepUpdate();

        std::size_t getMutationNum();
        double getLikelihood();
        double getDiffF();
        double getDiffQ();
        int getSampleNum();
        int getSeqSize();
        void printQ(const std::string& out_dir);        
        void printF(const std::string& out_dir);
        void printTheta(const std::string& out_dir);
};


#endif


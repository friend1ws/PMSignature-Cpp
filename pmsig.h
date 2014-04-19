#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include<map>
#include "mutinfo.h"
#include "random.h"

#ifndef __PMSIG_H__
#define __PMSIG_H__

class PMSig {

    private:
        std::size_t N_; // # of mutation
        int S_; // # of samples
        int L_; // # of bases to consider
    
        int K_; // # of cluster
        int K_max_; // the limit of # of cluster

        std::vector<MutInfo> mutInfo_;    // information of depth and variant read for each mutation

        std::map<std::string, int> sampleName2Index;

        int* z_;    // cluster assignment
        int* y_;    // direction from the mutation signature


        int*** A_;
        int*** B_;
        int** C_;
        int** D_;

        std::size_t repNum_; // # of repeats after burn-in

        double*** mean_Phi_;
        double*** mean_Psi_;
        double** mean_Theta_;

        double bayesianDeviance_;
        double* bayesianDevianceList_;

        double* p_pos;
        double* p_neg;

        double param_alpha_;  /// hyper parameter
        double param_beta_; /// auxially parameter
        double param_gamma_;

    public:
        PMSig(const std::string& file);
        void preparation(const int K, const long unsigned int repNum);
        void gibbsUpdate();
        double getLogLikelihood(); 
        void updateBayesianDeviance(const int num);
        void incrementParam();
        void printA();
        void printB();
        void printC();
        void printD();
        void printMean_Phi(const std::string& out_dir);
        void printMean_Psi(const std::string& out_dir);
        void printMean_Theta(const std::string& out_dir);
        void printBayesianDeviance();
        void printPenalizedMeanBayesianDeviance(const std::string& out_dir);
        std::size_t getMutationNum ();
        int getSampleNum();
        int getSeqSize();
};

inline int base2num (char base) {

    if (base == 'A') {
        return(0);
    } else if (base == 'C') {
        return(1);
    } else if (base == 'G') {
        return(2);
    } else if (base == 'T') {
        return(3);
    } else {
        std::cerr << "wrong input for the base2num function\n";
        return(-1);
    }
}

inline int compBase2num (char base) {

    if (base == 'A') {
        return(3);
    } else if (base == 'C') {
        return(2);
    } else if (base == 'G') {
        return(1);
    } else if (base == 'T') {
        return(0);
    } else {
        std::cerr << "wrong input for the base2num function\n";
        return(-1);
    }
}

#endif


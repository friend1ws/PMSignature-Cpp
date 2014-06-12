#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <map>
#include "mutinfo.h"
#include "random.h"

#ifndef __PMSIG_H__
#define __PMSIG_H__

class PMSig {

    protected:
        std::size_t N_; // # of mutation
        int S_; // # of samples
        int L_; // # of bases to consider
        int K_; // # of cluster

        std::vector<MutInfo> mutInfo_;    // information of depth and variant read for each mutation
        std::vector<MutInfo_simple> mutInfo_simple_;
        std::map<std::string, int> sampleName2Index;

        std::size_t repNum_; // # of repeats after burn-in

        int* z_;    // cluster assignment

        double bayesianDeviance_;
        double* bayesianDevianceList_;

        double* p_pos;
        double* p_neg;

    public:
        PMSig(const std::string& file);
        virtual void preparation(const int K, const long unsigned int repNum) = 0;
        virtual void gibbsUpdate() = 0;
        virtual void updateBayesianDeviance(const int num) = 0;
        virtual void incrementParam() = 0;
        void printBayesianDeviance();
        void printPenalizedMeanBayesianDeviance(const std::string& out_dir);
        std::size_t getMutationNum ();
        int getSampleNum();
        int getSeqSize();
};


class PMSig_independent : public PMSig {

    private:
        int* y_;    // direction from the mutation signature

        int*** A_;
        int*** B_;
        int** C_;
        int** D_;

        double*** mean_Phi_;
        double*** mean_Psi_;
        double** mean_Theta_;

        double param_alpha_;  /// hyper parameter
        double param_beta_; /// auxially parameter
        double param_gamma_;

    public:
        PMSig_independent(const std::string& file);
        void printA();
        void printB();
        void printC();
        void printD();
        void printMean_Phi(const std::string& out_dir);
        void printMean_Psi(const std::string& out_dir);
        void printMean_Theta(const std::string& out_dir);
        void preparation(const int K, const long unsigned int repNum);
        void gibbsUpdate();
        double getLogLikelihood();
        void updateBayesianDeviance(const int num);
        void incrementParam();

};


class PMSig_full : public PMSig {

    private:
        int patternNum_;

        int** mutCount_assign_pattern_;
        int* mutCount_assign_;
        int** mutCount_assign_sample_;

        double** mean_Phi_;
        double** mean_Theta_;

        double param_alpha_;  /// hyper parameter
        double param_gamma_;

    public:
        PMSig_full(const std::string& file);
        void preparation(const int K, const long unsigned int repNum);
        void gibbsUpdate();
        void incrementParam();
        void updateBayesianDeviance(const int num);
        void printMean_Phi(const std::string& out_dir);
        void printMean_Theta(const std::string& out_dir);
 
};


union PMSig_union {
    PMSig_independent independent;
    PMSig_full full;
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


inline int mut2patternID (std::string sequence, char altBase, int L_) {

    int patternID = 0;

    int mbase = 1;

    if (sequence.at((L_ - 1) / 2) == 'C' or sequence.at((L_ - 1) / 2) == 'T') {

        for (int l = 0; l < L_; l++) {

            if (l == (L_ - 1) / 2) {
                continue;
            }

            if (sequence.at(l) == 'A') {
                //
            } else if (sequence.at(l) == 'C') {
                patternID = patternID + 1 * mbase;
            } else if (sequence.at(l) == 'G') {
                patternID = patternID + 2 * mbase;
            } else if (sequence.at(l) == 'T') {
                patternID = patternID + 3 * mbase;
            } else {
                std::cerr << "wrong input for the mut2patternID function\n";
                return(-1);
            }
            mbase = mbase * 4;
        }

        if (sequence.at((L_ - 1) / 2) == 'C' and altBase == 'A') {
            //
        } else if (sequence.at((L_ - 1) / 2) == 'C' and altBase == 'G') {
            patternID = patternID + 1 * mbase;
        } else if (sequence.at((L_ - 1) / 2) == 'C' and altBase == 'T') {
            patternID = patternID + 2 * mbase;
        } else if (sequence.at((L_ - 1) / 2) == 'T' and altBase == 'A') {
            patternID = patternID + 3 * mbase;
        } else if (sequence.at((L_ - 1) / 2) == 'T' and altBase == 'C') {
            patternID = patternID + 4 * mbase;
        } else if (sequence.at((L_ - 1) / 2) == 'T' and altBase == 'G') {
            patternID = patternID + 5 * mbase;
        } else {
            std::cerr << "wrong input for the mut2patternID function\n";
            return(-1);
        }

    } else {

        for (int l = L_ - 1; l >=  0; l--) {
        
            if (l == (L_ - 1) / 2) {
                continue;
            }
        
            if (sequence.at(l) == 'T') {
                //
            } else if (sequence.at(l) == 'G') {
                patternID = patternID + 1 * mbase;
            } else if (sequence.at(l) == 'C') {
                patternID = patternID + 2 * mbase;
            } else if (sequence.at(l) == 'A') {
                patternID = patternID + 3 * mbase;
            } else {
                std::cerr << "wrong input for the mut2patternID function\n";
                return(-1);
            }
            mbase = mbase * 4;
        }

        if (sequence.at((L_ - 1) / 2) == 'G' and altBase == 'T') {
            //
        }    
        else if (sequence.at((L_ - 1) / 2) == 'G' and altBase == 'C') {
            patternID = patternID + 1 * mbase;
        } else if (sequence.at((L_ - 1) / 2) == 'G' and altBase == 'A') {
            patternID = patternID + 2 * mbase;
        } else if (sequence.at((L_ - 1) / 2) == 'A' and altBase == 'T') {
            patternID = patternID + 3 * mbase;
        } else if (sequence.at((L_ - 1) / 2) == 'A' and altBase == 'G') {
            patternID = patternID + 4 * mbase;
        } else if (sequence.at((L_ - 1) / 2) == 'A' and altBase == 'C') {
            patternID = patternID + 5 * mbase;
        } else {
            std::cerr << "wrong input for the mut2patternID function\n";
            return(-1);
        }

    }


    return(patternID);
}


#endif


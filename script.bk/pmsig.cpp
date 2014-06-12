#include "pmsig.h"
#include <cmath>


PMSig::PMSig(const std::string& file) {

    std::ifstream ifs(file.c_str());
    if (ifs.fail()) {
        std::cerr << "Error: Could not open the input file.\n";
        exit(8);
    }

    std::string line;

    // Get feature size
    int curNum = 0;
    int prevSeqSize = 0;
    int curSeqSize = 0;
    int sID = -1;
    while(std::getline(ifs, line)) {

        std::istringstream is(line);
        curNum = curNum + 1;

        std::string mutationID;
        std::string sampleName;
        std::string sequence;
        char altBase;
        is >> mutationID >> sampleName >> sequence >> altBase;
    

        std::map<std::string, int>::iterator item_loc;
        item_loc = sampleName2Index.find(sampleName);

        if (item_loc == sampleName2Index.end()) {
            sID = sID + 1;
            sampleName2Index.insert(std::pair<std::string, int>(sampleName, sID));
        }
 

        prevSeqSize = curSeqSize;
        curSeqSize = sequence.size();

        if (curNum == 1 and curSeqSize % 2 == 0) {
            std::cerr << "The size of sequence should be odd number.\n";
            exit(8);
        }

        if (curNum > 1 and curSeqSize != prevSeqSize) {
            std::cerr << "At the " << curNum << " th key, " << mutationID << " " << sampleName << " " << sequence << " " << altBase << ". The size of sequence is inconsistent.\n";
            exit(8);
        }

        for (std::size_t i = 0; i < sequence.size(); i++) {
            if (sequence.at(i) != 'A' and sequence.at(i) != 'C' and sequence.at(i) != 'G' and sequence.at(i) != 'T') {
                std::cerr << "At the " << curNum << " th key, " << mutationID << " " << sampleName << " " << sequence << " " << altBase << ". The sequence should consist of 'A', 'C', 'G', or 'T'.\n";
                exit(8);
            }
        }

        if (altBase != 'A' and altBase != 'C' and altBase != 'G' and altBase != 'T') {
            std::cerr << "At the " << curNum << " th key, " << mutationID << " " << sampleName << " " << sequence << " " << altBase << ". The changed base should be 'A', 'C', 'G', or 'T'.\n";
            exit(8);
        }

        if (sequence.at((curSeqSize - 1) / 2) == altBase) {
            std::cerr << "At the " << curNum << " th key, " << mutationID << " " << sampleName << " " << sequence << " " << altBase << ". The reference base and the changed base should be different.\n";
            exit(8);
        }


        mutInfo_.push_back(MutInfo(mutationID, sID, sequence, altBase));
    }


    N_ = curNum;
    L_ = curSeqSize;
    S_ = sID + 1;

}


void PMSig::printBayesianDeviance() {
    std::cout << bayesianDeviance_ << "\n";
}

void PMSig::printPenalizedMeanBayesianDeviance(const std::string& out_dir) {

    double meanBayesianDeviance = 0;
    for (std::size_t i = 0; i < repNum_; i++) {
        meanBayesianDeviance = meanBayesianDeviance + bayesianDevianceList_[i] / static_cast<double>(repNum_);
    }

    double varianceBayesianDeviance = 0;
    for (std::size_t i = 0; i < repNum_; i++) {
        varianceBayesianDeviance = varianceBayesianDeviance + (bayesianDevianceList_[i] - meanBayesianDeviance) * (bayesianDevianceList_[i] - meanBayesianDeviance) / static_cast<double>(repNum_);
    }

    std::string out_file = out_dir + ".PMBD.txt";
    std::ofstream ofs(out_file.c_str());
    if (ofs.fail()) {
        std::cerr << "Error: Could not open the output file.\n";
        exit(8);
    }

    ofs << meanBayesianDeviance + 0.25 * varianceBayesianDeviance << "\n";

}

std::size_t PMSig::getMutationNum () {
    return(N_);
}

int PMSig::getSampleNum() {
    return(S_);
}

int PMSig::getSeqSize() {
    return(L_);
}



#include "pmsig.h"
#include <cmath>

PMSig_full::PMSig_full(const std::string& file) : PMSig(file) {
}


void PMSig_full::preparation(const int K, const long unsigned int repNum) {

    K_ = K; 
    repNum_ = repNum;
    param_alpha_ = 1; 
    param_gamma_ = 1;

    patternNum_ = pow(4, L_ - 1) * 2 * 3;

    try {

        z_ = new int[N_];
        for (std::size_t i = 0; i < N_; i++) {
            z_[i] = gen_rand(K_);;
        }

        mutCount_assign_pattern_ = new int*[K_];
        for (int k = 0; k < K_; ++k) {
            mutCount_assign_pattern_[k] = new int[patternNum_];
            for (int p = 0; p < patternNum_; p++) {
                mutCount_assign_pattern_[k][p] = 0;
            }
        }

        mutCount_assign_ = new int[K_];
        for (int k = 0; k <= K_; ++k) {
            mutCount_assign_[k] = 0;
        }

        mutCount_assign_sample_ = new int*[K_];
        for (int k = 0; k < K_; ++k) {
            mutCount_assign_sample_[k] = new int[S_];
            for (int s = 0; s < S_; s++) {
                mutCount_assign_sample_[k][s] = 0;
            }
        }

        mean_Phi_ = new double*[K_];
        for (int k = 0; k < K_; ++k) {
            mean_Phi_[k] = new double[patternNum_];
            for (int p = 0; p < patternNum_; p++) {
                mean_Phi_[k][p] = 0;
            }
        }

        mean_Theta_ = new double*[K_];
        for (int k = 0; k < K_; ++k) {
            mean_Theta_[k] = new double[S_];
            for (int s = 0; s < S_; s++) {
                mean_Theta_[k][s] = 0;
            }
        }

        bayesianDevianceList_ = new double[repNum_];
        for (std::size_t i = 0; i < repNum_; ++i) {
            bayesianDevianceList_[i] = 0;
        }

        p_pos = new double[K_];
        for (int i = 0; i < K_; ++i) {
            p_pos[i] = 0.0;
        }
        
        for (std::size_t i = 0; i < N_; i++) {

            MutInfo MI = mutInfo_[i];
            int pID = mut2patternID(MI.sequence, MI.altBase, L_);
            mutInfo_simple_.push_back(MutInfo_simple(MI.sampleID, pID));

            mutCount_assign_pattern_[z_[i]][pID] = mutCount_assign_pattern_[z_[i]][pID] + 1;
            mutCount_assign_[z_[i]] = mutCount_assign_[z_[i]] + 1;
            mutCount_assign_sample_[z_[i]][MI.sampleID] = mutCount_assign_sample_[z_[i]][MI.sampleID] + 1;

        }

    }


    catch (...) {
        std::cerr << "init(): Out of memmory" << std::endl;
        exit(EXIT_FAILURE);
    }


}


void PMSig_full::gibbsUpdate() {

    for (std::size_t i = 0; i < N_; ++i) {

    MutInfo_simple MIS = mutInfo_simple_[i];
    int assign = z_[i];

    // ********************
    // remove the i-th sample and modify the sufficient statistics
    mutCount_assign_pattern_[assign][MIS.patternVec] = mutCount_assign_pattern_[assign][MIS.patternVec] - 1;
    mutCount_assign_[assign] = mutCount_assign_[assign] - 1;
    mutCount_assign_sample_[assign][MIS.sampleID] = mutCount_assign_sample_[assign][MIS.sampleID] - 1;
    // ********************


    // ********************
    // calculate the probability of each assignment and strand
    for (int k = 0; k < K_; k++) {
        p_pos[k] = 1.0;
    }

    for (int k = 0; k < K_; k++) {
        p_pos[k] = p_pos[k] * (mutCount_assign_pattern_[k][MIS.patternVec] + param_alpha_) / (mutCount_assign_[k] + patternNum_ * param_alpha_);
    }
        
    for (int k = 0; k < K_; k++) {
        p_pos[k] = p_pos[k] * (mutCount_assign_sample_[k][MIS.sampleID] + param_gamma_);
    }
    // ********************


    // ********************
    // generate the random variable and select the next assignment and strand
    double p_pos_sum = 0;
    for (int k = 0; k < K_; k++) {
        p_pos_sum = p_pos_sum + p_pos[k];
    }

    for (int k = 0; k < K_; k++) {
        p_pos[k] = p_pos[k] / (p_pos_sum);
    }

    double u = gen_rand(1.0);
    double cum_p_pos = 0;            
    for (int k = 0; k < K_; ++k) {
        if (u > cum_p_pos) {
            assign = k;
        }
        cum_p_pos = cum_p_pos + p_pos[k];
    }
    // ********************


    // ********************
        // update the i-th status and sufficient statistics
    mutCount_assign_pattern_[assign][MIS.patternVec] = mutCount_assign_pattern_[assign][MIS.patternVec] + 1;
    mutCount_assign_[assign] = mutCount_assign_[assign] + 1;
    mutCount_assign_sample_[assign][MIS.sampleID] = mutCount_assign_sample_[assign][MIS.sampleID] + 1;

    z_[i] = assign;
        // ********************

    }

}


void PMSig_full::updateBayesianDeviance(const int num) {

/* 
    bayesianDeviance_ = 0;
    bayesianDeviance_ = bayesianDeviance_ + K_ * L_ * (lgamma(4 * param_alpha_) - 4 * lgamma(param_alpha_));

    for (int k = 0; k < K_; k++) {
        for (int l = 0; l < L_; l++) {
            bayesianDeviance_ = bayesianDeviance_ + lgamma(A_[k][l][0] + param_alpha_) + lgamma(A_[k][l][1] + param_alpha_) + lgamma(A_[k][l][2] + param_alpha_) + lgamma(A_[k][l][3] + param_alpha_);
            bayesianDeviance_ = bayesianDeviance_ - lgamma(A_[k][l][0] + A_[k][l][1] + A_[k][l][2] + A_[k][l][3] + 4 * param_alpha_);
        }
    }

    bayesianDeviance_ = bayesianDeviance_ + 4 * K_ * (lgamma(4 * param_beta_) - 4 * lgamma(param_beta_));
    for (int k = 0; k < K_; k++) {
        for (int w = 0; w < 4; w++) {
            bayesianDeviance_ = bayesianDeviance_ + lgamma(B_[k][w][0] + param_beta_) + lgamma(B_[k][w][1] + param_beta_) + lgamma(B_[k][w][2] + param_beta_) + lgamma(B_[k][w][3] + param_beta_) - lgamma(B_[k][w][w] + param_beta_);
            bayesianDeviance_ = bayesianDeviance_ - lgamma(B_[k][w][0] + B_[k][w][1] + B_[k][w][2] + B_[k][w][3] - B_[k][w][w] + 3 * param_beta_);
        }
    }

    bayesianDeviance_ = bayesianDeviance_ + K_ * S_ * (lgamma(K_ * param_gamma_) - K_ * lgamma(param_gamma_));
    for (int s = 0; s < S_; s++) {
        int tempSumC = 0;
        for (int k = 0; k < K_; k++) {
            bayesianDeviance_ = bayesianDeviance_ + lgamma(C_[k][s] + param_gamma_);
            tempSumC = tempSumC + C_[k][s];
        }
        bayesianDeviance_ = bayesianDeviance_ - lgamma(tempSumC + K_ * param_gamma_);
    }

    bayesianDevianceList_[num] = -2 * bayesianDeviance_;
*/
}



void PMSig_full::incrementParam() {

    
    for (int k = 0; k < K_; ++k) {
        double sumPhi = 0;
        for (int p = 0; p < patternNum_; p++) {
            sumPhi = sumPhi + mutCount_assign_pattern_[k][p] + param_alpha_;
        }
        for (int p = 0; p < patternNum_; p++) {
            mean_Phi_[k][p] = mean_Phi_[k][p] + static_cast<double>(mutCount_assign_pattern_[k][p] + param_alpha_) / (static_cast<double>(repNum_) * sumPhi); 
        }
    }

    for (int s = 0; s < S_; ++s) {

        double sumC = 0;
        for (int k = 0; k < K_; k++) {
            sumC = sumC + static_cast<double>(mutCount_assign_sample_[k][s] + param_gamma_);
        }    

        for (int k = 0; k < K_; k++) {
            mean_Theta_[k][s] = mean_Theta_[k][s] + static_cast<double>(mutCount_assign_sample_[k][s] + param_gamma_) / (static_cast<double>(repNum_) * sumC);
        }

    }

}


void PMSig_full::printMean_Phi(const std::string& out_dir) {

    std::string out_file = out_dir + ".Phi.txt";
    std::ofstream ofs(out_file.c_str());
    if (ofs.fail()) {
        std::cerr << "Error: Could not open the output file.\n";
        exit(8);
    }

    for (int k = 0; k < K_; k++) {

        ofs << "signature " << k;
        for (int p = 0; p < patternNum_; p++) {
            ofs << "\t" << mean_Phi_[k][p];
        }
        if (k != K_ - 1) {
            ofs << "\n";
        }
    }

    ofs.close();
}


void PMSig_full::printMean_Theta(const std::string& out_dir) {

    std::string out_file = out_dir + ".Theta.txt";
    std::ofstream ofs(out_file.c_str());
    if (ofs.fail()) {
        std::cerr << "Error: Could not open the output file.\n";
        exit(8);
    }

    for (int k = 0; k < K_; k++) {
        ofs << "sample " << k << "\n";
        ofs << mean_Theta_[k][0];
        for (int s = 1; s < S_; s++) {
            ofs << "\t" << mean_Theta_[k][s];
        }
        ofs << "\n";
    }
    ofs.close();
}


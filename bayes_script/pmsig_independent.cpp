#include "pmsig.h"
#include <cmath>


PMSig_independent::PMSig_independent(const std::string& file) : PMSig(file) {
}


void PMSig_independent::preparation(const int K, const long unsigned int repNum) {

    K_ = K; 
    repNum_ = repNum;
    param_alpha_ = 1; 
    param_beta_ = 1; 
    param_gamma_ = 1;

    try {

        z_ = new int[N_];
        for (std::size_t i = 0; i < N_; i++) {
            z_[i] = gen_rand(K_);;
        }

        y_ = new int[N_];
        for (std::size_t i = 0; i < N_; i++) {
            MutInfo MI = mutInfo_[i];
            if (MI.sequence.at((L_ - 1) / 2) == 'C' or MI.sequence.at((L_ - 1) / 2) == 'T') {
                y_[i] = 1;
            } else {
                y_[i] = 0;
            }
        }

        A_ = new int**[K_];
        for (int i = 0; i < K_; ++i) {
            A_[i] = new int*[L_];
            for (int j = 0; j < L_; ++j) {
                A_[i][j] = new int[4];
                for (int w = 0; w < 4; w++) {
                    A_[i][j][w] = 0;
                }
            }
        }

        B_ = new int**[K_];
        for (int i = 0; i < K_; ++i) {
            B_[i] = new int*[4];
            for (int j = 0; j < 4; ++j) {
                B_[i][j] = new int[4];
                for (int w = 0; w < 4; w++) {
                    B_[i][j][w] = 0;
                }
            }
        }

        C_ = new int*[K_];
        for (int i = 0; i < K_; ++i) {
            C_[i] = new int[S_];
            for (int s = 0; s < S_; s++) {
                C_[i][s] = 0;
            }
        }

        D_ = new int*[K_];
        for (int i = 0; i < K_; ++i) {
            D_[i] = new int[2];
            D_[i][0] = 0;
            D_[i][1] = 0;
        }


        mean_Phi_ = new double**[K_];
        for (int i = 0; i < K_; ++i) {
            mean_Phi_[i] = new double*[L_];
            for (int j = 0; j < L_; ++j) {
                mean_Phi_[i][j] = new double[4];
                for (int w = 0; w < 4; w++) {
                    mean_Phi_[i][j][w] = 0;
                }
            }
        }

        mean_Psi_ = new double**[K_];
        for (int i = 0; i < K_; ++i) {
            mean_Psi_[i] = new double*[4];
            for (int j = 0; j < 4; ++j) {
                mean_Psi_[i][j] = new double[4];
                for (int w = 0; w < 4; w++) {
                    mean_Psi_[i][j][w] = 0;
                }
            }
        }

        mean_Theta_ = new double*[K_];
        for (int i = 0; i < K_; ++i) {
            mean_Theta_[i] = new double[S_];
            for (int s = 0; s < S_; s++) {
                mean_Theta_[i][s] = 0;
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
        p_neg = new double[K_];
        for (int i = 0; i < K_; ++i) {
            p_neg[i] = 0.0;
        }
        
        for (std::size_t i = 0; i < N_; i++) {

            MutInfo MI = mutInfo_[i];
            
            if (y_[i] == 1) {

                for (int l = 0; l < L_; l++) {
                    A_[z_[i]][l][base2num(MI.sequence.at(l))] = A_[z_[i]][l][base2num(MI.sequence.at(l))] + 1;
                }
                B_[z_[i]][base2num(MI.sequence.at((L_ - 1) / 2))][base2num(MI.altBase)] = B_[z_[i]][base2num(MI.sequence.at((L_ - 1) / 2))][base2num(MI.altBase)] + 1;
                D_[z_[i]][1] = D_[z_[i]][1] + 1;
 
            } else {

                for (int l = 0; l < L_; l++) {
                    A_[z_[i]][L_ - l - 1][compBase2num(MI.sequence.at(l))] = A_[z_[i]][L_ - l - 1][compBase2num(MI.sequence.at(l))] + 1;
                }
                B_[z_[i]][compBase2num(MI.sequence.at((L_ - 1) / 2))][compBase2num(MI.altBase)] = B_[z_[i]][compBase2num(MI.sequence.at((L_ - 1) / 2))][compBase2num(MI.altBase)] + 1;
                D_[z_[i]][0] = D_[z_[i]][0] + 1;
            }

            C_[z_[i]][MI.sampleID] = C_[z_[i]][MI.sampleID] + 1;

        }

    }


    catch (...) {
        std::cerr << "init(): Out of memmory" << std::endl;
        exit(EXIT_FAILURE);
    }


}


void PMSig_independent::gibbsUpdate() {

    for (std::size_t i = 0; i < N_; ++i) {

    MutInfo MI = mutInfo_[i];
    int assign = z_[i];
    int strand = y_[i];

    // ********************
    // remove the i-th sample and modify the sufficient statistics
    if (strand == 1) {

    for (int l = 0; l < L_; l++) {
        A_[assign][l][base2num(MI.sequence.at(l))] = A_[assign][l][base2num(MI.sequence.at(l))] - 1;
    }
        B_[assign][base2num(MI.sequence.at((L_ - 1) / 2))][base2num(MI.altBase)] = B_[assign][base2num(MI.sequence.at((L_ - 1) / 2))][base2num(MI.altBase)] - 1;

    } else {

    for (int l = 0; l < L_; l++) {
        A_[assign][L_ - l - 1][compBase2num(MI.sequence.at(l))] = A_[assign][L_ - l - 1][compBase2num(MI.sequence.at(l))] - 1;
    }
        B_[assign][compBase2num(MI.sequence.at((L_ - 1) / 2))][compBase2num(MI.altBase)] = B_[assign][compBase2num(MI.sequence.at((L_ - 1) / 2))][compBase2num(MI.altBase)] - 1;

    }
    C_[assign][MI.sampleID] = C_[assign][MI.sampleID] - 1;
    D_[assign][strand] = D_[assign][strand] - 1; 
    // ********************


    // ********************
    // calculate the probability of each assignment and strand
    for (int k = 0; k < K_; k++) {
        p_pos[k] = 1.0;
        p_neg[k] = 1.0;
    }

    for (int l = 0; l < L_; l++) {
        for (int k = 0; k < K_; k++) {
            p_pos[k] = p_pos[k] * (A_[k][l][base2num(MI.sequence.at(l))] + param_alpha_) / (A_[k][l][0] + A_[k][l][1] + A_[k][l][2] + A_[k][l][3] + 4 * param_alpha_);
            p_neg[k] = p_neg[k] * (A_[k][L_ - l - 1][compBase2num(MI.sequence.at(l))] + param_alpha_) / (A_[k][L_ - l - 1][0] + A_[k][L_ - l - 1][1] + A_[k][L_ - l - 1][2] + A_[k][L_ - l - 1][3] + 4 * param_alpha_);
        }
    }
        
    for (int k = 0; k < K_; k++) {
        p_pos[k] = p_pos[k] * (B_[k][base2num(MI.sequence.at((L_ - 1) / 2))][base2num(MI.altBase)] + param_beta_) / (B_[k][base2num(MI.sequence.at((L_ - 1) / 2))][0] + B_[k][base2num(MI.sequence.at((L_ - 1) / 2))][1] + B_[k][base2num(MI.sequence.at((L_ - 1) / 2))][2] + B_[k][base2num(MI.sequence.at((L_ - 1) / 2))][3] + 3 * param_beta_);
        p_neg[k] = p_neg[k] * (B_[k][compBase2num(MI.sequence.at((L_ - 1) / 2))][compBase2num(MI.altBase)] + param_beta_) / (B_[k][compBase2num(MI.sequence.at((L_ - 1) / 2))][0] + B_[k][compBase2num(MI.sequence.at((L_ - 1) / 2))][1] + B_[k][compBase2num(MI.sequence.at((L_ - 1) / 2))][2] + B_[k][compBase2num(MI.sequence.at((L_ - 1) / 2))][3] + 3 * param_beta_);
    }

    for (int k = 0; k < K_; k++) {
        p_pos[k] = p_pos[k] * (C_[k][MI.sampleID] + param_gamma_);
        p_neg[k] = p_neg[k] * (C_[k][MI.sampleID] + param_gamma_);
    }
    // ********************


    // ********************
    // generate the random variable and select the next assignment and strand
    double p_pos_sum = 0;
    double p_neg_sum = 0;
    for (int k = 0; k < K_; k++) {
        p_pos_sum = p_pos_sum + p_pos[k];
        p_neg_sum = p_neg_sum + p_neg[k];
    }

    for (int k = 0; k < K_; k++) {
        p_pos[k] = p_pos[k] / (p_pos_sum + p_neg_sum);
        p_neg[k] = p_neg[k] / (p_pos_sum + p_neg_sum);
    }

    p_pos_sum = p_pos_sum / (p_pos_sum + p_neg_sum);
    p_neg_sum = p_neg_sum / (p_pos_sum + p_neg_sum);


    double u = gen_rand(1.0);

    if (u > p_pos_sum) {
    
        strand = 0;
        u = u - p_pos_sum;
        double cum_p_neg = 0;            
        for (int k = 0; k < K_; ++k) {
            if (u > cum_p_neg) {
                assign = k;
            }
            cum_p_neg = cum_p_neg + p_neg[k];
        }

    } else {

        strand = 1;
        double cum_p_pos = 0;
        for (int k = 0; k < K_; ++k) {
            if (u > cum_p_pos) {
                assign = k;
            }
            cum_p_pos = cum_p_pos + p_pos[k];
        }
        
    }
    // ********************


    // ********************
        // update the i-th status and sufficient statistics
    if (strand == 1) {

        for (int l = 0; l < L_; l++) {
            A_[assign][l][base2num(MI.sequence.at(l))] = A_[assign][l][base2num(MI.sequence.at(l))] + 1;
        }
        B_[assign][base2num(MI.sequence.at((L_ - 1) / 2))][base2num(MI.altBase)] = B_[assign][base2num(MI.sequence.at((L_ - 1) / 2))][base2num(MI.altBase)] + 1;

    } else {

        for (int l = 0; l < L_; l++) {
            A_[assign][L_ - l - 1][compBase2num(MI.sequence.at(l))] = A_[assign][L_ - l - 1][compBase2num(MI.sequence.at(l))] + 1;
        }
        B_[assign][compBase2num(MI.sequence.at((L_ - 1) / 2))][compBase2num(MI.altBase)] = B_[assign][compBase2num(MI.sequence.at((L_ - 1) / 2))][compBase2num(MI.altBase)] + 1;

    }

    C_[assign][MI.sampleID] = C_[assign][MI.sampleID] + 1;
    D_[assign][strand] = D_[assign][strand] + 1;

    y_[i] = strand;
    z_[i] = assign;
        // ********************

    }

}


void PMSig_independent::updateBayesianDeviance(const int num) {

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
}


double PMSig_independent::getLogLikelihood() {

    double logLikelihood = 0;
    logLikelihood = logLikelihood + K_ * L_ * (lgamma(4 * param_alpha_) - 4 * lgamma(param_alpha_));

    for (int k = 0; k < K_; k++) {
        for (int l = 0; l < L_; l++) {
            logLikelihood = logLikelihood + lgamma(A_[k][l][0] + param_alpha_) + lgamma(A_[k][l][1] + param_alpha_) + lgamma(A_[k][l][2] + param_alpha_) + lgamma(A_[k][l][3] + param_alpha_);
            logLikelihood = logLikelihood - lgamma(A_[k][l][0] + A_[k][l][1] + A_[k][l][2] + A_[k][l][3] + 4 * param_alpha_);
        }
    }

    logLikelihood = logLikelihood + 4 * K_ * (lgamma(4 * param_beta_) - 4 * lgamma(param_beta_));
    for (int k = 0; k < K_; k++) {
        for (int w = 0; w < 4; w++) {
            logLikelihood = logLikelihood + lgamma(B_[k][w][0] + param_beta_) + lgamma(B_[k][w][1] + param_beta_) + lgamma(B_[k][w][2] + param_beta_) + lgamma(B_[k][w][3] + param_beta_) - lgamma(B_[k][w][w] + param_beta_);
            logLikelihood = logLikelihood - lgamma(B_[k][w][0] + B_[k][w][1] + B_[k][w][2] + B_[k][w][3] - B_[k][w][w] + 3 * param_beta_);
        }
    }

    logLikelihood = logLikelihood + K_ * S_ * (lgamma(K_ * param_gamma_) - K_ * lgamma(param_gamma_));
    for (int s = 0; s < S_; s++) {
        int tempSumC = 0;
        for (int k = 0; k < K_; k++) {
            logLikelihood = logLikelihood + lgamma(C_[k][s] + param_gamma_);
            tempSumC = tempSumC + C_[k][s];
        }
        logLikelihood = logLikelihood - lgamma(tempSumC + K_ * param_gamma_);
    }

    return(logLikelihood);
}


void PMSig_independent::incrementParam() {

    
    for (int k = 0; k < K_; ++k) {
        for (int l = 0; l < L_; ++l) {
            double sumA = A_[k][l][0] + A_[k][l][1] + A_[k][l][2] + A_[k][l][3] + 4 * param_alpha_;
            for (int w = 0; w < 4; w++) {
                mean_Phi_[k][l][w] = mean_Phi_[k][l][w] + static_cast<double>(A_[k][l][w] + param_alpha_) / (static_cast<double>(repNum_) * sumA);
            }
        }
    }

    for (int k = 0; k < K_; ++k) {
        for (int j = 0; j < 4; ++j) {
            double sumB = B_[k][j][0] + B_[k][j][1] + B_[k][j][2] + B_[k][j][3] + 3 * param_beta_;
            for (int w = 0; w < 4; w++) {
                mean_Psi_[k][j][w] = mean_Psi_[k][j][w] + static_cast<double>(B_[k][j][w] + param_beta_) / (static_cast<double>(repNum_) * sumB);
            }
            mean_Psi_[k][j][j] = 0;
        }
    }

    for (int s = 0; s < S_; ++s) {

        double sumC = 0;
        for (int k = 0; k < K_; k++) {
            sumC = sumC + static_cast<double>(C_[k][s] + param_gamma_);
        }    

        for (int k = 0; k < K_; k++) {
            mean_Theta_[k][s] = mean_Theta_[k][s] + static_cast<double>(C_[k][s] + param_gamma_) / (static_cast<double>(repNum_) * sumC);
        }

    }

}


void PMSig_independent::printA() {

    for (int k = 0; k < K_; k++) {
        std::cout << "A_[" << k << "][][]" << "\n";
        for (int l = 0; l < L_; l++) {
            double sumA = A_[k][l][0] + A_[k][l][1] + A_[k][l][2] + A_[k][l][3] + param_alpha_;
            std::cout << static_cast<double>(A_[k][l][0]) / sumA << "\t" << static_cast<double>(A_[k][l][1]) / sumA << "\t" << static_cast<double>(A_[k][l][2]) / sumA << "\t" << static_cast<double>(A_[k][l][3]) / sumA << "\n";
        }
        std::cout << "\n";
    }

}

void PMSig_independent::printMean_Phi(const std::string& out_dir) {

    std::string out_file = out_dir + ".Phi.txt";
    std::ofstream ofs(out_file.c_str());
    if (ofs.fail()) {
        std::cerr << "Error: Could not open the output file.\n";
        exit(8);
    }

    for (int k = 0; k < K_; k++) {

        ofs << "signature " << k << "\n";
        for (int l = 0; l < L_; l++) {
            ofs << mean_Phi_[k][l][0] << "\t" << mean_Phi_[k][l][1] << "\t" << mean_Phi_[k][l][2] << "\t" << mean_Phi_[k][l][3] << "\n";
        }
        if (k != K_ - 1) {
            ofs << "\n";
        }
    }

    ofs.close();
}

void PMSig_independent::printB() {

    for (int k = 0; k < K_; k++) {
        std::cout << "B_[" << k << "][][]" << "\n";
        for (int w = 0; w < 4; w++) {
            std::cout << B_[k][w][0] << "\t" << B_[k][w][1] << "\t" << B_[k][w][2] << "\t" << B_[k][w][3] << "\n";
        }
        std::cout << "\n";
    }

}

void PMSig_independent::printMean_Psi(const std::string& out_dir) {

    std::string out_file = out_dir + ".Psi.txt";
    std::ofstream ofs(out_file.c_str());
    if (ofs.fail()) {
        std::cerr << "Error: Could not open the output file.\n";
        exit(8);
    }


    for (int k = 0; k < K_; k++) {

        // ofs << "mean_Psi_[" << k << "][][]" << "\n";
        ofs << "signature " << k << "\n";
        for (int w = 0; w < 4; w++) {
            ofs << mean_Psi_[k][w][0] << "\t" << mean_Psi_[k][w][1] << "\t" << mean_Psi_[k][w][2] << "\t" << mean_Psi_[k][w][3] << "\n";
        }
        if (k != K_ - 1) {
            ofs << "\n";
        }
    }

    ofs.close();
}

void PMSig_independent::printC() {

    for (int k = 0; k < K_; k++) {
        std::cout << "C_[" << k << "][]" << "\n";
        std::cout << C_[k][0]; 
        for (int s = 1; s < S_; s++) {
            std::cout << "\t" << C_[k][s];
        }
        std::cout << "\n";
    }
    std::cout << "\n";

}

void PMSig_independent::printMean_Theta(const std::string& out_dir) {

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


void PMSig_independent::printD() {

    for (int k = 0; k < K_; k++) {
        std::cout << "D_[" << k << "][]" << "\n";
        std::cout << D_[k][0] << "\t" << D_[k][1] << "\n";
    }

}


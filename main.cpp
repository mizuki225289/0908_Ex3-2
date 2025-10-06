#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

void solve(int N, std::vector<double> a, std::vector<double> b, std::vector<double> c) {
    //Step7 ガウスの消去法
    std::vector<double> new_a(N);
    std::vector<double> new_b(N);
    std::vector<double> phi(N);
    
    new_a.at(0) = a.at(0);
    new_b.at(0) = b.at(0);
    for(int i=1; i < N; i++) {
        new_a.at(i) = a.at(i) - c.at(i-1) * c.at(i-1) / new_a.at(i-1);
    }
    for(int i=1; i < N; i++) {
        new_b.at(i) = b.at(i) - c.at(i-1) * new_b.at(i-1) / new_a.at(i-1);
    }

    phi.at(N-1) = new_b.at(N-1) / new_a.at(N-1);
    for(int i=N-2; i >= 0; i--) {
        phi.at(i) = (new_b.at(i) - c.at(i) * phi.at(i+1)) / new_a.at(i);
    }

    // for(int i=0; i < N; i++) {
    //     std::cout << "new_a " << new_a.at(i) <<std::endl;
    // }
    // for(int i=0; i < N; i++) {
    //     std::cout << "new_b " << new_b.at(i) <<std::endl;
    // }
    // for(int i=0; i < N; i++) {
    //     std::cout << "phi " << phi.at(i) <<std::endl;
    // }

    //結果ファイル出力
    std::string name;
    std::cout << "file name > \n";
    std::cin >> name;
    std::ofstream ofs(name);
    for(int i=0; i < N; i++) {
        ofs << 1.0 * i / (N-1) << " " << phi.at(i) << std::endl;
    }
    ofs.close();
}

int main(void) {
    // int M; //エレメントの個数
    // std::cout << "M > ";
    // std::cin >> M;
    // int N = M+1;
    // //double l = 1.0 / M; //エレメントのサイズ（今回は全部一定という想定）

    //ファイル読み込み
    std::string filename, line, term;
    std::cout << "input file > ";
    std::cin >> filename;
    std::ifstream file(filename);

    std::vector<double> alpha, beta, f, length;

    int li = 0, lj = 0;
    if(file.is_open()) {
        while(std::getline(file, line)) {
            //std::cout << line << std::endl;
            std::stringstream ss{line}; //入出力可能なsstreamに変換
            lj = 0;
            while(std::getline(ss, term, ',')) {
                if(lj == 1) {
                    alpha.push_back(std::stod(term));
                }
                if(lj == 2) {
                    beta.push_back(std::stod(term));
                }
                if(lj == 3) {
                    f.push_back(std::stod(term));
                }
                if(lj == 4) {
                    length.push_back(std::stod(term));
                }
                lj++;
            }
            li++;
        }
    } else {
        std::cout << "Error in file open" << std::endl;
        exit(1);
    }
    
    // if(M != std::size(alpha)) {
    //     std::cout << "Error in M size" << std::endl;
    //     exit(1);
    // }
    int M = std::size(alpha);
    int N = M+1;

    //境界条件設定
    int b_type1, b_type2; //境界条件のタイプ
    double p1, q1, gamma1, p2, q2, gamma2; //境界条件関係の定数
    std::cout << "Boundary Condition at x = 0 >" << std::endl;
    std::cout << "type 1 or 3" << std::endl;
    std::cin >> b_type1;
    if(b_type1 == 1) {
        std::cout << "p > " << std::endl;
        std::cin >> p1;
    } else if(b_type1 == 3) {
        std::cout << "q > " << std::endl;
        std::cin >> q1;
        std::cout << "gamma > " << std::endl;
        std::cin >> gamma1;
    } else {
        std::cout << "Error!!" << std::endl;
        exit(1);
    }
    std::cout << "Boundary Condition at x = L >" << std::endl;
    std::cout << "type 1 or 3" << std::endl;
    std::cin >> b_type2;
    if(b_type2 == 1) {
        std::cout << "p > " << std::endl;
        std::cin >> p2;
    } else if(b_type2 == 3) {
        std::cout << "q > " << std::endl;
        std::cin >> q2;
        std::cout << "gamma > " << std::endl;
        std::cin >> gamma2;
    } else {
        std::cout << "Error!!" << std::endl;
        exit(1);
    }

    //データ出力
    for(int i=0; i < M; i++) {
        std::cout << i << " : alpha = " << alpha[i] 
                       << ", beta = " << beta[i] 
                       << ", f = " << f[i] 
                       << ", l = " << length[i] << std::endl;
    }

    //K_iiなど生成
    std::vector<double> K_diag(N), K_offdiag(N-1);
    K_diag.at(0) = (alpha.at(0) / length.at(0)) + (beta.at(0) * length.at(0) / 3.0);
    K_diag.at(N-1) = (alpha.at(M-1) / length.at(M-1)) + (beta.at(M-1) * length.at(M-1) / 3.0);
    for(int i=2; i <= N-1; i++) {
        K_diag.at(i-1) = (alpha.at(i-2) / length.at(i-2) + beta.at(i-2) * length.at(i-2) / 3.0) 
                       + (alpha.at(i-1) / length.at(i-1) + beta.at(i-1) * length.at(i-1) / 3.0);
    }
    for(int i=1; i <= N-1; i++) {
        K_offdiag.at(i-1) = -alpha.at(i-1) / length.at(i-1) + beta.at(i-1) * length.at(i-1) / 6.0;
    }
    std::cout << "K_ii OK\n";

    //b生成
    std::vector<double> b(N);
    b.at(0) = f.at(0) * length.at(0) / 2.0;
    b.at(N-1) = f.at(M-1) * length.at(M-1) / 2.0;
    for(int i=2; i <= N-1; i++) {
        b.at(i-1) = f.at(i-2) * length.at(i-2) / 2.0 + f.at(i-1) * length.at(i-1) / 2.0;
    }
    std::cout << "b_i OK\n";

    //境界条件設定
    if(b_type1 == 1) {
        K_diag.at(0) = 1.0;
        b.at(0) = p1;
        b.at(1) = b.at(1) - K_offdiag.at(0) * p1;
        K_offdiag.at(0) = 0;
    } else if(b_type1 == 3) {
        K_diag.at(0) = K_diag.at(0) + gamma1;
        b.at(0) = b.at(0) + q1;
    }
    if(b_type2 == 1) {
        K_diag.at(N-1) = 1.0;
        b.at(N-1) = p2;
        b.at(N-2) = b.at(N-2) - K_offdiag.at(N-3) * p2;
        K_offdiag.at(N-2) = 0;
    } else if(b_type2 == 3) {
        K_diag.at(N-1) = K_diag.at(N-1) + gamma2;
        b.at(N-1) = b.at(N-1) + q2;
    }
    std::cout << "boundary OK\n";

    // for(int i=0; i < K_diag.size(); i++) {
    //     std::cout << "K_diag : " << i << " " << K_diag.at(i) << std::endl;
    // }
    // for(int i=0; i < K_offdiag.size(); i++) {
    //     std::cout << "K_offdiag : " << i << " " << K_offdiag.at(i) << std::endl;
    // }
    // for(int i=0; i < b.size(); i++) {
    //     std::cout << "b : " << i << " " << b.at(i) << std::endl;
    // }

    //解く
    solve(N, K_diag, b, K_offdiag);

    // std::ofstream ofs("result.dat");
    // for(int i=0; i < N; i++) {
    //     ofs << i / double(M) << " " << phi.at(i) << std::endl;
    // }
    // ofs.close();

    return 0;

}
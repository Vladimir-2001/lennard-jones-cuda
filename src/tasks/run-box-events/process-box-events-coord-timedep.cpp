#include "NumberStatistics.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <cmath> 

using namespace std;
using namespace SampleMoments;


const int tabsize = 15;

//std::string filenames_file = "run.N400.ust-0.0458.rhost0.3.seed1.filenames.dat";
std::vector<std::string> filenames_file = { "run.N400.ust-0.0458.rhost0.3.seed1.filenames.dat",
"run.N400.ust-0.0458.rhost0.3.seed2.filenames.dat",
"run.N400.ust-0.0458.rhost0.3.seed3.filenames.dat",
"run.N400.ust-0.0458.rhost0.3.seed4.filenames.dat",
"run.N400.ust-0.0458.rhost0.3.seed5.filenames.dat",
"run.N400.ust-0.0458.rhost0.3.seed6.filenames.dat",
"run.N400.ust-0.0458.rhost0.3.seed7.filenames.dat",
"run.N400.ust-0.0458.rhost0.3.seed8.filenames.dat",
"run.N400.ust-0.0458.rhost0.3.seed9.filenames.dat",
"run.N400.ust-0.0458.rhost0.3.seed10.filenames.dat",
"run.N400.ust-0.0458.rhost0.3.seed11.filenames.dat",
"run.N400.ust-0.0458.rhost0.3.seed12.filenames.dat",
"run.N400.ust-0.0458.rhost0.3.seed13.filenames.dat",
"run.N400.ust-0.0458.rhost0.3.seed14.filenames.dat" };
double dalpha = 0.02;
int Nalphas = 50;
double L = 1., dL = 0.1;
double t = 50.;
float Tst = 1.4;

double Tch = 0.150; // GeV
const double  mN = 0.938; // GeV

double ycm_from_ecm(double ecm) {
    return log((ecm + sqrt(ecm * ecm - 4. * mN * mN)) / 2. / mN);
}


double alpha_i(int i) {
    return dalpha * i;
}

int main(int argc, char* argv[]) {

 //   if (argc > 1) {
//        filenames_file = string(argv[1]);
 //   }

    if (argc > 2) {
        dalpha = atof(argv[2]);
    }

    Nalphas = round(1. / dalpha) + 1;
    if ((Nalphas - 1) * dalpha > 1.)
        Nalphas--;


    double ecm = 2. * mN;
    if (argc > 3) {
        ecm = atof(argv[3]);
    }
    double ycm = ycm_from_ecm(ecm);

    vector<double> times;
    vector<string> filenames;
    // Read the filenames for each time moment
    {
        ifstream fin(filenames_file[0]);
        if (fin.is_open()) {
            double t;
            string file;
            times.push_back(0);
            filenames.push_back("run.N400.ust-0.0458.rhost0.3.seed1.t0.000000.bin");
            while (fin >> t >> file) {
                times.push_back(t);
                filenames.push_back(file);
            }
            fin.close();
        }
        else {
            cerr << "Error: could not open file " << filenames_file[0] << endl;
            return 1;
        }
    }

    //ofstream fout(filenames_file + ".fluctuations-time-dep.dat");
    ofstream fout("factorialalphadep_randomprotons_rho03.dat");
    if (!fout.is_open()) {
        cerr << "Error: could not open output file " << filenames_file[0] << ".fluctuations-time-dep.dat" << endl;
        return 1;
    }

  //  fout << Nalphas << " # Number of subvolumes" << endl;
  //  fout << endl;

    cout << setw(tabsize) << "t" << " "
        << setw(tabsize) << "nevents" << " "
        << setw(tabsize) << "alpha" << " "
        << setw(tabsize) << "mean" << " "
        << setw(tabsize) << "error" << " "
        << setw(tabsize) << "wtil" << " "
        << setw(tabsize) << "error" << " "
        << setw(tabsize) << "Ss" << " "
        << setw(tabsize) << "s_error" << " "
        << setw(tabsize) << "Ks" << " "
        << setw(tabsize) << "k_error" << " ";
    cout << endl;


    float ycut[100] = {00};

    { 

        for (int ifile = 0; ifile < filenames.size(); ifile++) {
            long int ind1 = 0;
            double t = times[ifile];

            fout << t << " # Time" << endl;
            cout << t << " # Time" << endl;

            vector<NumberStatistics> stats(14 * Nalphas);

            //Mean shift for alpha*N
            for (int i=0;i<Nalphas;i++)
                stats[i].SetMeanShift(200 * i / Nalphas);

            long long nevents = 0;
            double rho = 0.1;
            int N = 400;

            for (int j = 0; j < 14; j++) {
                cout << "seed=" << j + 1 << endl;

                vector<string> filenames2;

                ifstream fin0(filenames_file[j]);
                if (fin0.is_open()) {
                    string file;
                    filenames2.push_back("run.N400.ust-0.0458.rhost0.3.seed"+to_string(j+1)+".t0.000000.bin");
                    while (fin0 >> t >> file) {

                        filenames2.push_back(file);
                    }
                    fin0.close();
                }
                else {
                    cerr << "Error: could not open file " << filenames_file[j] << endl;
                    return 1;
                }

                int flag[400] = {};
                const string& filename = filenames2[ifile];

                ifstream fin1(filename, std::ios::binary);
                string in1, in2;
                int in_int;
                double in_double;
                if (fin1.is_open()) {
                    float buff = 0;
                    int ibuff = 0;
                    int k = 0;
                    int i = 0;
                    int r = 0;
                    int s = 0;
                    int l1 = 0;
                    while (!fin1.eof()) {

                        //Random particles in acceptance
                        for (int num = 0; num < 400; num++)
                        {
                            if (num == 0)
                            {
                                srand(nevents);
                            }
                               
                            flag[num] = rand() % 2;
                        }
                        if ((k >= 2) & (k < 6))
                        {
                            k++;
                            fin1.read(reinterpret_cast<char*>(&buff), sizeof(float));
                            k++;
                            fin1.read(reinterpret_cast<char*>(&buff), sizeof(float));
                            if (k == 4)
                                Tst = buff;
                            k++;
                            fin1.read(reinterpret_cast<char*>(&buff), sizeof(float));
                            if (k == 5)
                            {
                                rho = buff;
                                L = pow(N / rho, 1. / 3.);
                            }

                        }
                        else
                        {
                            if (k < 3)
                            {
                                if (k == 0)
                                {
                                    k++;
                                    fin1.read(reinterpret_cast<char*>(&ibuff), sizeof(int));
                                    if ((nevents < ibuff) || ((j > 0) & (r < ibuff)))
                                    {
                                        nevents++;
                                        r = ibuff;
                                    }

                                }

                                if (k == 1)
                                {
                                    k++;
                                    fin1.read(reinterpret_cast<char*>(&ibuff), sizeof(int));
                                    N = ibuff;
                                }
                            }
                        }
                        if (k == 5)
                        {
                            k++;
                            fin1.read(reinterpret_cast<char*>(&ibuff), sizeof(int));
                            if (k == 6)
                            {
                                t = ibuff;
                             //   cout << "t = " << ibuff << endl;
                            }

                        }


                        // Initialize counters
                        vector<int> cnts(Nalphas);

                        //t=0 counters
                        vector<int> cnts0(Nalphas);

                           
                        // Process particles
                        double x, y, z, vx, vy, vz;

                        if (k > 5) {
                            for (int iN = 0; iN < N; iN++) {
                                k++;
                                fin1.read(reinterpret_cast<char*>(&buff), sizeof(float));

                                if (k == 7 + 6 * iN)
                                    x = buff;

                                k++;
                                fin1.read(reinterpret_cast<char*>(&buff), sizeof(float));

                                if (k == 8 + 6 * iN)
                                    y = buff;

                                k++;
                                fin1.read(reinterpret_cast<char*>(&buff), sizeof(float));
                                if (k == 9 + 6 * iN)
                                {
                                    z = buff;

                               //     int indz = static_cast<int>((z / L) / dalpha) + 1;
                                    if (flag[iN] == 1)
                                        for (int s = 0; s < Nalphas; s++) {
                                            if (abs(z / L) <= float(s) / Nalphas)
                                            {   
                                                cnts[s]++;
                                            }
                                        }
                                }

                                k++;
                                fin1.read(reinterpret_cast<char*>(&buff), sizeof(float));

                                if (k == 10 + 6 * iN)
                                    vx = buff;

                                k++;
                                fin1.read(reinterpret_cast<char*>(&buff), sizeof(float));

                                if (k == 11 + 6 * iN)
                                    vy = buff;

                                k++;
                                fin1.read(reinterpret_cast<char*>(&buff), sizeof(float));

                                if (k == 12 + 6 * iN)
                                {
                                    vz = buff;
                                }
                            }
                        }
                        // Compute the prefix sums
                  //      for (int i = 1; i < Nalphas; i++)
                  //          cnts[i] += cnts[i - 1];

                        // Add the counts to the statistics
                        if (cnts[1] > 0)
                            for (int i = 0; i < Nalphas; i++)
                            {

                                stats[i].AddObservation(cnts[i]);

                            }

                        if (k >= 2406)
                        {
                            k = 0;
                        }


                    }
                }
                else {
                    cout << "Cannot open file " << filename << " j=" << j+1 << endl;
                    return 1;
                }
            
            
            fout << nevents << " # Number of events" << endl;
            fout << setw(tabsize) << "alpha" << " "
                << setw(tabsize) << "mean" << " "
                << setw(tabsize) << "error" << " "
                << setw(tabsize) << "F2/F1" << " "
                << setw(tabsize) << "errorF21" << " "
                << setw(tabsize) << "F3/F1" << " "
                << setw(tabsize) << "errorF31" << " "
                << setw(tabsize) << "F4/F1" << " "
                << setw(tabsize) << "errorF41" << " "
                << setw(tabsize) << "F2/F1^2" << " "
                << setw(tabsize) << "errorC21" << " "
                << setw(tabsize) << "F3/F1^3" << " "
                << setw(tabsize) << "errorC31" << " "
                << setw(tabsize) << "F4/F1^4" << " "
                << setw(tabsize) << "errorC41" << " ";
            fout << endl;

            for (int i = 0; i < Nalphas; i++) {
            double alpha = stats[i].GetMean()/ stats[Nalphas-1].GetMean();
                fout << setw(tabsize) << alpha << " "
                    << setw(tabsize) << stats[i].GetMean() << " "
                    << setw(tabsize) << stats[i].GetMeanError() << " "
                    << setw(tabsize) << stats[i].GetFactorialCumulantRatio(2, 1) << " "
                    << setw(tabsize) << stats[i].GetFactorialCumulantRatioError(2, 1) << " "
                    << setw(tabsize) << stats[i].GetFactorialCumulantRatio(3, 1) << " "
                    << setw(tabsize) << stats[i].GetFactorialCumulantRatioError(3, 1) << " "
                    << setw(tabsize) << stats[i].GetFactorialCumulantRatio(4, 1) << " "
                    << setw(tabsize) << stats[i].GetFactorialCumulantRatioError(4, 1) << " "
                    << setw(tabsize) << stats[i].GetReducedFactorialCumulant(2) << " "
                    << setw(tabsize) << stats[i].GetReducedFactorialCumulantError(2) << " "
                    << setw(tabsize) << stats[i].GetReducedFactorialCumulant(3) << " "
                    << setw(tabsize) << stats[i].GetReducedFactorialCumulantError(3) << " "
                    << setw(tabsize) << stats[i].GetReducedFactorialCumulant(4) << " "
                    << setw(tabsize) << stats[i].GetReducedFactorialCumulantError(4) << " ";
                fout << endl;

            }
            fout << endl;
           
        }
    }
    fout.close();

    return 0;
}



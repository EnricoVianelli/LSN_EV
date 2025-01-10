#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../RndGen/random.h"
#include "GA.h"

using namespace std;

int main(int argc, char *argv[])
{
    int N_city;  // number of cities
    int DIM_POP; // number of travels
    int N_GEN;   // number of generations
    int SIM;     // type of simulation : 0 cir 1 sq
    string property;
    string line;

    array<string, 2> simu_type = {"Circumference", "Square"};

    if (argc != 2 || atoi(argv[1]) < 0 || atoi(argv[1]) > 1)
    {

        if (argc != 2)
            cerr << "Usage: " << argv[0] << "  <simulation type>" << endl;
        else
            cerr << " Invalid simulation type !" << endl;

        cerr << " Available type:" << endl;

        for (int i = 0; i < (int)std::size(simu_type); i++)
            cerr << "\t" << i << " = " << simu_type[i] << endl;

        cerr << endl;
        return 1;
    }

    SIM = atoi(argv[1]);
    ifstream inputdat("input.dat");
    while (getline(inputdat, line))
    {
        stringstream ss(line); // Usa stringstream per analizzare la riga
        ss >> property;        // Leggi la prima parola della riga

        if (property == "N_CITY:")
            ss >> N_city;
        else if (property == "DIM_POP:")

            ss >> DIM_POP;

        else if (property == "N_GEN:")

            ss >> N_GEN;

        // else if (property == "SIM_TYPE:")

        //    ss >> SIM;

        else
            cerr << line << " unknown input parameter" << endl;
    }
    inputdat.close();
    ////////////////////////////////////////////////////
    std::vector<city> cities;
    Random rnd;
    rnd.initialize();
    // parameters
    string directory{};
    for (int i = 0; i < N_city; i++)
    {
        if (SIM)
            cities.emplace_back(city{rnd.Rannyu(-1, 1), rnd.Rannyu(-1, 1)});
        else
        {
            double theta = rnd.Rannyu(0, 2 * M_PI);
            cities.emplace_back(city{cos(theta), sin(theta)});
        }
    }
    directory = "OUTPUT/" + simu_type[SIM] + "/";
    cout << "---------------------Map type == " << simu_type[SIM]
         << "---------------------" << endl;

    // Print starting configuration
    ofstream file_config(directory + "start_map.dat");
    file_config << setw(12) << "X" << setw(12) << "Y" << endl;
    for (auto c : cities)

        file_config << setw(12) << c["X"] << setw(12) << c["Y"] << endl;

    // final came back
    file_config << setw(12) << cities[0]["X"] << setw(12) << cities[0]["Y"]
                << endl;

    file_config.close();
    // start with the ordered travel
    vector<int> way(N_city - 1);
    for (int i = 1; i < N_city; i++)
        way[i - 1] = i;
    road route{way};
    vector<road> roads(DIM_POP);
    int i{}, j{};
    // generate DIM POP travels with random swaps
    for (int k = 0; k < DIM_POP; k++)
    {
        for (int l{}; l < 2 * N_city; l++)
        { // 2*N swaps
            i = rnd.Rannyu(0, N_city - 1);
            j = rnd.Rannyu(0, N_city - 1);
            while (i == j)
            {
                j = rnd.Rannyu(0, N_city - 1);
            }

            route.Swap_road(i % (N_city - 1), j % (N_city - 1)); // swap of two
                                                                 // cities
        }
        roads[k] = route;
    }
    GA ts{cities, roads}; // initialize

    // output a video
    std::cout << "Initial best road"
              << "\n";
    ts.Print_Best();

    ofstream file_averange{directory + "Average_L2.dat"};
    ofstream file_best{directory + "BestRoad_L2.dat"};

    file_best << "Gen:"
              << "\tL_2"
              << "\n";
    file_averange << "Gen:"
                  << "\t<L_2>"
                  << "\n";
    file_averange << 0 << "\t" << ts.Average_Len() << "\n";
    // print config of best travel along generetions
    ts.Best_road(directory + "/roads/BestRoad0.dat");

    for (int i = 1; i <= N_GEN; i++)
    {
        cerr << ".";
        ts.New_Gen();
        file_averange << i << "\t" << ts.Average_Len() << "\n";
        file_best << i << "\t" << ts.Best_Len() << "\n";
        // ts.Best_road(directory + "/roads/BestRoad" + to_string(i) +
        //             ".dat");
    }
    // last best travel
    std::cout << "\nLast best road\n";
    ts.Print_Best();
    // print the best travel to file
    ts.Best_road(directory + "BestRoad.dat");

    ofstream file_data{directory + "output.dat"};

    file_data << "---------------------Map type == " << simu_type[SIM]
              << "---------------------\n\n";
    file_data
        << "Tot Gen:" << N_GEN << "\nBest road Lenght:" << ts.Best_Len()
        << "\n\n"
        << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"
           "\\";

    file_data.close();
    file_averange.close();
    file_best.close();
    return 0;
}

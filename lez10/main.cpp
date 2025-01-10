#include "../RndGen/random.h"
#include "GA.h"
#include "mpi.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main(int argc, char *argv[])
{
    int N_MIGR;  // number of cities
    int DIM_POP; // number of roads
    int N_GEN;   // number of generations
    int SIM;     // type of simulation : 0 no migr 1 migr
    string property;

    ifstream inputdat("INPUT/input.dat");
    while (!inputdat.eof())
    {
        inputdat >> property;
        if (property == "N_MIGR:") { inputdat >> N_MIGR; }
        else if (property == "DIM_POP:")
        {
            inputdat >> DIM_POP;
        }
        else if (property == "N_GEN:")
        {
            inputdat >> N_GEN;
        }
        else if (property == "SIM_TYPE:")
        {
            inputdat >> SIM;
        }
        else
            cerr << "PROBLEM: unknown input" << endl;
    }
    inputdat.close();
    ////////////////////////////////////////////////////
    int size{}, rank{};
    // Start MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size); // world size, how many continents
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // address
    MPI_Status stat;                      // status of the communication
    ////////////////////////////////////////////////////

    Random rnd;
    rnd.initialize(rank);
    cout << rank << setw(12) << size << endl;
    // parameters
    string directory{};
    switch (SIM)
    {
    case 0:
    {
        directory = "OUTPUT/Indipendent/";
        cout << "---------------------Indipendent---------------------\n";
        break;
    }
    case 1:
        directory = "OUTPUT/Migration/";
        cout << "---------------------Migration---------------------\n";
        break;
    default:
        cerr << "Invalid map type\n";
        cerr << "0: Indipendent\n1: Migration\n";
        exit(-1);
    }

    directory += to_string(rank);

    // read input ciyties position from file
    std::vector<city> cities;
    ifstream Input("INPUT/cap_prov_ita.dat");
    if (!Input)
    {
        cerr << "Error opening the input file"
             << "\n";
        return -1;
    }
    int N_city = 0;
    double x, y;
    while (Input >> x >> y)
    {
        cities.emplace_back(city(x, y));
        N_city++;
    }
    cities.reserve(N_city);
    cities.shrink_to_fit();
    Input.close();

    // Print cities configuration
    if (rank == 0)
    {
        ofstream config(directory + "start_map.dat");
        config << setw(12) << "X" << setw(12) << "Y"
               << "\n";
        for (auto c : cities)
        {
            config << setw(12) << c.getX() << setw(12) << c.getY() << "\n";
        }

        config << setw(12) << cities[0].getX() << setw(12) << cities[0].getY()
               << "\n";

        config.close();
    }
    // start with the ordered road
    vector<int> way(N_city - 1);
    for (int i = 1; i < N_city; i++)
    {
        way[i - 1] = i;
    }
    road route{way};
    int i, j;
    vector<road> roads(DIM_POP);
    // generate DIM POP roads with random swaps
    for (int k{}; k < DIM_POP; k++)
    {
        for (int l{}; l < 2 * N_city; l++)
        {

            i = rnd.Rannyu(0, N_city - 1);
            j = rnd.Rannyu(0, N_city - 1);
            while (i == j)
            {
                j = rnd.Rannyu(0, N_city - 1);
            }

            route.swap(i % (N_city - 1), j % (N_city - 1));
        }
        roads[k] = route;
    }

    GA ts{cities, roads, rnd}; // initialize

    // std::cout << "\nRank:" << rank << "\tInitial best road lenght = " <<
    // ts.BestLength() << "\n";

    ofstream ave{directory + "AverageL_2.dat"};
    ofstream best{directory + "BestTourL_2.dat"};

    best << "Gen:"
         << "\tL_2"
         << "\n";
    ave << "Gen:"
        << "\t<L_2>"
        << "\n";
    ave << 0 << "\t" << ts.AverageLength() << "\n";
    // print config of best road along generetions
    // ts.Besttour(directory + "/tours_coord/BestTour0.dat");
    //	cout << "+" << endl;

    for (int i{1}; i <= N_GEN; i++)
    {

        cout << ".";
        if (SIM == 0)
        {
            ts.NewGeneration(); // indipendent
        }
        else
        {
            // ofstream Migrout("OUTPUT/Migration/peregrinage.dat");
            if (i % N_MIGR == 0)
            {
		cout << endl;
                // first migrate then generate new generation
                int giver{};
                int receiver{};
                vector<int> migrator(N_city -
                                     1); // int vect to store de index order of
                                         // a road aka cities order
                int itag = 1;
                // only the master will choose the individual to migrate
                if (rank == 0)
                {
                    giver = static_cast<int>(
                        rnd.Rannyu(0, size)); // choose the starting core
                    receiver = static_cast<int>(
                        rnd.Rannyu(0, size)); // choose the destination core
                    while (receiver == giver)
                    {
                        receiver = static_cast<int>(rnd.Rannyu(0, size));
                    }
                    // Migrout << "\nMigration at Gen: " << i
                    //<< "\tGiver" << giver << "\tReceiver:" << receiver ;
                }
                // master broadcasts the migration data
                MPI_Bcast(&giver, 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Bcast(&receiver, 1, MPI_INT, 0, MPI_COMM_WORLD);

                // Giver gives best individual
                if (rank == giver)
                {
                    for (int i = 0; i < N_city - 1; i++)
                    {
                        migrator[i] = ts[0][i];
                    }
                    MPI_Send(migrator.data(), N_city - 1, MPI_INT, receiver,
                             itag, MPI_COMM_WORLD);
                }
                //  receiver ubstitute the worst individual
                if (rank == receiver)
                {
                    MPI_Recv(migrator.data(), N_city - 1, MPI_INT, giver, itag,
                             MPI_COMM_WORLD, &stat);
                    ts[DIM_POP - 1] = road(migrator);
                }
            }
            // Migrout.close();
            ts.NewGeneration();
        }

        ave << i << "\t" << ts.AverageLength() << "\n";
        best << i << "\t" << ts.BestLength() << "\n";
        /*if(i%10==0 && rank == 0){
            ts.Besttour(directory + "/tours_coord/BestTour" + to_string(i)
        +".dat");
        }*/
    }

    ave.close();
    best.close();
    std::cout << "\nRank:" << rank
              << "\tLast best road lenght = " << ts.BestLength() << "\n";

    ofstream data{directory + "output.dat"};

    if (SIM == 0)
    { data << "---------------------Indipendent---------------------\n\n"; }
    else
        data << "---------------------Migration---------------------\n\n";

    data << "Tot Gen:" << N_GEN << "\nBest road Lenght:" << ts.BestLength()
         << "\n\n"
         << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"
            "\\";

    data.close();

    // print the final best road to file
    ts.Besttour(directory + "Besttour.dat");
    MPI_Finalize();
    return 0;
}

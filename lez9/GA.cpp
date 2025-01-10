#include "GA.h"

#include <iomanip>

// Constructor of the road class
road::road(std::vector<int> road) : _road(std::move(road))
{
    if (!Check_road(_road))
    {
        std::cerr << "error: the road is not valid\n";
        exit(-1);
    }
}

// check if all IDs are unique
bool road::Check_road(const std::vector<int> &vec)
{
    std::unordered_set<int> labels;

    // Controlla se il numero è già presente nel set
    if (!std::all_of(vec.begin(), vec.end(),
                     [&](int num) { return labels.insert(num).second; }))
    {
        std::cerr << "error: the road contains duplicate cities\n";
        return false;
    }

    return true;
}

// Calculate the distance between two cities using the L1 norm
double road::Dist_L1(city c, city d)
{
    return fabs(c["X"] - d["X"]) + fabs(c["Y"] - d["Y"]);
}
double road::Len_L1(const std::vector<city> &cities)
{
    if (_road.empty()) return 0.0;
    // initialize the length with the distance between the first and the last
    // city and the first and second city the first city is the one with ID 0
    // and the last one is the one with ID _road[_road.size()-1]
    double len = {Dist_L1(cities[0], cities[_road[_road.size() - 1]]) +
                  Dist_L1(cities[0], cities[_road[0]])};
    // add the distance between the current city and the next one
    for (unsigned int i = 0; i < _road.size() - 1; i++)
        len += Dist_L1(cities[_road[i]], cities[_road[i + 1]]);

    return len;
}

// Calculate the distance between two cities using the L2 norm
double road::Dist_L2(city c, city d)
{
    return sqrt(pow(c["X"] - d["X"], 2) + pow(c["Y"] - d["Y"], 2));
}
double road::Len_L2(const std::vector<city> &cities)
{
    if (_road.empty()) return 0.0;
    // initialize the length with the distance between the first and the last
    // city and the first and second city
    double length = {Dist_L2(cities[0], cities[_road[_road.size() - 1]]) +
                     Dist_L2(cities[0], cities[_road[0]])};
    for (unsigned int i = 0; i < _road.size() - 1; i++)
        length += Dist_L2(cities[_road[i]], cities[_road[i + 1]]);

    return length;
}

/**
Print the road
**/
void road::Print_road()
{
    std::cout << "[0, ";
    std::for_each(_road.begin(), _road.end(),
                  [](int i) { std::cout << i << ", "; });
    std::cout << "\b\b]\n";
}

/**
Shift of m contiguous cities of n positions starting from the i-th city
**/
void road::Shift_road(int i, int m, int n)
{
    // apply swap m times
    for (int j = 0; j < m; j++)
    {
        // swap the i+j-th city with the (i+j+n)-th city
        std::swap(_road[(i + j) % _road.size()],
                  _road[(i + j + n) % _road.size()]);
    }
}

/**
Permutation of m cities with other different m cities
**/
void road::Permutation_road(int m, int start, int sep)
{
    m %= _road.size() / 2; // m must be less than the half of the road
    int size = _road.size();
    int diff = size - m;
    sep %= diff;
    for (int j = 0; j < m; j++)
    {
        // swap the (start+j)-th city with the (start+m+sep+j)-th city
        std::swap(_road[(start + j) % _road.size()],
                  _road[(start + m + sep + j) % _road.size()]);
    }
}

/**
Inversion of the order of m contiguous cities starting from the i-th city
**/
void road::Inversion_road(int m, int i)
{
    m %= _road.size(); // m must be less than the length of the road
    int size = _road.size();
    for (int j = 0; j < m / 2; j++)
    {
        // swap the (i+j)-th city with the (i+m-j-1)-th city
        //  i, i+1, i+2, ..., i+m-1 (are m cities)
        std::swap(_road[(i + j) % size], _road[(i + m - 1 - j) % size]);
    }
}

// Constructor
GA::GA() { _rnd.initialize(); }

GA::GA(std::vector<city> cities) : _cities(std::move(cities))
{
    _rnd.initialize();
}

GA::GA(std::vector<city> cities, std::vector<road> roads)
    : _cities(std::move(cities)), _roads(std::move(roads))
{
    _rnd.initialize();
}
/**
Selection
**/
int GA::Select_road()
{
    // use M*rand^p to select the road
    // lowering p the prob of pick the parent in best half is descendent
    return _roads.size() * pow(_rnd.Rannyu(), 4);
}

/**
New generation of roads
**/
void GA::New_Gen()
{
    std::vector<road> new_gen(_roads.size());
    // CrossOver

    sort(_roads.begin(), _roads.end(), [&](road &a, road &b) {
        return a.Len_L2(_cities) < b.Len_L2(_cities);
    });

    unsigned int start{0u};
    if (_roads.size() % 2 != 0)
    {
        new_gen[0] = _roads[0];
        start++; // start from the second element
    }
    else
    {
        new_gen[0] = _roads[0];
        new_gen[1] = _roads[1];
        start += 2u; // start from the third element
    }

    int cut{}, i{}, j{};
    int N = _cities.size();
    road path;
    for (unsigned int k{start}; k < _roads.size(); k += 2)
    {
        // choose random road
        i = Select_road(); // first parent
        j = Select_road(); // second parent
        while (i == j)
            j = Select_road(); // ensure that the two parents are different
        std::vector<int> child1(N - 1);
        std::vector<int> child2(N - 1);

        // Crossover
        if (_rnd.Rannyu() < 0.53)
        {
            cut = _rnd.Rannyu(0, N - 1);

            std::copy_n(_roads[i].getroad().begin(), cut, child1.begin());
            std::copy_n(_roads[j].getroad().begin(), cut, child2.begin());

            for (int l{}, m1{cut}, m2{cut}; l < N - 1; l++)
            {
                // if the value is not already present in the child, add it
                if (std::count(child1.begin(), child1.end(),
                               _roads[j].getroad()[l]) == 0)
                { child1[m1++] = _roads[j].getroad()[l]; }
                if (std::count(child2.begin(), child2.end(),
                               _roads[i].getroad()[l]) == 0)
                { child2[m2++] = _roads[i].getroad()[l]; }
            }

            new_gen[k] = road(std::move(child1));
            new_gen[k + 1] = road(std::move(child2));
        }
        else
        {
            // if the crossover is not performed, copy the parents to the new
            // generation
            new_gen[k] = _roads[i];
            new_gen[k + 1] = _roads[j];
        }
    }

    // Mutations
    for (auto &path : new_gen)
        Mutations(path);

    // Assign the new generation of roads
    _roads = new_gen;

    // Sort the roads
    sort(_roads.begin(), _roads.end(), [&](road &a, road &b) {
        return a.Len_L2(_cities) < b.Len_L2(_cities);
    });
}

// Mutations
road &GA::Mutations(road &path)
{

    if (_rnd.Rannyu() < 0.1)

        path.Shift_road(_rnd.Rannyu(0, _cities.size() - 1),
                        _rnd.Rannyu(1, _cities.size() - 1),
                        _rnd.Rannyu(0, _cities.size() - 1));
    if (_rnd.Rannyu() < 0.1)

        path.Inversion_road(_rnd.Rannyu(0, _cities.size() - 1),
                            _rnd.Rannyu(0, _cities.size() - 1));
    if (_rnd.Rannyu() < 0.15)

        path.Permutation_road(_rnd.Rannyu(1, _cities.size() - 1),
                              _rnd.Rannyu(0, _cities.size() - 1),
                              _rnd.Rannyu(0, _cities.size() - 1));
    if (_rnd.Rannyu() < 0.2)
        path.Swap_road(_rnd.Rannyu(0, _cities.size() - 1),
                       _rnd.Rannyu(0, _cities.size() - 1));

    return path;
}

// Print the road with length
void GA::Print()
{
    for (auto &p : _roads)
    {
        p.Print_road();
        std::cout << "Length L2: " << p.Len_L2(_cities) << "\n";
    }
}

// Print the best road
void GA::Print_Best()
{
    _roads[0].Print_road();
    std::cout << "Length L2: " << _roads[0].Len_L2(_cities) << "\n";
}

// Average of the best half of the roads
double GA::Average_Len()
{
    double sum{};
    for (unsigned int i = 0; i < _roads.size() / 2; i++)
    {
        sum += _roads[i].Len_L2(_cities);
    }
    return sum / double(_roads.size() / 2);
}
// return the best length of the current population
double GA::Best_Len() { return _roads[0].Len_L2(_cities); }

// Print the best road on file
void GA::Best_road(const std::string &filename)
{
    std::ofstream out(filename);
    if (!out)
    {
        std::cerr << "error: open the file" << filename << "\n";
        exit(-1);
    }
    // Position of the first city
    out << std::setw(12) << "X" << std::setw(12) << "Y"
        << "\n";
    out << std::setw(12) << _cities[0]["X"] << std::setw(12) << _cities[0]["Y"]
        << "\n";
    for (auto i : _roads[0].getroad())
    {
        out << std::setw(12) << _cities[i]["X"] << std::setw(12)
            << _cities[i]["Y"] << "\n";
    }
    out << std::setw(12) << _cities[0]["X"] << std::setw(12) << _cities[0]["Y"]
        << "\n";

    out.close();
}

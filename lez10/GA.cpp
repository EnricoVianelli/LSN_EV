#include "GA.h"

// Constructor of the road class
road::road(std::vector<int> road) : _road(std::move(road))
{
    if (!Checktour(_road))
    {
        std::cerr << "error: the road is not valid\n";
        exit(-1);
    }
}

// check if all IDs are unique
bool road::Checktour(const std::vector<int> &vec)
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
double road::distanceL1(city c, city d)
{
    return fabs(c.getX() - d.getX()) + fabs(c.getY() - d.getY());
}

// Calculate the distance between two cities using the L2 norm
double road::distanceL2(city c, city d)
{
    return sqrt(pow(c.getX() - d.getX(), 2) + pow(c.getY() - d.getY(), 2));
}

double road::LengthL1(const std::vector<city> &cities)
{
    if (_road.empty()) return 0.0;
    // initialize the length with the distance between the first and the last
    // city and the first and second city the first city is the one with ID 0
    // and the last one is the one with ID _road[_road.size()-1]
    double length = {distanceL1(cities[0], cities[_road[_road.size() - 1]]) +
                     distanceL1(cities[0], cities[_road[0]])};
    // add the distance between the current city and the next one
    for (unsigned int i{}; i < _road.size() - 1; i++)
    {
        length += distanceL1(cities[_road[i]], cities[_road[i + 1]]);
    }
    return length;
}

double road::LengthL2(const std::vector<city> &cities)
{
    if (_road.empty()) return 0.0;
    // initialize the length with the distance between the first and the last
    // city and the first and second city
    double length = {distanceL2(cities[0], cities[_road[_road.size() - 1]]) +
                     distanceL2(cities[0], cities[_road[0]])};
    for (unsigned int i{}; i < _road.size() - 1; i++)
    {
        length += distanceL2(cities[_road[i]], cities[_road[i + 1]]);
    }
    return length;
}

// Print the road
void road::Printtour()
{
    std::cout << "[0, ";
    std::for_each(_road.begin(), _road.end(),
                  [](int i) { std::cout << i << ", "; });
    std::cout << "\b\b]\n";
}

// shift of m contiguous cities of n positions starting from the i-th city
void road::shift(int i, int m, int n)
{
    // apply swap m times
    for (int j = 0; j < m; j++)
    {
        // swap the i+j-th city with the (i+j+n)-th city
        std::swap(_road[(i + j) % _road.size()],
                  _road[(i + j + n) % _road.size()]);
    }
}

// Permutation of m cities with other different m cities
void road::Permutation(int m, int start, int sep)
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

// inversion of the order of m contiguous cities starting from the i-th city
void road::inversion(int m, int i)
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

GA::GA(std::vector<city> cities, std::vector<road> roads, Random rnd)
    : _cities(std::move(cities)), _roads(std::move(roads)), _rnd(rnd)
{
}

// Selection
int GA::Selection()
{
    // use M*rand^p to select the road
    // lowering p the prob of pick the parent in best half is descendent
    return _roads.size() * pow(_rnd.Rannyu(), 4);
}

// New generation of roads
void GA::NewGeneration()
{
    std::vector<road> new_gen(_roads.size());
    // CrossOver

    sort(_roads.begin(), _roads.end(), [&](road &a, road &b) {
        return a.LengthL2(_cities) < b.LengthL2(_cities);
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
        i = Selection(); // first parent
        j = Selection(); // second parent
        while (i == j)
        {
            j = Selection(); // ensure that the two parents are different
        }
        std::vector<int> child1(N - 1);
        std::vector<int> child2(N - 1);

        // Crossover
        if (_rnd.Rannyu() < 0.53)
        {
            cut = _rnd.Rannyu(0, N - 1);

            std::copy_n(_roads[i].gettour().begin(), cut, child1.begin());
            std::copy_n(_roads[j].gettour().begin(), cut, child2.begin());

            for (int l{}, m1{cut}, m2{cut}; l < N - 1; l++)
            {
                // if the value is not already present in the child, add it
                if (std::count(child1.begin(), child1.end(),
                               _roads[j].gettour()[l]) == 0)
                { child1[m1++] = _roads[j].gettour()[l]; }
                if (std::count(child2.begin(), child2.end(),
                               _roads[i].gettour()[l]) == 0)
                { child2[m2++] = _roads[i].gettour()[l]; }
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
    {
        Mutate(path);
    }

    // Assign the new generation of roads
    _roads = new_gen;

    // Sort the roads
    sort(_roads.begin(), _roads.end(), [&](road &a, road &b) {
        return a.LengthL2(_cities) < b.LengthL2(_cities);
    });
}

// Mutations
road &GA::Mutate(road &path)
{
    if (_rnd.Rannyu() < 0.2)
    {
        path.swap(_rnd.Rannyu(0, _cities.size() - 1),
                  _rnd.Rannyu(0, _cities.size() - 1));
    }
    if (_rnd.Rannyu() < 0.1)
    {
        path.shift(_rnd.Rannyu(0, _cities.size() - 1),
                   _rnd.Rannyu(1, _cities.size() - 1),
                   _rnd.Rannyu(0, _cities.size() - 1));
    }
    if (_rnd.Rannyu() < 0.15)
    {
        path.Permutation(_rnd.Rannyu(1, _cities.size() - 1),
                         _rnd.Rannyu(0, _cities.size() - 1),
                         _rnd.Rannyu(0, _cities.size() - 1));
    }
    if (_rnd.Rannyu() < 0.1)
    {
        path.inversion(_rnd.Rannyu(0, _cities.size() - 1),
                       _rnd.Rannyu(0, _cities.size() - 1));
    }
    return path;
}

// Print the road with length
void GA::Print()
{
    for (auto &p : _roads)
    {
        p.Printtour();
        std::cout << "Length L2: " << p.LengthL2(_cities) << "\n";
    }
}

// Print the best road
void GA::PrintBest()
{
    _roads[0].Printtour();
    std::cout << "Length L2: " << _roads[0].LengthL2(_cities) << "\n";
}

// Average of the best half of the roads
double GA::AverageLength()
{
    double sum{};
    for (unsigned int i{}; i < _roads.size() / 2; i++)
    {
        sum += _roads[i].LengthL2(_cities);
    }
    return sum / double(_roads.size() / 2);
}
// return the best length of the current population
double GA::BestLength() { return _roads[0].LengthL2(_cities); }

// Print the best road on file
void GA::Besttour(const std::string &filename)
{
    std::ofstream out(filename);
    if (!out)
    {
        std::cerr << "error: unable to open the file" << filename << "\n";
        exit(-1);
    }
    // Position of the first city
    out << std::setw(12) << "X" << std::setw(12) << "Y"
        << "\n";
    out << std::setw(12) << _cities[0].getX() << std::setw(12)
        << _cities[0].getY() << "\n";
    for (auto i : _roads[0].gettour())
    {
        out << std::setw(12) << _cities[i].getX() << std::setw(12)
            << _cities[i].getY() << "\n";
    }
    out << std::setw(12) << _cities[0].getX() << std::setw(12)
        << _cities[0].getY() << "\n";

    out.close();
}

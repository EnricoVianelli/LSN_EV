#pragma once

#include "../RndGen/random.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <unordered_set>
#include <vector>

class city
{
  private:
    double _x, _y;

  public:
    city() : _x(0), _y(0) {}
    city(double x, double y) : _x(x), _y(y) {}

    // Access operator var
    double &operator[](const std::string &var)
    {
        if (var == "X")
            return _x; // return the value of x
        else if (var == "Y")
            return _y; // return the value of y
        throw std::invalid_argument(var + ": Invalid variable name");
    }
};

class road
{
  private:
    std::vector<int> _road;

  public:
    road() = default;
    explicit road(const road &p) : _road(p._road) {}
    road(road &&p) noexcept : _road(std::move(p._road)) {}
    road(std::vector<int> road);

    road &operator=(const road &p)
    {
        if (this != &p) { _road = p._road; }
        return *this;
    }
    road &operator=(road &&p) noexcept
    {
        if (this != &p) { _road = std::move(p._road); }
        return *this;
    }

    // check if all IDs are unique
    bool Check_road(const std::vector<int> &vec);

    // swap of two cities
    void Swap_road(int i, int j) { std::swap(_road[i], _road[j]); }
    // shift of m contiguous cities of n positions starting from the i-th city
    void Shift_road(int i, int m, int n);
    // permutation among m contiguous cities starting from the i-th city
    void Permutation_road(int m, int i, int sep);
    // inversion of the order of m contiguous cities starting from the i-th city
    void Inversion_road(int m, int i);

    // Distance between two cities using the L1 norm
    double Dist_L1(city c, city d);
    // Total distance of the road using the L1 norm
    double Len_L1(const std::vector<city> &cities);
    // Distance between two cities using the L2 norm
    double Dist_L2(city c, city d);
    // Total distance of the road using the L2 norm
    double Len_L2(const std::vector<city> &cities);
    // Show the road
    void Print_road();

    // access to the road
    std::vector<int> getroad() const { return _road; }
};

// Traveling Salesman Problem
class GA
{
  private:
    // vector of cities ordered from 0 to _cities.size()-1
    std::vector<city> _cities;
    std::vector<road> _roads;
    Random _rnd;

  public:
    GA();
    GA(std::vector<city> cities);
    GA(std::vector<city> cities, std::vector<road> roads);

    // Selections
    int Select_road();

    // New generation of roads, update the population, old population are
    // deleted
    void New_Gen();

    // Mutations
    road &Mutations(road &way);

    // Print the road with length
    void Print();
    void Print_Best();

    // Average length of best half of the population
    double Average_Len();
    // Best length of the current population
    double Best_Len();
    // Print to file the Best road in Cartesian coordinates
    void Best_road(const std::string &filename);
};
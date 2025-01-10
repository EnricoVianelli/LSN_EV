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
    double getX() const { return _x; }
    double getY() const { return _y; }
    void setX(double x) { _x = x; }
    void setY(double y) { _y = y; }
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
    // Define operator[] for accessing the road elements
    int &operator[](size_t index) { return _road[index]; }

    const int &operator[](size_t index) const { return _road[index]; }
    // check if all IDs are unique
    bool Checktour(const std::vector<int> &vec);

    // Calculate the distance between two cities using the L1 norm
    double distanceL1(city c, city d);
    // Calculate the distance between two cities using the L2 norm
    double distanceL2(city c, city d);
    // Calculate the total distance of the road using the L1 norm
    double LengthL1(const std::vector<city> &cities);
    // Calculate the total distance of the road using the L2 norm
    double LengthL2(const std::vector<city> &cities);
    // Print the road
    void Printtour();
    // swap of two cities
    void swap(int i, int j) { std::swap(_road[i], _road[j]); }
    // shift of m contiguous cities of n positions starting from the i-th city
    void shift(int i, int m, int n);
    // permutation among m contiguous cities starting from the i-th city
    void Permutation(int m, int i, int sep);
    // inversion of the order of m contiguous cities starting from the i-th city
    void inversion(int m, int i);

    // access to the road
    std::vector<int> gettour() const { return _road; }
};

// GA alg
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
    GA(std::vector<city> cities, std::vector<road> tours);
    GA(std::vector<city> cities, std::vector<road> paths, Random rnd);

    inline std::vector<city> getCities() const { return _cities; }
    inline std::vector<road> getPaths() const { return _roads; }

    // access to the tourr
    road &operator[](int i)
    {
        if (i < 0 || (unsigned)i >= _roads.size())
        { throw std::out_of_range("Index out of range"); }
        return _roads[i];
    }

    // Selections
    int Selection();

    // New generation of tours, update the population, old population are
    // deleted
    void NewGeneration();

    // Mutations
    road &Mutate(road &way);

    // Print the road with lengt
    void Print();
    void PrintBest();

    // Average length of best half of the population
    double AverageLength();
    // Best length of the current population
    double BestLength();
    // Print to file the Best road in Cartesian coordinates
    void Besttour(const std::string &filename);
};

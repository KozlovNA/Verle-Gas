#include <iostream>
#include <array>
#include <vector>
#include <random>
#include <cmath>
using namespace std;


array<double, 3> get_random_v(double mi, double ma)
{   
    random_device r;
    default_random_engine e1(r());
    std::uniform_real_distribution<> uniform_dist(mi, ma);
    double v_x = uniform_dist(e1);
    double v_y = uniform_dist(e1);
    double v_z = uniform_dist(e1);
    array<double, 3> v{v_x, v_y, v_z};
    return v;
}

template<int N>
class Lattice
{
public:
    Lattice(double m_delta_x): _delta_x(m_delta_x)
    {    
    }

    array<double, 3> lat_pos()
    {   
        default_random_engine e1(42);
        std::uniform_int_distribution<int> uniform_dist(0, N-1);
        int x = uniform_dist(e1);
        int y = uniform_dist(e1);
        int z = uniform_dist(e1);
        while (_lat[x][y][z]!=0)
        {
            x = uniform_dist(e1);
            y = uniform_dist(e1);
            z = uniform_dist(e1);
        }
        _lat[x][y][z] = 1;
        array<double, 3> coordinates{x*_delta_x, y*_delta_x, z*_delta_x};
        return coordinates;
    }

    bool get_lat(int i, int j, int k) const
    {
        return _lat[i][j][k];
    }

protected:
    array<array<array<bool, N>, N>, N> _lat{};
    double _delta_x;
};


class Init
{
public:

    Init(double dt, int npart, double temp) : _npart(npart), _x(npart), _v(npart), _xm(npart), _temp(temp), _dt(dt)
    {
        _sumv2 = 0;
        Lattice<20> lat(0.1);
        for(int i = 0; i < _npart; i++)
        {
            _x[i] = lat.lat_pos();
            array<double, 3> vel = get_random_v(-1, 1);
            for (int j = 0; j < 3; j++)
            {
            _v[i][j] = vel[j];
            _sumv[j] += _v[i][j];
            }
            _sumv2 += _v[i][0]*_v[i][0] + _v[i][1]*_v[i][1] + _v[i][2]*_v[i][2];
        }
        for(int i = 0; i < 3; i++)
        {
            _sumv[i] = _sumv[i]/_npart;
        }
        _sumv2 = _sumv2/_npart;
        double fs = std::sqrt(_temp/_sumv2);
        for(int i = 0; i < _npart; i++)
        {
            for (int j = 0; j < 3; j++)
            {
            _v[i][j] = (_v[i][j] - _sumv[j])*fs;
            _xm[i][j] = _x[i][j] - _v[i][j]*_dt;
            }
        }
    }
    
    vector<array<double, 3>> get_x() const
    {
        return _x;
    }
    vector<array<double, 3>> get_xm() const
    {
        return _xm;
    }
    vector<array<double, 3>> get_v() const
    {
        return _v;
    }
    int get_npart() const
    {
        return _npart;
    } 

protected:     
    array<double, 3> _sumv{};
    double _sumv2;
    const int _npart;
    vector<array<double, 3>> _x;
    vector<array<double, 3>> _xm;
    vector<array<double, 3>> _v;
    double _temp;
    double _dt;
};

//void integrate();
//void sample();
//void force();


int main()
{   
    Init init(0.01, 10, 3);
    vector<array<double, 3>> x = init.get_x();
    vector<array<double, 3>> xm = init.get_xm();
    vector<array<double, 3>> v = init.get_v();
 
    /*
    double t = 0;
    double t_max = 1;
    double delt = 0.1;
    while(t < t_max)
    {
        force(f, en);
        integrate(f, en);
        t = t+delt;
        sample();
    }
    */
    return 0;
}
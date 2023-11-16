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

    Init(int npart, double delta_x, double dt, double temp) : _npart(npart), _x(npart), _v(npart), _xm(npart), _temp(temp), _dt(dt)
    {
        _sumv2 = 0;
        Lattice<20> lat(delta_x);
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
    
    vector<array<double, 3>>* get_x()
    {
        return &_x;
    }
    vector<array<double, 3>>* get_xm()
    {
        return &_xm;
    }
    vector<array<double, 3>>* get_v()
    {
        return &_v;
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


class Force
{
public:
    Force(int npart, double delta_x, double rc): _npart(npart), _rc(rc), _f(npart), _box(npart*delta_x) 
    {
        _rc2 = rc*rc;
        _en = 0;
        _ecut = 4*(1/pow(_rc2, 6)-1/pow(_rc2, 3));
    }
    void calculate_f(vector<array<double, 3>>* x)
    {
        for(int i = 0; i < _npart-1; i++)
        {
            for(int j = i+1; j < _npart; j++)
            {
                array<double, 3> xr;
                for (int k = 0; k < 3; k++)
                {
                   xr[k] = x->at(i)[k]-x->at(j)[k];
                   xr[k] = xr[k] - _box*round(xr[k]/_box); 
                }
                double r2 = xr[0]*xr[0] + xr[1]*xr[1] + xr[2]*xr[2];
                if (r2 < _rc2)
                {
                    double r2i = 1/r2;
                    double r6i = r2i*r2i*r2i;
                    double ff = 48*r2i*r6i*(r6i-0.5);
                    for (int k = 0; k < 3; k++)
                    {
                        _f[i][k] = _f[i][k] + ff*xr[k];
                        _f[j][k] = _f[j][k] - ff*xr[k];
                        _en = _en + 4*r6i*(r6i - 1) - _ecut;
                    } 
                } 
            }
        }
    }
    vector<array<double, 3>>* get_f()
    {
        return &_f;
    }
    double get_en()
    {
        return _en;
    }
protected:
    double _en;
    int _npart;
    vector<array<double, 3>> _f;
    double _box;
    double _rc;
    double _rc2;
    double _ecut;     
};

class Integrate_Verle{
public:
    Integrate_Verle(int npart, double dt): _sumv2(0), _npart(npart), _dt(dt), _temp(0)
    {
    }

    void integrate(vector<array<double, 3>>* f, vector<array<double, 3>>* x, vector<array<double, 3>>* xm, double en)
    {   
        _sumv2 = 0;
        for (int i = 0; i < _npart; i++)
        {   
            vector<array<double, 3>> v(_npart);
            double xx;
            for(int j = 0; j < 3; j++)
            {
                xx = 2*(x->at(i)[j]) - xm->at(i)[j] + _dt*_dt*f->at(i)[j];
                v[i][j] = (xx - xm->at(i)[j])/(2*_dt);
                _sumv[j] = _sumv[j] + v[i][j];
                _sumv2 = _sumv2 + v[i][j]*v[i][j];
                xm->at(i)[j] = x->at(i)[j];
                x->at(i)[j] = xx;
            }    
        }
        _temp = _sumv2/(3*_npart);
        _etot = (en + 0.5*_sumv2)/_npart;
    }
protected:
    array<double, 3> _sumv{};
    double _sumv2;
    int _npart;
    double _dt;
    double _temp;
    double _etot;
};
//void integrate();
//void sample();



int main()
{   
    int npart = 10;
    double dt = 0.01;
    double dx = 1;
    Init init(npart, dx, dt, 2);
    for (int i = 0; i < npart; i++)
    {
        vector<array<double, 3>>* x = init.get_x();
        cout << x->at(i)[0] << " " << x->at(i)[1] << " " << x->at(i)[2] << "\n";
    }
    cout << "\n\n\n";
    Force force(npart, dx, 6*dx);
    force.calculate_f(init.get_x());
    for (int i = 0; i < npart; i++)
    {
        vector<array<double, 3>>* f = force.get_f();
        cout << f->at(i)[0] << " " << f->at(i)[1] << " " << f->at(i)[2] << "\n";
    }
    cout << "\n\n\n";
    Integrate_Verle verle(npart, dt);
    verle.integrate(force.get_f(), init.get_x(), init.get_xm(), force.get_en());
    for (int i = 0; i < npart; i++)
    {
        vector<array<double, 3>>* x = init.get_x();
        cout << x->at(i)[0] << " " << x->at(i)[1] << " " << x->at(i)[2] << "\n";
    }
    cout << "\n\n\n";

    
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
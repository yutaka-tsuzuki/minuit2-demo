#include "Minuit2/FCNBase.h"
#include <vector>
#include <cmath>

class Fcn : public ROOT::Minuit2::FCNBase {
private:
	// array of data points
	std::vector<double> a_x;
	std::vector<double> a_y;
	// error
	double error_def;
public:
	// inherited methods
	void SetErrorDef(double edef) { error_def = edef; }
	virtual double Up() const { return error_def; }
	// const, dest
	Fcn() {}
	virtual ~Fcn() {}
	// add data points
	void add_point(double, double);
	// call
	virtual double operator()(const std::vector<double>&) const;
};

void Fcn::add_point(double x, double y) {
	a_x.push_back(x);
	a_y.push_back(y);
}

double Fcn::operator()(const std::vector<double>& a_par) const {
	double err = 0.0;
	int n_point = a_x.size();
	for(int i_point = 0; i_point < n_point; ++i_point) {
		double val_fcn = 0.0;
		int n_par = a_par.size();
		for(int i_par = 0; i_par < n_par; ++i_par) {
			val_fcn += a_par[i_par] * std::pow(a_x[i_point], i_par);
		}
		err += (val_fcn - a_y[i_point]) * (val_fcn - a_y[i_point]);
	}
	return err;
}


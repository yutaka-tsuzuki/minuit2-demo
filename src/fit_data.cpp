#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "Fcn.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnMigrad.h"


void read_csv(std::string fname, std::vector<double> &a_x, std::vector<double> &a_y) {
	// file stream
	std::ifstream ifs(fname);
	// data buffer
	std::string buf_str, bufbuf_str;
	// read each line
	while(std::getline(ifs, buf_str)) {
		// something ridiculous
		std::istringstream strstr(buf_str);
		std::vector<double> a_var;
		while(std::getline(strstr, bufbuf_str, ',')) {
			a_var.push_back(std::stod(bufbuf_str));
		}
		a_x.push_back(a_var[0]);
		a_y.push_back(a_var[1]);
	}
}


int fit_data(std::string fname) {
	// read csv file
	std::vector<double> a_x, a_y;
	read_csv(fname, a_x, a_y);
	
	// function object
	Fcn fcn;

	// add data points
	int n_point = a_x.size();
	for(int i_point = 0; i_point < n_point; ++i_point) {
		fcn.add_point(a_x[i_point], a_y[i_point]);
	}

	// fit
	ROOT::Minuit2::MnUserParameters par;
	par.Add("p0", 0.0, 0.1);
	par.Add("p1", 1.0, 0.1);
	par.Add("p2", 0.0, 0.1);
	par.Add("p3", 0.0, 0.1);
	ROOT::Minuit2::MnMigrad migrad(fcn, par);
	ROOT::Minuit2::FunctionMinimum min = migrad();

	// print parameters
	std::cout << min << std::endl;

	return 0;
}

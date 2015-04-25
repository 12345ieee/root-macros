#include <TLegend.h>
#include <TRandom3.h>
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include <string>
#include <cstring>
#include <fstream>
#include <iostream>

const double xmax= 76.1;
const double ymax= 168;

inline bool is_inside_x(double x)
{
	return 0<x && x<xmax;
}

inline bool is_inside_y(double y)
{
	return 0<y && y<ymax;
}

const int ydiv = 27;
const double ystep = ymax/ydiv;

const double eff1[24][ydiv]={{0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431}, {0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431}, {0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431}, {0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431}, {0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431}, {0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.833333, 0.722222, 0.75, 0.533333, 0.307692, 0.9, 0.625, 0.6, 0.625, 0.25, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431}, {0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.333333, 0.486431, 0.486486, 0.544118, 0.407407, 0.538462, 0.298507, 0.428571, 0.428571, 0.405405, 0.34375, 0.190476, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431}, {0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.75, 0.6, 0.759259, 0.722892, 0.745098, 0.744681, 0.925, 0.5, 0.583333, 0.558824, 0.583333, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.5, 0.486431, 0.71875, 0.602941, 0.887324, 0.454545, 0.764706, 0.758065, 0.604651, 0.649123, 0.65625, 0.6, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431}, {0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.5, 0.486431, 0.529412, 0.666667, 0.615385, 0.521739, 0.516667, 0.560976, 0.533333, 0.317073, 0.65, 0.130435, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431}, {0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.2, 1, 0.954545, 0.4875, 0.365385, 0.472527, 0.307692, 0.5, 0.217949, 0.301887, 0.275, 0.0833333, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431}, {0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.6, 0.491228, 0.61039, 0.425287, 0.414286, 0.46, 0.338235, 0.375, 0.32, 0.333333, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431}, {0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 1, 0.8, 0.450549, 0.477477, 0.405941, 0.461538, 0.488636, 0.263889, 0.396825, 0.226415, 0.121212, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431}, {0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 1, 0.486431, 0.486431, 0.6, 0.854167, 0.402439, 0.403509, 0.608696, 0.78125, 0.35, 0.555556, 0.607143, 0.333333, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431}, {0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.653846, 0.37931, 0.37931, 0.376623, 0.28169, 0.321429, 0.611111, 0.235294, 0.4375, 0.047619, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431}, {0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.6, 0.666667, 0.428571, 0.368421, 0.333333, 0.266667, 0.7, 0.6, 0.25, 0.3, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431}, {0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.153846, 0.176471, 0.3, 0.222222, 0.117647, 0.0625, 0.0909091, 0.0769231, 0.454545, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431}, {0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431}, {0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431, 0.486431}};
const double eff2[24][ydiv]={{0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}, {0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}, {0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}, {0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}, {0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}, {0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 1, 0.875, 0.6, 0.75, 0.6, 0.571429, 1, 0.6, 1, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}, {0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 1, 0.769231, 0.769231, 0.885714, 0.72093, 0.911111, 0.6, 0.772727, 0.909091, 1, 1, 1, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}, {0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.5, 0.875, 0.634921, 0.6, 0.770833, 0.6, 0.6, 0.8125, 0.857143, 0.166667, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}, {0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.833333, 0.74359, 0.793103, 0.811765, 0.606742, 0.712329, 0.661017, 0.653846, 0.604651, 0.55, 0.666667, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.372093, 0.835443, 0.777778, 0.971698, 0.795181, 0.857143, 0.894737, 0.6, 0.925926, 0.6, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}, {0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.6, 1, 0.644068, 0.823529, 0.661017, 0.806723, 0.851485, 0.754902, 0.602564, 0.571429, 0.85, 0.571429, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}, {0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.5, 0.142857, 0.333333, 0.657143, 0.846154, 0.674157, 0.85, 0.821918, 0.84375, 0.769231, 0.6, 0.6, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}, {0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.2, 0.689655, 0.79661, 0.731183, 0.757576, 0.706522, 0.848101, 0.847458, 0.805556, 0.5, 1, 0.666667, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}, {0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.333333, 0.53125, 0.59322, 0.573333, 0.863636, 0.779412, 0.723404, 0.535714, 0.740741, 0.6, 0.307692, 0.6, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}, {0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.263158, 0.288462, 0.403226, 0.366667, 0.5, 0.595238, 0.730769, 0.45, 0.6, 0.8, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}, {0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.6, 0.866667, 0.733333, 0.814815, 0.64, 0.3125, 1, 0.666667, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}, {0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.5, 0.411765, 0.142857, 0.266667, 0.166667, 0.454545, 0.2, 0.333333, 0.25, 0.25, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}, {0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}, {0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}, {0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}, {0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}, {0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}, {0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974, 0.661974}};
const double eff3[24][ydiv]={{0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984}, {0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984}, {0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984}, {0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984}, {0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984}, {0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.5, 0.6, 0.37931, 0.6, 0.5, 0.470588, 0.4, 0.333333, 0.428571, 0.6, 0.583333, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984}, {0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.25, 0.6, 0.6, 0.6, 0.923077, 0.6, 0.6, 0.916667, 0.6, 0.6, 0.6, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984}, {0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.2, 0.7, 0.611111, 0.569444, 0.722222, 0.701754, 0.465116, 0.571429, 0.782609, 0.447368, 0.8, 0.25, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984}, {0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.6, 0.6, 0.864865, 0.957447, 0.478261, 0.6, 0.782609, 0.814815, 0.45, 0.647059, 0.6, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984}, {0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.6, 0.6, 0.807692, 0.650794, 0.561644, 0.528302, 0.470588, 0.533333, 0.571429, 0.695652, 0.714286, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984}, {0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.6, 1, 0.696429, 0.520548, 0.644068, 0.826923, 0.681818, 0.514286, 0.5, 0.59375, 1, 0.25, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984}, {0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.25, 0.2, 0.618984, 0.5, 0.6, 0.769231, 0.844444, 0.636364, 0.697674, 0.55814, 0.5, 0.740741, 0.384615, 0.692308, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984}, {0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.6, 0.2, 0.618984, 1, 0.6, 0.945946, 0.6, 0.659574, 0.85, 0.851852, 0.6, 1, 0.6, 0.631579, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984}, {0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.4, 0.6, 0.689655, 0.730159, 0.493976, 0.558824, 0.548387, 0.659091, 0.677419, 0.65, 0.419355, 0.166667, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984}, {0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.6, 0.545455, 0.587302, 0.54717, 0.704918, 0.509804, 0.724138, 0.62963, 1, 0.4, 0.142857, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984}, {0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.6, 0.6, 0.928571, 0.877193, 0.583333, 0.717949, 0.510638, 0.607143, 0.76, 0.75, 0.222222, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984,}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 1, 0.111111, 0.73913, 1, 0.74359, 0.842105, 0.833333, 0.954545, 0.6, 0.458333, 1, 0.777778, 0.5, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984}, {0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984}, {0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984, 0.618984}};

const double geff1 = 0.486431;
const double geff2 = 0.661974;
const double geff3 = 0.618984;

using namespace std;

double Traccia(int N_volte, int mode){

	TRandom *rand= new TRandom ();
	TCanvas *c1 = new TCanvas("c1","multipads",800,600);
	TCanvas *c2 = new TCanvas("c2","multipads",800,600);
	TF1 *f1= new TF1 ("f1", "sin(x)*pow(cos(x),2)", 0, TMath::Pi()/2);
	TH1D *h1= new TH1D ("Phi", "Distribuzione in Phi", 78, 0, 360);
	TH1D *h2= new TH1D ("Theta", "Distribuzione in Theta",100 , 0, 90);

	int N;
	int Na=0;
	const double z12 = 149   - 95.5;
	const double z23 =  95.5 - 49;
	const double xstep = 3.2;
	const double zt = z12+z23;

	for(N=0; N<N_volte; N++){
		
		Double_t x1=rand->Uniform(0,xmax);
		Double_t y1=rand->Uniform(0,ymax);
		Double_t theta1=f1->GetRandom();					//Coordinate della prima camera generate
		Double_t phi1=rand->Uniform(0, 2*TMath::Pi());
		
		Double_t d2=z12*tan(theta1);		//Coordinate polari della seconda camera
		Double_t x2=d2*cos(phi1)+x1;
		if (is_inside_x(x2)) {
			
			Double_t y2=d2*sin(phi1)+y1;
			if (is_inside_y(y2)) {
				
				Double_t d3=z23*tan(theta1);		//Coordinate polari della terza camera
				Double_t x3=d3*cos(phi1)+x2;
				if (is_inside_x(x3)) {
					
					Double_t y3=d3*sin(phi1)+y2;
					if (is_inside_y(y3)) {
						
						if (mode==1) {
							int x_1=floor(x1/xstep);
							if (x_1!=9 || x_1!=19 || x_1!=23) {
								int x_2=floor(x2/xstep);
								if (x_2!=10) {
									int x_3=floor(x3/xstep);
									if (x_3!=17 || x_3!=19 || x_3!=21) {
										Na++;
										h2->Fill(theta1*180/TMath::Pi());
										h1->Fill(phi1*180/TMath::Pi());
									}
								}
							}
						}
						else if (mode==0) {
							Na++;
							h2->Fill(theta1*180/TMath::Pi());
							h1->Fill(phi1*180/TMath::Pi());
						}
						else if (mode==2) {
							int x_1=floor(x1/xstep);
							int y_1=floor(y1/ystep);
							//cout << x_1 << " - " << y_1 << endl;
							if (rand->Uniform(0, 1) < eff1[x_1][y_1]) {
								int x_2=floor(x2/xstep);
								int y_2=floor(y2/ystep);
								if (rand->Uniform(0, 1) < eff2[x_2][y_2]) {
									int x_3=floor(x3/xstep);
									int y_3=floor(y3/ystep);
									if (rand->Uniform(0, 1) < eff3[x_3][y_3]) {
										Na++;
										h2->Fill(theta1*180/TMath::Pi());
										h1->Fill(phi1*180/TMath::Pi());
									}
								}
							}
						}
						else if (mode==3) {
							int x_1=floor(x1/xstep);
							int y_1=floor(y1/ystep);
							//cout << x_1 << " - " << y_1 << endl;
							if (rand->Uniform(0, 1) < eff1[x_1][y_1]) {
								int x_2=floor(x2/xstep);
								int y_2=floor(y2/ystep);
								if (rand->Uniform(0, 1) < eff2[x_2][y_2]) {
									int x_3=floor(x3/xstep);
									int y_3=floor(y3/ystep);
									if (rand->Uniform(0, 1) < eff3[x_3][y_3]) {
										Na++;
										
										double yt = sqrt((x_1-x_3)*(x_1-x_3)*xstep*xstep + (y1-y3)*(y1-y3));
										h2->Fill(atan2(yt, zt)*180/TMath::Pi());
										
										double yp = (y1-y3);
										double xp = (x_1-x_3)*xstep;
										h1->Fill(atan2(yp, xp)*180/TMath::Pi()+180);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	double acc=(double)Na/N;
	printf("N:\t%d\nNa:\t%d\nAcc:\t%f\n", N, Na, acc);

	c1->cd();
	h1->SetXTitle("Phi[deg]");
	h1->SetMinimum(0);
	h1->SetYTitle("Eventi");
	h1->Draw();


	c2->cd();
	h2->SetXTitle("Theta[deg]");
	h2->SetYTitle("Eventi");
	h2->Draw();

	return acc*1.23648*2*TMath::Pi();

}


		//printf("%f %f %f\n", x1, x2, x3);

		/*
			int pos_x=x/xstep;
			pos_y=;
			double eff_cella=eff[pos_x][pos_y];
			if (eff==1) {x}
			else if (eff==0) {y}
			else {
			double e = rand->Uniform(0,1);
			if (e>eff_cella) {x}
			else {y}
			}
		*/


void montecarlo(int mode=0)
{
	int N_volte=10200000;
	std::cout << std::endl << "Mode: " << mode << std::endl;
	printf("L'accettanza è %f\n", Traccia(N_volte, mode));
}

int main(int, char** argv)
{
	int mode=atoi(argv[1]); // normale - cmorti - full efficienza - ricostr
	montecarlo(mode);
}

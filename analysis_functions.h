#pragma once

#include "LLT/LiftingLine.h"

#include <matplot/matplot.h>
#include <Eigen/Dense>

#include <cmath>
#include <string>



typedef struct Wing {
    float span;
    std::function<float (const float)> theta, chord, alpha, alpha0, m;
    std::vector<float> y;
} Wing;

matplot::figure_handle plot_wing_planform(const Wing& wing);

typedef struct airfoil_parameters {
    float m, alpha0;
} airfoil_parameters;

airfoil_parameters get_airfoil_parameters(const std::string& airfoil_name);



Eigen::VectorXf solve_system(const Wing& testing_wing);

std::vector<float> calculate_gamma0(const Eigen::VectorXf& A, const Wing& testing_wing);

typedef struct Results {
    float wing_area, aspect_ratio, delta, E, Cl, Cd;
} Results;

Results calculate_results(const Wing& testing_wing, const Eigen::VectorXf& A);

void airfoil_select(const Wing& straightWing, float lift_required, float velocity);

void aspect_ratio_select(float span, const airfoil_parameters& airfoil, float lift_required, float velocity);

void planform_select(float span, const airfoil_parameters& airfoil, float lift_required, float velocity);

void twist_select(float span, const airfoil_parameters& airfoil, float lift_required, float velocity);

std::function<float (const float)> get_required_alpha(const Wing& current_wing, float lift_required, float velocity);
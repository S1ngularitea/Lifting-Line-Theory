#pragma once

#include <functional>
#include <vector>

#include <Eigen/Dense>

namespace LLT {

const float pi = 3.14159265;

float wing_area(float wing_span, float root_chord, float tip_chord);

float wing_area(float wing_span, const std::function<float (const float)>& chord_function);

float aspect_ratio(float wing_span, float wing_area);

std::function<float (const float)> theta(float wing_span);

std::function<float (const float)> linear_function(/*float y, */float root_value, float tip_value, float wing_span);

std::function<float (const float)> const_function(const float value);

std::function<float (const float)> elliptical_function(float root_value, float wing_span);

std::vector<float> uniform_grid(float wing_span, int number_of_point);

float C_yn(float y, int n, float wing_span, const std::function<float (const float)>& Theta, const std::function<float (const float)>& chord, const std::function<float (const float)>& airfoil_lift_slope);

typedef struct {
    Eigen::MatrixXf C;
    Eigen::VectorXf D;
} linear_system;

linear_system build_linear_system(const std::vector<float>& y, float wing_span, const std::function<float (const float)>& theta,  const std::function<float (const float)>& chord, const std::function<float (const float)>& airfoil_lift_slope, const std::function<float (const float)>& alpha, const std::function<float (const float)>& zero_lift_alpha);

float degrees_to_radians(float angle);

std::vector<float> gamma_theta(const std::vector<float>& y, const Eigen::VectorXf& A, const std::function<float (const float)>& theta);

float delta(const Eigen::VectorXf& A);

float efficiency_factor(float delta);

float lift_coefficient(const Eigen::VectorXf& A, float aspect_ratio);

float drag_coefficient(float Cl, float aspect_ratio, float efficency_factor);

float maximum(const std::vector<float>& input);

Eigen::VectorXf solve_system(const LLT::linear_system& system);

}
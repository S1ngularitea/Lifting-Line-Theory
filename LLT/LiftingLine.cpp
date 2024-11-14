#include "LiftingLine.h"

#include <cmath>

float LLT::wing_area(float wing_span, float root_chord, float tip_chord) {
    float area = wing_span *(root_chord+tip_chord)/2;
    return area;
}

float LLT::wing_area(float wing_span, const std::function<float(const float)> &chord_function) {
    int number_of_points = 50;
    float area = 0;
    auto y = LLT::uniform_grid(wing_span, number_of_points);
    for (int i = 0; i < number_of_points-1; i++) {
        area += std::abs(y[i+1]-y[i])*(chord_function(y[i])+chord_function(y[i+1]))/2;
    }
    area += y[0]*(chord_function(0)+chord_function(y[0]))/2;
    area += y[0]*(chord_function(y[number_of_points-1])+chord_function(wing_span/2))/2;

    return area*2;
}

float LLT::aspect_ratio(float wing_span, float wing_area) {
    return pow(wing_span, 2)/wing_area;
}

std::function<float (const float)> LLT::theta(float wing_span) {
    auto theta_function= [wing_span](const float y) {return std::acos(2 * y/wing_span);};

    return theta_function;
}

std::function<float (const float)> LLT::linear_function(const float root_value, const float tip_value, const float wing_span) {
    // set up a lambda function to return
    auto function = [root_value, tip_value, wing_span] (const float y) {
        float m = (root_value - tip_value) / (wing_span/2);

        return root_value - std::abs(y) * m;
    };

    return function;
}

std::function<float (const float)> LLT::const_function(const float value) {
    auto function = [value] (const float y ) {return value;};

    return function;
}

std::function<float(const float)> LLT::elliptical_function(float root_value, float wing_span) {
    auto function = [root_value, wing_span](const float y) {
        const float a = wing_span/2;
        const float b = 2*a*root_value/wing_span;

        return b*std::sqrt(1-(y*y)/(a*a));
    };

    return function;
}

std::vector<float> LLT::uniform_grid(float wing_span, int number_of_points) {
    std::vector<float> points(number_of_points);

    float half_span = wing_span/2;
    float delta = half_span/number_of_points;

    for (int i = 0; i < number_of_points; i++) {
        points[i] = delta/2 + delta * i;
    }

    return points;
}

float LLT::C_yn(float y, int n, float wing_span, const std::function<float (const float)>& Theta, const std::function<float (const float)>& chord, const std::function<float (const float)>& airfoil_lift_slope) {
    float n_odd = 2*n+1;
    float t1 = (4*wing_span)/(airfoil_lift_slope(y) * chord(y));
    float t2 = n_odd/sin(Theta(y));
    float t3 = t1+t2;
    return t3*std::sin(n_odd*Theta(y));
}

float LLT::degrees_to_radians(float angle) {
    return angle * LLT::pi/180;
}

LLT::linear_system LLT::build_linear_system(const std::vector<float>& y, float wing_span, const std::function<float (const float)>& theta,  const std::function<float (const float)>& chord, const std::function<float (const float)>& airfoil_lift_slope, const std::function<float (const float)>& alpha, const std::function<float (const float)>& zero_lift_alpha) {
    int nPoints = y.size();

    Eigen::MatrixXf C(nPoints, nPoints);
    Eigen::VectorXf D(nPoints);

    for (int j=0; j < nPoints; j++) {
        for (int i=0; i < nPoints; i++) {
            C(i,j) = LLT::C_yn(y[i], j, wing_span, theta, chord, airfoil_lift_slope);
            
        }

        D(j) = LLT::degrees_to_radians(alpha(y[j])-zero_lift_alpha(y[j]));
    }

    LLT::linear_system return_struct;
    return_struct.C = C;
    return_struct.D = D;

    return return_struct;
}

std::vector<float> LLT::gamma_theta(const std::vector<float>& y, const Eigen::VectorXf& A, const std::function<float (const float)>& theta) {
    int nPoints = y.size();
    std::vector<float> G(nPoints);

    for (int i = 0; i < nPoints; i++) {
        for (int n = 0; n < nPoints; n++) {
            G[i] += A[n]*sin((2*n+1)*theta(y[i]));
        }
    }

    return G;
}

float LLT::delta(const Eigen::VectorXf& A) {
    float sum = 0;
    int nCoeffs = A.size();

    for (int n = 1; n < nCoeffs; n++) {
        sum += (2*n+1)*pow((A[n]/A[0]),2);
    }

    return sum;
}

float LLT::efficiency_factor(float delta) {
    return 1/(1+delta);
}

float LLT::lift_coefficient(const Eigen::VectorXf& A, float aspect_ratio) {
    return LLT::pi * A(0) * aspect_ratio;
}

float LLT::drag_coefficient(float Cl, float aspect_ratio, float efficency_factor) {
    return pow(Cl,2)/(LLT::pi * efficency_factor *aspect_ratio);
}

float LLT::maximum(const std::vector<float>& input) {
    float max = input[0];

    for (int i = 0; i < input.size(); i++) {
        if (input[i] > max)
            max = input[i];
    }

    return max;
}

Eigen::VectorXf LLT::solve_system(const LLT::linear_system& system) {
    return (system.C).fullPivLu().solve(system.D);
}

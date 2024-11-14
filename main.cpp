#include "LLT/LiftingLine.h"
#include "analysis_functions.h"

#include <iostream>
#include <cmath>
#include <matplot/matplot.h>
#include <Eigen/Dense>

// definition from cmake to remove directory errors
#ifdef Airfoil_Folder 
std::string airfoil_folder = Airfoil_Folder;
#endif




int main() {

    // Testing Airfoils with straight wing

    float U = 40;
    float span = 4.5f;
    float chord = 1.5f;
    float alpha = 3.5f;

    Wing straightWing = {};
    straightWing.theta = LLT::theta(span);
    straightWing.span = span;
    straightWing.chord = LLT::const_function(chord);
    straightWing.alpha = LLT::const_function(alpha);
    straightWing.y = LLT::uniform_grid(span,10);

    float mass = 100;
    float air_density = 1.225;
    float g = 9.81;
    float weight = mass*g;
    float lift = weight;
    float wing_area = LLT::wing_area(straightWing.span ,straightWing.chord);
    float cl_required = (2*lift)/(air_density * U * U * wing_area);

    // pick airfoil
    //std::cout << "required cl = " << cl_required << '\n';
    airfoil_select(straightWing, lift, U);

    
    
    float wing_area_2 = 4.5*4.5/7;
    float cl_required_2 = (2*lift)/(air_density*U*U*wing_area_2);
    
    std::cout << "required cl = " << cl_required_2 << '\n';
    auto airfoil = get_airfoil_parameters("EPPLER_715");

    aspect_ratio_select(span, airfoil, lift, U);
    
    planform_select(span, airfoil, lift, U);

    twist_select(span, airfoil, lift, U);

    wing_area = span*span/7;

    Wing final_design = {};
    final_design.theta = LLT::theta(span);
    final_design.span = span;
    final_design.y = LLT::uniform_grid(span,10);
    final_design.alpha0 = LLT::const_function(airfoil.alpha0);
    final_design.m = LLT::const_function(airfoil.m);

    final_design.chord = LLT::const_function(wing_area/span);
    final_design.alpha = LLT::const_function(0);

    float prop_efficiency = 0.8;
    float fuel_consumpion = 0.523; // kg/(kwh)

    float dry_mass = 60;

    std::vector<float> test_masses = {90, 100, 110, 120, 125, 130, 140, 150};
    for (int i = 0; i < test_masses.size(); i++) {
        float weight = test_masses[i]*g;
        final_design.alpha = get_required_alpha(final_design, weight, U);

        Eigen::VectorXf avs_A = solve_system(final_design);
        Results fd_results = calculate_results(final_design, avs_A);

        std::cout << "mass = " << test_masses[i] <<" cl = " << fd_results.Cl << " E = " << fd_results.E << "\n";
        float fd_glide_slope = fd_results.Cl/fd_results.Cd;
        std::cout << "cl/cd = " << fd_glide_slope << '\n';
        std::cout << "alpha = " << final_design.alpha(0) << '\n';

        float drag = 0.5*1.225*U*U*fd_results.Cd*wing_area;

        float fuel_mass = test_masses[i]-dry_mass; // assume fuel mass is half of total mass
        float range = fuel_mass/(fuel_consumpion*prop_efficiency*drag) * 3600 /*convert from hours to s*/; // range in km
        std::cout << "range = " << range << '\n';
    }

    plot_wing_planform(final_design);
    
    matplot::show();

    return EXIT_SUCCESS;
}
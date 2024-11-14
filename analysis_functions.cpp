#include "analysis_functions.h"

#include <fstream>
// definition from cmake to remove directory errors
#ifdef Airfoil_Folder 
const std::string airfoil_folder = Airfoil_Folder;
#endif

matplot::figure_handle plot_wing_planform(const Wing& wing) {
    
    float span = wing.span;
    std::function<float (const float)> chord_function = wing.chord;

    std::vector<double> x = matplot::linspace(-span/2, span/2);
    auto leading_edge = [chord_function](const double y){ return (double)chord_function(y)/2; }; // has to be cast to double for matplot++ compiler throws warnings otherwise
    std::vector<double> leading_edge_line = matplot::transform(x, leading_edge);
    auto trailing_edge = [chord_function](const double y){ return -(double)chord_function(y)/2; };
    std::vector<double> trailing_edge_line = matplot::transform(x, trailing_edge);
    std::vector<double> wing_tips = {chord_function(span/2)/2, -chord_function(span/2)/2};
    matplot::hold(matplot::on);
    matplot::plot(x, leading_edge_line);
    matplot::plot(x, trailing_edge_line);
    matplot::plot({-span/2,-span/2}, wing_tips);
    matplot::plot({span/2, span/2}, wing_tips);

    matplot::title("Wing Planform");
    matplot::xlabel("spanwise coordinate (m)");
    matplot::ylabel("streamwise coordinate (m)");

    matplot::axis({-(span/2+0.5f),(span/2+0.5f), -(chord_function(0)/2+1), chord_function(0)/2+1});

    return matplot::figure();
}

typedef struct airfoil_data {
    std::vector<float> alpha, CL;
} airfoil_data;

const airfoil_data read_airfoil_data(const std::string& airfoil_file_name) {
#ifndef Airfoil_Folder
    return {std::vector<float>(),std::vector<float>()};// returns no data when not defined
#else
    std::string airfoil_path = airfoil_folder + "/" + airfoil_file_name;
    std::fstream airfoil_data_file(airfoil_path);

    airfoil_data raw_data;
    std::string line;
    int line_count = 0;
    int i = 0;
    
    while (std::getline(airfoil_data_file, line)) {
        if (i >= 11 && line.length() > 1) {
            line_count++;
        }
        i++;
    }

    airfoil_data_file.clear();
    airfoil_data_file.seekg(0);

    raw_data.alpha.resize(line_count);
    raw_data.CL.resize(line_count);
    std::string alpha_as_string, CL_as_string;
    std::string::size_type sz;

    i = 0;
    while (std::getline(airfoil_data_file, line)) {
        if(i >= 11 && line.length() > 1) {
            alpha_as_string = line.substr(1,7);
            CL_as_string = line.substr(10,7);
            raw_data.alpha[i-11]= std::stof(alpha_as_string);
            raw_data.CL[i-11] = std::stof(CL_as_string);
        }
        i++;
    }
    airfoil_data_file.close();

    return raw_data;

#endif
}

const std::string replace_underscores(const std::string& string) {
    std::string buffer_string(string.length(), ' ');

    for (int i = 0; i < string.length(); i++) {
        if (string[i] != '_') {
            buffer_string[i] = string[i];
        }
    }

    return buffer_string;
}

const std::array<char[6],6> reynolds_numbers = {
    "0.500",
    "1.000",
    "1.500",
    "2.000",
    "2.500",
    "3.000"
};

float find_alpha0(const airfoil_data& data) {
    float alpha0_lower_bound, alpha0_upper_bound;
    float cl_upper_bound, cl_lower_bound;
    for(int i = 0; i < data.alpha.size()-1; i++) {
        if (data.CL[i+1] > 0 && data.CL[i] < 0) {
            cl_upper_bound = data.CL[i+1];
            cl_lower_bound = data.CL[i];
            alpha0_upper_bound = data.alpha[i+1];
            alpha0_lower_bound = data.alpha[i];

            float easing_value = (0-cl_lower_bound)/(cl_upper_bound-cl_lower_bound);

            return /*LLT::degrees_to_radians*/(alpha0_lower_bound + easing_value*(alpha0_upper_bound-alpha0_lower_bound));
        }
    }

    return 0;
}

float find_m(const airfoil_data& data) {
    float Cl1, Cl2;
    for (int i = 0; i < data.CL.size(); i++) {
        if (data.alpha[i] <= -3)
            Cl1 = data.CL[i];
        
        if (data.alpha[i] <= 3)
            Cl2 = data.CL[i];
    }

    return (Cl2-Cl1)/LLT::degrees_to_radians(3-(-3));
}


airfoil_parameters get_airfoil_parameters(const std::string& airfoil_name) {
    
    std::string airfoil_name_reformatted = replace_underscores(airfoil_name);

    float alpha0_sum=0, m_sum=0;

    int i = 0;
    for (const auto& reynolds_number: reynolds_numbers) {
        airfoil_data data = read_airfoil_data(airfoil_name + '/' + airfoil_name_reformatted + "_T1_Re" + reynolds_number + "_M0.00_N9.0.txt");
        if (data.CL[0] == 0 && data.alpha[0] == 0)
            return {0,0};
        alpha0_sum += find_alpha0(data);
        m_sum += find_m(data);

        i++;
    }


    airfoil_parameters return_parameters;
    return_parameters.alpha0 = alpha0_sum/6;
    return_parameters.m = m_sum/6;
    return return_parameters;
}

Eigen::VectorXf solve_system(const Wing &testing_wing) {

    LLT::linear_system system = LLT::build_linear_system(testing_wing.y, testing_wing.span, testing_wing.theta, testing_wing.chord, testing_wing.m, testing_wing.alpha, testing_wing.alpha0);

    return LLT::solve_system(system);
}

std::vector<float> calculate_gamma0(const Eigen::VectorXf& A, const Wing& testing_wing) {
    std::vector<float> G0 = LLT::gamma_theta(testing_wing.y, A, testing_wing.theta);
    std::vector<float> G_G0(testing_wing.y.size());

    for (int i = 0; i < testing_wing.y.size(); i++) {
        G_G0[i] = G0[i]/LLT::maximum(G0);
    }

    return G_G0;
}

Results calculate_results(const Wing &testing_wing, const Eigen::VectorXf& A) {
    Results results = {};
    float area = LLT::wing_area(testing_wing.span, testing_wing.chord);
    float aspect_ratio = LLT::aspect_ratio(testing_wing.span,area);
    float delta = LLT::delta(A);
    float E = LLT::efficiency_factor(delta);
    float Cl = LLT::lift_coefficient(A, aspect_ratio);
    float Cd = LLT::drag_coefficient(Cl, aspect_ratio, E);

    return {area, aspect_ratio, delta, E, Cl, Cd};
}



std::function<float (const float)> get_required_alpha(const Wing& current_wing, float lift_required, float velocity) {
    Wing test_wing = Wing(current_wing);
    float area = LLT::wing_area(current_wing.span, current_wing.chord);
    float Cl = 0;
    float lift_generated = 0;
    std::function<float (const float)> current_wing_alpha = nullptr;
    if (current_wing.alpha) {
        current_wing_alpha = current_wing.alpha;
    } else {
        current_wing_alpha = LLT::const_function(0);
    }
    float alpha = 0;
    while (std::abs(lift_required-lift_generated)>0.1) {
        auto alpha_func = [alpha, current_wing_alpha](const float y) {
            return current_wing_alpha(y) + alpha;
        };
        test_wing.alpha = alpha_func;

        Eigen::VectorXf A = solve_system(test_wing);
        Results results = calculate_results(test_wing, A);

        Cl = results.Cl;
        lift_generated = 0.5 * 1.225 * area * Cl * velocity*velocity;
        
        alpha += (lift_required-lift_generated)/lift_required;
        if (alpha >=45) {
            alpha = 45;
            break;
        }
    }
    auto return_func = [alpha, current_wing_alpha](const float y) {return current_wing_alpha(y) + alpha; };
    return return_func;
}


void airfoil_select(const Wing& straightWing, float lift_required, float velocity) {

    // Testing AVISTAR
    auto avistar_parameters = get_airfoil_parameters("AVISTAR");

    float avs_alpha0_const = avistar_parameters.alpha0; // functions of airfoil
    float avs_m_const = avistar_parameters.m; // function of airfoil

    Wing AVISTAR = Wing(straightWing);
    AVISTAR.alpha0 = LLT::const_function(avs_alpha0_const);
    AVISTAR.m = LLT::const_function(avs_m_const);
    AVISTAR.alpha = get_required_alpha(AVISTAR, lift_required, velocity);

    Eigen::VectorXf avs_A = solve_system(AVISTAR);
    Results avistar_results = calculate_results(AVISTAR, avs_A);

    std::cout << "avistar cl = " << avistar_results.Cl << " avistar E = " << avistar_results.E << "\n";
    float avs_glide_slope = avistar_results.Cl/avistar_results.Cd;
    std::cout << "cl/cd = " << avs_glide_slope << '\n';
    std::cout << "alpha = " << AVISTAR.alpha(0) << '\n';


    // Testing E209
    auto E209_parameters = get_airfoil_parameters("E209");

    float E209_alpha0_const = E209_parameters.alpha0;
    float E209_m_const= E209_parameters.m;

    Wing E209 = Wing(straightWing);
    E209.alpha0 = LLT::const_function(E209_alpha0_const);
    E209.m = LLT::const_function(E209_m_const);
    E209.alpha = get_required_alpha(E209, lift_required, velocity);

    Eigen::VectorXf E209_A = solve_system(E209);
    Results E209_results = calculate_results(E209, E209_A);

    std::cout << "E209 cl = " << E209_results.Cl << " E209 E = " << E209_results.E << '\n';
    float E209_glide_slope = E209_results.Cl/E209_results.Cd;
    std::cout << "cl/cd = " << E209_glide_slope << '\n';

    // Testing EPPLER 715
    auto EP715_parameters = get_airfoil_parameters("EPPLER_715");

    float EP715_alpha0_const = EP715_parameters.alpha0;
    float EP715_m_const = EP715_parameters.m;

    Wing EP715 = Wing(straightWing);
    EP715.alpha0 = LLT::const_function(EP715_alpha0_const);
    EP715.m = LLT::const_function(EP715_m_const);
    EP715.alpha = get_required_alpha(EP715, lift_required, velocity);

    Eigen::VectorXf EP715_A = solve_system(EP715);
    Results EP715_results = calculate_results(EP715, EP715_A);

    std::cout << "EPPLER 715 cl = " << EP715_results.Cl << " E = " << EP715_results.E << '\n';
    float EP715_glide_slope = EP715_results.Cl/EP715_results.Cd;
    std::cout << "cl/cd = " << EP715_glide_slope << '\n';


    // Testing NACA 2414
    auto NACA_2414_parameters = get_airfoil_parameters("NACA_2414");

    float N2414_alpha0_const = NACA_2414_parameters.alpha0;
    float N2414_m_const = NACA_2414_parameters.m;

    Wing N2414 = Wing(straightWing);
    N2414.alpha0 = LLT::const_function(N2414_alpha0_const);
    N2414.m = LLT::const_function(N2414_m_const);
    N2414.alpha = get_required_alpha(N2414, lift_required, velocity);

    Eigen::VectorXf N2414_A = solve_system(N2414);
    Results N2414_results = calculate_results(N2414, N2414_A);

    std::cout << "NACA 2414 cl = " << N2414_results.Cl << " E = " << N2414_results.E << '\n';
    float N2414_glide_slope = N2414_results.Cl/N2414_results.Cd;
    std::cout << "cl/cd = " << N2414_glide_slope << '\n';


    // Testing NACA 2415
    auto NACA_2415_parameters = get_airfoil_parameters("NACA_2415");

    float N2415_alpha0_const = NACA_2415_parameters.alpha0;
    float N2415_m_const = NACA_2415_parameters.m;

    Wing N2415 = Wing(straightWing);
    N2415.alpha0 = LLT::const_function(N2415_alpha0_const);
    N2415.m = LLT::const_function(N2415_m_const);
    N2415.alpha = get_required_alpha(N2415, lift_required, velocity);

    Eigen::VectorXf N2415_A = solve_system(N2415);
    Results N2415_results = calculate_results(N2415, N2415_A);

    std::cout << "NACA 2415 cl = " << N2415_results.Cl << " E = " << N2415_results.E << '\n';
    float N2415_glide_slope = N2415_results.Cl/N2415_results.Cd;
    std::cout << "cl/cd = " << N2415_glide_slope << '\n';


    // Testing S2027
    auto S2027_parameters = get_airfoil_parameters("S2027");

    float S2027_alpha0_const = S2027_parameters.alpha0;
    float S2027_m_const = S2027_parameters.m;

    Wing S2027 = Wing(straightWing);
    S2027.alpha0 = LLT::const_function(S2027_alpha0_const);
    S2027.m = LLT::const_function(S2027_m_const);
    S2027.alpha = get_required_alpha(S2027, lift_required, velocity);

    Eigen::VectorXf S2027_A = solve_system(S2027);
    Results S2027_results = calculate_results(S2027, S2027_A);

    std::cout << "S2027 cl = " << S2027_results.Cl << " E = " << S2027_results.E << '\n';
    float S2027_glide_slope = S2027_results.Cl/S2027_results.Cd;
    std::cout << "cl/cd = " << S2027_glide_slope << '\n';


    // Testing SA7036
    auto SA7036_parameters = get_airfoil_parameters("SA7036");

    float SA7036_alpha0_const = SA7036_parameters.alpha0;
    float SA7036_m_const = SA7036_parameters.m;

    Wing SA7036 = Wing(straightWing);
    SA7036.alpha0 = LLT::const_function(SA7036_alpha0_const);
    SA7036.m = LLT::const_function(SA7036_m_const);
    SA7036.alpha = get_required_alpha(SA7036, lift_required, velocity);

    Eigen::VectorXf SA7036_A = solve_system(SA7036);
    Results SA7036_results = calculate_results(SA7036, SA7036_A);

    std::cout << "SA7036 cl = " << SA7036_results.Cl << " E = " << SA7036_results.E << '\n';
    float SA7036_glide_ratio = SA7036_results.Cl/SA7036_results.Cd;
    std::cout << "cl/cd = " << SA7036_glide_ratio << '\n';


}

void aspect_ratio_select(float span, const airfoil_parameters& airfoil, float lift_required, float velocity) {
    
    float AR_testing_wing_area = 4.5*4.5/7;

    Wing AR_test_wing_base = {};
    AR_test_wing_base.theta = LLT::theta(span);
    AR_test_wing_base.y = LLT::uniform_grid(span,10);
    AR_test_wing_base.alpha0 = LLT::const_function(airfoil.alpha0);
    AR_test_wing_base.m = LLT::const_function(airfoil.m);

    Wing AR_3 = Wing(AR_test_wing_base);
    float AR_3_span = std::sqrt(AR_testing_wing_area*3);
    float AR_3_chord = std::sqrt(AR_testing_wing_area/3);
    AR_3.span = AR_3_span;
    AR_3.chord = LLT::const_function(AR_3_chord);
    AR_3.alpha = get_required_alpha(AR_3, lift_required, velocity);


    Eigen::VectorXf AR_3_A = solve_system(AR_3);
    Results AR_3_results = calculate_results(AR_3, AR_3_A);

    std::cout << "AR_3 cl = " << AR_3_results.Cl << " E = " << AR_3_results.E << '\n';
    std::cout << "Area = " << LLT::wing_area(AR_3.span ,AR_3.chord) << '\n';
    float AR_3_glide_slope = AR_3_results.Cl/AR_3_results.Cd;
    std::cout << "cl/cd = " << AR_3_glide_slope << '\n';

    
    Wing AR_4 = Wing(AR_test_wing_base);
    float AR_4_span = std::sqrt(AR_testing_wing_area*4);
    float AR_4_chord = std::sqrt(AR_testing_wing_area/4);
    AR_4.span = AR_4_span;
    AR_4.chord = LLT::const_function(AR_4_chord);
    AR_4.alpha = get_required_alpha(AR_4, lift_required, velocity);

    Eigen::VectorXf AR_4_A = solve_system(AR_4);
    Results AR_4_results = calculate_results(AR_4, AR_4_A);

    std::cout << "AR_4 cl = " << AR_4_results.Cl << " E = " << AR_4_results.E << '\n';
    std::cout << "Area = " << LLT::wing_area(AR_4.span ,AR_4.chord) << '\n';
    float AR_4_glide_slope = AR_4_results.Cl/AR_4_results.Cd;
    std::cout << "cl/cd = " << AR_4_glide_slope << '\n';


    Wing AR_5 = Wing(AR_test_wing_base);
    float AR_5_span = std::sqrt(AR_testing_wing_area*5);
    float AR_5_chord = std::sqrt(AR_testing_wing_area/5);
    AR_5.span = AR_5_span;
    AR_5.chord = LLT::const_function(AR_5_chord);
    AR_5.alpha = get_required_alpha(AR_5, lift_required, velocity);

    Eigen::VectorXf AR_5_A = solve_system(AR_5);
    Results AR_5_results = calculate_results(AR_5, AR_5_A);

    std::cout << "AR_5 cl = " << AR_5_results.Cl << " E = " << AR_5_results.E << '\n';
    std::cout << "Area = " << LLT::wing_area(AR_5.span ,AR_5.chord) << '\n';
    float AR_5_glide_slope = AR_5_results.Cl/AR_5_results.Cd;
    std::cout << "cl/cd = " << AR_5_glide_slope << '\n';


    Wing AR_6 = Wing(AR_test_wing_base);
    float AR_6_span = std::sqrt(AR_testing_wing_area*6);
    float AR_6_chord = std::sqrt(AR_testing_wing_area/6);
    AR_6.span = AR_6_span;
    AR_6.chord = LLT::const_function(AR_6_chord);
    AR_6.alpha = get_required_alpha(AR_6, lift_required, velocity);

    Eigen::VectorXf AR_6_A = solve_system(AR_6);
    Results AR_6_results = calculate_results(AR_6, AR_6_A);

    std::cout << "AR_6 cl = " << AR_6_results.Cl << " E = " << AR_6_results.E << '\n';
    std::cout << "Area = " << LLT::wing_area(AR_6.span ,AR_6.chord) << '\n';
    float AR_6_glide_slope = AR_6_results.Cl/AR_6_results.Cd;
    std::cout << "cl/cd = " << AR_6_glide_slope << '\n';


    Wing AR_7 = Wing(AR_test_wing_base);
    float AR_7_span = std::sqrt(AR_testing_wing_area*7);
    float AR_7_chord = std::sqrt(AR_testing_wing_area/7);
    AR_7.span = AR_7_span;
    AR_7.chord = LLT::const_function(AR_7_chord);
    AR_7.alpha = get_required_alpha(AR_7, lift_required, velocity);

    Eigen::VectorXf AR_7_A = solve_system(AR_7);
    Results AR_7_results = calculate_results(AR_7, AR_7_A);

    std::cout << "AR_7 cl = " << AR_7_results.Cl << " E = " << AR_7_results.E << '\n';
    std::cout << "Area = " << LLT::wing_area(AR_7.span ,AR_7.chord) << '\n';
    float AR_7_glide_slope = AR_7_results.Cl/AR_7_results.Cd;
    std::cout << "cl/cd = " << AR_7_glide_slope << '\n';


    std::vector<float> glide_slopes = {
        AR_3_glide_slope,
        AR_4_glide_slope,
        AR_5_glide_slope,
        AR_6_glide_slope,
        AR_7_glide_slope
    };

    std::vector<float> aspect_ratios = {
        3,4,5,6,7
    };

    matplot::figure();
    matplot::plot(aspect_ratios, glide_slopes);
    matplot::title("Glide Ratio/Aspec Ratio");
    matplot::xlabel("Aspect Ratio");
    matplot::ylabel("Glide Slope");

}

void planform_select(float span, const airfoil_parameters& airfoil, float lift_required, float velocity) {
    
    float wing_area = span*span/7; // using aspect ratio of 7
    
    // Taper ratio 1 (rectangular wing)
    Wing TR_test_wing_base = {};
    TR_test_wing_base.theta = LLT::theta(span);
    TR_test_wing_base.span = span;
    TR_test_wing_base.y = LLT::uniform_grid(span,10);
    TR_test_wing_base.alpha0 = LLT::const_function(airfoil.alpha0);
    TR_test_wing_base.m = LLT::const_function(airfoil.m);
    

    Wing TR_1 = Wing(TR_test_wing_base);
    float TR_1_root_chord = 2*wing_area/(span*(1+1)); // (2*area)/(s*(1+taper_ratio))
    float TR_1_tip_chord = 2*wing_area*1/(span*(1+1)); // (2*area*taper_ratio)/(s*(1+taper_ratio))
    TR_1.chord = LLT::linear_function(TR_1_root_chord, TR_1_tip_chord, span);
    TR_1.alpha = get_required_alpha(TR_1, lift_required, velocity);

    Eigen::VectorXf TR_1_A = solve_system(TR_1);
    Results TR_1_results = calculate_results(TR_1, TR_1_A);

    std::cout << "TR_1 cl = " << TR_1_results.Cl << " E = " << TR_1_results.E << '\n';
    std::cout << "Area = " << LLT::wing_area(TR_1.span ,TR_1.chord) << '\n';
    float TR_1_glide_slope = TR_1_results.Cl/TR_1_results.Cd;
    std::cout << "cl/cd = " << TR_1_glide_slope << '\n';

    //Taper Ratio 0.8
    Wing TR_0_8 = Wing(TR_test_wing_base);
    float TR_0_8_root_chord = 2*wing_area/(span*(1+0.8)); // (2*area)/(s*(1+taper_ratio))
    float TR_0_8_tip_chord = 2*wing_area*0.8/(span*(1+0.8)); // (2*area*taper_ratio)/(s*(1+taper_ratio))
    TR_0_8.chord = LLT::linear_function(TR_0_8_root_chord, TR_0_8_tip_chord, span);
    TR_0_8.alpha = get_required_alpha(TR_0_8, lift_required, velocity);

    Eigen::VectorXf TR_0_8_A = solve_system(TR_0_8);
    Results TR_0_8_results = calculate_results(TR_0_8, TR_0_8_A);

    std::cout << "TR_0_8 cl = " << TR_0_8_results.Cl << " E = " << TR_0_8_results.E << '\n';
    std::cout << "Area = " << LLT::wing_area(TR_0_8.span ,TR_0_8.chord) << '\n';
    float TR_0_8_glide_slope = TR_0_8_results.Cl/TR_0_8_results.Cd;
    std::cout << "cl/cd = " << TR_0_8_glide_slope << '\n';


    //Taper Ratio 0.6
    Wing TR_0_6 = Wing(TR_test_wing_base);
    float TR_0_6_root_chord = 2*wing_area/(span*(1+0.6)); // (2*area)/(s*(1+taper_ratio))
    float TR_0_6_tip_chord = 2*wing_area*0.6/(span*(1+0.6)); // (2*area*taper_ratio)/(s*(1+taper_ratio))
    TR_0_6.chord = LLT::linear_function(TR_0_6_root_chord, TR_0_6_tip_chord, span);
    TR_0_6.alpha = get_required_alpha(TR_0_6, lift_required, velocity);

    Eigen::VectorXf TR_0_6_A = solve_system(TR_0_6);
    Results TR_0_6_results = calculate_results(TR_0_6, TR_0_6_A);

    std::cout << "TR_0_6 cl = " << TR_0_6_results.Cl << " E = " << TR_0_6_results.E << '\n';
    std::cout << "Area = " << LLT::wing_area(TR_0_6.span ,TR_0_6.chord) << '\n';
    float TR_0_6_glide_slope = TR_0_6_results.Cl/TR_0_6_results.Cd;
    std::cout << "cl/cd = " << TR_0_6_glide_slope << '\n';


    //Taper Ratio 0.5
    Wing TR_0_5 = Wing(TR_test_wing_base);
    float TR_0_5_root_chord = 2*wing_area/(span*(1+0.5)); // (2*area)/(s*(1+taper_ratio))
    float TR_0_5_tip_chord = 2*wing_area*0.5/(span*(1+0.5)); // (2*area*taper_ratio)/(s*(1+taper_ratio))
    TR_0_5.chord = LLT::linear_function(TR_0_5_root_chord, TR_0_5_tip_chord, span);
    TR_0_5.alpha = get_required_alpha(TR_0_5, lift_required, velocity);

    Eigen::VectorXf TR_0_5_A = solve_system(TR_0_5);
    Results TR_0_5_results = calculate_results(TR_0_5, TR_0_5_A);

    std::cout << "TR_0_5 cl = " << TR_0_5_results.Cl << " E = " << TR_0_5_results.E << '\n';
    std::cout << "Area = " << LLT::wing_area(TR_0_5.span ,TR_0_5.chord) << '\n';
    float TR_0_5_glide_slope = TR_0_5_results.Cl/TR_0_5_results.Cd;
    std::cout << "cl/cd = " << TR_0_5_glide_slope << '\n';

    //Taper Ratio 0.2
    Wing TR_0_2 = Wing(TR_test_wing_base);
    float TR_0_2_root_chord = 2*wing_area/(span*(1+0.2)); // (2*area)/(s*(1+taper_ratio))
    float TR_0_2_tip_chord = 2*wing_area*0.2/(span*(1+0.2)); // (2*area*taper_ratio)/(s*(1+taper_ratio))
    TR_0_2.chord = LLT::linear_function(TR_0_2_root_chord, TR_0_2_tip_chord, span);
    TR_0_2.alpha = get_required_alpha(TR_0_2, lift_required, velocity);

    Eigen::VectorXf TR_0_2_A = solve_system(TR_0_2);
    Results TR_0_2_results = calculate_results(TR_0_2, TR_0_2_A);

    std::cout << "TR_0_2 cl = " << TR_0_2_results.Cl << " E = " << TR_0_2_results.E << '\n';
    std::cout << "Area = " << LLT::wing_area(TR_0_2.span ,TR_0_2.chord) << '\n';
    float TR_0_2_glide_slope = TR_0_2_results.Cl/TR_0_2_results.Cd;
    std::cout << "cl/cd = " << TR_0_2_glide_slope << '\n';


    //Taper Ratio 0
    Wing TR_0 = Wing(TR_test_wing_base);
    float TR_0_root_chord = 2*wing_area/(span*(1+0)); // (2*area)/(s*(1+taper_ratio))
    float TR_0_tip_chord = 2*wing_area*0/(span*(1+0)); // (2*area*taper_ratio)/(s*(1+taper_ratio))
    TR_0.chord = LLT::linear_function(TR_0_root_chord, TR_0_tip_chord, span);
    TR_0.alpha = get_required_alpha(TR_0, lift_required, velocity);

    Eigen::VectorXf TR_0_A = solve_system(TR_0);
    Results TR_0_results = calculate_results(TR_0, TR_0_A);

    std::cout << "TR_0 cl = " << TR_0_results.Cl << " E = " << TR_0_results.E << '\n';
    std::cout << "Area = " << LLT::wing_area(TR_0.span ,TR_0.chord) << '\n';
    float TR_0_glide_slope = TR_0_results.Cl/TR_0_results.Cd;
    std::cout << "cl/cd = " << TR_0_glide_slope << '\n';

    //Elliptical
    Wing EL = Wing(TR_test_wing_base);
    float EL_root_chord = 2*wing_area/(LLT::pi*span/2);
    EL.chord = LLT::elliptical_function(EL_root_chord, span);
    EL.alpha = get_required_alpha(EL, lift_required, velocity);

    Eigen::VectorXf EL_A = solve_system(EL);
    Results EL_results = calculate_results(EL, EL_A);

    std::cout << "EL cl = " << EL_results.Cl << " E = " << EL_results.E << '\n';
    std::cout << "Area = " << LLT::wing_area(EL.span, EL.chord) << '\n';
    float EL_glide_slope = EL_results.Cl/EL_results.Cd;
    std::cout << "cl/cd = " << EL_glide_slope << '\n';

    std::vector<float> glide_slopes = {
        TR_0_glide_slope,
        TR_0_2_glide_slope,
        TR_0_5_glide_slope,
        TR_0_6_glide_slope,
        TR_0_8_glide_slope,
        TR_1_glide_slope
    };

    std::vector<float> taper_ratios = {
        0.0, 0.2, 0.5, 0.6, 0.8, 1.0
    };

    matplot::figure();
    matplot::plot(taper_ratios, glide_slopes);
    matplot::title("Glide Slope/Taper Ratio");
    matplot::xlabel("Taper Ratio");
    matplot::ylabel("Glide Slope");
}

void twist_select(float span, const airfoil_parameters &airfoil, float lift_required, float velocity) {

    float wing_area = span*span/7;

    Wing TW_test_wing_base = {};
    TW_test_wing_base.theta = LLT::theta(span);
    TW_test_wing_base.span = span;
    TW_test_wing_base.y = LLT::uniform_grid(span,10);
    TW_test_wing_base.alpha0 = LLT::const_function(airfoil.alpha0);
    TW_test_wing_base.m = LLT::const_function(airfoil.m);
    TW_test_wing_base.chord = LLT::const_function(wing_area/span);

    Wing TW_0 = Wing(TW_test_wing_base);
    TW_0.alpha = LLT::const_function(0);
    TW_0.alpha = get_required_alpha(TW_0, lift_required, velocity);

    Eigen::VectorXf TW_0_A = solve_system(TW_0);
    Results TW_0_results = calculate_results(TW_0, TW_0_A);

    std::cout << "Twist 0 cl = " << TW_0_results.Cl << " E = " << TW_0_results.E << '\n';
    float TW_0_glide_slope = TW_0_results.Cl/TW_0_results.Cd;
    std::cout << "cl/cd = " << TW_0_glide_slope << '\n';


    Wing TW_1 = Wing(TW_test_wing_base);
    TW_1.alpha = LLT::linear_function(1,0,span);
    TW_1.alpha = get_required_alpha(TW_1, lift_required, velocity);

    Eigen::VectorXf TW_1_A = solve_system(TW_1);
    Results TW_1_results = calculate_results(TW_1, TW_1_A);

    std::cout << "Twist 1 cl = " << TW_1_results.Cl << " E = " << TW_1_results.E << '\n';
    float TW_1_glide_slope = TW_1_results.Cl/TW_1_results.Cd;
    std::cout << "cl/cd = " << TW_1_glide_slope << '\n';

    
    Wing TW_2 = Wing(TW_test_wing_base);
    TW_2.alpha = LLT::linear_function(2,0,span);
    TW_2.alpha = get_required_alpha(TW_2, lift_required, velocity);

    Eigen::VectorXf TW_2_A = solve_system(TW_2);
    Results TW_2_results = calculate_results(TW_2, TW_2_A);

    std::cout << "Twist 2 cl = " << TW_2_results.Cl << " E = " << TW_2_results.E << '\n';
    float TW_2_glide_slope = TW_2_results.Cl/TW_2_results.Cd;
    std::cout << "cl/cd = " << TW_2_glide_slope << '\n';

    Wing TW_neg_1 = Wing(TW_test_wing_base);
    TW_neg_1.alpha = LLT::linear_function(-1,0,span);
    TW_neg_1.alpha = get_required_alpha(TW_neg_1, lift_required, velocity);

    Eigen::VectorXf TW_neg_1_A = solve_system(TW_neg_1);
    Results TW_neg_1_results = calculate_results(TW_neg_1, TW_neg_1_A);

    std::cout << "Twist -1 cl = " << TW_neg_1_results.Cl << " E = " << TW_neg_1_results.E << '\n';
    float TW_neg_1_glide_slope = TW_neg_1_results.Cl/TW_neg_1_results.Cd;
    std::cout << "cl/cd = " << TW_neg_1_glide_slope << '\n';


    Wing TW_neg_2 = Wing(TW_test_wing_base);
    TW_neg_2.alpha = LLT::linear_function(-2,0,span);
    TW_neg_2.alpha = get_required_alpha(TW_neg_2, lift_required, velocity);

    Eigen::VectorXf TW_neg_2_A = solve_system(TW_neg_2);
    Results TW_neg_2_results = calculate_results(TW_neg_2, TW_neg_2_A);

    std::cout << "Twist -2 cl = " << TW_neg_2_results.Cl << " E = " << TW_neg_2_results.E << '\n';
    float TW_neg_2_glide_slope = TW_neg_2_results.Cl/TW_neg_2_results.Cd;
    std::cout << "cl/cd = " << TW_neg_2_glide_slope << '\n';


    // Elliptical twist?
    Wing EL_TW = Wing(TW_test_wing_base);
    EL_TW.alpha = LLT::elliptical_function(1, span);
    EL_TW.alpha = get_required_alpha(EL_TW, lift_required, velocity);

    Eigen::VectorXf EL_TW_A = solve_system(EL_TW);
    Results EL_TW_results = calculate_results(EL_TW, EL_TW_A);

    std::cout << "Twist Elyptical cl = " << EL_TW_results.Cl << " E = " << EL_TW_results.E << '\n';
    float EL_TW_glide_slope = EL_TW_results.Cl/EL_TW_results.Cd;
    std::cout << "cl/cd = " << EL_TW_glide_slope << '\n';

    //Elliptical
    Wing EL = Wing(TW_test_wing_base);
    float EL_root_chord = 2*wing_area/(LLT::pi*span/2);
    EL.chord = LLT::elliptical_function(EL_root_chord, span);
    EL.alpha = get_required_alpha(EL, lift_required, velocity);

    Eigen::VectorXf EL_A = solve_system(EL);
    Results EL_results = calculate_results(EL, EL_A);

    std::cout << "EL cl = " << EL_results.Cl << " E = " << EL_results.E << '\n';
    std::cout << "Area = " << LLT::wing_area(EL.span, EL.chord) << '\n';
    float EL_glide_slope = EL_results.Cl/EL_results.Cd;
    std::cout << "cl/cd = " << EL_glide_slope << '\n';


    std::vector<float> glide_ratio = {
        TW_neg_2_glide_slope,
        TW_neg_1_glide_slope,
        TW_0_glide_slope,
        TW_1_glide_slope,
        TW_2_glide_slope
    };

    std::vector<float> Twist = {
        -2, -1, 0, 1, 2
    };

    matplot::figure();
    matplot::plot(Twist, glide_ratio);
    matplot::title("Glide Ratio/Twist");
    matplot::xlabel("Twist (Root angle - Tip angle) (degrees)");
    matplot::ylabel("Glide Ratio");

    /*matplot::figure();
    matplot::hold(matplot::on);
    auto gamma0 = calculate_gamma0(TW_0_A, TW_0);
    matplot::plot(TW_test_wing_base.y, gamma0);
    gamma0 = calculate_gamma0(TW_1_A, TW_1);
    //matplot::plot(TW_test_wing_base.y, gamma0);
    gamma0 = calculate_gamma0(TW_2_A, TW_2);
    //matplot::plot(TW_test_wing_base.y, gamma0);
    gamma0 = calculate_gamma0(TW_neg_2_A, TW_neg_2);
    //matplot::plot(TW_test_wing_base.y, gamma0);
    gamma0 = calculate_gamma0(EL_TW_A, EL_TW);
    matplot::plot(TW_test_wing_base.y, gamma0);
    gamma0 = calculate_gamma0(EL_A, EL);
    matplot::plot(TW_test_wing_base.y, gamma0);
    matplot::show();*/
    
}

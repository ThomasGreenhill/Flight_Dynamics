/*  Test for aircraft_eom.cpp
 *  
 *      Contents:
 *  
 *      History:
 *      
 */

#include "../Utilities/rk4.cpp"
#include "../Figure_Formatting/plotting_utils.cpp"
#include "aircraft_eom.cpp"

int main()
{
    auto t_start = std::chrono::high_resolution_clock::now();
    // Jiffy Jerboa in hover
    std::vector<double> tspan = utils::linspace(0, 1, 101);
    std::vector<double> x0 = {0.000001, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    std::vector<double> x;
    x = rk4(JJ_hover_nc, tspan, x0);

    std::vector<double> Vt(tspan.size()), alp(tspan.size()), bet(tspan.size()), Phi(tspan.size()), The(tspan.size()), Psi(tspan.size());
    std::vector<double> p(tspan.size()), q(tspan.size()), r(tspan.size()), xE(tspan.size()), yE(tspan.size()), hE(tspan.size());

    // Unpack the states

    for (int ii = 0; ii < tspan.size(); ii++)
    {
        Vt[ii] = x[ii * 12];
        alp[ii] = x[ii * 12 + 1];
        bet[ii] = x[ii * 12 + 2];
        Phi[ii] = x[ii * 12 + 3];
        The[ii] = x[ii * 12 + 4];
        Psi[ii] = x[ii * 12 + 5];
        p[ii] = x[ii * 12 + 6];
        q[ii] = x[ii * 12 + 7];
        r[ii] = x[ii * 12 + 8];
        xE[ii] = x[ii * 12 + 9];
        yE[ii] = x[ii * 12 + 10];
        hE[ii] = x[ii * 12 + 11];
    }

    auto t_stop = std::chrono::high_resolution_clock::now();
    utils::print_exectime(t_start, t_stop);

    int counter = 1;
    plt_utils::formatfigures();
    counter = plt_utils::named_plotls("Phi", tspan, Phi, "", counter);
    counter = plt_utils::named_plotls("The", tspan, The, "", counter);
    counter = plt_utils::named_plotls("Psi", tspan, Psi, "", counter);
    plt::legend();
    plt::figure();
    counter = plt_utils::plotls(tspan, xE, "", 1);

    plt::figure();
    counter = 1;
    counter = plt_utils::plotls(tspan, yE, "", counter);
    counter = plt_utils::plotls(tspan, hE, "", counter);

    plt::show();

    return 0;
}
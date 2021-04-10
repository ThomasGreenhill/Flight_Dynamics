/*  aircraft_eom.cpp
 *      A collection of equations of motion for aircraft for simulation and testing of control design.
 *      
 *      Contents:
 * 
 * 
 *      Notes:
 * 
 * 
 *      History: 
 *          
 */

#include <vector>
#include <iostream>
#include <math.h>
#include "../Utilities/utils.cpp"

// Function signatures
std::vector<double> JJ_hover_nc(double t, std::vector<double> x);
std::vector<double> JJ_hover_dynamics(std::vector<double> x, std::vector<double> ctrl);
std::vector<double> eom_6DOFnonlin(std::vector<double> x, double Cx, double Cy, double Cz, double Cl, double Cm, double Cn);

std::vector<double> JJ_hover_nc(double t, std::vector<double> x)
/*  JJ_hover_nc (nc indicates no controller)
     *      Nonlinear 6-DOF equations of motion for the Jiffy Jerboa in the hover configuration
     *      
     *      Notes:
     *          States vector is organized as follows:
     *              std::vector<double> x ({Vt, alp, bet, Phi, The, Psi, p, q, r, xE, yE, zE});
     * 
     *          Control vector is organized as follows (where each control is the % thrust of each motor)
     *              std::vector<double> ctrl ({d1, d2, d3, d4, d5, d6, d7, d8});
     * 
     *          All units are base SI
     * 
     *      History: 
     *          04.05.2021 Created, TVG
     */
{
    std::vector<double> ctrl;
    if (t < 25)
    {
        ctrl = utils::timesones(8, 0.755814);
        // ctrl[0] = 1;
        // ctrl[7] = 1;
    }
    else
    {
        ctrl = utils::timesones(8, -2);
        // ctrl = utils::timesones(8, 0.755814);
        // ctrl[0] = 1.2;
        // ctrl[7] = 1.2;
        // ctrl = utils::timesones(8, 0.755814);
    }

    std::vector<double> xdot = JJ_hover_dynamics(x, ctrl);

    return (xdot);
}

std::vector<double> JJ_hover_dynamics(std::vector<double> x, std::vector<double> ctrl)
/* JJ_hover_dynamics
        *      Hover dynamics of Jiffy Jerboa (plant dynamics), uses eom_6DOFnonlin
     *
     *      Inputs:
     *          x: state vector
     *          ctrl: control vector 
     * 
     *      States vector is organized as follows:
     *          std::vector<double> x ({Vt, alp, bet, Phi, The, Psi, p, q, r, xE, yE, zE});
     * 
     *      Control vector is organized as follows (where each control is the % thrust of each motor)
     *          std::vector<double> ctrl ({d1, d2, d3, d4, d5, d6, d7, d8});
     */
{
    // atmospheric and enviromental properties:
    double qbar = 0.5 * 1.225 * pow(x[0], 2);

    double qbarx = 0.5 * 1.225 * pow(x[0] * cos(x[1]) * cos(x[2]), 2);
    double qbary = 0.5 * 1.225 * pow(x[0] * sin(x[2]), 2);
    double qbarz = 0.5 * 1.225 * pow(x[0] * sin(x[1]) * cos(x[2]), 2);
    double g = 9.81;

    // geometric and inertial properties of aircraft
    double S = 16;
    double bspan = 12.65;
    double cbar = 1.2697;
    double m = 13000 / g; // kg
    // CG based on CG_locator.py
    std::vector<double> CG_xyz = {-3.5544829, 0, 0};

    // locations of motors relative to AVL origin (origin as nose tip)
    std::vector<double> d1_xyz = {-0.8, 1.3, 0};
    std::vector<double> d2_xyz = {-4.959, 5.915, -0.590};
    std::vector<double> d3_xyz = {-4.903, 3.770, -0.590};
    std::vector<double> d4_xyz = {-4.848, 1.620, -0.590};
    std::vector<double> d5_xyz = {-4.848, -1.620, -0.590};
    std::vector<double> d6_xyz = {-4.903, -3.770, -0.590};
    std::vector<double> d7_xyz = {-4.959, -5.915, -0.590};
    std::vector<double> d8_xyz = {-0.8, -1.3, 0};

    double d1 = ctrl[0];
    double d2 = ctrl[1];
    double d3 = ctrl[2];
    double d4 = ctrl[3];
    double d5 = ctrl[4];
    double d6 = ctrl[5];
    double d7 = ctrl[6];
    double d8 = ctrl[7];

    std::vector<double> xdot;

    // Motor max thrust and torque
    double Fmax = 2150; // N (from 'final_i/variable_pitch_analysis_trunc.png')
    double Tmax = 155;  // N.m

    // x-body axis force is only due to drag (CDmin assumed due to no wing lift induced drag)
    double CDx = 0; // 0.01991
    double Fx = -std::abs(qbarx * S * CDx / m);

    // y-body axis force is only due to drag
    double CDy = 0; // 0.4 Approximation for now, use FS or something later to get a better approx
    double Fy = -std::abs(qbary * S * CDy / m);

    // z-body axis force is just the sum of the motor thrusts minus the drag with V_inf in the -z direction
    double CDz = 0; // 0.7 Approximation for now, use FS or something later to get a better approx
    double Fz = -(d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8) * Fmax - std::abs(qbarz * S * CDz);

    // In moment calculations drag forces are assumed to act through CG and that y_cg = 0 (no offset from CG)
    double L = -(d1_xyz[1] * d1 * Fmax + d2_xyz[1] * d2 * Fmax + d3_xyz[1] * d3 * Fmax + d4_xyz[1] * d4 * Fmax + d5_xyz[1] * d5 * Fmax + d6_xyz[1] * d6 * Fmax + d7_xyz[1] * d7 * Fmax + d8_xyz[1] * d8 * Fmax);

    double M = (d1_xyz[0] - CG_xyz[0]) * d1 * Fmax + (d2_xyz[0] - CG_xyz[0]) * d2 * Fmax + (d3_xyz[0] - CG_xyz[0]) * d3 * Fmax + (d4_xyz[0] - CG_xyz[0]) * d4 * Fmax + (d5_xyz[0] - CG_xyz[0]) * d5 * Fmax + (d6_xyz[0] - CG_xyz[0]) * d6 * Fmax + (d7_xyz[0] - CG_xyz[0]) * d7 * Fmax + (d8_xyz[0] - CG_xyz[0]) * d8 * Fmax;

    // Assumed linear relationship between torque and thrust, BAD ASSUMPTION
    double N = -d1 * Tmax + d2 * Tmax - d3 * Tmax + d4 * Tmax - d5 * Tmax + d6 * Tmax - d7 * Tmax + d8 * Tmax;

    xdot = eom_6DOFnonlin(x, Fx, Fy, Fz, L, M, N);

    return xdot;
}

std::vector<double> eom_6DOFnonlin(std::vector<double> x, double Fx, double Fy, double Fz, double Ml, double Mm, double Mn)
/* eom_6DOFnonlin
     *      General 6DOF eom based on current state vetor and applied forces and moments 
     *      States vector is organized as follows:
     *          std::vector<double> x ({Vt, alp, bet, Phi, The, Psi, p, q, r, xE, yE, zE});
     *
     *      Inputs:
     *          x: state vector
     *          C{}: force coefficient applied along {}-body axis
     *          Cl: moment coefficient applied along x-body axis
     *          Cm: moment coefficient applied along y-body axis
     *          Cn: moment coefficient applied along z-body axis
     */
{
    // atmospheric and enviromental properties:
    double qbar = 0.5 * 1.225 * pow(x[0], 2);
    double g = 9.81;

    // geometric and inertial properties of aircraft
    double S = 16;
    double bspan = 12.65;
    double cbar = 1.2697;
    double m = 13000 / g; // kg
    double Ixx = 2353;    // kg.m^2
    double Iyy = 2327;    // kg.m^2
    double Izz = 3974;    // kg.m^2
    double Ixz = 1;       // Assumed (not yet calculated)

    // Unpack the state input vector
    double Vt = x[0];
    double alp = x[1];
    double bet = x[2];
    double Phi = x[3];
    double The = x[4];
    double Psi = x[5];
    double p = x[6];
    double q = x[7];
    double r = x[8];
    double xE = x[9];
    double yE = x[10];
    double zE = x[11];

    std::vector<double> xdot(12);

    // Body-axis velocities
    double u = Vt * cos(alp) * cos(bet);
    double v = Vt * sin(bet);
    double w = Vt * sin(alp) * cos(bet);

    // Body-axis accelerations
    double udot = r * v - q * w - g * sin(The) + Fx / m;
    double vdot = p * w - r * u + g * cos(The) * sin(Phi) + Fy / m;
    double wdot = q * u - p * v + g * cos(The) * cos(Phi) + Fz / m;
    std::cout << g * cos(The) * cos(Phi) + Fz / m << std::endl;

    // Total velocity and wind angles
    double Vtdot = (u * udot + v * vdot + w * wdot) / Vt;
    double alpdot = (u * wdot - w * udot) / (pow(u, 2) + pow(w, 2));
    double betdot = (Vt * vdot - v * Vtdot) / (pow(Vt, 2) * pow(1 - pow((v / Vt), 2), 0.5));

    // Some constants
    double c1 = ((Iyy - Izz) * Izz - pow(Ixz, 2)) / (Ixx * Izz - pow(Ixz, 2));
    double c2 = ((Ixx - Iyy + Izz) * Ixz) / (Ixx * Izz - pow(Ixz, 2));
    double c3 = Izz / (Ixx * Izz - pow(Ixz, 2));
    double c4 = Ixz / (Ixx * Izz - pow(Ixz, 2));
    double c5 = (Izz - Ixx) / Iyy;
    double c6 = Ixz / Iyy;
    double c7 = 1 / Iyy;
    double c8 = ((Ixx - Iyy) * Ixx - pow(Ixz, 2)) / (Ixx * Izz - pow(Ixz, 2));
    double c9 = Ixx / (Ixx * Izz - pow(Ixz, 2));

    // Angular accelerations (angular momentum of rotating masses are ignored)
    double pdot = (c1 * r + c2 * p) * q + c3 * Ml + c4 * Mn;
    double qdot = (c5 * p) * r - c6 * (pow(p, 2) - pow(r, 2)) + c7 * Mm;
    double rdot = (c8 * p - c2 * r) * q + c4 * Ml + c9 * Mn;

    // Navigation equations
    double Phidot = p + tan(The) * (q * sin(Phi) + r * cos(Phi));
    double Thedot = q * cos(Phi) - r * sin(Phi);
    double Psidot = (q * sin(Phi) + r * cos(Phi)) / cos(The);

    double xEdot = u * cos(Psi) * cos(The) + v * (cos(Psi) * sin(The) * sin(Phi) - sin(Psi) * cos(Phi)) + w * (cos(Psi) * sin(The) * cos(Phi) + sin(Psi) * sin(Phi));
    double yEdot = u * sin(Psi) * cos(The) + v * (sin(Psi) * sin(The) * sin(Phi) + cos(Psi) * cos(Phi)) + w * (sin(Psi) * sin(The) * cos(Phi) - cos(Psi) * sin(Phi));
    double hEdot = u * sin(The) - v * cos(The) * sin(Phi) - w * cos(The) * cos(Phi);

    xdot[0] = std::move(Vtdot);
    xdot[1] = std::move(alpdot);
    xdot[2] = std::move(betdot);
    xdot[3] = std::move(Phidot);
    xdot[4] = std::move(Thedot);
    xdot[5] = std::move(Psidot);
    xdot[6] = std::move(pdot);
    xdot[7] = std::move(qdot);
    xdot[8] = std::move(rdot);
    xdot[9] = std::move(xEdot);
    xdot[10] = std::move(yEdot);
    xdot[11] = std::move(hEdot);

    return (xdot);
}

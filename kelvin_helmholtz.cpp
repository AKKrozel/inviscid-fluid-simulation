#define _USE_MATH_DEFINES
#include "fluid_solver_2d.h"

int main() {
    int Nx = 1082; // number of cells in x-direction. two of the cells are reserved for ghost cells, the image x-resolution is Nx-2. 
    int Ny = 1082; // number of cells in y-direction. two of the cells are reserved for ghost cells, the image y-resolution is Ny-2. 
    double Lx = 1.0; // domain x-length
    double Ly = 1.0; // domain y-length
    double gamma = 5.0 / 3.0; // adiabatic index corresponding to ideal gas with no rotational or vibrational degrees of freedom
    double output_period = 0.01; // time period for data output
    double final_time = 10.0;
    double initial_time = 0.0;

    // parameter values used for initial conditions
    double a = 0.45;
    double b = 0.55;
    double rho_d = 4.0;
    double rho_0 = 1.0;
    double v_0 = 0.5;
    double P_0 = 2.5;
    double v_small = 0.02;
    double k = 6.0 * M_PI;
    double sigma = 0.05;

    fluid_solver_2d solver(Lx, Ly, Nx, Ny, gamma, output_period);

    std::vector<double> rho(Nx * Ny);
    std::vector<double> vx(Nx * Ny);
    std::vector<double> vy(Nx * Ny);
    std::vector<double> P(Nx * Ny);

    double x_step = Lx / (Nx - 2);
    double y_step = Ly / (Ny - 2);

    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            int idx = i + j * Nx;
            double x = (i - 1) * x_step;
            double y = (j - 1) * y_step;

            // initial values for primitive variables chosen to display kelvin-helmholtz instability patterns
            if (y >= a && y <= b) {
                rho[idx] = rho_d;
                vx[idx] = v_0;
            }
            else {
                rho[idx] = rho_0;
                vx[idx] = -v_0;
            }
            vy[idx] = 0.0;
            int layer_thickness = 8;
            double mult = layer_thickness / 2.0;
            if ((y > (a - mult * y_step) && y < (a + mult * y_step)) || (y > (b - mult * y_step) && y < (b + mult * y_step))) {
                double exp_0 = -(x - Lx/2.0) * (x - Lx/2.0) / (0.25);
                double coeff = v_small * std::sin(k * x) * std::exp(exp_0);
                double sigmaSQ = sigma * sigma;
                double exp_1 = -(y - a) * (y - a) / (sigmaSQ);
                double exp_2 = -(y - b) * (y - b) / (sigmaSQ);
                vy[idx] = coeff * (std::exp(exp_1) + std::exp(exp_2));
            }

            P[idx] = P_0;
        }
    }

    solver.init(rho, vx, vy, P);
    solver.solve(initial_time, final_time);

    return 0;

}
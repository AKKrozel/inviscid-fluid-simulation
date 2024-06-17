
#pragma once

#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>


double square(double x) { return x * x; }

class fluid_solver_2d {

public:

    fluid_solver_2d(double Lx0, double Ly0, int Nx0, int Ny0, double gamma0,
        double output_period0)
        : gamma(gamma0), Lx(Lx0), Ly(Ly0), Nx(Nx0), Ny(Ny0) {

        output_period = output_period0;
        dx = Lx / (Nx - 2);
        dy = Ly / (Ny - 2);

        rho.resize(Nx * Ny);
        vx.resize(Nx * Ny);
        vy.resize(Nx * Ny);
        P.resize(Nx * Ny);

        mass.resize(Nx * Ny);
        mom_x.resize(Nx * Ny);
        mom_y.resize(Nx * Ny);
        energy.resize(Nx * Ny);

        rho_tmp.resize(Nx * Ny);
        vx_tmp.resize(Nx * Ny);
        vy_tmp.resize(Nx * Ny);
        P_tmp.resize(Nx * Ny);

        rho_Lx.resize(Nx * Ny);
        rho_Rx.resize(Nx * Ny);
        rho_Ly.resize(Nx * Ny);
        rho_Ry.resize(Nx * Ny);

        vx_Lx.resize(Nx * Ny);
        vx_Rx.resize(Nx * Ny);
        vx_Ly.resize(Nx * Ny);
        vx_Ry.resize(Nx * Ny);
        vy_Lx.resize(Nx * Ny);
        vy_Rx.resize(Nx * Ny);
        vy_Ly.resize(Nx * Ny);
        vy_Ry.resize(Nx * Ny);
        P_Lx.resize(Nx * Ny);
        P_Rx.resize(Nx * Ny);
        P_Ly.resize(Nx * Ny);
        P_Ry.resize(Nx * Ny);

        mass_flux_x.resize(Nx * Ny);
        mass_flux_y.resize(Nx * Ny);
        momx_flux_x.resize(Nx * Ny);
        momx_flux_y.resize(Nx * Ny);
        momy_flux_x.resize(Nx * Ny);
        momy_flux_y.resize(Nx * Ny);
        energy_flux_x.resize(Nx * Ny);
        energy_flux_y.resize(Nx * Ny);
    }

    ~fluid_solver_2d() {}

    // HELPER FUNCTIONS ADDED //

    // helper function for finding the array index of a given cell
    int get_index(int i, int j) {
        return i + j * Nx;
    }

    // returns U based on primitive variables
    double U_n(double rho_n, double vx_n, double vy_n, double P_n) {
        double v_sq = square(vx_n) + square(vy_n);
        double u = P_n / ((gamma - 1) * rho_n);
        return rho_n * (0.5 * v_sq + u);
    }

    // calculates the special v term used in some functions
    double v_term(int n) {
        double Q1 = gamma * P[n] / rho[n];
        double Q2 = square(vx[n]) + square(vy[n]);
        return std::sqrt(Q1) + std::sqrt(Q2);
    }

    // END HELPER FUNCTIONS //

    void primitive_to_conserved() {
        // Compute conserved variables from primitive ones
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                int n = get_index(i, j);
                mass[n] = rho[n] * dx * dy;
                mom_x[n] = rho[n] * vx[n] * dx * dy;
                mom_y[n] = rho[n] * vy[n] * dx * dy;
                energy[n] = U_n(rho[n], vx[n], vy[n], P[n]) * dx * dy;
            }
        }
        // apply boundary conditions
        periodic_boundary(mass);
        periodic_boundary(mom_x);
        periodic_boundary(mom_y);
        periodic_boundary(energy);
    }

    void conserved_to_primitive() {
        // Compute primitive variables from conserved ones
        double u;
        double v_sq;
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                int n = get_index(i, j);
                rho[n] = mass[n] / (dx * dy);
                vx[n] = mom_x[n] / mass[n];
                vy[n] = mom_y[n] / mass[n];
                double U = energy[n] / (dx * dy);
                double v_sq = square(vx[n]) + square(vy[n]);
                double u = U / rho[n] - 0.5 * v_sq;
                P[n] = (gamma - 1) * rho[n] * u;
            }
        }
        periodic_boundary(rho);
        periodic_boundary(vx);
        periodic_boundary(vy);
        periodic_boundary(P);
    }

    void init(const std::vector<double>& rho0, const std::vector<double>& vx0,
        const std::vector<double>& vy0, const std::vector<double>& P0) {
        // Initialize the primitive variables using the given initial condition
        rho = rho0;
        vx = vx0;
        vy = vy0;
        P = P0;
        // apply boundary conditions
        periodic_boundary(rho);
        periodic_boundary(vx);
        periodic_boundary(vy);
        periodic_boundary(P);
        // update conserved variables
        primitive_to_conserved();
    }

    double find_dt() {
        // find the optimal dt that satisfies the CFL condition, and return its value
        double CFL = 0.1; // if CFL is too large, the simulation becomes unstable, smaller CFL increases runtime however. finding the best value of CFL is non-trivial

        // find max value of v
        double v_max = 0.0;
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                int n = get_index(i, j);
                v_max = std::max(v_max, v_term(n));
            }
        }
        double dt = CFL * std::min(dx, dy) / v_max;
        return dt;
    }

    void solve(double t0, double t_end) {
        // solve the fluid equations, starting from t0 and stoping at t_end
        double t = t0;
        int n = 0; // n labels the output file
        while (t < t_end) {
            if (t >= output_period * n) {
                output(n);
                n += 1;
            }
            // ensure data output files are not spaced more than output_period in time
            double dt = find_dt();
            if (dt > output_period) {
                dt = output_period;
            }
            step(dt);
            t += dt;
        }
    }

    void step(double dt) {
        // extrapolate a half step in time using primitive equations
        primitive_update(0.5 * dt);

        // compute fluxes
        compute_fluxes();

        // update solultion from fluxes
        update_conserved(dt);

        // update primitive variables
        conserved_to_primitive();
    }

    void periodic_boundary(std::vector<double>& f) {
        // apply periodic boundary conditions to array f
        for (int i = 1; i < Nx - 1; i++) {
            f[get_index(i, 0)] = f[get_index(i, Ny - 2)]; // top
            f[get_index(i, Ny - 1)] = f[get_index(i, 1)]; // bottom
        }
        for (int j = 1; j < Ny - 1; j++) {
            f[get_index(0, j)] = f[get_index(Nx - 2, j)]; // left
            f[get_index(Nx - 1, j)] = f[get_index(1, j)]; // right
        }
    }

    void primitive_update(double dt) {
        // update the primitive variables using Euler equations in primitive
        // form using an FTCS scheme
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                // get relevant indices
                int n = get_index(i, j);
                int nx_minus = get_index(i - 1, j);
                int nx_plus = get_index(i + 1, j);
                int ny_minus = get_index(i, j - 1);
                int ny_plus = get_index(i, j + 1);
                // compute partial derivatives
                double drho_dx = (rho[nx_plus] - rho[nx_minus]) / (2.0 * dx);
                double drho_dy = (rho[ny_plus] - rho[ny_minus]) / (2.0 * dy);
                double dvx_dx = (vx[nx_plus] - vx[nx_minus]) / (2.0 * dx);
                double dvx_dy = (vx[ny_plus] - vx[ny_minus]) / (2.0 * dy);
                double dvy_dx = (vy[nx_plus] - vy[nx_minus]) / (2.0 * dx);
                double dvy_dy = (vy[ny_plus] - vy[ny_minus]) / (2.0 * dy);
                double dP_dx = (P[nx_plus] - P[nx_minus]) / (2.0 * dx);
                double dP_dy = (P[ny_plus] - P[ny_minus]) / (2.0 * dy);
                // update temporary arrays
                rho_tmp[n] = rho[n] - dt * (vx[n] * drho_dx + vy[n] * drho_dy + rho[n] * (dvx_dx + dvy_dy));
                vx_tmp[n] = vx[n] - dt * (vx[n] * dvx_dx + vy[n] * dvx_dy + dP_dx) / rho[n];
                vy_tmp[n] = vy[n] - dt * (vx[n] * dvy_dx + vy[n] * dvy_dy + dP_dy) / rho[n];
                P_tmp[n] = P[n] - dt * (vx[n] * dP_dx + vy[n] * dP_dy + gamma * P[n] * (dvx_dx + dvy_dy));
            }
        }
        // assign updated values to actual arrays
        rho = rho_tmp;
        vx = vx_tmp;
        vy = vy_tmp;
        P = P_tmp;
        // apply boundary conditions
        periodic_boundary(rho);
        periodic_boundary(vx);
        periodic_boundary(vy);
        periodic_boundary(P);
        // convert
    }

    void extrapolate_to_interface() {
        // compute rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, P_R 
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                // get relevant indices
                int n = get_index(i, j);
                int nx_minus = get_index(i - 1, j);
                int nx_plus = get_index(i + 1, j);
                int ny_minus = get_index(i, j - 1);
                int ny_plus = get_index(i, j + 1);
                // extrapolate
                rho_Lx[n] = rho[n] - 0.25 * (rho[nx_plus] - rho[nx_minus]); 
                rho_Rx[n] = rho[n] + 0.25 * (rho[nx_plus] - rho[nx_minus]);
                rho_Ly[n] = rho[n] - 0.25 * (rho[ny_plus] - rho[ny_minus]);
                rho_Ry[n] = rho[n] + 0.25 * (rho[ny_plus] - rho[ny_minus]);
                vx_Lx[n] = vx[n] - 0.25 * (vx[nx_plus] - vx[nx_minus]); 
                vx_Rx[n] = vx[n] + 0.25 * (vx[nx_plus] - vx[nx_minus]);
                vx_Ly[n] = vx[n] - 0.25 * (vx[ny_plus] - vx[ny_minus]);
                vx_Ry[n] = vx[n] + 0.25 * (vx[ny_plus] - vx[ny_minus]);
                vy_Lx[n] = vy[n] - 0.25 * (vy[nx_plus] - vy[nx_minus]); 
                vy_Rx[n] = vy[n] + 0.25 * (vy[nx_plus] - vy[nx_minus]);
                vy_Ly[n] = vy[n] - 0.25 * (vy[ny_plus] - vy[ny_minus]);
                vy_Ry[n] = vy[n] + 0.25 * (vy[ny_plus] - vy[ny_minus]);
                P_Lx[n] = P[n] - 0.25 * (P[nx_plus] - P[nx_minus]); 
                P_Rx[n] = P[n] + 0.25 * (P[nx_plus] - P[nx_minus]);
                P_Ly[n] = P[n] - 0.25 * (P[ny_plus] - P[ny_minus]);
                P_Ry[n] = P[n] + 0.25 * (P[ny_plus] - P[ny_minus]);
            }
        }
        // apply boundary conditions
        periodic_boundary(rho_Lx);
        periodic_boundary(rho_Ly);
        periodic_boundary(rho_Rx);
        periodic_boundary(rho_Ry);
        periodic_boundary(vx_Lx); 
        periodic_boundary(vx_Ly);
        periodic_boundary(vx_Rx);
        periodic_boundary(vx_Ry);
        periodic_boundary(vy_Lx); 
        periodic_boundary(vy_Ly);
        periodic_boundary(vy_Rx);
        periodic_boundary(vy_Ry);
        periodic_boundary(P_Lx); 
        periodic_boundary(P_Ly);
        periodic_boundary(P_Rx);
        periodic_boundary(P_Ry);
    }

    void compute_fluxes() {
        // update cell boundaries
        extrapolate_to_interface();

        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                // get relevant indices
                int n = get_index(i, j);
                int nx_minus = get_index(i - 1, j);
                int nx_plus = get_index(i + 1, j);
                int ny_minus = get_index(i, j - 1);
                int ny_plus = get_index(i, j + 1);

                // calculate v_max
                double v_n = v_term(n);
                double v_nx_plus = v_term(nx_plus);
                double v_ny_plus = v_term(ny_plus);
                double v_max_x = std::max(v_n, v_nx_plus);
                double v_max_y = std::max(v_n, v_ny_plus);

                // mass x flux
                double F_L = rho_Rx[n] * vx_Rx[n];
                double F_R = rho_Lx[nx_plus] * vx_Lx[nx_plus];
                double Q_L = rho_Rx[n];
                double Q_R = rho_Lx[nx_plus];
                mass_flux_x[n] = 0.5 * (F_L + F_R) - (v_max_x / 2.0) * (Q_R - Q_L);

                // mass y flux
                F_L = rho_Ry[n] * vy_Ry[n];
                F_R = rho_Ly[ny_plus] * vy_Ly[ny_plus];
                Q_L = rho_Ry[n];
                Q_R = rho_Ly[ny_plus];
                mass_flux_y[n] = 0.5 * (F_L + F_R) - (v_max_y / 2.0) * (Q_R - Q_L);

                // momx x flux
                F_L = rho_Rx[n] * square(vx_Rx[n]) + P_Rx[n];
                F_R = rho_Lx[nx_plus] * square(vx_Lx[nx_plus]) + P_Lx[nx_plus];
                Q_L = rho_Rx[n] * vx_Rx[n];
                Q_R = rho_Lx[nx_plus] * vx_Lx[nx_plus];
                momx_flux_x[n] = 0.5 * (F_L + F_R) - (v_max_x / 2.0) * (Q_R - Q_L);

                // momx y flux
                F_L = rho_Ry[n] * vx_Ry[n] * vy_Ry[n];
                F_R = rho_Ly[ny_plus] * vx_Ly[ny_plus] * vy_Ly[ny_plus];
                Q_L = rho_Ry[n] * vx_Ry[n];
                Q_R = rho_Ly[ny_plus] * vx_Ly[ny_plus];
                momx_flux_y[n] = 0.5 * (F_L + F_R) - (v_max_y / 2.0) * (Q_R - Q_L);

                // momy x flux
                F_L = rho_Rx[n] * vx_Rx[n] * vy_Rx[n];
                F_R = rho_Lx[nx_plus] * vx_Lx[nx_plus] * vy_Lx[nx_plus];
                Q_L = rho_Rx[n] * vy_Rx[n];
                Q_R = rho_Lx[nx_plus] * vy_Lx[nx_plus];
                momy_flux_x[n] = 0.5 * (F_L + F_R) - (v_max_x / 2.0) * (Q_R - Q_L);

                // momy y flux
                F_L = rho_Ry[n] * square(vy_Ry[n]) + P_Ry[n];
                F_R = rho_Ly[ny_plus] * square(vy_Ly[ny_plus]) + P_Ly[ny_plus];
                Q_L = rho_Ry[n] * vy_Ry[n];
                Q_R = rho_Ly[ny_plus] * vy_Ly[ny_plus];
                momy_flux_y[n] = 0.5 * (F_L + F_R) - (v_max_y / 2.0) * (Q_R - Q_L);

                // energy x flux
                F_L = (U_n(rho_Rx[n], vx_Rx[n], vy_Rx[n], P_Rx[n]) + P_Rx[n]) * vx_Rx[n];
                F_R = (U_n(rho_Lx[nx_plus], vx_Lx[nx_plus], vy_Lx[nx_plus], P_Lx[nx_plus]) + P_Lx[nx_plus]) * vx_Lx[nx_plus];
                Q_L = U_n(rho_Rx[n], vx_Rx[n], vy_Rx[n], P_Rx[n]);
                Q_R = U_n(rho_Lx[nx_plus], vx_Lx[nx_plus], vy_Lx[nx_plus], P_Lx[nx_plus]);
                energy_flux_x[n] = 0.5 * (F_L + F_R) - (v_max_x / 2.0) * (Q_R - Q_L);

                // energy y flux
                F_L = (U_n(rho_Ry[n], vx_Ry[n], vy_Ry[n], P_Ry[n]) + P_Ry[n]) * vy_Ry[n];
                F_R = (U_n(rho_Ly[ny_plus], vx_Ly[ny_plus], vy_Ly[ny_plus], P_Ly[ny_plus]) + P_Ly[ny_plus]) * vy_Ly[ny_plus];
                Q_L = U_n(rho_Ry[n], vx_Ry[n], vy_Ry[n], P_Ry[n]);
                Q_R = U_n(rho_Ly[ny_plus], vx_Ly[ny_plus], vy_Ly[ny_plus], P_Ly[ny_plus]);
                energy_flux_y[n] = 0.5 * (F_L + F_R) - (v_max_y / 2.0) * (Q_R - Q_L);
            }
        }
        // apply boundary conditions
        periodic_boundary(mass_flux_x); 
        periodic_boundary(mass_flux_y);
        periodic_boundary(momx_flux_x); 
        periodic_boundary(momx_flux_y);
        periodic_boundary(momy_flux_x); 
        periodic_boundary(momy_flux_y);
        periodic_boundary(energy_flux_x); 
        periodic_boundary(energy_flux_y);
    }

    void update_conserved(double dt) {
        // update the conserved variables using the fluxes
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                // get relevant indices
                int n = get_index(i, j);
                int nx_minus = get_index(i - 1, j);
                int ny_minus = get_index(i, j - 1);
                // update
                mass[n] = mass[n] - (mass_flux_x[n] - mass_flux_x[nx_minus]) * dy * dt - (mass_flux_y[n] - mass_flux_y[ny_minus]) * dx * dt;
                mom_x[n] = mom_x[n] - (momx_flux_x[n] - momx_flux_x[nx_minus]) * dy * dt - (momx_flux_y[n] - momx_flux_y[ny_minus]) * dx * dt;
                mom_y[n] = mom_y[n] - (momy_flux_x[n] - momy_flux_x[nx_minus]) * dy * dt - (momy_flux_y[n] - momy_flux_y[ny_minus]) * dx * dt;
                energy[n] = energy[n] - (energy_flux_x[n] - energy_flux_x[nx_minus]) * dy * dt - (energy_flux_y[n] - energy_flux_y[ny_minus]) * dx * dt;
            }
        }
        // apply boundary conditions
        periodic_boundary(mass);
        periodic_boundary(mom_x);
        periodic_boundary(mom_y);
        periodic_boundary(energy);
    }

    void output(int n) {
        std::ofstream outfile(std::to_string(n) + ".csv");
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                int idx = i + j * Nx;
                outfile << rho[idx];
                if (i != Nx - 2)
                    outfile << ", ";
                else
                    outfile << std::endl;
            }
        }
        outfile.close();
    }

    int Nx, Ny;
    double Lx, Ly;
    double dx, dy;
    double gamma, output_period;
    std::vector<double> rho, vx, vy, P;                 // primitive variables
    std::vector<double> mass, mom_x, mom_y, energy;     // conserved variables
    // arrays to hold the results during primitive_update
    std::vector<double> rho_tmp, vx_tmp, vy_tmp, P_tmp;
    // arrays of fluxes for each conserved variable:
    std::vector<double> mass_flux_x, mass_flux_y;
    std::vector<double> momx_flux_x, momx_flux_y;
    std::vector<double> momy_flux_x, momy_flux_y;
    std::vector<double> energy_flux_x, energy_flux_y;
    // arrays for extrapolating to cell interfaces:
    std::vector<double> rho_Lx, rho_Ly, rho_Rx, rho_Ry;
    std::vector<double> vx_Lx, vx_Ly, vx_Rx, vx_Ry;
    std::vector<double> vy_Lx, vy_Ly, vy_Rx, vy_Ry;
    std::vector<double> P_Lx, P_Ly, P_Rx, P_Ry;
};

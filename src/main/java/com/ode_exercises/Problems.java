package com.ode_exercises;
import java.util.ArrayList;

import java.lang.Math;

public class Problems {
    // Computational Exercises for the U of C class M273: Basic Theory of Ordinary Differential Equations

    public static ArrayList<ArrayList<Double>> Problem1_10(int N[], String methodType) {
        // First Order Methods of Numerical Integration
        // System: x' = x w/ initial condition x(0) = 1

        double E = Math.E;
        ArrayList<ArrayList<Double>> result = new ArrayList<ArrayList<Double>>(3);
        // Col 1: The first column should contain the step size,
        // specifically the integer N such that ∆ = 1/N. Take N=2k for k ∈ [10,20].
        result.add(new ArrayList<Double>());
        // Col 2: The second column should contain the numerically computed error of x(1).
        result.add(new ArrayList<Double>());
        // Col 3: The third column should contain the size of the error to the actual value of e = 2.71828
        result.add(new ArrayList<Double>());
        
        for(int i = 0; i < N.length; i++) {
            double dt = 1.0 / N[i];
            double values[] = new double[N[i]+1];
            // initial condition
            values[0] = 1.0;
            
            for (int j = 1; j <= N[i]; j++) {
                double t_n = dt * j;
                if(methodType.toLowerCase().equals("euler")) {
                    values[j] = values[j-1] + dt * values[j-1];
                } else if(methodType.toLowerCase().equals("midpoint")) {
                    double k_temp = values[j-1] + (dt/2) * values[j-1];
                    values[j] = values[j-1] + dt * k_temp;
                } else if(methodType.toLowerCase().equals("runge-kutta")) {
                    double k_1 = values[j-1];
                    double k_2 = values[j-1] + (dt/2) * k_1;
                    double k_3 = values[j-1] + (dt/2) * k_2;
                    double k_4 = values[j-1] + dt * k_3;
                    values[j] = values[j-1] + (dt/6) * (k_1 + 2*k_2 + 2*k_3 + k_4);
                } else {
                    throw new Error("Not a recognized integration method.");
                }
            }
            result.get(0).add((double)N[i]);
            result.get(1).add(Math.abs(E - values[values.length-1]));
            result.get(2).add(Math.abs(E - values[values.length-1]) / E);
        }
        return result;
    }

    public static double[][] Problem7_8(
        double sigma, double beta, double rho,
        double dt, double T,
        double x0, double y0, double z0) {
        // Simulating the Lorenz System
        // x' = σ(y − x)
        // y' = x(ρ−z)−y
        // z' = xy − βz

        int steps = (int) (T/dt);

        // Coordinates x,y,z
        double[] x_seq = new double[steps+1];
        double[] y_seq = new double[steps+1];
        double[] z_seq = new double[steps+1];

        x_seq[0] = x0;
        y_seq[0] = y0;
        z_seq[0] = z0;

        // Use Runge-Kutta 4
        for (int t = 1; t <= steps; t++){
            double[] k_1 = f_lorenz(x_seq[t-1],y_seq[t-1],z_seq[t-1],sigma,rho,beta);
            double[] k_2 = f_lorenz(
                x_seq[t-1]+(dt/2)*k_1[0],
                y_seq[t-1]+(dt/2)*k_1[1],
                z_seq[t-1]+(dt/2)*k_1[2],
                sigma,rho,beta);
            double[] k_3 = f_lorenz(
                x_seq[t-1]+(dt/2)*k_2[0],
                y_seq[t-1]+(dt/2)*k_2[1],
                z_seq[t-1]+(dt/2)*k_2[2],
                sigma,rho,beta);
            double[] k_4 = f_lorenz(
                x_seq[t-1]+dt*k_3[0],
                y_seq[t-1]+dt*k_3[1],
                z_seq[t-1]+dt*k_3[2],
                sigma,rho,beta);
            x_seq[t] = x_seq[t-1] + (dt/6)*(k_1[0] + 2*k_2[0] + 2*k_3[0] + k_4[0]);
            y_seq[t] = y_seq[t-1] + (dt/6)*(k_1[1] + 2*k_2[1] + 2*k_3[1] + k_4[1]);
            z_seq[t] = z_seq[t-1] + (dt/6)*(k_1[2] + 2*k_2[2] + 2*k_3[2] + k_4[2]);
        }
        double[][] result = {x_seq,y_seq,z_seq};
        return (result);
    }
    public static double[][] Problem7_8(double sigma, double beta, double rho)
    {return Problem7_8(sigma,beta,rho,Math.pow(10,-6),10,0,0,0);}
    private static double[] f_lorenz(Double x, Double y, Double z, double sigma, double rho, double beta) {
        double[] result = {sigma*(y-x), (x*(rho-z)-y), x*y - beta*z};
        return(result);
    }


    public static void main(String[] args) {
        // Problem 1.10:
        int[] N = new int[11];
        for(int k = 10; k <= 20; k++) {
            N[k-10] = (int) Math.pow(2, k);
        }
        System.out.println(Problem1_10(N, "runge-kutta"));


        // double[][] Lorenz = Problem7_8(
        //     10.0, 8.0/3, 28.0,
        //     Math.pow(10,-2), 50.0,
        //     50, 0, 0);
        // // Rescale
        // for (int i = 0; i < 3; i++){
        //     for (int j = 0; j < Lorenz[0].length; j++){
        //         Lorenz[i][j] = Lorenz[i][j] / 20;
        //     }
        // }

        // // Visualize the trajectory
        // ODEVisualizer.visualize(Lorenz[0], Lorenz[1], Lorenz[2]);

        double[][] init_conditions = {{52,0,0},{50,0,0},{48,0,0},{52,2,0},{50,2,0},{48,2,0},{52,-2,0},{50,-2,0},{48,-2,0}};
        ArrayList<ArrayList<double[]>> trajectories = new ArrayList<>();
        for (int i = 0; i < init_conditions.length; i++){
            double[][] result = Problem7_8(
                10.0, 8.0/3, 28.0,
                Math.pow(10,-3), 50.0,
                init_conditions[i][0], init_conditions[i][1], init_conditions[i][2]);

            for (int a = 0; a < 3; a++){
                for (int b = 0; b < result[0].length; b++){
                    result[a][b] = result[a][b] / 20;
                }
            }

            ArrayList<double[]> particle = new ArrayList<>();
            particle.add(result[0]);
            particle.add(result[1]);
            particle.add(result[2]);

            trajectories.add(particle);
        }
        MultiParticleVisualizer.visualize(trajectories);
    }
}
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
    public static double[][] Problem2_8(
        double k, double l, double a, double b,
        double dt, double T,
        double x0, double y0,
        String methodType
        ) {
        // Problem 2.8 Lotka-Volterra System
        // x' = kx - axy
        // y' = -ly + bxy

        int steps = (int) (T/dt);
        // Coordinates x,y,z
        double[] x_seq = new double[steps+1];
        double[] y_seq = new double[steps+1];
 
        x_seq[0] = x0;
        y_seq[0] = y0;
        
        for(int i = 1; i <= steps; i++) {
            if(methodType.toLowerCase().equals("euler")) {
                double[] f = f_lotka(x_seq[i-1],y_seq[i-1],k,a,l,b);
                x_seq[i] = x_seq[i-1] + dt * f[0];
                y_seq[i] = y_seq[i-1] + dt * f[1];
            } else if(methodType.toLowerCase().equals("runge-kutta")) {
                double[] k1 = f_lotka(x_seq[i-1],y_seq[i-1],k,a,l,b);
                double[] k2 = f_lotka(x_seq[i-1] + (dt/2)*k1[0], y_seq[i-1] + (dt/2)*k1[1],k,a,l,b);
                double[] k3 = f_lotka(x_seq[i-1] + (dt/2)*k2[0], y_seq[i-1] + (dt/2)*k2[1],k,a,l,b);
                double[] k4 = f_lotka(x_seq[i-1] + dt*k3[0], y_seq[i-1] + dt*k3[1],k,a,l,b);
                x_seq[i] = x_seq[i-1] + (dt/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]);
                y_seq[i] = y_seq[i-1] + (dt/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
            } else {
                throw new Error("Not a recognized integration method.");
            }
        }

        double[][] result = {x_seq,y_seq};
        return (result);
    }
    public static double[] f_lotka(double x, double y, double k, double a, double l, double b) {
        double [] result = {k*x - a*x*y, b*x*y-l*y};
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

    public static double[][] Problem3_9(double x0, double v0, double mu, double dt, double T){
        // Problem 3.9: Van der Pol Oscillator

        // v' = mu(1-x^2)v - x
        // x' = v

        int steps = (int) (T/dt);

        // Coordinates x,v
        double[] x_seq = new double[steps+1];
        double[] v_seq = new double[steps+1];

        x_seq[0] = x0;
        v_seq[0] = v0;
        
        for (int t = 1; t <= steps; t++) {
            double[] f = f_vanderpol(x_seq[t-1], v_seq[t-1], mu);
            x_seq[t] = x_seq[t-1] + dt*f[0];
            v_seq[t] = v_seq[t-1] + dt*f[1];
        }

        double [][] result = {x_seq, v_seq};
        return result;
    }

    public static double[] f_vanderpol(double x, double v, double mu) {
        double [] result = {v,mu*(1-Math.pow(x,2))*v-x};
        return result;
    }


    public static void main(String[] args) {
        // Problem 1.10:
        int[] N = new int[11];
        for(int k = 10; k <= 20; k++) {
            N[k-10] = (int) Math.pow(2, k);
        }
        System.out.println(Problem1_10(N, "midpoint"));
        
        // // Problem 2.8: Lotka-Volterra System
        // double k = 1.0, l = 1.0, a = 1.0, b = 1.0;  // parameters
        // double T = 20.0;  // time span
        // double x0 = 2.0, y0 = 2.0;  // single initial condition
        
        // // Different step sizes to show numerical stability
        // double[] stepSizes = {0.1, 0.05, 0.01, 0.005, 0.001};

        // // Create collection of trajectories for visualization
        // ArrayList<ArrayList<double[]>> trajectories = new ArrayList<>();
        
        // // Compute trajectories for all step sizes with both methods
        // for (double dt : stepSizes) {
        //     // Euler method
        //     double[][] eulerResult = Problem2_8(k, l, a, b, dt, T, x0, y0, "euler");
        //     ArrayList<double[]> eulerTraj = new ArrayList<>();
        //     eulerTraj.add(eulerResult[0]);  // x coordinates
        //     eulerTraj.add(eulerResult[1]);  // y coordinates
        //     eulerTraj.add(new double[eulerResult[0].length]);  // z = 0 for 2D
        //     trajectories.add(eulerTraj);

        //     // Runge-Kutta method
        //     double[][] rkResult = Problem2_8(k, l, a, b, dt, T, x0, y0, "runge-kutta");
        //     ArrayList<double[]> rkTraj = new ArrayList<>();
        //     rkTraj.add(rkResult[0]);
        //     rkTraj.add(rkResult[1]);
        //     rkTraj.add(new double[rkResult[0].length]);
        //     trajectories.add(rkTraj);
        // }

        // // Visualize all trajectories
        // System.out.println("Lotka-Volterra Phase Space");
        // System.out.println("Initial condition: (x₀,y₀) = (2.0,2.0)");
        // System.out.println("Step sizes: 0.1, 0.05, 0.01, 0.005, 0.001");
        // System.out.println("Even indices: Euler method");
        // System.out.println("Odd indices: Runge-Kutta method");
        // MultiParticleVisualizer.visualize(trajectories);

        // // Problem 3.9: Van der Pol Oscillator

        // // Van der Pol oscillator with different μ values
        // double T = 50.0;  // longer time span to see limit cycles
        // double dt = 0.001;
        // double x0 = 1.0, v0 = 1.0;  // initial conditions
        
        // // Different μ values to compare
        // double[] muValues = {0.1, 0.5, 1.0, 2.0, 4.0};
        
        // // Create collection of trajectories for visualization
        // ArrayList<ArrayList<double[]>> trajectories = new ArrayList<>();
        
        // // Compute trajectories for different μ values
        // for (double mu : muValues) {
        //     double[][] result = Problem3_9(x0, v0, mu, dt, T);
            
        //     ArrayList<double[]> traj = new ArrayList<>();
        //     traj.add(result[0]);  // x coordinates
        //     traj.add(result[1]);  // v (dx/dt) coordinates
        //     traj.add(new double[result[0].length]);  // z = 0 for 2D phase space
        //     trajectories.add(traj);
        // }

        // // Visualize phase space
        // System.out.println("Van der Pol Oscillator Phase Space");
        // System.out.println("Initial condition: (x₀,v₀) = (1.0,1.0)");
        // System.out.println("μ values: 0.1, 0.5, 1.0, 2.0, 4.0");
        // System.out.println("Each color represents a different μ value");
        // System.out.println("Observe how the limit cycle changes with μ");
        // MultiParticleVisualizer.visualize(trajectories);


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

        // double[][] init_conditions = {{52,0,0},{50,0,0},{48,0,0},{52,2,0},{50,2,0},{48,2,0},{52,-2,0},{50,-2,0},{48,-2,0}};
        // ArrayList<ArrayList<double[]>> trajectories = new ArrayList<>();
        // for (int i = 0; i < init_conditions.length; i++){
        //     double[][] result = Problem7_8(
        //         10.0, 8.0/3, 28.0,
        //         Math.pow(10,-3), 50.0,
        //         init_conditions[i][0], init_conditions[i][1], init_conditions[i][2]);

        //     for (int a = 0; a < 3; a++){
        //         for (int b = 0; b < result[0].length; b++){
        //             result[a][b] = result[a][b] / 20;
        //         }
        //     }

        //     ArrayList<double[]> particle = new ArrayList<>();
        //     particle.add(result[0]);
        //     particle.add(result[1]);
        //     particle.add(result[2]);

        //     trajectories.add(particle);
        // }
        // // MultiParticleVisualizer.visualize(trajectories);
        

        // // Plot x-coordinates over time
        // ArrayList<ArrayList<double[]>> timeSeriesPlot = new ArrayList<>();
        // double dt = Math.pow(10,-3);
        // double T = 50.0;
        // int steps = (int) (T/dt);
        // double[] timePoints = new double[steps+1];
        // for (int i = 0; i <= steps; i++) {
        //     timePoints[i] = i * dt;
        // }

        // // Create time series plot for each trajectory
        // for (int i = 0; i < init_conditions.length; i++) {
        //     ArrayList<double[]> timeSeries = new ArrayList<>();
        //     timeSeries.add(timePoints);  // x-axis: time
        //     timeSeries.add(trajectories.get(i).get(0));  // y-axis: x-coordinate
        //     timeSeries.add(new double[timePoints.length]); // z-axis: 0 (for 2D visualization)
        //     timeSeriesPlot.add(timeSeries);
        // }

        // // Show the time series plot in a separate window
        // MultiParticleVisualizer.visualize(timeSeriesPlot);

        // // Points with large coordinates
        // double[][] large_init = {{100,100,100}, {-100,-100,100}, {100,-100,-100}};
        // ArrayList<ArrayList<double[]>> largeTrajectories = new ArrayList<>();

        // for (int i = 0; i < large_init.length; i++){
        //     double[][] result = Problem7_8(
        //         10.0, 8.0/3, 28.0,
        //         Math.pow(10,-3), 50.0,
        //         large_init[i][0], large_init[i][1], large_init[i][2]);

        //     // Scale down for visualization
        //     for (int a = 0; a < 3; a++){
        //         for (int b = 0; b < result[0].length; b++){
        //             result[a][b] = result[a][b] / 20;
        //         }
        //     }

        //     ArrayList<double[]> particle = new ArrayList<>();
        //     particle.add(result[0]);
        //     particle.add(result[1]);
        //     particle.add(result[2]);

        //     largeTrajectories.add(particle);
        // }
        // MultiParticleVisualizer.visualize(largeTrajectories);
    }
}
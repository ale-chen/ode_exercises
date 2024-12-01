package com.ode_exercises;
import java.util.ArrayList;
import java.lang.Math; 

public class Problems{
    // Computational Exercises for the U of C class M273: Basic Theory of Ordinary Differential Equations
    public static ArrayList<ArrayList<Double>> Problem1_10(int N[], String methodType){
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

        for(int i = 0; i < N.length; i++){
            double dt = 1.0 / N[i];
            double values[] = new double[N[i]+1];
            // initial condition
            values[0] = 1.0;
            for (int j = 1; j <= N[i]; j++){
                double t_n = dt * j;
                if(methodType.toLowerCase().equals("euler")){
                    // 1.10.1: Euler Method
                    values[j] = values[j-1] + dt * values[j-1];
                }else if(methodType.toLowerCase().equals("midpoint")){
                    // 1.10.2: Midpoint Method i.e. Verlet Integration Method
                    double k_temp = values[j-1] + (dt/2) * values[j-1];
                    values[j] = values[j-1] + dt * k_temp;
                }else if(methodType.toLowerCase().equals("runge-kutta")){
                    // 1.10.3: Runge-Kutta Method
                    double k_1 = values[j-1];
                    double k_2 = values[j-1] + (dt/2) * k_1;
                    double k_3 = values[j-1] + (dt/2) * k_2;
                    double k_4 = values[j-1] + dt * k_3;
                    values[j] = values[j-1] + (dt/6) * (k_1 + 2*k_2 + 2*k_3 + k_4);
                }else{
                    throw new Error("Not a recognized integration method.");
                }
            }
            result.get(0).add((double)N[i]);
            result.get(1).add(Math.abs(E - values[values.length-1]));
            result.get(2).add(Math.abs(E - values[values.length-1]) / E);
        }
        return result;
    }

    public static void main(String[] args) {
        // Problem 1.10:
        int[] N = new int[11];
        for(int k = 10; k <= 20; k++){
            N[k-10] = (int) Math.pow(2, k);
        }

        System.out.println(Problem1_10(N, "runge-kutta"));
        
    }
}
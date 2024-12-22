#include <cstdio>
#include <cstdlib>
#include <cmath>
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

int main(int argc, char* argv[]){
    const static double Epsilon = 1e-10;
    const static int N_of_Divisions = 1000;
    const static int N_max = 5;

    FILE *Wavefunction_file, *Energy_file, *Potential_file;
    
    Wavefunction_file = fopen("psi-osc.dat","w");
    Energy_file = fopen("E_n_Oszillator.dat","w");
    Potential_file = fopen("HarmonicPotentialNoDim.dat", "w");
    
    if (!(Wavefunction_file && Energy_file && Potential_file)) {
        printf("Problems to create files output.\n");
        exit(2);
    }

    const static double xRange = 12;
    const static double h_0 = xRange / N_of_Divisions;
    double* E_pot = new double[N_of_Divisions+1];
    double dist;

    for (int i = 0; i <= N_of_Divisions; ++i){
        dist = i*h_0 - 0.5*xRange;
        E_pot[i] = 0.5*dist*dist;
        fprintf(Potential_file, "%16.12e \t\t %16.12e\n", dist, E_pot[i]);
    }
    fclose(Potential_file);

    const static double Psi_left = 1.0e-3;
    const static double Psi_right = 0.0;

    double *Psi, *EigenEnergies;
    Psi = new double[N_of_Divisions+1];
    EigenEnergies = new double[N_max+1];
    Psi[0] = Psi_left;
    Psi[1] = Psi_left + 1.0e-3;

    int N_quantum;
    int Nodes_plus;

    double K_square;

    double E_lowerLimit = 0.0;
    double E_upperLimit = 10.0;
    int End_sign = -1;
    bool Limits_are_defined = false;
    double Normalization_coefficient;
    double E_trial;

    for(N_quantum=1;N_quantum <= N_max; ++N_quantum){
        Limits_are_defined = false;
        while (Limits_are_defined == false){
            Nodes_plus = 0;
            E_trial = E_upperLimit;
            for (int i=2; i <= N_of_Divisions; ++i){
                K_square = 2.0*(E_trial - E_pot[i]);
                Psi[i] = 2.0*Psi[i-1]*(1.0-(5.0*h_0*h_0*K_square / 12.0)) / (1.0 + (h_0*h_0*K_square/12.0))-Psi[i-2];
                if (Psi[i]*Psi[i-1] < 0.0){
                    ++Nodes_plus;
                }
            }
            if (E_upperLimit < E_lowerLimit){
                E_upperLimit = MAX(2*E_upperLimit, -2*E_upperLimit);
            }

            if (Nodes_plus > N_quantum){
                E_upperLimit *= 0.7;
            } else if (Nodes_plus < N_quantum){
                E_upperLimit *= 2.0;
            } else {
                Limits_are_defined = true;
            }
        }

        End_sign = -End_sign;
        while((E_upperLimit - E_lowerLimit) > Epsilon){
            E_trial = (E_upperLimit + E_lowerLimit) / 2.0;
            for (int i=2; i <= N_of_Divisions; i++){
                K_square = 2.0*(E_trial - E_pot[i]);
                Psi[i] = 2.0*Psi[i-1]*(1.0-(5.0*h_0*h_0*K_square / 12.0)) / (1.0 + (h_0*h_0*K_square/12.0))-Psi[i-2];
            }
            if (End_sign*Psi[N_of_Divisions] > Psi_right){
                E_lowerLimit = E_trial;
            } else {
                E_upperLimit = E_trial;
            }
        }
        E_trial = (E_upperLimit + E_lowerLimit)/2.0;
        EigenEnergies[N_quantum] = E_trial;
        E_upperLimit = E_trial;
        E_lowerLimit = E_trial;

        double Integral = 0.0;
        for (int i=1; i<= N_of_Divisions; ++i){
            Integral += 0.5*h_0*(Psi[i-1]*Psi[i-1]+Psi[i]*Psi[i]);
        }
        Normalization_coefficient = sqrt(1.0/Integral);

        for (int i=0; i <= N_of_Divisions; ++i){
            fprintf(Wavefunction_file, "%16.12e \t\t %16.12e\n", i*h_0-0.5*xRange, Normalization_coefficient*Psi[i]);
        }
        fprintf(Wavefunction_file, "\n");
    }
    fclose(Wavefunction_file);

    const static double Energyscale = 3.5;
    printf("Quantum Harmonic Oscillator, program osci.cpp\n");
    printf("Energies in MeV:\n");
    printf("n \t\t E_n\n");
    for (N_quantum=1; N_quantum <= N_max; ++N_quantum){
        fprintf(Energy_file, "%d \t\t %16.12e\n", N_quantum-1, Energyscale*EigenEnergies[N_quantum]);
        printf("%d \t\t %16.12e \n", N_quantum-1, Energyscale*EigenEnergies[N_quantum]);
    }
    
}
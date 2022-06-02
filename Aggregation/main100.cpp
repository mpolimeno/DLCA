#include "common100.h"
#include <cassert>
#include <random>
#include <iostream>
#include <ctime>

// This program simulates the formation mechanism of aggregates,
// accounting for both rotational and translational random motion
int main(int argc, char *argv[])
{
    /**************************************************************/
    time_t start, end;
    // First, let us check that parameter values are sensible
    // Currently the code only works for a 3D implementation -> We hard-coded the Orientation Matrix Q
    if (DIM!=3)
    {
        printf("ERROR: DIM = 3 expected\n");
        exit(0);
    }

    assert(argc==2);

    // SEED
    int s = atoi(argv[1]);
    srand(s);

    // Define static variables
    int nparcl[TOTNPARCL];                  // Array to hold the number of parcels present in a given aggregate
    int nagg = TOTNPARCL;                   // Totall number of aggregates in the system
    int Nleaders = TOTNPARCL;               // Integer to see how many aggregates there are at a given moment
    double Q_n[TOTNPARCL][DIM][DIM];        // Orientation matrix for Rotational Dynamics: one per aggregate
    double com_agg[TOTNPARCL][DIM];         // Center of Mass of each aggregate


    // At t=0 aggregates only have one particle each
    for (int n=0;n<nagg;n++) nparcl[n] = 1;

    // Initial position (lab frame) for each parcel in the aggregates
    double* xlab = new double[TOTNPARCL*MAXNPARCL*DIM]; // Position of aggregates in LAB frame of reference (dynamic allocation)
    init_aggs(xlab,nparcl,nagg);

    for (int n=0;n<nagg;n++)
    {
        double com_t0[DIM];
        int n_agg = n*nagg; // need this to extract the correct data from the pointer that holds the 3d array
        compute_com(xlab,n_agg,nparcl[n],com_t0);

        for (int d=0;d<DIM;d++) com_agg[n][d] = com_t0[d];
    }

    // Initialize Orientation Matrix Q
    for (int n=0;n<nagg;n++)
    {
        for(int i=0;i<DIM;i++)
        {
            // Initialize Q to the identity
            for(int j=0;j<DIM;j++)
            {
                Q_n[n][i][j] = (i==j) ? 1. : 0.;
            }
        }
    }

    // Initialize position in body frame (xrel) to zero vector
    double* xrel = new double[TOTNPARCL*MAXNPARCL*DIM]; // Position of aggregates in BODY frame of reference (dynamic allocation)
    for (int n=0;n<nagg;n++)
    {
        for(int i=0;i<nparcl[n];i++)
        {
            for(int j=0;j<DIM;j++)
            {
                *(xrel+n*nagg*DIM+i*DIM+j) = 0.;
            }
        }
    }

    // Arrays to store the radius of gyration and the maximum radius of each aggregate
    double max_out[MAXNPARCL];
    double gyr_out[MAXNPARCL];
    double maxrad[MAXNPARCL];
    double Grad[MAXNPARCL];


    // Store data at t = 0
    for (int n=0;n<nagg;n++)
    {
        if (nparcl[n]==0) continue;
        int n_agg = n*nagg;
        max_out[n] = find_max_rad(maxrad,n_agg,xlab,com_agg[n],nparcl[n]);
        gyr_out[n] = find_gyr_rad(Grad,n_agg,xlab,com_agg[n],nparcl[n]);
    }
    // random numbers to perform random translation and random rotation
    std::default_random_engine generator(s);
    std::normal_distribution<double> distribution(0.,1.);

    /**************************************************************/
    // Driver
    int countsteps = 0;
    int t_stop = 0;
    time(&start);
    while(Nleaders>NAGG && t_stop<=STOP)
    //while (Nleaders>NAGG)
    {
        /**************************************************************/
        // Update the data collected at each time-step
        for (int n=0;n<nagg;n++)
        {
            if (nparcl[n]==0) continue;
            // Save data for post-processing
            char filename[MAXFNLEN];
            sprintf(filename,"res%d_rotation_N_%d_Dt_%g",s,TOTNPARCL,Dt); // one output per seed
            FILE *fp = fopen(filename,"a");
            if ((countsteps<50000) && (countsteps%1000)==0)
            {
                fprintf(fp,"%d\t%d\t%g\t%g\t%d\n",n,nparcl[n],max_out[n],gyr_out[n],countsteps);
            }
            if ((countsteps>=50000) && (countsteps<=200000) && (countsteps%10000)==0)
            {
                fprintf(fp,"%d\t%d\t%g\t%g\t%d\n",n,nparcl[n],max_out[n],gyr_out[n],countsteps);
            }
            if ((countsteps>200000) && (countsteps%50000)==0)
            {
                fprintf(fp,"%d\t%d\t%g\t%g\t%d\n",n,nparcl[n],max_out[n],gyr_out[n],countsteps);
            }
            fclose(fp);
        }
        /**************************************************************/
        // DYNAMICS
        // FOR SANITY compute all the translation first and then all the rotation

        // TRANSLATION
        for (int n=0;n<nagg;n++)
        {
            if (nparcl[n]==0) continue;

            // Variable diffusion coefficient for translation
            double K_T = (nparcl[n]==1) ? D_T : D_T/gyr_out[n]; // Diff_Translation ~ 1/size

            // Array to hold "translational step"
            double dx[DIM];
            dx[0] = sqrt(2.*K_T*Dt)*distribution(generator);
            dx[1] = sqrt(2.*K_T*Dt)*distribution(generator);
            // Bias
            if (nparcl[n]>1)
            {
                double Udt = (C_u*(nparcl[n]/gyr_out[n]) - U_0)*Dt;
                dx[2] = sqrt(2.*K_T*Dt)*distribution(generator) - Udt; // variance stays the same. It is the mean that has a shift by U*t (advection diffusion)
            }
            else dx[2] = sqrt(2.*K_T*Dt)*distribution(generator);
            
            // Translate center of mass by dx accounting for periodic boundaries
            for (int d=0;d<DIM;d++)
            {
                com_agg[n][d] += dx[d];
                com_agg[n][d] -= BOXL*rint(com_agg[n][d]/BOXL);
            }
        }
        // ROTATION
        for (int n=0;n<nagg;n++)
        {
            // Skip rotation routine for aggregates less than 2 in size
            if (nparcl[n]<2) continue;

            // Variable diffusion coefficient for rotation
            double Rgcb = (gyr_out[n]*gyr_out[n]*gyr_out[n]);
            double K_R = D_R/Rgcb; // Diff_rotation ~ 1/size^3

            // Random principal angle for "rotational step"
            double dtheta[DIM];
            for(int d=0;d<DIM;d++)
            {
                dtheta[d] = sqrt(2.*K_R*Dt)*distribution(generator);
            }

            // Initialize antisymmetric matrix
            double A[DIM][DIM]=
            {
                {0., 0., 0.},
                {0., 0., 0.}
            };

            // Build antisymmetric matrix
            A[0][1] = dtheta[2];
            A[0][2] = -1.*dtheta[1];
            A[1][2] = dtheta[0];
            A[1][0] = -1.*A[0][1];
            A[2][0] = -1.*A[0][2];
            A[2][1] = -1.*A[1][2];
            // Diagonal elements are all zero
            A[0][0] = 0.;
            A[1][1] = 0.;
            A[2][2] = 0.;

            // To compute Frobenius Norm to be used for Rodrigues's formula
            double t0sq = dtheta[0]*dtheta[0];
            double t1sq = dtheta[1]*dtheta[1];
            double t2sq = dtheta[2]*dtheta[2];
            // Frobenius norm for antisymmetric matrix with real coeffs
            double v = sqrt(t0sq+t1sq+t2sq);

            // Identity matrix for Rodrigues's Formula
            double I[DIM][DIM];
            for(int i=0;i<DIM;i++)
            {
                for(int j=0;j<DIM;j++)  I[i][j] = (i==j) ? 1. : 0.;
            }

            // Build rotation matrix exp(A) = RotMat
            double Asq[DIM][DIM]=
            {
                {0., 0., 0.},
                {0., 0., 0.}
            };
            // For Rodrigues's Formula
            for (int i=0;i<DIM;i++)
            {
                for (int j=0;j<DIM;j++)
                {
                    for (int k=0;k<DIM;k++)
                    {
                        Asq[i][j] += A[i][k]*A[k][j];
                    }
                }
            }

            // Now Rodrigues's Formula
            double expA[DIM][DIM]=
            {
                {0., 0., 0.},
                {0., 0., 0.}
            };
            // Rodrigues formula for expm of antisymmetric matrix
            for (int i=0;i<DIM;i++)
            {
                for (int j=0;j<DIM;j++)
                {
                    expA[i][j] = I[i][j] + (sin(v)/v)*A[i][j] + ((1-cos(v))/(v*v))*Asq[i][j];
                }
            }

            // Update Q for each aggregate
            double Q_new[DIM][DIM]=
            {
                {0., 0., 0.},
                {0., 0., 0.}
            };
            double Q_old[DIM][DIM]=
            {
                {0., 0., 0.},
                {0., 0., 0.}
            };
            // Assign Q_n to Q_old = Q(t)
            for (int i=0;i<DIM;i++)
            {
                for (int j=0;j<DIM;j++)
                {
                    Q_old[i][j] = Q_n[n][i][j];
                }
            }
            // Compute Q_new = Q(t+Dt)
            for (int i=0;i<DIM;i++)
            {
                for (int j=0;j<DIM;j++)
                {
                    for (int k=0;k<DIM;k++)
                    {
                        Q_new[i][j] += expA[i][k]*Q_old[k][j];
                    }
                }
            }
            // Assign Q_new to Q_n
            for (int i=0;i<DIM;i++)
            {
                for (int j=0;j<DIM;j++) Q_n[n][i][j] = Q_new[i][j];
            }
        }
        // update step count
        countsteps++;
        t_stop++;
        /**************************************************************/
        // NOW particles have both translated and rotated within a Dt

        // Check if they got within attaching range and, if so, attach them
        // Compute distance between two aggregates (lab frame)
        for (int m=0;m<nagg;m++)
        {
            if (nparcl[m]==0) continue;
                int m_agg = m*nagg;
                compute_xlab(xrel,m_agg,com_agg[m],nparcl[m],Q_n[m],xlab);
                for (int n=m+1;n<nagg;n++)
                {
                    if (nparcl[n]==0) continue;
                        int n_agg = n*nagg;
                        compute_xlab(xrel,n_agg,com_agg[n],nparcl[n],Q_n[n],xlab);
                        double distsq = compute_distsq(xlab,m_agg,nparcl[m],n_agg,nparcl[n]);
                        // Combine aggregates that are within attaching range
                        if (distsq<=ATTACH*ATTACH)
                        {
                            combine_aggs(xlab,m_agg,nparcl[m],n_agg,nparcl[n]);
                            // Update aggregate count
                            Nleaders = Nleaders - 1;
                            // Reset Q_m to be the identity
                            for(int i=0;i<DIM;i++)
                            {
                                for(int j=0;j<DIM;j++)  Q_n[m][i][j] = (i==j) ? 1. : 0.;
                            }
                            // Update center of mass for new m-th aggregate
                            double com_mn[DIM];

                            compute_com(xlab,m_agg,nparcl[m],com_mn);
                            for (int d=0;d<DIM;d++)
                            {
                                com_agg[m][d] = com_mn[d];
                            }
                            // Update relative position of parcels in m-th aggregate
                            compute_xrel(xlab,m_agg,nparcl[m],com_agg[m],xrel);

                            // Update Radius of Gyration and Maximum Radius for the aggregate with lower index
                            max_out[m] = find_max_rad(maxrad,m_agg,xlab,com_agg[m],nparcl[m]);
                            gyr_out[m] = find_gyr_rad(Grad,m_agg,xlab,com_agg[m],nparcl[m]);
                        }
                }
        }
    }
    time(&end);

    /**************************************************************/
    // Dynamics is done
    // Compute data at the final time and append to file
    if (Nleaders==NAGG)
    {
        for(int n=0;n<nagg;n++)
        {
            if(nparcl[n]==0) continue;
            int n_agg = n*nagg;
            // final center of mass
            double com_tf[DIM];
            compute_com(xlab,n_agg,nparcl[n],com_tf);
            // save data for post processing
            double maxrad[MAXNPARCL];
            double Grad[MAXNPARCL];
            double max_radius = (nparcl[n]==1) ? R : find_max_rad(maxrad,n_agg,xlab,com_tf,nparcl[n]);
            double gyr_radius = (nparcl[n]==1) ? (sqrt(3./5.))*R : find_gyr_rad(Grad,n_agg,xlab,com_tf,nparcl[n]);
            printf("%d\t%d\t%g\t%g\t%d\n",n,nparcl[n],max_radius,gyr_radius,countsteps);
            char filename2[MAXFNLEN];
            sprintf(filename2,"res%d_rotation_N_%d_Dt_%g",s,TOTNPARCL,Dt); // one output per seed
            FILE *fp2 = fopen(filename2,"a");
            fprintf(fp2,"%d\t%d\t%g\t%g\t%d\n",n,nparcl[n],max_radius,gyr_radius,countsteps);
            fclose(fp2);
        }
    }

    /*
    // print the position of the aggregates left when simulation stops
    char filename_pos[MAXFNLEN];
    sprintf(filename_pos,"pos%d_rotation_N_%d_Dt_%g",s,TOTNPARCL,Dt);
    FILE *fp_pos = fopen(filename_pos,"w");
    print_aggs_data(fp_pos,xlab,nparcl,nagg);
    fclose(fp_pos);*/

    double time_taken = double(end - start);

    std::cout << "Execution time is: " << time_taken << " sec" << std::endl;
    std::cout << "Nsteps = " << countsteps << std::endl;
    std::cout << "L = " << BOXL << std::endl;
    delete[] xlab;
    delete[] xrel;

    return 0;

}

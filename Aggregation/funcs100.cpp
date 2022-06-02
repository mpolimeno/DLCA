#include "common100.h"

// initialize position of the cubes in the aggregates
// at first each aggregate has only one cube
void init_aggs(double* aggs,int nparcl[],int nagg)
{

    // initialize first aggregate
    for (int d=0;d<DIM;d++)
    {
        double unif = rand()/(RAND_MAX+1.);
        *(aggs+d) = BOXL*(unif-0.5);
    }
    nparcl[0] = 1;

    // initialize the other aggregates
    for (int i=1;i<nagg;i++)
    {
        double tmp[DIM];

        bool keep_init = true;
        while(keep_init)
        {
            for (int d=0;d<DIM;d++)
            {
                double unif = rand()/(RAND_MAX+1.);
                tmp[d] = BOXL*(unif-0.5);
            }
            double distsq;
            for (int j=0;j<i;j++)
            {
                double tmp2;
                double sum = 0.;
                for (int d=0;d<DIM;d++)
                {
                    tmp2 = tmp[d] - (*(aggs+j*nagg*DIM+d));
                    tmp2 -= BOXL*rint(tmp2/BOXL);

                    sum += tmp2*tmp2;
                }
                // check if two aggregates overlap
                distsq = sum;
                if (distsq<=SEPDISTSQ) // I prefer particles not being within attaching range upon generation
                {
                    keep_init = true;
                    break; //break j loop, b/c aggregates are too close upon generation
                }
            }

            if (distsq>SEPDISTSQ)
            {
                keep_init = false; // end while loop
            }
        }
        for (int d=0;d<DIM;d++) (*(aggs+i*nagg*DIM+d)) = tmp[d]; // assign position to i aggregate
        nparcl[i] = 1;
    }
    return;
}

void compute_com(double *agg,int idx, int nparcl,double com[DIM])
{

    // for periodic boundary conditions, in order to compute the center of Mass
    // we map each cartesian coordinate to an angle

    int nagg = TOTNPARCL;
    // each angle for each coordinate of each particle in each aggregate
    // first index is the aggregate "leader", then we have the number of particles in the aggregate and their x,y,z pos
    double theta[TOTNPARCL][DIM];
    for (int i=0;i<nparcl;i++)
    {
        for (int d=0;d<DIM;d++)
        {
            // theta for each particle in each aggregate, in x y and z
            theta[i][d] = *(agg+idx*DIM+i*DIM+d);
            theta[i][d] = 2.*Pi*(theta[i][d]/BOXL);
        }
    }


    // define two new points xi and zeta from this angle theta
    double xi[TOTNPARCL][DIM];
    double zeta[TOTNPARCL][DIM];

    for (int i=0;i<nparcl;i++)
    {
        for (int d=0;d<DIM;d++)
        {
            // xi and zeta for each particle in each aggregate
            xi[i][d] = cos(theta[i][d]);
            zeta[i][d] = sin(theta[i][d]);
        }
    }

    // compute mean values of xi and zeta for each aggregate
    double xibar[DIM];
    double zetabar[DIM];

    for (int d=0;d<DIM;d++)
    {
        xibar[d] = 0.;
        zetabar[d] = 0.;
    }

    for (int i=0;i<nparcl;i++)
    {
        xibar[0] += xi[i][0];
        zetabar[0]+= zeta[i][0];
        xibar[1] += xi[i][1];
        zetabar[1] += zeta[i][1];
        xibar[2] += xi[i][2];
        zetabar[2] +=zeta[i][2];
    }

    for (int d=0;d<DIM;d++)
    {
        xibar[d] /= nparcl;
        zetabar[d] /= nparcl;
    }


    // now we map back onto new thetabar and use it to compute the center of mass for each aggregate

    double thetabarx = atan2(-zetabar[0],-xibar[0])+Pi;
    com[0] = BOXL*(thetabarx/(2.*Pi));

    double thetabary = atan2(-zetabar[1],-xibar[1])+Pi;
    com[1] = BOXL*(thetabary/(2.*Pi));

    double thetabarz = atan2(-zetabar[2],-xibar[2])+Pi;
    com[2] = BOXL*(thetabarz/(2.*Pi));

    // ad hoc way to deal with the fact that our domain is [-BOXL/2,BOXL)
    for (int d=0;d<DIM;d++)
    {
        com[d] = com[d] - BOXL*rint(com[d]/BOXL);
    }

}

// compute position in the body frame of reference
void compute_xrel(double *agg,int idx,int nparcl,double com[DIM],double *xrel)
{

    for (int i=0;i<nparcl;i++)
    {
        double tmp[DIM];
        for (int d=0;d<DIM;d++)
        {
            tmp[d] = *(agg+idx*DIM+i*DIM+d);
            tmp[d]-= com[d];
            tmp[d] -= BOXL*rint(tmp[d]/BOXL);
            *(xrel+idx*DIM+i*DIM+d) = tmp[d];
        }
    }
    return;
}

// compute position in the lab frame of reference
void compute_xlab(double *xrel,int idx,double com[DIM],int nparcl,double Q[][DIM],double *xlab)
{

    double xrot_mat[nparcl][DIM];
    for (int n=0;n<nparcl;n++)
    {
        for (int d=0;d<DIM;d++) xrot_mat[n][d] = 0.;
    }

    for (int n=0;n<nparcl;n++)
    {
        double xtemp[DIM];
        for (int d=0;d<DIM;d++) xtemp[d] = *(xrel+idx*DIM+n*DIM+d);
        double xrot[DIM];
        for (int i=0;i<DIM;i++)
        {
            xrot[i] = 0.;
            for (int j=0;j<DIM;j++)
            {
                xrot[i] += Q[i][j]*xtemp[j];
            }
        }
        for (int d=0;d<DIM;d++) xrot_mat[n][d] = xrot[d];
    }

    for (int n=0;n<nparcl;n++)
    {
        for (int d=0;d<DIM;d++) *(xlab+idx*DIM+n*DIM+d) = 0.;
    }

    for (int n=0; n<nparcl; n++)
    {
        for (int d=0;d<DIM;d++)
        {
            *(xlab+idx*DIM+n*DIM+d) = com[d] + xrot_mat[n][d];
            *(xlab+idx*DIM+n*DIM+d) -= BOXL*rint(*(xlab+idx*DIM+n*DIM+d)/BOXL);
        }
    }

    return;
}

double find_max_rad(double maxrad[],int idx,double *agg,double com[DIM],int nparcl)
{
    double radmax = 0.;

    for (int i=0; i<nparcl; i++)
    {
        double tmpx = *(agg+idx*DIM+i*DIM+0);
        tmpx -= com[0];
        tmpx -= BOXL*rint(tmpx/BOXL);

        double tmpy = *(agg+idx*DIM+i*DIM+1);
        tmpy -= com[1];
        tmpy -= BOXL*rint(tmpy/BOXL);

        double tmpz = *(agg+idx*DIM+i*DIM+2);
        tmpz -= com[2];
        tmpz -= BOXL*rint(tmpz/BOXL);

        maxrad[i] = sqrt((tmpx*tmpx)+(tmpy*tmpy)+(tmpz*tmpz));

        if ( maxrad[i] > radmax)
        {
            radmax = maxrad[i];
        }

    }
    return radmax+R;
}

double find_gyr_rad(double Grad[],int idx,double *agg,double com[DIM],int nparcl)
{
    double sum = 0.;
    for (int i=0; i<nparcl; i++)
    {
        double tmpx = *(agg+idx*DIM+i*DIM+0);
        tmpx -= com[0];
        tmpx -= BOXL*rint(tmpx/BOXL);

        double tmpy = *(agg+idx*DIM+i*DIM+1);
        tmpy -= com[1];
        tmpy -= BOXL*rint(tmpy/BOXL);

        double tmpz = *(agg+idx*DIM+i*DIM+2);
        tmpz -= com[2];
        tmpz -= BOXL*rint(tmpz/BOXL);

        Grad[i] = ((tmpx*tmpx)+(tmpy*tmpy)+(tmpz*tmpz));

        sum = sum + Grad[i];

    }
    return sqrt((3./5.)*R*R + (sum/nparcl));
}

// compute distance between two aggregates
double compute_distsq(double *agg,int idx_1,int nparcl1,int idx_2,int nparcl2)
{
    double distmin = 0.5*BOXL*sqrt(DIM);
    double distminsq = distmin*distmin;

    for (int i=0;i<nparcl1;i++)
    {
        for (int j=0;j<nparcl2;j++)
        {
            double tmp;
            double sum = 0.;
            for (int d=0;d<DIM;d++)
            {
                tmp = *(agg+idx_1*DIM+i*DIM+d) - *(agg+idx_2*DIM+j*DIM+d);//agg1[i][d]-agg2[j][d];
                tmp -= BOXL*rint(tmp/BOXL);

                sum += tmp*tmp;
            }
            distminsq = (sum<distminsq) ? sum : distminsq;
        }
    }
    return distminsq;
}

// combine two aggregates: agg1 + agg2 -> agg1 (extended) + nothing
void combine_aggs(double *agg,int idx_1,int &nparcl1,int idx_2,int &nparcl2)
{
    if (nparcl1+nparcl2>MAXNPARCL)
    {
        printf("ERROR: maximum number of cubes reached\n");
        exit(0);
    }

    for (int i=0;i<nparcl2;i++)
        for (int d=0; d<DIM; d++)
            *(agg+idx_1*DIM+(i+nparcl1)*DIM+d) = *(agg+idx_2*DIM+i*DIM+d);

    nparcl1 = nparcl1 + nparcl2;
    nparcl2 = 0;

    return;
}

// print position of the aggregate
void print_aggs_data(FILE *fp,double *aggs,int nparcl[],int nagg)
{
    for (int n=0;n<nagg;n++)
    {
        if (nparcl[n]==0) continue;
        for(int i=0;i<nparcl[n];i++)
        {
            for(int j=0;j<DIM;j++)
            {
                fprintf(fp,"%g\t",*(aggs+n*nagg*DIM+i*DIM+j));
            }
            fprintf(fp,"\n");
        }
    }

    return;
}

#include "ising_model.h"

void initialize()
{
    allo_aux();
    // set the observables to 2: the energy and magnetization
    Obser = new Obs[2];
    def_latt();
    init_cnf();
    set_Obs();
    def_prob();
    return;
}

void allo_aux()
{
    if (D == 2)
    {
        Ly = Lx;
        Vol = Lx * Ly;
    }
    else
    {
        std::cerr << "D should be 2" << std::endl;
        exit(1);
    }

    // there are 2* (4D+1) possibilities (2 stands for 2 initial state -1 and +1)
    // (4D+1) stands for -2D to 2D overall spin states
    p_flip = new double[(4 * D + 1) * 2];
    dE_flip = new double[(4 * D + 1) * 2];

    if (D == 2)
    {
        nnsite = new int[4 * Vol]; // nnsite[4][Vol]
        backdir = new int[4];
    }
    else
    {
        std::cerr << "D should be 2" << std::endl;
        exit(1);
    }
    return;
}

////////////////////////////////////////////////////////////////////////////
// initialize configuration
////////////////////////////////////////////////////////////////////////////
void init_cnf()
{
    std::cout << "init_cnf" << std::endl;
    sp = new int[Vol];
    // fill sp with 1
    std::fill_n(sp, Vol, 1);
    // for (int i = 0; i< Vol; i++)
    // {
    //     std::cout << sp[i] << std::endl;
    // }
    return;
}

////////////////////////////////////////////////////////////////////////////
// Core simulation module
////////////////////////////////////////////////////////////////////////////

void markov(void)
{
    // std::cout << "markov" << std::endl;
    for (int istep_ = 0; istep_ < N_measure; istep_++)
    {
        worm();
    }
    return;
}

void worm(void)
{
    // I is a random integer from 0 to dir - 1
    int I = floor(Vol * ((double)rand() / (RAND_MAX))); // TODO: I = round(Vol * rand()); is incorrect
    int s = 0;
    for (int dir_ = 0; dir_ < nnb; dir_++) // nnb=4
    {
        // print*, sp(nnsite(I,dir))
        s += sp[nnsite[dir_ * Vol + I]]; // nnsite[4][Vol] //WRONG
    }

    // TODO: when converting from Fortran to cpp, the two dimential array should have it's x and y axis flipped
    int ind_ = spin2index(s, sp[I]);
    double p = p_flip[ind_];
    dE = dE_flip[ind_];
    if (p > 1 || ((double)rand() / (RAND_MAX)) < p)
    {
        sp[I] = -1 * sp[I];
        MP = MP + 2 * sp[I];
        E = E + dE;
    }
    return;
}

int spin2index(int d_spin, int s)
{
    int index = 0;
    index = (d_spin + 2 * D) + (s + 1) / 2 * (4 * D + 1);
    return index;
}
////////////////////////////////////////////////////////////////////////////
// define flipping probability
////////////////////////////////////////////////////////////////////////////
void def_prob()
{
    std::cout << "def_prob" << std::endl;
    for (int i_ = -1; i_ <= 1; i_ = i_ + 2) // for (int i_ = -1; i_ <= 1; i_++)
    {
        for (int j_ = -2 * D; j_ <= 2 * D; j_++)
        {
            /*define the energy change for flipping a spin here*/
            int ind_ = spin2index(j_, i_);
            dE = 2 * i_ * (t * j_ + h);
            dE_flip[ind_] = dE;
            p_flip[ind_] = exp(-1 * beta * dE);
        }
    }
    return;
}

////////////////////////////////////////////////////////////////////////////
// definite Lattice
////////////////////////////////////////////////////////////////////////////
void def_latt()
{
    /* Lx ganna be an even number */
    if (Lx % 2 != 0)
    {
        std::cerr << "L be even?" << std::endl;
        exit(1);
    }
    if (D == 2)
    {
        nnb = 4;
        for (int i_ = 0; i_ < nnb; i_++)
        {                                   //  13 14 15 16
            backdir[i_] = (i_ + 1) % 4; //  9 10 11 12
        }                                   // 5  6  7  8
        for (int j_ = 0; j_ < Vol; j_++)    //   1  2  3  4
        {
            nnsite[0 * Vol + j_] = j_ / Lx * Lx + (j_ + 1) % Lx;
            nnsite[1 * Vol + j_] = (j_ - Lx + Vol) % Vol;
            nnsite[2 * Vol + j_] = j_ / Lx * Lx + (j_ - 1 + Lx) % Lx;
            nnsite[3 * Vol + j_] = (j_ + Lx) % Vol;
        }
        // Lx = 4: Vol = 16
        // backdir = {2: down, 3: left, 4: up, 1: right}
        // 1 5 9  13
        // 2 6 10 14
        // 3 7 11 15
        // 4 8 12 16
        // nnsite: 
        // 2  3  4  1  6  7  8  5  10 11 12 9  14 15 16 13
        // 13 14 15 16 1  2  3  4  5  6  7  8  9  10 11 12
        // 4  1  2  3  8  5  6  7  12 9  10 11 16 13 14 15
        // 5  6  7  8  9  10 11 12 13 14 15 16 1  2  3  4


        // print out for checking
        // for (int i_ = 0; i_ < nnb; i_++) {
        //     std::cout<< backdir[i_] << " ";
        // }
        // std::cout<<std::endl;
        //  for (int j_ = 0; j_ < Vol; j_++)    //   transposed output
        // {
        //     std::cout<< nnsite[0 * Vol + j_]  << " ";
        //     std::cout<< nnsite[1 * Vol + j_]  << " ";
        //     std::cout<< nnsite[2 * Vol + j_]  << " ";
        //     std::cout<< nnsite[3 * Vol + j_]  << " ";
        //     std::cout<<std::endl;
        // }
        // exit(0);
    }
    else
    {
        std::cerr << "wrong! D must be 2" << std::endl;
        exit(1);
    }
    nnb1 = nnb + 1;
    nnb2 = nnb / 2; // for nnb = 4, nnb2 = 2
    return;
}

////////////////////////////////////////////////////////////////////////////
// module for getting the energy of the system (in a global way)
////////////////////////////////////////////////////////////////////////////

double getE()
{
    double getE_ = 0;

    for (int i = 0; i < Vol; i++)
    {
        getE_ -= h * sp[i];
        for (int dir_ = 0; dir_ < nnb2; dir_++)
        {
            getE_ -= t * sp[i] * sp[nnsite[dir_ * Vol + i]];
        }
    }
    return getE_;
}

double getMP()
{
    double getMP_ = 0.0;

    for (int Vc_ = 0; Vc_ < Vol; Vc_++)
    {
        getMP_ += sp[Vc_];
    }
    return getMP_;
}

/*initialize the spin configuration*/
void init_stat()
{
    std::cout << "init_stat" << std::endl;
    /*get the initial energy*/
    E = getE();
    std::cout << "init energy" << " "  << E  <<  std::endl;
    /*get the initial magnetization*/
    MP = getMP();
    std::cout << "init magnetization" << " " << MP <<  std::endl;
    wormstep = 0;
    return;
}

void measure()
{
    Obser[0].vnr = E;
    Obser[1].vnr = MP;

    if (prtflg)
    {
        coll_data();
    }

    wormstep++;

    return;
}

void set_Obs()
{
    std::string obs_name1 = "energy  ";
    std::string obs_name2 = "M_z     ";
    init_Obs(0, obs_name1, 1.0 / Vol, 0.0);
    init_Obs(1, obs_name2, 1.0 / Vol, 0.0);
    Npb = 1;
    Nbk = 0;
    Np = 0;
    return;
}

void init_Obs(int i_, std::string nam0_, double a0, double b0)
{
    Obser[i_].nam = nam0_;
    Obser[i_].a = a0;
    Obser[i_].b = b0;
    Obser[i_].vnr = 0.0;
    Obser[i_].val = 0.0;
    Obser[i_].cor = 1.00;
    Obser[i_].erb = 1.00;
    std::fill_n(Obser[i_].blk, Max_block, 0.0);

    return;
}

////////////////////////////////////////////////////////////////////////////
// Core module for calculating the error and correlations
////////////////////////////////////////////////////////////////////////////
void coll_data()
{
    // std::cout << "coll_data" << std::endl;
    if (Np == Npb)
    {
        if (flg_cor)
        {
            if (Nbk == Max_block)
            {
                merge_blk();
            }
        }
        else
        {
            if (Nbk >= N_block && Nbk % 2 * 2 == Nbk)
            {
                if (cal_cor_dev())
                {
                    flg_cor = true;
                }
                else
                {
                    merge_blk();
                }
            }
        }
        Nbk++;
        Np = 0;
    }
    Np++;
    for (int i = 0; i < NObs; i++)
    {
        Obser[i].blk[Nbk] = Obser[i].blk[Nbk] + Obser[i].vnr;
    }
    return;
}

void merge_blk()
{
    Nbk = Nbk / 2;
    Npb = 2 * Npb;
    for (int i = 0; i < NObs; i++)
    {
        for (int j = 0; j < Nbk; j++)
        {
            Obser[i].blk[j] = Obser[i].blk[2 * j - 1] + Obser[i].blk[2 * j];
        }
        std::fill(std::begin(Obser[i].blk)+Max_block/2, std::end(Obser[i].blk), 0.0);
    }
    return;
}

bool cal_cor_dev()
{
    int Nbk0;
    double cor, dev, devp, devn;
    bool cal_cor_dev_ = true;

    if (Np != Npb)
    {
        Nbk0 = Nbk - 1;
    }
    else
    {
        Nbk0 = Nbk;
    }
    for (int i = 0; i < NObs; i++)
    {
        std::accumulate(std::begin(Obser[i].blk), std::begin(Obser[i].blk) + Nbk0, Obser[i].val);

        Obser[i].val /= Nbk0;
        cor = 0.0;
        dev = 0.0;
        devp = 0.0;
        for (int iblck = 0; iblck < Nbk0; iblck++)
        {
            devn = Obser[i].blk[iblck] - Obser[i].val;
            dev += devn * devn;
            cor += devn * devp;
            devp = devn;
        }
        Obser[i].cor = cor / dev;
        Obser[i].val = Obser[i].val / Npb;
        Obser[i].erb = sqrt(dev / Nbk0 / (Nbk0 - 1.0)) / Npb;
        if (abs(Obser[i].cor) > Max_cor)
        {
            cal_cor_dev_ = false;
        }
    }
    return cal_cor_dev_;
}

void printing()
{

    std::cout << "========================================" << std::endl;
    std::cout << "Now	:" << std::endl;
    std::cout << "Z=" << Npb * (Nbk - 1) + Np << std::endl;
    for (int i = 0; i < NObs; i++)
    {
        std::cout << i << " " << Obser[i].nam << std::endl;
    }
    if (prtflg)
    {
        std::cout << "Average:" << std::endl;
        // Total step: total number of flips
        std::cout << "Total steps= " << Npb * (Nbk - 1) << std::endl;
        // <energy>     average      error      correlation
        // <M_z>        average      error      correlation
        for (int i = 0; i < NObs; i++)
        {
            std::cout << i + 1 << " " << '<' << Obser[i].nam << '>' << " " << Obser[i].val * Obser[i].a + Obser[i].b << " " << Obser[i].erb * Obser[i].a << " " << Obser[i].cor << std::endl;
        }
        std::cout << "flg_cor= " << flg_cor << " "
                  << "Npb= " << Npb << " "
                  << "Nbk= " << Nbk << " "
                  << "Np=  " << Np << std::endl;
        std::cout << "wormstep= " << wormstep << std::endl;
        // std::cout << "simulation time :" << " " << t_simu <<std::endl
    }
    else
    {
        std::cout << "therm step: " << itoss << std::endl;
        // std::cout << "therm time:"<< " " << t_toss << std::endl;
    }
    std::cout << "========================================" << std::endl;
    return;
}

bool check(void)
{

    int MP0;
    double E0;

    std::cout << "checking ..." << std::endl;
    E0 = getE();
    MP0 = getMP();
    if (std::abs(E - E0) > 1e-7)
    {
        std::cout << "wrong!"
                  << " \n"
                  << "Current:" << E << " "
                  << "Correct:"
                  << " " << E0 << std::endl;
        return false;
    }
    if (MP0 != MP)
    {
        std::cout << "wrong!"
                  << " \n"
                  << "Current:" << MP << " "
                  << "Correct:"
                  << " " << MP0 << std::endl;
        return false;
    }

    std::cout << "You are good!"
              << " \n"
              << "Current:" << E << " "
              << "Correct:"
              << " " << E0 << std::endl;
    std::cout << "You are good!"
              << " \n"
              << "Current:" << MP << " "
              << "Correct:"
              << " " << MP0 << std::endl;
    std::cout << "check done" << std::endl;
    return true;
}

int main(void)
{
    int i2_ = 0;
    int itoss = 0;

    /* input the system size */
    std::cout << "please input the system size" << std::endl;
    std::cin >> Lx;
    std::cout << "The simulation scale: " << Lx << std::endl;

    /* input the system parameter */
    double t, h, beta;
    std::cout << "please input the system t, h and beta" << std::endl;
    std::cin >> t >> h >> beta;
    std::cout << "The simulation parameter: " << t << " " << h << " " << beta << std::endl;

    /* initialize the system */
    initialize();

    std::cout << "create a new simulation:" << std::endl;
    init_stat();

    ////////////////////////////////////////////////////////////////////////////
    // Thermialization
    ////////////////////////////////////////////////////////////////////////////
    prtflg = false;
    for (int isamp = 1; isamp <= Ntoss; isamp++)
    {
        markov();
        measure();
        i2_++;
        if (i2_ >= N_print)
        {
            itoss = itoss + N_print;
            check();
            printing();
            i2_ = 0;
        }
    }
    wormstep = 0;

    ////////////////////////////////////////////////////////////////////////////
    // Simulation
    ////////////////////////////////////////////////////////////////////////////
    prtflg = true;
    for (int isamp = 1; isamp <= Nsamp; isamp++)
    {
        markov();
        measure();

        i2_++;
        if (i2_ >= N_print)
        {
            printing();
            check();
            flg_cor = cal_cor_dev();
            i2_ = 0;
        }
    }

    exit(0);
}

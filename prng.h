/*------------------------Random number generator-------------------------------*/


/*
void PRNG_random_seed()
{
    // Seed with external entropy
    
    pcg128_t seeds[2];
    entropy_getbytes((void*)seeds, sizeof(seeds));
    pcg64_srandom_r(&rng, seeds[0], seeds[1]);
}

static inline double PRNG_double01()
{
    return (double)pcg64_random_r(&rng)/(double)UINT64_MAX;
}

static inline void PRNG_gaussian01(double *random)
{
    double u1, u2, r = 0.;
    
    while(r==0 || r >= 1)
    {
        u1 = 2.*PRNG_double01() - 1.;
        u2 = 2.*PRNG_double01() - 1.;
        r = u1*u1 + u2*u2;
    }
    
    r = sqrt(-2.*log(r)/r);
    random[0] = u1*r;
    random[1] = u2*r;
}
*/

static inline double randNum(){
/*
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.,1.);
    return dis(gen);*/

    pcg_extras::seed_seq_from<std::random_device> seed_source; // seeding
    pcg32 rng(seed_source); // rng

    std::uniform_real_distribution<float> uni_dist(0.,1.);
    return uni_dist(rng);
}

static inline int randSite(const int N){
/*
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0,N-1);
    return dist(gen);*/

    pcg_extras::seed_seq_from<std::random_device> seed_source;
    pcg32 rng(seed_source);

    std::uniform_int_distribution<int> uniform_dist(0,N-1);
    return uniform_dist(rng);
}

static inline int randDir(){
    
    double coinToss = randNum();
    
    if (coinToss < 0.5)
        return 1;
    else
        return -1;
}

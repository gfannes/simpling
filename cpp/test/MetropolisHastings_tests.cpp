#include "simpling/MetropolisHastings.hpp"
#include <iostream>
#include <algorithm>
#include <array>

namespace  { 

    class Target
    {
        public:
            typedef double result_type;
            double radius = 2.0;

            //A ^-shaped distribution centered around 0 with domain [-radius, radius]
            bool log_prob(result_type &lp, result_type x) const
            {
                const auto d = std::abs(x);
                if (d >= radius)
                    return false;
                lp = std::log(radius - d);
                return true;
            }
    };

    //Initialization and generation of potential new states for the markov chain
    class Proposal
    {
        public:
            typedef double result_type;
            typedef std::normal_distribution<> Noise;

            //This might generate a few initial states with zero probability
            Noise initialization_noise = Noise(0.0, 100.0);
            template <typename G>
                bool initialize(result_type &x, G &g)
                {
                    x = initialization_noise(g);
                    std::cout << "Initialization: " << x << std::endl;
                    return true;
                }

            //Rule of thumb is to aim for acceptance rate in [0.5, 0.85]
            Noise generation_noise = Noise(0.1, 0.5);
            template <typename G>
                bool generate(result_type &x, double &ratio, const result_type &prev, G &g)
                {
                    //Generate a new proposal state starting from the previous
                    x = prev + generation_noise(g);
                    //Normally, the proposal distribution is chosen to be symmetric (generation_noise.mean() == 0.0),
                    //in which case ratio should be set to 1.0 always.
                    //If this is _not_ the case, make sure you adjust for it correctly.
                    ratio = generation_density_unnorm(prev, x)/generation_density_unnorm(x, prev);
                    return true;
                }
            //Density of the new proposal state x being generated starting from given
            double generation_density_unnorm(double x, double given) const
            {
                const double d = (x - (given+generation_noise.mean()))/generation_noise.stddev();
                return std::exp(-0.5*d*d);
            }
    };

    //We combine the wanted Target distribution together with the Proposal distribution into our Metropolis-Hastings markov chain
    typedef simpling::MetropolisHastings<Target, Proposal> MyMHMC;

    template <size_t NrBins>
        class Histogram
        {
            public:
                Histogram(double min, double max):
					min_(min), max_(max)
			{
				for (auto &cnt: bins_)
					if (cnt != 0)
						//This is a compiler error, std::array should be default initialized, which is 0 for unsigned long
						cnt = 0;
			}
                bool add(double v)
                {
                    if (v < min_)
                        return false;
                    if (v > max_)
                        return false;

                    size_t ix = (v-min_)/width_;
                    if (ix >= bins_.size())
                        ix = bins_.size()-1;
                    ++bins_[ix];

                    return true;
                }
                void stream(std::ostream &os, unsigned long maxNrCols) const
                {
                    unsigned long maxBinCnt = *std::max_element(bins_.begin(), bins_.end());
                    const double scale = double(maxNrCols)/maxBinCnt;
                    unsigned long total = 0;
                    for (auto cnt: bins_)
                    {
                        total += cnt;
                        os << std::string((size_t)(scale*cnt), '*') << std::endl;
                    }
                    os << "Domain: [" << min_ << ", " << max_ << "], " << total << " observations, " << NrBins << " bins" << std::endl;
                }
            private:
                const double min_;
                const double max_;
                const double width_ = (max_-min_)/NrBins;
                std::array<unsigned long, NrBins> bins_;
        };

    class Acceptance
    {
        public:
            void add(bool accepted)
            {
                if (accepted)
                    ++nrAccepted_;
                ++nrTotal_;
            }
            double rate() const {return double(nrAccepted_)/nrTotal_;}
            void stream(std::ostream &os) const
            {
                os << "Acceptance rate: " << rate() << std::endl;
            }
        private:
            unsigned long nrAccepted_ = 0;
            unsigned long nrTotal_ = 0;
    };
} 

int main()
{
    //Some well-known RNG
    std::mt19937 engine;

    //Our markov chain
    MyMHMC mhmc;

    //Burn-in
    for (int i = 0; i < 1000000; ++i)
    {
        if (!mhmc(engine))
            std::cout << "Problem in iteration " << i << ", probably with the initialization" << std::endl;
    }
    assert(mhmc.isInitialized());

    //Generation, collection into histogram and computation of acceptance rate
    Histogram<40> histogram(-mhmc.target().radius, mhmc.target().radius);
    Acceptance acceptance;
    for (int i = 0; i < 1000000; ++i)
    {
        if (!mhmc(engine))
            //Something went wrong
            continue;
        //This is a good value, record it into the histogram and the acceptance rate
        histogram.add(mhmc.value());
        acceptance.add(mhmc.isNew());
    }

    //Print the histogram and acceptance rate
    histogram.stream(std::cout, 80);
    acceptance.stream(std::cout);

    return 0;
}

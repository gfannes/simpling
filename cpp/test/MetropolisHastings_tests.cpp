#include "simpling/MetropolisHastings.hpp"
#include <iostream>
#include <map>
#include <algorithm>

namespace  { 

    class Target
    {
        public:
            typedef double result_type;
            double radius = 2.0;

            bool log_prob(result_type &lp, result_type x) const
            {
                const auto d = std::abs(x);
                if (d >= radius)
                    return false;
                lp = std::log(radius - d);
                return true;
            }
        private:
    };
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
                    std::cout << x << std::endl;
                    return true;
                }

            //Rule of thumb is to aim for acceptance rate in [0.5, 0.85]
            Noise generation_noise = Noise(0.0, 0.5);
            template <typename G>
                bool generate(result_type &x, double &ratio, const result_type &prev, G &g)
                {
                    //Our simple proposal distribution is symmetric, so ratio == 1.0
                    ratio = 1.0;
                    x = prev + generation_noise(g);
                    return true;
                }
        private:
    };

    typedef simpling::MetropolisHastings<Target, Proposal> MyMHMC;

    template <size_t NrBins>
        class Histogram
        {
            public:
                Histogram(double min, double max): min_(min), max_(max) {}
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
                void stream(std::ostream &os, unsigned long max) const
                {
                    unsigned long m = *std::max_element(bins_.begin(), bins_.end());
                    const double scale = double(max)/m;
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
            unsigned long nrAccepted_;
            unsigned long nrTotal_;
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
            std::cout << "Problem in iteration " << i << std::endl;
    }
    assert(mhmc.isInitialized());

    //Generation, collection into histogram and computation of acceptance rate
    Acceptance acceptance;
    Histogram<40> histogram(-mhmc.target().radius, mhmc.target().radius);
    for (int i = 0; i < 1000000; ++i)
    {
        if (!mhmc(engine))
            //Something went wrong
            continue;
        histogram.add(mhmc.value());
        acceptance.add(mhmc.isNew());
    }

    histogram.stream(std::cout, 80);
    acceptance.stream(std::cout);

    return 0;
}

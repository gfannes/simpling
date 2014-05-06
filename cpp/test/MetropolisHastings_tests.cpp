#include "simpling/MetropolisHastings.hpp"
#include <iostream>

namespace  { 

    class Target
    {
        public:
            typedef double result_type;

            double log_prob(result_type x) const
            {
                return 0;
            }
        private:
    };
    class Proposal
    {
        public:
            typedef double result_type;

            template <typename G>
                bool initialize(result_type &x, G &g)
                {
                    x = 0;
                    return true;
                }
            template <typename G>
                bool generate(result_type &x, double &ratio, const result_type &prev, G &g)
                {
                    ratio = 1.0;
                    x = prev;
                    return true;
                }
        private:
    };

    typedef simpling::MetropolisHastings<Target, Proposal> MyMHMC;
} 

int main()
{
    std::mt19937 engine;
    MyMHMC mhmc;
    mhmc(engine);
    return 0;
}

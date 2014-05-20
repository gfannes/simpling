#ifndef HEADER_simpling_MetropolisHastings_hpp_ALREADY_INCLUDED
#define HEADER_simpling_MetropolisHastings_hpp_ALREADY_INCLUDED

#include <random>
#include <cassert>

namespace simpling { 

    //User should ensure that Proposal::initialize() and Proposal::generate()
    //generate a value compatible with the layout of Target
    //Target: The target distribution you want to sample from
    //        Provide log_prob(logprob, location), the logarithm of the unnormalized target density
    //Proposal: The distribution that generates the initial state, and candidate states based on the
    //previous state
    template <typename Target, typename Proposal>
        class MetropolisHastings: public Target, public Proposal
        {
            public:
                typedef typename Target::result_type result_type;

                //Access to the target distribution
                Target &target() {return static_cast<Target&>(*this);}
                const Target &target() const {return static_cast<const Target&>(*this);}

                //Access to the proposal distribution
                Proposal &proposal() {return static_cast<Proposal&>(*this);}
                const Proposal &proposal() const {return static_cast<const Proposal&>(*this);}

                bool isInitialized() const {return isInitialized_;}
                bool isNewState() const {return isNewState_;}

                //The current state of the chain
                const result_type &value() const
                {
                    assert(isInitialized_);
                    return x_;
                }

                //Performs one iteration of the Metropolis-Hastings algorithm
                //This will generate a new state, which might be the same as the previous one
                //Returns false if something went wrong (e.g., failure to generate the initial state, or a new candidate state)
                //The return value does not reflect whether the generated candidate state was accepted as the new state
                template <typename G>
                    bool operator()(G &g)
                    {
                        //Initialize the chain, if not already done so
                        if (!initialize_(g))
                            //Could not initialize the chain
                            return false;

                        //Generate a new candidate, together with the ratio of the proposal probs (p(x_|cand)/p(cand|x_))
                        result_type candidate;
                        double ratio = -1;
                        if (!Proposal::generate(candidate, ratio, const_cast<const result_type &>(x_), g))
                            //Could not generate a candidate
                            return false;
                        assert(ratio >= 0);

                        //We have a candidate, so nothing can fail anymore. By default, we indicate that
                        //this new state is rejected
                        isNewState_ = false;

                        double candidate_lp;
                        if (!Target::log_prob(candidate_lp, candidate))
                        {
                            //This candidate has zero probability, we always reject it
                            //We still return true, this is not an error.
                            return true;
                        }

                        const double accept_prob = ratio*std::exp(candidate_lp - log_prob_);
                        if (accept_prob >= 1.0 or with_prob_(accept_prob, g))
                        {
                            //We accept the newly generated candidate: move it, together with its prob
                            x_ = std::move(candidate);
                            log_prob_ = candidate_lp;
                            isNewState_ = true;
                        }

                        return true;
                    }

            private:
                //Returns true with prob
                template <typename G>
                    bool with_prob_(const double prob, G &g)
                    {
                        return urd(g) < prob;
                    }
                //Initializes the internal state, if not already done so
                template <typename G>
                    bool initialize_(G &g)
                    {
                        if (isInitialized_)
                            return true;
                        if (!Proposal::initialize(x_, g))
                            return false;
                        isNewState_ = isInitialized_ = Target::log_prob(log_prob_, x_);
                        return isInitialized_;
                    }

                //Tracks if the chain is already initialized
                bool isInitialized_ = false;
                bool isNewState_ = false;
                //The state itself
                result_type x_;
                //Cache of Target::log_prob(lp, x_)
                double log_prob_;
                //A uniform distribution used to decide if the candidate state should be accepted or not
                std::uniform_real_distribution<double> urd;
        };

} 

#endif

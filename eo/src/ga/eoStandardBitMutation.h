
#ifndef _eoStandardBitMutation_h_
#define _eoStandardBitMutation_h_

#include "../utils/eoRNG.h"

/** Standard bit mutation with mutation rate p:
 * choose k from the binomial distribution Bin(n,p) and apply flip_k(x).
 *
 * @ingroup Bitstrings
 * @ingroup Variators
 */
template<class EOT>
class eoStandardBitMutation : public eoMonOp<EOT>
{
    public:
        eoStandardBitMutation(double rate = 0) :
            _rate(rate),
            _nb(1),
            _bitflip(_nb)
        {}

        virtual bool operator()(EOT& chrom)
        {
            if (this->_rate == 0) {this->_rate = (double) 1/chrom.size();}
            this->_nb = eo::rng.binomial(chrom.size(),_rate);
            this->_bitflip.number_bits(this->_nb);
            return this->_bitflip(chrom);
        }

        virtual std::string className() const {return "eoStandardBitMutation";}

    protected:
        double _rate;
        unsigned _nb;
        eoDetSingleBitFlip<EOT> _bitflip;
};

/** Uniform bit mutation with mutation rate p:
 * choose k from the uniform distribution U(0,n) and apply flip_k(x).
 *
 * @ingroup Bitstrings
 * @ingroup Variators
 */
template<class EOT>
class eoUniformBitMutation : public eoMonOp<EOT>
{
    public:
        eoUniformBitMutation() :
            _nb(1),
            _bitflip(_nb)
        {}

        virtual bool operator()(EOT& chrom)
        {
            this->_nb = eo::rng.random(chrom.size());
            this->_bitflip.number_bits(this->_nb);

            return this->_bitflip(chrom);
        }

        virtual std::string className() const {return "eoUniformBitMutation";}

    protected:
        double _rate;
        unsigned _nb;
        eoDetSingleBitFlip<EOT> _bitflip;
};


/** Conditional standard bit mutation with mutation rate p:
 * choose k from the binomial distribution Bin(n,p) until k >0
 * and apply flip_k(x).
 *
 * This is identical to sampling k from the conditional binomial
 * distribution Bin>0(n,p) which re-assigns the probability to sample
 * a 0 proportionally to all values i ∈ [1..n]. 
 *
 * @ingroup Bitstrings
 * @ingroup Variators
 */
template<class EOT>
class eoConditionalBitMutation : public eoStandardBitMutation<EOT>
{
    public:
        eoConditionalBitMutation(double rate = 0) :
            eoStandardBitMutation<EOT>(rate)
        {}

        virtual bool operator()(EOT& chrom)
        {
            if (this->_rate == 0) {this->_rate = (double) 1/chrom.size();}
            this->_nb = 0;
            while (this->_nb < 1)
            {
                this->_nb = eo::rng.binomial(chrom.size(),this->_rate);
            }
            this->_bitflip.number_bits(this->_nb);

            return this->_bitflip(chrom);
        }

        virtual std::string className() const {return "eoConditionalBitMutation";}
};

/** Shifted standard bit mutation with mutation rate p:
 * choose k from the binomial distribution Bin(n,p).
 * When k= 0, set k= 1. Apply flip_k(x).
 *
 * This is identical to sampling k from the conditional binomial
 * distribution Bin0→1(n,p) which re-assigns the probability to
 * sample a 0 to sampling k= 1.
 *
 * @ingroup Bitstrings
 * @ingroup Variators
 */
template<class EOT>
class eoShiftedBitMutation : public eoStandardBitMutation<EOT>
{
    public:
        eoShiftedBitMutation(double rate = 0.5) :
            eoStandardBitMutation<EOT>(rate)
        {}

        virtual bool operator()(EOT& chrom)
        {
            if (this->_rate == 0) {this->_rate = (double) 1/chrom.size();}
            this->_nb = eo::rng.binomial(chrom.size(),this->_rate);
            if(this->_nb == 0) {
                this->_nb = 1;
            }

            this->_bitflip.number_bits(this->_nb);

            return this->_bitflip(chrom);
        }

        virtual std::string className() const {return "eoShiftedBitMutation";}
};

/** Mutation which size is sample in a gaussian.
 *
 * sample k from the normal distribution N(pn,σ^2)
 * and apply flip_k(x).
 *
 * From:
 * Furong Ye, Carola Doerr, and Thomas Back.
 * Interpolating local and global search by controllingthe variance of standard bit mutation.
 * In 2019 IEEE Congress on Evolutionary Computation(CEC), pages 2292–2299.
 *
 * In contrast to standard bit mutation, this operators allows to scale
 * the variance of the mutation strength independently of the mean.
 *
 * @ingroup Bitstrings
 * @ingroup Variators
 */
template<class EOT>
class eoNormalBitMutation : public eoMonOp<EOT>
{
    public:
        eoNormalBitMutation(double expected = 0, double variance = 0) :
            _expected(expected),
            _variance(variance)
        {
            assert(rate <= 1);
        }

        virtual bool operator()(EOT& chrom)
        {
            if (this->_expected == 0) {this->_expected = (double) 1/chrom.size();}
            if (this->_variance == 0) {this->_variance = log(chrom.size());}

            this->_nb = eo::rng.normal(this->_expected, this->_variance);
            if(this->_nb >= chrom.size()) {
                this->_nb = eo::rng.random(chrom.size());
            }
            this->_bitflip.number_bits(this->_nb);

            return this->_bitflip(chrom);
        }

        virtual std::string className() const {return "eoNormalBitMutation";}

    protected:
        double _expected;
        unsigned _nb;
        eoDetSingleBitFlip<EOT> _bitflip;
        double _variance;
};

/** Fast mutation which size is sampled from an adaptive power law.
 *
 * From:
 * Benjamin Doerr, Huu Phuoc Le, Régis Makhmara, and Ta Duy Nguyen.
 * Fast genetic algorithms.
 * In Proc. of Genetic and Evolutionary Computation Conference (GECCO’17), pages 777–784.ACM, 2017.
 *
 * @ingroup Bitstrings
 * @ingroup Variators
 */
template<class EOT>
class eoFastBitMutation : public eoMonOp<EOT>
{
    public:
        eoFastBitMutation(double beta = 1.5) :
            _beta(beta)
        {
            assert(beta > 1);
        }

        virtual bool operator()(EOT& chrom)
        {
            this->_nb = powerlaw(chrom.size(),_beta);
            // BitFlip operator is bound to the _nb reference,
            // thus one don't need to re-instantiate.
            this->_bitflip.number_bits(this->_nb);
            return this->_bitflip(chrom);
        }

        virtual std::string className() const {return "eoFastBitMutation";}

    protected:

        double powerlaw(unsigned n, double beta)
        {
            double cnb = 0;
            for(unsigned i=1; i<=n/2; ++i) {
                cnb += std::pow(i,-beta);
            }
            double trigger = ((double) rand() / (RAND_MAX));
            double cursor = 0;
            double rate = 1;
            for(unsigned i=1; i<=n/2; ++i) {
                cursor += std::pow(i,-beta)/cnb;
                if (cursor >= trigger)
                {
                    rate = (double) i/n;
                    break;
                }
            }
            return eo::rng.binomial(n,rate);
        }

        double _beta;
        unsigned _nb;
        eoDetSingleBitFlip<EOT> _bitflip;
};

/** Bucket mutation which assign probability for each bucket
 *
 * From:
 * Carola Doerr, Johann Dréo, Alexis Robbes
 */
template<class EOT>
class eoBucketBitMutation : public eoMonOp<EOT>
{
    public:
        eoBucketBitMutation(std::vector<std::vector<int>> buckets, std::vector<double> bucketsValues) :
            _buckets(buckets),
            _bucketsValues(bucketsValues)
    
        {
            assert(map_bucket_bit.size() == bucketValues.size());
        }

        virtual bool operator()(EOT& chrom)
        {

            this->_nb = customlaw(chrom.size(), _buckets, _bucketsValues);
            // BitFlip operator is bound to the _nb reference,
            // thus one don't need to re-instantiate.
            //std::cout << "nb flip: " << this->_nb << std::endl;
            this->_bitflip.number_bits(this->_nb);
            return this->_bitflip(chrom);
        }

        virtual std::string className() const {return "eoBucketBitMutation";}

    protected:

        double customlaw(unsigned n, std::vector<std::vector<int>> buckets, std::vector<double> bucketsValues)
        {
            int bucketIndex = eo::rng.roulette_wheel(bucketsValues);
            int startBit = buckets[bucketIndex][0];
            int endBit = buckets[bucketIndex][1];
            int gapBit = endBit - startBit;
            int nbBits;
            if (gapBit > 0) { nbBits = rand() % gapBit + startBit;}
            else { nbBits = endBit;}
            //std::cout << nbBits << std::endl;
            return nbBits;
        }

        std::vector<double> _bucketsValues;
        std::vector<std::vector<int>> _buckets;

        unsigned _nb;
        eoDetSingleBitFlip<EOT> _bitflip;
};
#endif // _eoStandardBitMutation_h_

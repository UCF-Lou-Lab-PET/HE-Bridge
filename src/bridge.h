#ifndef BRIDGE_H
#define BRIDGE_H

#include <helib/helib.h>
#include <helib/Ptxt.h>
#include <helib/norms.h>
#include <NTL/mat_ZZ.h>

using namespace std;
using namespace NTL;
using namespace helib;

namespace he_bridge{

// the type of interpolation polynomial in beFV
enum CircuitType{UNI, BI, TAN};

class Bridge{
    const Context& m_context;
    unsigned long m_slotDeg;
    unsigned long m_expansionLen;
    vector<DoubleCRT> m_mulMasks;
    vector<double> m_mulMasksSize;
    // slot generator
    ZZX m_slot_gen;
    // secret key
    SecKey m_sk;
    // public key
    PubKey m_pk;
    // Polynomial evaluation type
    CircuitType m_type;
    
    // univariate comparison polynomial of the less-than function
    ZZX m_univar_less_poly;
    // univariate comparison polynomial of the less-than function
    ZZX m_univar_min_max_poly;
    // bivariate comparison polynomial coefficients of the less-than function
    mat_ZZ m_bivar_less_coefs;

    // polynomial evaluation parameters of the Patterson-Stockmeyer algorithm
    // number of baby steps
    long m_bs_num_comp;
    long m_bs_num_min;
    // number of giant steps
    long m_gs_num_comp;
    long m_gs_num_min;
    // leading coefficient
    ZZ m_top_coef_comp;
    ZZ m_top_coef_min;
    // extra coefficient
    ZZ m_extra_coef_comp;
    ZZ m_extra_coef_min;

    // indexes to compute x^{p-1}
    long m_baby_index;
    long m_giant_index;

    // Using finite extension field
    // this is basically the 2D matrix of the \rho in Tan's paper and XA=I, the A matrix
    vector<vector<DoubleCRT>> m_extraction_const;
    vector<vector<double>> m_extraction_const_size;

    // print/hide flag for debugging
  	bool m_verbose;

    // Define functions for aggregation
    // create multiplicative masks for shifts
  	DoubleCRT create_shift_mask(double& size, long shift);
  	void create_all_shift_masks();
    // shifts ciphertext slots to the left by shift within batches of size m_expansionLen starting at start. Slots shifted outside their respective batches are zeroized.
    void batch_shift(Ctxt& ctxt, long start, long shift) const;
    // shifts ciphertext slots to the left by shift within batches of size m_expansionLen starting at start. Slots shifted outside their respective batches filled with 1.
    void batch_shift_for_mul(Ctxt& ctxt, long start, long shift) const;
    // running sums of slot batches
    void shift_and_add(Ctxt& x, long start, long shift_direction = false) const;
    // running products of slot batches
    void shift_and_mul(Ctxt& x, long start, long shift_direction = false) const;

    // compute Patterson-Stockmeyer parameters to evaluate the comparison polynomial
    void compute_poly_params();
    // create the comparison polynomial
    void create_poly();
    // univariate comparison polynomial evaluation
    void evaluate_univar_less_poly(Ctxt& ret, Ctxt& ctxt_p_1, const Ctxt& x) const;

    // send non-zero elements of a field F_{p^d} to 1 and zero to 0
    // if pow = 1, this map operates on elements of the prime field F_p
    void mapTo01_subfield(Ctxt& ctxt, long pow) const;
    // exact equality 
    void is_zero(Ctxt& ctxt_res, const Ctxt& ctxt_z, long pow = 1) const;

    public:
    // constructor
	Bridge(const Context& context, CircuitType type, unsigned long d, unsigned long expansion_len, const SecKey& sk, bool verbose);

    const DoubleCRT& get_mask(double& size, long index) const;
    const ZZX& get_less_than_poly() const;

    // comparison x>0?
    void compare(Ctxt& ctxt_res, const Ctxt& ctxt_x) const;
    void reduce(std::vector<Ctxt>& digits, const Ctxt& c, long r) const;
    void lift(Ctxt& res, const Ctxt& c, long r) const;
    void print_decrypted(const Ctxt& ctxt) const;

    // input is from FV, evalute a relu in beFV and switch the result back to FV
    // the relu is based on the comparison function
    // test the arith-relu function with HEBridge
    void test_bridge(long runs) const;

    void m_test_bridge(long runs) const;
};
}

#endif // #ifndef COMPARATOR_H
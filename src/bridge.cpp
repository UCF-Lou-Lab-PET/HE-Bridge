#include "bridge.h"
#include "tools.h"
#include <helib/debugging.h>
#include <helib/polyEval.h>
#include <random>
#include <map> 
#include <NTL/ZZ_pE.h>
#include <NTL/mat_ZZ_pE.h>
#include <helib/Ptxt.h>

using namespace he_bridge;

DoubleCRT Bridge::create_shift_mask(double& size, long shift)
{
	cout << "Mask for shift " << shift << " is being created" << endl;
	// get EncryptedArray
  	const EncryptedArray& ea = m_context.getEA();
  	//extract slots
	long nSlots = ea.size();
	//number of batches in one slot
	// number of integers in one ciphertext
	long batch_size = nSlots / m_expansionLen;
	// create a mask vector
	vector<long> mask_vec(nSlots,1);

	//starting position of all batches
	long start = 0;
	// set zeros in the unused slots
	long nEndZeros = nSlots - batch_size * m_expansionLen;
	for (int i = 1; i <= nEndZeros; i++)
	{
	long indx = (start + nSlots - i) % nSlots;
	mask_vec[indx] = 0;
	}

	// masking values rotated outside their batches
	for (long i = 0; i < batch_size; i++)
	{
	if (shift < 0)
	{
	  for (long j = 0;  j < -shift; j++)
	  {
	    long indx = (start + (i + 1) * m_expansionLen - j - 1) % nSlots;
	    mask_vec[indx] = 0;
	  }
	}
	else if (shift > 0)
	{
	  for (long j = 0;  j < shift; j++)
	  {
	    long indx = (start + i * m_expansionLen + j) % nSlots;
	    mask_vec[indx] = 0;
	  }
	}
	}
	ZZX mask_zzx;
	ea.encode(mask_zzx, mask_vec);

	size = conv<double>(embeddingLargestCoeff(mask_zzx, m_context.getZMStar()));

	DoubleCRT mask_crt = DoubleCRT(mask_zzx, m_context, m_context.allPrimes());
	return mask_crt;
}

void Bridge::create_all_shift_masks()
{
	long shift = 1;
	while (shift < m_expansionLen)
	{
		double size;
	    DoubleCRT mask_ptxt = create_shift_mask(size, -shift);
	    m_mulMasks.push_back(mask_ptxt);
	    m_mulMasksSize.push_back(size);

	    shift <<=1;
	}
	cout << "All masks are created" << endl;
}

void Bridge::batch_shift(Ctxt& ctxt, long start, long shift) const
{
	HELIB_NTIMER_START(BatchShift);
	// get EncryptedArray
	const EncryptedArray& ea = m_context.getEA();
	
	// if shift is zero, do nothing
	if(shift == 0)
		return;

	// left cyclic rotation
	ea.rotate(ctxt, shift);

	// masking elements shifted out of batch
	long index = static_cast<long>(intlog(2, -shift));
	//cout << "Mask index: " << index << endl;
	double size;
	DoubleCRT mask = get_mask(size, index);
	ctxt.multByConstant(mask, size);
	HELIB_NTIMER_STOP(BatchShift);
}

void Bridge::batch_shift_for_mul(Ctxt& ctxt, long start, long shift) const
{
	HELIB_NTIMER_START(BatchShiftForMul);
	// get EncryptedArray
	const EncryptedArray& ea = m_context.getEA();
	
	// if shift is zero, do nothing
	if(shift == 0)
		return;
	// left cyclic rotation
	ea.rotate(ctxt, shift);
	
	long index = static_cast<long>(intlog(2, -shift));
	//cout << "Mask index: " << index << endl;
	double mask_size;
	DoubleCRT mask = get_mask(mask_size, index);
	ctxt.multByConstant(mask, mask_size);

	// add 1 to masked slots
	ctxt.addConstant(ZZ(1));
	mask.Negate();
	ctxt.addConstant(mask, mask_size);

	HELIB_NTIMER_STOP(BatchShiftForMul);
}

void Bridge::shift_and_add(Ctxt& x, long start, long shift_direction) const
{
  HELIB_NTIMER_START(ShiftAdd);
  long shift_sign = -1;
  if(shift_direction)
    shift_sign = 1;

  long e = 1;

  // shift and add
  while (e < m_expansionLen){
    Ctxt tmp = x;
    batch_shift(tmp, start, e * shift_sign);
    x += tmp;
    e <<=1;
  }
  HELIB_NTIMER_STOP(ShiftAdd);
}

void Bridge::shift_and_mul(Ctxt& x, long start, long shift_direction) const
{
  HELIB_NTIMER_START(ShiftMul);
    // const EncryptedArray& ea = m_context.getEA();
	// long nslots = ea.size();
	// cout<< nslots<<endl;
  long shift_sign = -1;
  if(shift_direction)
    shift_sign = 1;

  long e = 1;

  // shift and add
  while (e < m_expansionLen){
// while (e < nslots){
    Ctxt tmp = x;
    batch_shift_for_mul(tmp, start, e * shift_sign);
    x.multiplyBy(tmp);
    e <<=1;
  }
  HELIB_NTIMER_STOP(ShiftMul);
}

void Bridge::mapTo01_subfield(Ctxt& ctxt, long pow) const
{
// if pow = 1, this map operates on elements of the prime field F_p
	// pow is set to 1 by defualt
	HELIB_NTIMER_START(MapTo01);
	// get EncryptedArray
	const EncryptedArray& ea = m_context.getEA();
	// cout << "ctxt.getPtxtSpace(): "<< ctxt.getPtxtSpace() <<endl;
	// cout << "ea.getPAlgebra().getP(): "<< ea.getPAlgebra().getP() <<endl;
	cout << "getP(): "<< m_context.getP() <<endl;
	cout << "getR(): "<< m_context.getR() <<endl;
	cout << "getP2R(): "<< m_context.getPPowR() <<endl;

	// get p
	long p = ctxt.getPtxtSpace();
	long r = m_context.getR();
	long p2r = m_context.getPPowR();

	if (p % pow != 0)
		throw helib::LogicError("Exponent must divide p");

	HELIB_NTIMER_START(FERMAT);

	if (p != ea.getPAlgebra().getP()){
		cout << "[modified r!=1] map to 01, F_p^r" << endl;
		ctxt.frobeniusAutomorph(p2r-1);
		// throw helib::LogicError("mapTo01 not implemented for r>1");
	} // ptxt space is p^r for r>1 (p!=p^r)
	else if (p > 2){
		cout << "[modified r==1] map to 01, F_p" << endl;
		ctxt.power((p - 1) / pow); // set y = x^{p-1}
	}
	HELIB_NTIMER_STOP(FERMAT);
	HELIB_NTIMER_STOP(MapTo01);
}


void Bridge::is_zero(Ctxt& ctxt_res, const Ctxt& ctxt_z, long pow) const
{
  HELIB_NTIMER_START(EqualityCircuit);

	ctxt_res = ctxt_z;

	//compute mapTo01: (z_i)^{p^d-1}
	cout << "Mapping to 0 and 1" << endl;
	cout << "pow: " << pow << endl;
	print_decrypted(ctxt_res);
	mapTo01_subfield(ctxt_res, pow);

	if(m_verbose)
	{
		cout << "Fermat little (if res!=0, map res to 1)" << endl;
		print_decrypted(ctxt_res);
		cout << endl;
	}

	//cout << "Computing NOT" << endl;
	//compute 1 - mapTo01(z_i)
	ctxt_res.negate();
	ctxt_res.addConstant(ZZ(1));

	if(m_verbose)
	{
		cout << "return is_zero()" << endl;
		print_decrypted(ctxt_res);
		cout << endl;
	}

  HELIB_NTIMER_STOP(EqualityCircuit);
}


// interpolation related function
// taken from https://eprint.iacr.org/2021/315
void Bridge::compute_poly_params()
{
	// get p
	ZZ p = ZZ(m_context.getP());
	long p_long = conv<long>(p);

	// hardcoded babysteps sizes
	map<unsigned long, unsigned long> bs_nums
	{
		{5, 1},
		{7, 2}, // 4
		{11, 3}, // 3 (6), 1..4
  		{13, 3}, // 3 (6), 1..5
  		{17, 4}, // 4 (7), 1..5
  		{19, 3}, // 3 (8), 1..4
  		{23, 5}, // 5 (9), 3..6
  		{29, 5}, // 5 (10), 1..6
  		{31, 5}, // 5 (10), 4..6
  		{37, 5}, // 5 (12)
  		{47, 5}, // 5 (13), 2..11 
  		{61, 6}, // 6 (14), 4..8 
  		{67, 5}, // 5 (15), 4..8
  		{71, 4}, // 4 (15), 3..7
  		{101, 7}, // 7 (16), 4..8
  		{109, 7}, // 7 (19)
  		{131, 8}, // 8 (19), 4..11
  		{167, 10}, // 10 (21), 8..12
  		{173, 10},  // 10 (21), 8..12
  		{271, 9},  // 9 (26), 9..10
  		{401, 12},  // 12 (28), 9..14
  		{659, 11}	// 11 (41), 11..12
  	};

  	m_bs_num_comp = -1;
  	m_bs_num_min = -1;
  	if(bs_nums.count(p_long) > 0)
  	{
  		m_bs_num_comp = bs_nums[p_long];
  		m_bs_num_min = bs_nums[p_long];
  	}

  	// if p > 3, d = (p-3)/2
  	long d_comp = deg(m_univar_less_poly);
  	// if p > 3, d = (p-1)/2
  	long d_min = deg(m_univar_min_max_poly);

  	// How many baby steps: set sqrt(d/2), rounded up/down to a power of two

	// FIXME: There may be some room for optimization here: it may be possible to choose this number as something other than a power of two and still maintain optimal depth, in principle we can try all possible values of m_babystep_num between two consecutive powers of two and choose the one that gives the least number of multiplies, conditioned on minimum depth.

  	if (m_bs_num_comp <= 0) 
	{
		long kk = static_cast<long>(sqrt(d_comp/2.0)); //sqrt(d/2)
		m_bs_num_comp = 1L << NextPowerOfTwo(kk);

    	// heuristic: if #baby_steps >> kk then use a smaler power of two
    	if ((m_bs_num_comp==16 && d_comp>167) || (m_bs_num_comp>16 && m_bs_num_comp>(1.44*kk)))
      		m_bs_num_comp /= 2;
  	}
  	if (m_bs_num_min <= 0) 
	{
		long kk = static_cast<long>(sqrt(d_min/2.0)); //sqrt(d/2)
		m_bs_num_min = 1L << NextPowerOfTwo(kk);

    	// heuristic: if #baby_steps >> kk then use a smaler power of two
    	if ((m_bs_num_min==16 && d_min>167) || (m_bs_num_min>16 && m_bs_num_min>(1.44*kk)))
      		m_bs_num_min /= 2;
  	}

	if(m_verbose)
	{
		cout << "Number of baby steps for comparison: " << m_bs_num_comp << endl;
		cout << "Number of baby steps for min/max: " << m_bs_num_min << endl;
	}

	// #giant_steps = ceil(d/#baby_steps), d >= #giant_steps * #baby_steps
	m_gs_num_comp = divc(d_comp,m_bs_num_comp);
	m_gs_num_min = divc(d_min,m_bs_num_min);

	if(m_verbose)
	{
		cout << "Number of giant steps for comparison: " << m_bs_num_comp << endl;
		cout << "Number of giant steps for min/max: " << m_bs_num_min << endl;
	}      

	// If #giant_steps is not a power of two, ensure that poly is monic and that
	// its degree is divisible by #baby_steps, then call the recursive procedure

	// top coefficient is equal to (p^2 - 1)/8 mod p
	// its inverse is equal to -8 mod p
	m_top_coef_comp = LeadCoeff(m_univar_less_poly);
	m_top_coef_min = LeadCoeff(m_univar_min_max_poly);
	ZZ topInv_comp = ZZ(-8) % p; // the inverse mod p of the top coefficient of poly (if any)
	ZZ topInv_min = ZZ(-8) % p; // the inverse mod p of the top coefficient of poly (if any)
	bool divisible_comp = (m_gs_num_comp * m_bs_num_comp == d_comp); // is the degree divisible by #baby_steps?
	bool divisible_min = (m_gs_num_min * m_bs_num_min == d_min); // is the degree divisible by #baby_steps?

	// FIXME: There may be some room for optimization below: instead of
	// adding a term X^{n*k} we can add X^{n'*k} for some n'>n, so long
	// as n' is smaller than the next power of two. We could save a few
	// multiplications since giantStep[n'] may be easier to compute than
	// giantStep[n] when n' has fewer 1's than n in its binary expansion.

	m_extra_coef_comp = ZZ::zero();    // extra!=0 denotes an added term extra*X^{#giant_steps * #baby_steps}
	m_extra_coef_min = ZZ::zero();    // extra!=0 denotes an added term extra*X^{#giant_steps * #baby_steps}

	if (m_gs_num_comp != (1L << NextPowerOfTwo(m_gs_num_comp)))
	{
		if (!divisible_comp) 
		{  // need to add a term
	    	m_top_coef_comp = NTL::to_ZZ(1);  // new top coefficient is one
	    	topInv_comp = m_top_coef_comp;    // also the new inverse is one
	    	// set extra = 1 - current-coeff-of-X^{n*k}
	    	m_extra_coef_comp = SubMod(m_top_coef_comp, coeff(m_univar_less_poly, m_gs_num_comp * m_bs_num_comp), p);
	    	SetCoeff(m_univar_less_poly, m_gs_num_comp * m_bs_num_comp); // set the top coefficient of X^{n*k} to one
		}

		if (!IsOne(m_top_coef_comp)) 
		{
	    	m_univar_less_poly *= topInv_comp; // Multiply by topInv to make into a monic polynomial
	    	for (long i = 0; i <= m_gs_num_comp * m_bs_num_comp; i++) rem(m_univar_less_poly[i], m_univar_less_poly[i], p);
	    	m_univar_less_poly.normalize();
		}
	}

	/*
	cout << "Less-than poly: ";
	printZZX(cout, m_univar_less_poly, conv<long>(p));
	cout << endl;
	*/

	if (m_gs_num_min != (1L << NextPowerOfTwo(m_gs_num_min)))
	{
		if (!divisible_min) 
		{  // need to add a term
	    	m_top_coef_min = NTL::to_ZZ(1);  // new top coefficient is one
	    	topInv_min = m_top_coef_min;    // also the new inverse is one
	    	// set extra = 1 - current-coeff-of-X^{n*k}
	    	m_extra_coef_min = SubMod(m_top_coef_min, coeff(m_univar_min_max_poly, m_gs_num_min * m_bs_num_min), p);
	    	SetCoeff(m_univar_min_max_poly, m_gs_num_min * m_bs_num_min); // set the top coefficient of X^{n*k} to one
		}

		if (!IsOne(m_top_coef_min)) 
		{
	    	m_univar_min_max_poly *= topInv_min; // Multiply by topInv to make into a monic polynomial
	    	for (long i = 0; i <= m_gs_num_min * m_bs_num_min; i++) rem(m_univar_min_max_poly[i], m_univar_min_max_poly[i], p);
	    	m_univar_min_max_poly.normalize();
		}
	}

	/*
	cout << "Min-max poly: ";
	printZZX(cout, m_univar_min_max_poly, conv<long>(p));
	cout << endl;
	*/

	long top_deg = conv<long>(p-1) >> 1;
	m_baby_index = top_deg % m_bs_num_comp;
	m_giant_index = top_deg / m_bs_num_comp;
	if(m_baby_index == 0)
	{
		m_baby_index = m_bs_num_comp;
		m_giant_index -= 1;
	}
}

void Bridge::create_poly()
{
	cout << "Creating comparison polynomial" << endl;
	// get p
	unsigned long p = m_context.getP();;

	if(m_type == UNI)
	{
		// polynomial coefficient
		ZZ_p coef;
		coef.init(ZZ(p));

		// field element
		ZZ_p field_elem;
		field_elem.init(ZZ(p));

		// initialization of the univariate comparison polynomial
		m_univar_less_poly = ZZX(INIT_MONO, 0, 0);

		// loop over all odd coefficient indices
		for (long indx = 1; indx < p - 1; indx+=2)
		{ 
			// coefficient f_i = sum_a a^{p-1-indx} where a runs over [1,...,(p-1)/2]
			coef = 1;
			for(long a = 2; a <= ((p-1) >> 1); a++)
			{
			  field_elem = a;
			  coef += power(field_elem, p - 1 - indx);
			}

			m_univar_less_poly += ZZX(INIT_MONO, (indx-1) >> 1, rep(coef));
		}

		/*
		cout << "Less-than poly: ";
		printZZX(cout, m_univar_less_poly, p);
		cout << endl;
		*/

		m_univar_min_max_poly = m_univar_less_poly * ZZX(INIT_MONO, 1, 1);

		/*
		cout << "Min-max poly: ";
		printZZX(cout, m_univar_min_max_poly, p);
		cout << endl;
		*/

		compute_poly_params();
	}
	else if (m_type == TAN)
	{
		// computing the coefficients of the bivariate polynomial of Tan et al.
		m_bivar_less_coefs.SetDims(p,p);

		// y^{p-1}
		m_bivar_less_coefs[0][p-1] = ZZ(1);

		// (p+1)/2 * x^{(p-1)/2} * y^{(p-1)/2}
		m_bivar_less_coefs[(p-1) >> 1][(p-1) >> 1] = ZZ((p+1) >> 1);

		// iterator
		ZZ_p field_elem;
		field_elem.init(ZZ(p));

		// inner sum
		ZZ_p inner_sum;
		inner_sum.init(ZZ(p));
		
		// outer sum
		ZZ_p outer_sum;
		outer_sum.init(ZZ(p));

		for (long i = 1; i < p; i++)
		{
			for (long j = 1; j < p; j++)
			{
				// x^i * y^i have the zero coefficient except for i = (p-1)/2
				if (i == j)
					continue;

				outer_sum = 0;
				// sum_{a=1}^{p-1} a^{p-1-i} sum_{b = a+1}^{p-1} b^{p-1-j} 
				for (long a = 1; a < p; a++)
				{
					inner_sum = 0;
					// sum_{b = a+1}^{p-1} b^{p-1-j} 
					for (long b = a+1; b < p; b++)
					{
						// b^{p-1-j}
						field_elem = b;
						field_elem = power(field_elem, p - 1 - j);

						inner_sum += field_elem;
					}
					// a^{p-1-i}
					field_elem = a;
					field_elem = power(field_elem, p - 1 - i);

					inner_sum *= field_elem;
					outer_sum += inner_sum;
				}
				m_bivar_less_coefs[i][j] = rep(outer_sum);
			}
		}

		cout << "Bivariate coefficients" << endl << m_bivar_less_coefs << endl;

		if (m_verbose)
		{
			cout << "Comparison polynomial: " << endl;
			printZZX(cout, m_univar_less_poly, (p-1)>>1);
			cout << endl;
		}
	}

	cout << "Comparison polynomial is created" << endl;
}

void Bridge::evaluate_univar_less_poly(Ctxt& ret, Ctxt& ctxt_p_1, const Ctxt& x) const
{
	HELIB_NTIMER_START(ComparisonCircuitUnivar);
	// get p
	ZZ p = ZZ(m_context.getP());

	if (p > ZZ(3)) //if p > 3, use the generic Paterson-Stockmeyer strategy
	{
	  // z^2
	  	Ctxt x2 = x;
	  	x2.square();

		DynamicCtxtPowers babyStep(x2, m_bs_num_comp);
		const Ctxt& x2k = babyStep.getPower(m_bs_num_comp);

		DynamicCtxtPowers giantStep(x2k, m_gs_num_comp);

		// Special case when #giant_steps is a power of two
		if (m_gs_num_comp == (1L << NextPowerOfTwo(m_gs_num_comp))) 
		{
			//cout << "I'm computing degPowerOfTwo" << endl;
	    	degPowerOfTwo(ret, m_univar_less_poly, m_bs_num_comp, babyStep, giantStep);
	    }
	    else
	    {
		  	recursivePolyEval(ret, m_univar_less_poly, m_bs_num_comp, babyStep, giantStep);

		  	if (!IsOne(m_top_coef_comp)) 
		  	{
		    	ret.multByConstant(m_top_coef_comp);
			}

			if (!IsZero(m_extra_coef_comp)) 
			{ // if we added a term, now is the time to subtract back
		    	Ctxt topTerm = giantStep.getPower(m_gs_num_comp);
		    	topTerm.multByConstant(m_extra_coef_comp);
		    	ret -= topTerm;
			}
		}
		ret.multiplyBy(x);

		// TODO: depth here is not optimal
		Ctxt top_term = babyStep.getPower(m_baby_index);
		top_term.multiplyBy(giantStep.getPower(m_giant_index));

		ctxt_p_1 = top_term; 

		top_term.multByConstant(ZZ((p+1)>> 1));

		ret += top_term;

		/*
		cout << "Computed baby steps" << endl;
		for(int i = 0; i < babyStep.size(); i++)
		{
			cout << i + 1 << ' ' << babyStep.isPowerComputed(i+1) << endl; 
		}

		cout << "Computed giant steps" << endl;
		for(int i = 0; i < giantStep.size(); i++)
		{
			cout << i + 1 << ' ' << giantStep.isPowerComputed(i+1) << endl;
		}
		*/
	}
	else //circuit for p=3
	{
		ret = x;

		ctxt_p_1 = x;
		ctxt_p_1.square();

		Ctxt top_term = ctxt_p_1;
		top_term.multByConstant(ZZ(2));

		ret += top_term;
	}
	HELIB_NTIMER_STOP(ComparisonCircuitUnivar);
}

Bridge::Bridge(const Context& context, CircuitType type, unsigned long d, unsigned long expansion_len, const SecKey& sk, bool verbose):
	m_context(context), m_type(type), m_slotDeg(d), m_expansionLen(expansion_len), m_sk(sk), m_pk(sk), m_verbose(verbose)
{
	std::cout <<"[construct] gen mask" <<  std::endl;
	create_all_shift_masks();
	std::cout <<"[construct] gen interpolation poly" <<  std::endl;
	create_poly();
	std::cout <<"[construct] done" <<  std::endl;
}

const DoubleCRT& Bridge::get_mask(double& size, long index) const
{
	size = m_mulMasksSize[index];
	return m_mulMasks[index];
}

const ZZX& Bridge::get_less_than_poly() const
{
	return m_univar_less_poly;
}

void Bridge::compare(Ctxt& ctxt_res, const Ctxt& ctxt_x) const{
	// vector of comparison result over F_p
	vector<Ctxt> ctxt_less_p;
	// vector of comparison result over F_p
	vector<Ctxt> ctxt_eq_p;

	unsigned long p = m_context.getP();
	unsigned long ord_p = m_context.getOrdP();
	unsigned long r = m_context.getR();
	unsigned long p2r = m_context.getPPowR();

	Ctxt ctxt_z = ctxt_x;

	// decompose z to mod p digits
	vector<Ctxt> ctxt_z_p;
	HELIB_NTIMER_START(Reduction);
	reduce(ctxt_z_p, ctxt_z, r);
	HELIB_NTIMER_STOP(Reduction);

	std::cout << "[Reduction] Reduced to: " <<ctxt_z_p.size() << " digits"  << std::endl;
	CheckCtxt(ctxt_z_p[0], "[Reduction] Reduced digits (in beFV)");

	cout << "[beFV] Interpolation: compute the less-than and equality functions modulo p" << endl;
	for (long iCoef = 0; iCoef < r; iCoef++){
		Ctxt ctxt_tmp = Ctxt(ctxt_z.getPubKey());
		Ctxt ctxt_tmp_eq = Ctxt(ctxt_z.getPubKey());

		// compute polynomial function for 'z < 0'
		// cout << "Compute univariate comparison polynomial" << endl;
		evaluate_univar_less_poly(ctxt_tmp, ctxt_tmp_eq, ctxt_z_p[iCoef]);

		if(m_verbose)
		{
			cout << "[beFV] Result of the less-than function" << endl;
			print_decrypted(ctxt_tmp);
			cout << endl;
		}
		ctxt_less_p.push_back(ctxt_tmp);

		//cout << "Computing NOT" << endl;
		//compute 1 - mapTo01(r_i*(x_i - y_i))
		ctxt_tmp_eq.negate();
		ctxt_tmp_eq.addConstant(ZZ(1));
		if(m_verbose)
		{
			cout << "[beFV] Result of the equality function" << endl;
			print_decrypted(ctxt_tmp_eq);
			cout << endl;
		}

		ctxt_eq_p.push_back(ctxt_tmp_eq);
	}	

	// HELIB_NTIMER_STOP(Comparison);
	
	HELIB_NTIMER_START(Aggregation);
	// digits result -> integer result
	cout << "[beFV] Aggregation" <<endl;
	Ctxt ctxt_less = ctxt_less_p[m_slotDeg-1];
	Ctxt ctxt_eq = ctxt_eq_p[m_slotDeg-1];

	for (long iCoef = m_slotDeg-2; iCoef >= 0; iCoef--)
	{
		Ctxt tmp = ctxt_eq;
		tmp.multiplyBy(ctxt_less_p[iCoef]);
		ctxt_less += tmp;

		ctxt_eq.multiplyBy(ctxt_eq_p[iCoef]);
	}

	if(m_expansionLen == 1)
	{
		ctxt_res = ctxt_less;
		return;
	}
	HELIB_NTIMER_STOP(Aggregation);
}


void Bridge::print_decrypted(const Ctxt& ctxt) const
{
	// get EncryptedArray
	const EncryptedArray& ea = m_context.getEA();

	// get order of p
	unsigned long ord_p = m_context.getOrdP();

    long nSlots = ea.size();
    vector<ZZX> decrypted(nSlots);
    ea.decrypt(ctxt, m_sk, decrypted);

    for(int i = 0; i < nSlots; i++)
    {
      printZZX(cout, decrypted[i], ord_p);
    }
	cout << endl;
	// helib::Ptxt<helib::BGV> new_plaintext_result(m_context);
	// m_sk.Decrypt(new_plaintext_result, ctxt);
	// std::cout << "[DEC] ct: " << new_plaintext_result << std::endl;

}

// Compute a degree-p polynomial poly(x) s.t. for any t<e and integer z of the
// form z = z0 + p^t*z1 (with 0<=z0<p), we have poly(z) = z0 (mod p^{t+1}).
//
// We get poly(x) by interpolating a degree-(p-1) polynomial poly'(x)
// s.t. poly'(z0)=z0 - z0^p (mod p^e) for all 0<=z0<p, and then setting
// poly(x) = x^p + poly'(x).
static void buildDigitPolynomial(NTL::ZZX& result, long p, long e)
{
  if (p < 2 || e <= 1)
    return; // nothing to do
  HELIB_TIMER_START;
  long p2e = NTL::power_long(p, e); // the integer p^e

  // Compute x - x^p (mod p^e), for x=0,1,...,p-1
  NTL::vec_long x(NTL::INIT_SIZE, p);
  NTL::vec_long y(NTL::INIT_SIZE, p);
  long bottom = -(p / 2);
  for (long j = 0; j < p; j++) {
    long z = bottom + j;
    x[j] = z;
    y[j] =
        z - NTL::PowerMod((z < 0 ? z + p2e : z), p, p2e); // x - x^p (mod p^e)

    while (y[j] > p2e / 2)
      y[j] -= p2e;
    while (y[j] < -(p2e / 2))
      y[j] += p2e;
  }
  interpolateMod(result, x, y, p, e);
  // interpolating p points, should get deg<=p-1
  assertTrue(deg(result) < p, "Interpolation error.  Degree too high.");
  SetCoeff(result, p); // return result = x^p + poly'(x)
  //  cerr << "# digitExt mod "<<p<<"^"<<e<<"="<<result<<endl;
  HELIB_TIMER_STOP;
}

void Bridge::reduce(std::vector<Ctxt>& digits, const Ctxt& c, long r) const
{
	std::cout<<"[Reduction] FV to beFV"<<std::endl;
	const Context& context = c.getContext();
	long rr = c.effectiveR();
	if (r <= 0 || r > rr)
		r = rr; // how many digits to extract

	long p = context.getP();
	NTL::ZZX x2p;
	if (p > 3) {
		buildDigitPolynomial(x2p, p, r);
	}
	std::cout<<"[Reduction] p=" << p<<std::endl;
	std::cout<<"[Reduction] r=" << rr<<std::endl;
	std::cout<<"[Reduction] ptxtSpace=" << c.getPtxtSpace() <<std::endl;

	//test G_e
	// NTL::ZZX magic;
	// helib::compute_magic_poly(magic, p, r);
	// std::cout << "# magic polynomial "<<p<<"^"<<r<<"="<<magic<< std::endl;

	Ctxt tmp(c.getPubKey(), c.getPtxtSpace());
	digits.resize(r, tmp); // allocate space

	#ifdef HELIB_DEBUG
	fprintf(stderr, "***\n");
	#endif
	for (long i = 0; i < r; i++) {
		tmp = c;
		for (long j = 0; j < i; j++) {
			if (p == 2)
				digits[j].square();
			else if (p == 3)
				digits[j].cube();
			else
				polyEval(digits[j], x2p, digits[j]);
				// "in spirit" digits[j] = digits[j]^p
			tmp -= digits[j];
			tmp.divideByP();
		}
		digits[i] = tmp; // needed in the next round
	}

	// mod p digits
	// the digits will be mod p
	for (long i = 0; i < r-1; i++) {
		for (long j=i; j<r-1; j++){
			digits[i].divideModByP();
		}
	}

}

static void compute_a_vals(NTL::Vec<NTL::ZZ>& a, long p, long e)
// computes a[m] = a(m)/m! for m = p..(e-1)(p-1)+1,
// as defined by Chen and Han.
// a.length() is set to (e-1)(p-1)+2

{
  NTL::ZZ p_to_e = NTL::power_ZZ(p, e);
  NTL::ZZ p_to_2e = NTL::power_ZZ(p, 2 * e);

  long len = (e - 1) * (p - 1) + 2;

  NTL::ZZ_pPush push(p_to_2e);

  NTL::ZZ_pX x_plus_1_to_p = power(NTL::ZZ_pX(NTL::INIT_MONO, 1) + 1, p);
  NTL::ZZ_pX denom =
      InvTrunc(x_plus_1_to_p - NTL::ZZ_pX(NTL::INIT_MONO, p), len);
  NTL::ZZ_pX poly = MulTrunc(x_plus_1_to_p, denom, len);
  poly *= p;

  a.SetLength(len);

  NTL::ZZ m_fac(1);
  for (long m = 2; m < p; m++) {
    m_fac = MulMod(m_fac, m, p_to_2e);
  }

  for (long m = p; m < len; m++) {
    m_fac = MulMod(m_fac, m, p_to_2e);
    NTL::ZZ c = rep(coeff(poly, m));
    NTL::ZZ d = GCD(m_fac, p_to_2e);
    if (d == 0 || d > p_to_e || c % d != 0)
      throw RuntimeError("cannot divide");
    NTL::ZZ m_fac_deflated = (m_fac / d) % p_to_e;
    NTL::ZZ c_deflated = (c / d) % p_to_e;
    a[m] = MulMod(c_deflated, InvMod(m_fac_deflated, p_to_e), p_to_e);
  }
}


static void compute_magic_poly(NTL::ZZX& poly1, long p, long e)
{
  HELIB_TIMER_START;

  NTL::Vec<NTL::ZZ> a;

  compute_a_vals(a, p, e);

  NTL::ZZ p_to_e = NTL::power_ZZ(p, e);
  long len = (e - 1) * (p - 1) + 2;

  NTL::ZZ_pPush push(p_to_e);

  NTL::ZZ_pX poly(0);
  NTL::ZZ_pX term(1);
  NTL::ZZ_pX X(NTL::INIT_MONO, 1);

  poly = 0;
  term = 1;

  for (long m = 0; m < p; m++) {
    term *= (X - m);
  }

  for (long m = p; m < len; m++) {
    poly += term * NTL::conv<NTL::ZZ_p>(a[m]);
    term *= (X - m);
  }

  // replace poly by poly(X+(p-1)/2) for odd p
  if (p % 2 == 1) {
    NTL::ZZ_pX poly2(0);

    for (long i = deg(poly); i >= 0; i--)
      poly2 = poly2 * (X + (p - 1) / 2) + poly[i];

    poly = poly2;
  }

  poly = X - poly;
  poly1 = NTL::conv<NTL::ZZX>(poly);
}

void Bridge::lift(Ctxt& res, const Ctxt& c, long r) const{

	std::cout<<"[Lifting] beFV to FV"<<std::endl;

	const Context& context = c.getContext();
	long p = context.getP();
	NTL::ZZX Ge;
	compute_magic_poly(Ge, p, r+1);
	polyEval(res, Ge, c);
	CheckCtxt(res, "[Lifting] Logic result after lifting (in FV)");

}

// test compare on Z_pr
// compare with 0
void Bridge::test_bridge(long runs) const{
	//reset timers
  	setTimersOn();

	//generate random data
	random_device rd;
  	mt19937 eng(rd());
  	uniform_int_distribution<unsigned long> distr_u;
  	uniform_int_distribution<long> distr_i;

	// get EncryptedArray
	const EncryptedArray& ea = m_context.getEA();
	//extract number of slots
	long nslots = ea.size();
	unsigned long p = m_context.getP();
	unsigned long ord_p = m_context.getOrdP();
	unsigned long r = m_context.getR();
	unsigned long p2r = m_context.getPPowR();
	// input digit is from [-p/2, p/2]
	unsigned long enc_base = p;
	unsigned long input_range = p2r;

	long min_capacity = 1000;
  	long capacity;
  	for (int run = 0; run < runs; run++){
        printf("\n");
    	printf("Run %d started\n", run);

		vector<ZZX> expected_result(nslots);
    	vector<ZZX> decrypted(nslots);
		vector<int> expected_sign(nslots);

		// check x > 0 ?
        // The comparison logic depends on how integers / fixed-point numbers are encoded in the plaintext space
        // For x \in [p/2, p], we return 1. (In spirit, we encode 0 as p/2) 
		long input_x;
		std::vector<long> x_vec(nslots);

		for (int i = 0; i < nslots; i++){
			input_x = distr_u(eng) % input_range;
			x_vec[i] = input_x;

			if (input_x > (p2r/2)){
				expected_result[i] = ZZX(INIT_MONO, 0, 1);
				expected_sign[i] = 1;
			}
			else{
				expected_result[i] = ZZX(INIT_MONO, 0, 0);
				expected_sign[i] = 0;
			}
		}

		// encrypt
		Ctxt ctxt_x(m_pk);
		ea.encrypt(ctxt_x, m_pk, x_vec);
		CheckCtxt(ctxt_x, "[Initial] Airthmetic value (in FV)");

        // res = x > 0?
        // x in FV -> reduce -> interpolation in beFV -> result in beFV -> lift -> result in FV
		Ctxt ctxt_res(m_pk);
		compare(ctxt_res, ctxt_x);
		CheckCtxt(ctxt_res, "[beFV] Logic result (in beFV)");
		// check logic operation in beFV
        // whether the comparison succeeded
		ea.decrypt(ctxt_res, m_sk, decrypted);
        for(int i = 0; i < nslots; i++)
        { 
			if (decrypted[i] != expected_result[i])
			{
				printf("Slot %ld: ", i);
				cout << endl;
				cout << "Failure" << endl;
				return;
			}
        }
        cout << "[beFV] Success" << endl;

		HELIB_NTIMER_START(Lifting);
		ctxt_res.multiplyModByP2R();
		Ctxt ctxt_res_fv(m_pk);
		lift(ctxt_res_fv, ctxt_res, r);
		HELIB_NTIMER_STOP(Lifting);

		// compute relu(x) = x \times (x>0?)
		Ctxt ctxt_res_relu(m_pk);
		ctxt_res_relu = ctxt_x;
		ctxt_res_relu.addConstant(-long(p2r)/2);
		ctxt_res_relu.multiplyBy(ctxt_res_fv);

        cout << endl;
        printNamedTimer(cout, "Reduction");
        printNamedTimer(cout, "ComparisonCircuitUnivar");
		printNamedTimer(cout, "Aggregation");
		printNamedTimer(cout, "Lifting");
        printNamedTimer(cout, "ReLU");
		cout << endl;

		helib::Ptxt<helib::BGV> result(m_context);
		m_sk.Decrypt(result, ctxt_res_relu);
		helib::Ptxt<helib::BGV> result_sign(m_context);
		m_sk.Decrypt(result_sign, ctxt_res_fv);
		auto decrypted_res = result.getSlotRepr();
		auto decrypted_res_sign = result_sign.getSlotRepr();

		cout << "Input: " ;
		for(int i = 0; i < 10; i++){
			cout << int(x_vec[i])-int(p2r)/2 << " ";
		}
		cout<<endl;

		cout << "Decrypted ReLU: " ;
		for(int i = 0; i < 10; i++){
			cout << decrypted_res[i] << " ";
		}
		cout<<endl;

		cout << "Expected Sign: " ;
		for(int i = 0; i < 10; i++){
			cout << expected_sign[i] << " ";
		}
		cout<<endl;

		cout << "Decrypted Sign: " ;
		for(int i = 0; i < 10; i++){
			cout << decrypted_res_sign[i] << " ";
		}
		cout<<endl;

	}



}

// test compare on Z_pr
// compare with 0
void Bridge::m_test_bridge(long runs) const{
	//reset timers
  	setTimersOn();

	//generate random data
	random_device rd;
  	mt19937 eng(rd());
  	uniform_int_distribution<unsigned long> distr_u;
  	uniform_int_distribution<long> distr_i;

	// get EncryptedArray
	const EncryptedArray& ea = m_context.getEA();
	//extract number of slots
	long nslots = ea.size();
	unsigned long p = m_context.getP();
	unsigned long ord_p = m_context.getOrdP();
	unsigned long r = m_context.getR();
	unsigned long p2r = m_context.getPPowR();
	// input digit is from [-p/2, p/2]
	unsigned long enc_base = p;
	unsigned long input_range = p2r;

	long min_capacity = 1000;
  	long capacity;
  	for (int run = 0; run < runs; run++){
        printf("\n");
    	printf("Run %d started\n", run);

		vector<ZZX> expected_result(nslots);
    	vector<ZZX> decrypted(nslots);
		vector<int> expected_sign(nslots);

		// check x > 0 ?
        // The comparison logic depends on how integers / fixed-point numbers are encoded in the plaintext space
        // For x \in [p/2, p], we return 1. (In spirit, we encode 0 as p/2) 
		long input_x;
		std::vector<long> x_vec(nslots);

		for (int i = 0; i < nslots; i++){
			input_x = distr_u(eng) % input_range;
			x_vec[i] = input_x;

			if (input_x > (p2r/2)){
				expected_result[i] = ZZX(INIT_MONO, 0, 1);
				expected_sign[i] = 1;
			}
			else{
				expected_result[i] = ZZX(INIT_MONO, 0, 0);
				expected_sign[i] = 0;
			}
		}

		// encrypt
		Ctxt ctxt_x(m_pk);
		Ctxt ctxt_res(m_pk);
		ea.encrypt(ctxt_x, m_pk, x_vec);
		CheckCtxt(ctxt_x, "[Initial] Airthmetic value (in FV)");
		HELIB_NTIMER_START(Linear);
		ctxt_res = ctxt_x;
		long scale = 2;
		// ctxt_l.multiplyBy(tmp);
		for (long i =0; i<128; i++){
			ctxt_res.multByConstant(scale);
		}
		for (long i=0; i<20; i++){
			long shift = 2;
			if(shift == 0)
				break;
			ea.rotate(ctxt_res, shift);
		}
		HELIB_NTIMER_STOP(Linear);

        // res = x > 0?
        // x in FV -> reduce -> interpolation in beFV -> result in beFV -> lift -> result in FV
		compare(ctxt_res, ctxt_x);
		CheckCtxt(ctxt_res, "[beFV] Logic result (in beFV)");
		// check logic operation in beFV
        // whether the comparison succeeded
		ea.decrypt(ctxt_res, m_sk, decrypted);
        for(int i = 0; i < nslots; i++)
        { 
			if (decrypted[i] != expected_result[i])
			{
				printf("Slot %ld: ", i);
				cout << endl;
				cout << "Failure" << endl;
				return;
			}
        }
        cout << "[beFV] Success" << endl;

		HELIB_NTIMER_START(Lifting);
		ctxt_x.multiplyModByP2R();
		Ctxt ctxt_res_fv(m_pk);
		lift(ctxt_res_fv, ctxt_x, r-1);
		HELIB_NTIMER_STOP(Lifting);

        cout << endl;
		printNamedTimer(cout, "Linear");
        printNamedTimer(cout, "Reduction");
        printNamedTimer(cout, "ComparisonCircuitUnivar");
		printNamedTimer(cout, "Aggregation");
		printNamedTimer(cout, "Lifting");
		cout << endl;
	}
}
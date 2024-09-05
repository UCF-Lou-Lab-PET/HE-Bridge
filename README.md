# HEBridge: Connecting arithmetic and logic operations in FV-style HE schemes

This is a demo code for HEBridge: Connecting arithmetic and logic operations in FV-style
HE schemes. Our implementaion is based on HElib and polynomial interpolation-based methods for FV-style HE schemes, such as [Faster Comparison](https://eprint.iacr.org/2021/315). Please note that the codes are still under development and should not be used in any security-sensitive environments. Future version of this code will be released soon.

## Installation
[HElib](https://github.com/homenc/HElib) should be first installed. Additionally, add the path of HElib to the makefile

    set(helib_DIR ".../helib_install/helib_pack/share/cmake/helib")

and then, tun

    cmake .

finally, build the project by

    make

## Testing
### Switch from FV to beFV and evaluate logic operaions
A simple test includes (1) the reduction function, which map the input from FV to beFV. (2) the evaluation in beFV, i.e., an interpolation polynomial of degree p. (3) Lifting from beFV to FV. To run a test, set the parameter as arguments:
  
    ./bin/m_test p=17 r=2 m=13201 b=256 t=64
    
where:
- **p**: a prime number for the base plaintext modulus in beFV.
- **r**: a positive integer denoting the power of prime. The plaintext modulus in FV is \(p^r\).
- **q**: the initial ciphertext modulus.
- **t**: the Hamming weight of sk.
- **m**: the cyclotomic order of the polynomial ring.
- **n**: the degree of the polynomial ring, where \(n = \Phi(m)\).
- **d**: the multiplicative order of \(p\) modulo \(m\). \(d\) is the smallest positive integer such that \(p^d = 1 mod m\).
- **ℓ**: the number of SIMD slots. ℓ satisfies the following relationship with \(n\) and \(d\): \(n=ℓd\).

We use uni-variate interpolation polynomials only. Parameter can be chosen differently according to the input bit-width.
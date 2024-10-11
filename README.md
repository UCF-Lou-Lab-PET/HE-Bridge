# HEBridge: Connecting arithmetic and logic operations in FV-style HE schemes

This is a demo code for HEBridge: Connecting Arithmetic and Logic Operations in
FV-style HE Schemes (in CCS WAHC'24).

Our implementaion is based on HElib and polynomial interpolation-based methods for FV-style HE schemes, such as [Faster Comparison](https://eprint.iacr.org/2021/315). Please note that the codes are still under development and should not be used in any security-sensitive environments. Future version of this code will be released soon.

## Build and installation
Before building and installing HElib, please add the functions in ```src/Ctxt.h``` to ```include/helib/Ctxt.h```. Then, follow the instruciton in [HElib](https://github.com/homenc/HElib) for installation. After installation, add the path of HElib to the makefile as:

    set(helib_DIR ".../helib_install/helib_pack/share/cmake/helib")

Once HElib is installed, you can build this project. Under the  ```./src``` directory, run

    mkdir build && cd build

and then, run:

    cmake .

finally, build the project by:

    make

## Testing
### Switch from FV to beFV and evaluate logic operaions
A simple test includes (1) the reduction function, which map the input from FV to beFV. (2) the evaluation in beFV, i.e., an interpolation polynomial of degree p. (3) Lifting from beFV to FV. We provide some parameter settings as follows:
  
    ./bin/bridge_circuit p=3 r=4 m=14401 b=300 t=64

    ./bin/bridge_circuit p=5 r=4 m=25351 b=512 t=64
    
    ./bin/bridge_circuit p=67 r=2 m=31159 b=690 t=64
    
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

## Running HEBridge
To make sure you are in the good position, we provide the logs from a successful running. By executing the command above, you are excpeted to see the following logs:

Encryption paramaters. This step indicates HElib is propoerly installed.

    m=31159, p=67, r=2, bits=690, c=2, skHwt=64
    Initialising context object...
    m = 31159, p = 67, phi(m) = 31158
    ord(p) = 54
    normBnd = 1.27324
    polyNormBnd = 1.27324
    factors = [31159]
    generator 21 has order (== Z_m^*) of 577
    [param] getP(): 67
    [param] getR(): 2
    [param] getP2R(): 4489
    [param] order(p) = 54

Initialization. This step generates the interpolation polynomials given the parameters.

    Creating secret key...
    Generating key-switching matrices...
    [test] comparator begin
    [construct] gen mask
    All masks are created
    [construct] gen interpolation poly
    Creating comparison polynomial
    Comparison polynomial is created
    [construct] done

Running HEBridge test.

    Run 0 started
    [Initial] Airthmetic value (in FV), capacity=655.515, p^r=4489
    [Reduction] FV to beFV
    [Reduction] p=67
    [Reduction] r=2
    [Reduction] ptxtSpace=4489
    [Reduction] Reduced to: 2 digits
    [Reduction] Reduced digits (in beFV), capacity=476.823, p^r=67
    [beFV] Interpolation: compute the less-than and equality functions modulo p
    [beFV] Aggregation
    [beFV] Logic result (in beFV), capacity=289.268, p^r=67
    [beFV] Success
    [Lifting] beFV to FV
    [Lifting] Logic result after lifting (in FV), capacity=61.332, p^r=4489

Runtime breakdown. We use the built-in timer in HElib to log the runtime of each component.
 
    Reduction: 7.25465 / 1 = 7.25465 
    ComparisonCircuitUnivar: 10.8985 / 2 = 5.44924 
    Aggregation: 0.562972 / 1 = 0.562972 
    Lifting: 5.72996 / 1 = 5.72996 
    ReLU: 18.7374 / 1 = 18.7374

Result verification. We generate ramdom numbers for testing. Under SIMD, ℓ ReLU can be computed in paralell. Here we show the first 10 random inputs and the decrypted results.

    Input: 1643 1203 -1478 1675 -1772 -75 976 2113 1592 189
    Decrypted ReLU: [1643] [1203] [0] [1675] [0] [0] [976] [2113] [1592] [189]
    Expected Sign: 1 1 0 1 0 0 1 1 1 1
    Decrypted Sign: [1] [1] [0] [1] [0] [0] [1] [1] [1] [1]
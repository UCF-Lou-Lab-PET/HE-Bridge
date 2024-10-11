#include <iostream>
#include <time.h>
#include <random>
#include <helib/helib.h>
#include <helib/debugging.h>
#include <helib/Context.h>
#include <helib/polyEval.h>
#include "tools.h"
#include "bridge.h"
#include "ArgMapping.h"

using namespace std;
using namespace NTL;
using namespace helib;
using namespace he_bridge;

int main(int argc, char *argv[]) {
    // Plaintext prime modulus
    unsigned long p = 3;
    // Cyclotomic polynomial - defines phi(m)
    unsigned long m = 31 * 41;
    // Hensel lifting (default = 1)
    unsigned long r = 1;
    // Number of bits of the modulus chain (budget)
    unsigned long bits = 256;
    // Number of columns of Key-Switching matrix (default = 2 or 3)
    unsigned long c = 2;
    // int c_m = 100;
    unsigned long t = 64;

    /************* Argument Parsing  ************/
    /********************************************/
    ArgMapping amap;
    amap.arg("p", p, "the plaintext modulus");
    amap.arg("r", r, "the lifting parameter for plaintext space p^r");
    amap.arg("m", m, "the order of the cyclotomic ring");
    // b represents L, the level
    amap.arg("b", bits, "the bitsize of the ciphertext modulus in ciphertexts");
    amap.arg("c", c, "Number of columns of Key-Switching matrix");
    amap.arg("t", t, "The hamming weight of sk");
    amap.parse(argc, argv);
    unsigned long d = r;
    std::cout << "m=" << m
            << ", p=" << p
            << ", r=" << r
            << ", bits=" << bits
            << ", c=" << c
            << ", skHwt=" << t
            << std::endl;
    std::cout << "Initialising context object..." << std::endl;
    helib::Context context = helib::ContextBuilder<helib::BGV>()
                                .m(m)
                                .p(p)
                                .r(r)
                                .bits(bits)
                                .c(c)
                                .skHwt(t)
                                .build();

    // Print the context
    context.getZMStar().printout();
    cout << "[param] getP(): "<< context.getP() <<endl;
	cout << "[param] getR(): "<< context.getR() <<endl;
	cout << "[param] getP2R(): "<< context.getPPowR() <<endl;
    cout << "[param] order(p) = " << context.getOrdP() <<endl;
    cout << endl;

    // Secret key management
    cout << "Creating secret key..." << endl;
    // Create a secret key associated with the context
    SecKey secret_key(context);
    // Generate the secret key
    secret_key.GenSecKey();
    // Compute key-switching matrices that we need
    cout << "Generating key-switching matrices..." << endl;
    addSome1DMatrices(secret_key);
    addFrbMatrices(secret_key);
    if (d > 1)
        addFrbMatrices(secret_key); //might be useful only when d > 1
    helib::PubKey& public_key = secret_key;
    const helib::EncryptedArray& ea = context.getEA();
    unsigned long expansion_len = 1;
    bool verbose = false;
    CircuitType type = UNI;


    cout << "[test] comparator begin" << endl; 

    Bridge bridge(context, type, r, expansion_len, secret_key, verbose);

    int runs = 1;
    bridge.test_bridge(runs);

    return 0;
}

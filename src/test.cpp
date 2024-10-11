#include <iostream>
#include <helib/helib.h>
#include <helib/matmul.h>
#include <ctime>
#include "ArgMapping.h"
#include <random>

using namespace std;

// Function to generate a vector of random numbers
std::vector<long> generateRandomVector(long n, unsigned long p) {
    std::vector<long> result(n);
    std::random_device rd;  // Obtain a random number from hardware
    std::mt19937 gen(rd()); // Seed the generator
    std::uniform_int_distribution<> distr(0, p/2-1); // Define the range

    for (int i = 0; i < n; i++) {
        result[i] = distr(gen); // Generate random numbers and store them in the vector
    }
    return result;
}

int main(int argc, char* argv[])
{
    /*  Example of BGV scheme  */

    // Plaintext prime modulus
    unsigned long p = 3;
    // Cyclotomic polynomial - defines phi(m)
    unsigned long m = 31 * 41;
    // Hensel lifting (default = 1)
    unsigned long r = 3;
    // Number of bits of the modulus chain (budget)
    unsigned long bits = 16;
    // Number of columns of Key-Switching matrix (default = 2 or 3)
    unsigned long c = 2;
    int c_m = 100;
    unsigned long t = 64;
    std::vector<long> mvec = {31, 41};
    std::vector<long> gens = {1026, 249};
    std::vector<long> ords = {30, -2};


    /************* Argument Parsing  ************/
    /********************************************/
    ArgMapping amap;
    amap.arg("p", p, "the plaintext modulus");
    amap.arg("r", r, "the dimension of a vector space over a finite field");
    amap.arg("m", m, "the order of the cyclotomic ring");
    // b represents L, the level
    amap.arg("b", bits, "the bitsize of the ciphertext modulus in ciphertexts");
    amap.arg("c", c, "Number of columns of Key-Switching matrix");
    amap.parse(argc, argv);

    m = helib::FindM(128, bits, c, p, 0, 1, 0);
    // clang-format off
    std::cout << "m=" << m
            << ", p=" << p
            << ", r=" << r
            << ", bits=" << bits
            << ", c=" << c
            << ", skHwt=" << t
            << ", c_m=" << c_m
            // << ", mvec=" << helib::vecToStr(mvec)
            // << ", gens=" << helib::vecToStr(gens)
            // << ", ords=" << helib::vecToStr(ords)
            << std::endl;
    
    // Initialize context
    // This object will hold information about the algebra created from the
    // previously set parameters
    std::cout << "Initialising context object..." << std::endl;
    helib::Context context = helib::ContextBuilder<helib::BGV>()
                                .m(m)
                                .p(p)
                                .r(r)
                                // .gens(gens)
                                // .ords(ords)
                                .bits(bits)
                                .c(c)
                                // .bootstrappable(true)
                                // .skHwt(t)
                                // .mvec(mvec)
                                // .thickboot()
                                .build();


    // Print the context
    context.printout();
    std::cout << "Creating secret key..." << std::endl;
    helib::SecKey secret_key(context);
    secret_key.GenSecKey();
    
    std::cout << "Generating key-switching matrices..." << std::endl;
    addSome1DMatrices(secret_key);
    addFrbMatrices(secret_key);
    // Generate bootstrapping data
    // secret_key.genRecryptData();
    helib::PubKey& public_key = secret_key;
    const helib::EncryptedArray& ea = context.getEA();
    unsigned long d = r;
    long nslots = ea.size();
    std::cout << "Number of slots: " << nslots << std::endl;
    helib::Ctxt ctxt1(public_key);
    helib::Ctxt ctxt2(public_key);
    cout << "[param] ptxtspace " << ctxt1.getPtxtSpace() <<endl;
    cout << "[param] getP(): "<< context.getP() <<endl;
	cout << "[param] getR(): "<< context.getR() <<endl;
	cout << "[param] getP2R(): "<< context.getPPowR() <<endl;
    // when checking security level, skHwt(t) should not be set in the build
    std::cout << "Security: " << context.securityLevel() << std::endl;


    std::vector<long> ptxt1 = generateRandomVector(nslots, p);


    // std::cout << "Initial Message: " << ptxt1 << std::endl;

    ea.encrypt(ctxt1, public_key, ptxt1);
    long capacity;
    capacity = ctxt1.bitCapacity();
    cout << "Initial capacity: " << capacity << endl;
    cout << "Initial size: " << ctxt1.logOfPrimeSet()/log(2.0) << endl;

    // ss is the string to hold the ciphertext
    std::stringstream ss_bin;
    // bit serialization
    ctxt1.writeTo(ss_bin);

    // Open a file in write mode
    // std::ofstream file1("bin.txt");
    // if (file1.is_open()) {
    //     // Write the stringstream content to the file
    //     file1 << ss_bin.str();
    //     file1.close();
    // } else {
    //     std::cerr << "Unable to open file.";
    //     return 1;
    // }

    // Open a file in write mode (in binary mode)
    std::ofstream file1("bin.txt", std::ios::binary);
    if (file1.is_open()) {
        // Write the stringstream content to the file
        file1 << ss_bin.rdbuf();  // Use rdbuf() to write binary data correctly
        file1.close();
    } else {
        std::cerr << "Unable to open file.";
        return 1;
    }

    // Open the file in input mode to get the file size
    std::ifstream fileIn1("bin.txt", std::ios::binary | std::ios::ate);
    if (fileIn1.is_open()) {
        // The position of the get pointer, which is at the end due to ios::ate, gives the size of the file
        std::streamsize size = fileIn1.tellg();
        fileIn1.close();
        std::cout << "Binary Ciphertext size: " << size << " bytes." << std::endl;
    } else {
        std::cerr << "Unable to open file for size checking.";
        return 1;
    }

    // Open the file in input mode to get the file size
    // std::ifstream fileIn1("ct_bin.txt",std::ios::binary);
    // if (fileIn1.is_open()) {
    //     // The position of the get pointer, which is at the end due to ios::ate, gives the size of the file
    //     std::streamsize size = fileIn1.tellg();
    //     fileIn1.close();
    //     std::cout << "Binary Ciphertext size: " << size << " bytes." << std::endl;
    // } else {
    //     std::cerr << "Unable to open file for size checking.";
    //     return 1;
    // }

    ea.encrypt(ctxt2, public_key, ptxt1);
    std::stringstream ss_js;
    ss_js << ctxt2;

    // Open a file in write mode
    std::ofstream file2("ct_js.txt");
    if (file2.is_open()) {
        // Write the stringstream content to the file
        file2 << ss_js.str();
        file2.close();
    } else {
        std::cerr << "Unable to open file.";
        return 1;
    }

    // Open the file in input mode to get the file size
    std::ifstream fileIn2("ct_js.txt", std::ios::ate | std::ios::binary);
    if (fileIn2.is_open()) {
        // The position of the get pointer, which is at the end due to ios::ate, gives the size of the file
        std::streamsize size = fileIn2.tellg();
        fileIn2.close();
        std::cout << "Json Ciphertext size: " << size << " bytes." << std::endl;
    } else {
        std::cerr << "Unable to open file for size checking.";
        return 1;
    }


    helib::Ctxt ctxt_res(public_key);
    ctxt1 += ctxt2;
    ctxt_res = ctxt1;

    helib::Ptxt<helib::BGV> result(context);
    secret_key.Decrypt(result, ctxt_res);

    capacity = ctxt_res.bitCapacity();
    cout << "Add capacity: " << capacity << endl;
    cout << "Add size: " << ctxt_res.logOfPrimeSet()/log(2.0) << endl;


    // std::cout << "Decrypt to: " << result.getSlotRepr() << std::endl;

    return 0;
}
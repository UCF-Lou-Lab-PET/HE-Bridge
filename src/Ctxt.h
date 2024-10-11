// add to Ctxt.h
void divideModByP();

// add to Ctxt.h
void multiplyModByP2R()
  {
    long p2r = getContext().getPPowR();
    long p = context.getP();
    ptxtSpace = p2r*p;
    intFactor %= ptxtSpace; // adjust intFactor
  }


// add to Ctxt.cpp
// decrease the mod only
void Ctxt::divideModByP()
{
  // Special case: if *this is empty then do nothing
  if (this->isEmpty())
    return;

  long p = getContext().getP();
  assertEq(ptxtSpace % p, 0l, "p must divide ptxtSpace");
  assertTrue(ptxtSpace > p, "ptxtSpace must be strictly greater than p");

  noiseBound /= p;        // noise is reduced by a p factor
  ptxtSpace /= p;         // and so is the plaintext space
  intFactor %= ptxtSpace; // adjust intFactor
}
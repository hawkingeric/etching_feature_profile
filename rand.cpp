#include "rand.h"
#include "etching.h"
    long iRandTag = 0;
    int RandomNumberCount = 0;
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

/*
Long period (> 2 × 1018) random number generator of L’Ecuyer with Bays-Durham shuffle
and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter
idum between successive deviates in a sequence. RNMX should approximate the largest floating
value that is less than 1.
*/

float ran2(long *idum)
{
int j;
long k;
static long idum2=123456789;
static long iy=0;
static long iv[NTAB];
float temp;

if (*idum <= 0) {                   //Initialize.
if (-(*idum) < 1) *idum=1;          //Be sure to prevent idum = 0.
else *idum = -(*idum);
idum2=(*idum);
for (j=NTAB+7;j>=0;j--) {           //Load the shuffle table (after 8 warm-ups).
k=(*idum)/IQ1;
*idum=IA1*(*idum-k*IQ1)-k*IR1;
if (*idum < 0) *idum += IM1;
if (j < NTAB) iv[j] = *idum;
}
iy=iv[0];
}
k=(*idum)/IQ1;                         //Start here when not initializing.
*idum=IA1*(*idum-k*IQ1)-k*IR1;         //Compute idum=(IA1*idum) % IM1 without
if (*idum < 0) *idum += IM1;           //overflows by Schrage’s method.
k=idum2/IQ2;
idum2=IA2*(idum2-k*IQ2)-k*IR2;         //Compute idum2=(IA2*idum) % IM2 likewise.
if (idum2 < 0) idum2 += IM2;
j=iy/NDIV;                             //Will be in the range 0..NTAB-1.
iy=iv[j]-idum2;                        //Here idum is shuffled, idum and idum2 are
iv[j] = *idum;                         //combined to generate output.
if (iy < 1) iy += IMM1;
if ((temp=AM*iy) > RNMX) return RNMX;  //Because users don’t expect endpoint values.
else return temp;
}

#include <stdlib.h>                     //Change to math.h in K&R C.
#define MBIG 1000000000
//#define MSEED 161803398
static long MSEED = 161803398;
#define MZ 0
#define FAC (1.0/MBIG)
//According to Knuth, any large MBIG, and any smaller (but still large) MSEED can be substituted
//for the above values.

float ran3(long *idum)

//Returns a uniform random deviate between 0.0 and 1.0. Set idum to any negative value to
//initialize or reinitialize the sequence.
{
static int inext,inextp;
static long ma[56];                         //The value 56 (range ma[1..55]) is special and
static int iff=0;                           //should not be modified; see Knuth.
long mj,mk;
int i,ii,k;

RandomNumberCount++;

if (*idum < 0 || iff == 0) {                //Initialization.
    iff=1;
    mj=labs(MSEED-labs(*idum));                 //Initialize ma[55] using the seed idum and the
    mj %= MBIG;                                 //large number MSEED.
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++) {                       //Now initialize the rest of the table,
        ii=(21*i) % 55;                             //in a slightly random order,
        ma[ii]=mk;                                  //with numbers that are not especially random.
        mk=mj-mk;
        if (mk < MZ) mk += MBIG;
        mj=ma[ii];
    }
    for (k=1;k<=4;k++)                          //We randomize them by “warming up the generator.”
        for (i=1;i<=55;i++) {
            ma[i] -= ma[1+(i+30) % 55];
            if (ma[i] < MZ) ma[i] += MBIG;
        }
        inext=0;                                    //Prepare indices for our first generated number.
        inextp=31;                                  //The constant 31 is special; see Knuth.
        *idum=1;
}
                                            //Here is where we start, except on initialization.
if (++inext == 56) inext=1;                 //Increment inext and inextp, wrapping around
if (++inextp == 56) inextp=1;               //56 to 1.
mj=ma[inext]-ma[inextp];                    //Generate a new random number subtractively.
if (mj < MZ) mj += MBIG;                    //Be sure that it is in range.
ma[inext]=mj;                               //Store it,

return mj*FAC;                              //and output the derived uniform deviate.

}

void ran3_ini(int seed){
    MSEED +=seed;
}

#define NITER 4
void psdes(unsigned long *lword, unsigned long *irword)
//“Pseudo-DES” hashing of the 64-bit word (lword,irword). Both 32-bit arguments are returned
//hashed on all bits.
{
unsigned long i,ia,ib,iswap,itmph=0,itmpl=0;
static unsigned long c1[NITER]={0xbaa96887L, 0x1e17d32cL, 0x03bcdc3cL, 0x0f33d1b2L};
static unsigned long c2[NITER]={0x4b0f3b58L, 0xe874f0c3L, 0x6955c5a6L, 0x55a7ca46L};
for (i=0;i<NITER;i++) {
//Perform niter iterations of DES logic,us ing a simpler (non-cryptographic) nonlinear function
//instead of DES’s.
    ia=(iswap=(*irword)) ^ c1[i];                       //The bit-rich constants c1 and (below)
                                                        //c2 guarantee lots of nonlinear mixing.
    itmpl = ia & 0xffff;
    itmph = ia >> 16;
    ib=itmpl*itmpl+ ~(itmph*itmph);
    *irword=(*lword) ^ (((ia = (ib >> 16) |((ib & 0xffff) << 16)) ^ c2[i])+itmpl*itmph);
    *lword=iswap;
}
}

float ran4(long *idum)
//Returns a uniform random deviate in the range 0.0 to 1.0,gener ated by pseudo-DES (DESlike)
//hashing of the 64-bit word (idums,idum),where idums was set by a previous call with
//negative idum. Also increments idum. Routine can be used to generate a random sequence
//by successive calls,le aving idum unaltered between calls; or it can randomly access the nth
//deviate in a sequence by calling with idum = n. Different sequences are initialized by calls with
//differing negative values of idum.
{
    void psdes(unsigned long *lword, unsigned long *irword);
    unsigned long irword,itemp,lword;
    static long idums = 0;
//The hexadecimal constants jflone and jflmsk below are used to produce a floating number
//between 1. and 2. by bitwise masking. They are machine-dependent. See text.

#if defined(vax) || defined(_vax_) || defined(__vax__) || defined(VAX)
    static unsigned long jflone = 0x00004080;
    static unsigned long jflmsk = 0xffff007f;
#else
    static unsigned long jflone = 0x3f800000;
    static unsigned long jflmsk = 0x007fffff;
#endif

if (*idum < 0) {                        //Reset idums and prepare to return the first
    idums = -(*idum);                   //deviate in its sequence.
    *idum=1;
}
irword=(*idum);
lword=idums;
psdes(&lword,&irword);                  //“Pseudo-DES” encode the words.
itemp=jflone | (jflmsk & irword);       //Mask to a floating number between 1 and 2.
++(*idum);
return (*(float *)&itemp)-1.0;          // Subtraction moves range to 0. to 1.
}







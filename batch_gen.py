
# For each d in diffs, 2**31 - d is a prime
#diffs = [1,19,61,69,85,99,105,151,159,171]
diffs = [1,19,61,69,85]
moduli = list(map(lambda d: 2**31 - d, diffs)) + [67]

# Generate primes_list.h

output = ("#ifndef __IDIGM_PRIMES_LIST_H__\n"
         +"#define __IDIGM_PRIMES_LIST_H__\n"
         +"#include \"util.cpp\"\n"
         +"vec_t<uint_t> IDIGM_PRIMES {")
for modulus in moduli:
    output += str(modulus) + ", "

output += "};\n#endif  // __IDIGM_PRIMES_LIST_H__\n"

with open("primes_list.h", "w") as f:
    f.write(output)


# Generate demux_mod.cpp

output = (
    "str_t demuxMod(vec_t<str_t> workFiles, str_t cmd, Options opts=Options())\n"
    + "{\n"
    + "  FileHandler fh(workFiles[0]);\n"
    + "  uint_t modulus;\n"
    + "  assert ( loadParamFromFile<uint_t>(fh, \"modulus\", modulus) );\n"
    + "  if ("
)

for i in range(len(moduli)):
    output += (
        "modulus == " + str(moduli[i]) + ") {\n"
        + "    return demuxCmd<MP"+str(moduli[i])+">(cmd, workFiles, opts);\n"
    )
    if i < len(moduli) - 1:
        output += "  } else if ("
    else:
        output += ("  } else {\n"
            + "    std::cout << \"Error: modulus \" << modulus << \" "
            + "not available\\n\";\n"
            + "    return \"\";\n"
            + "  }\n}\n")

with open("demux_mod.cpp", "w") as f:
    f.write(output)


# Helper functions for generating montgomeryDefinitions.cpp

def normalize(x, modulus):
    r = x % modulus
    if r < 0: r += m
    return r

def firstBezoutFactor(a, b):
    old_s = 1
    old_t = 0
    old_r = a
    s = 0
    t = 1
    r = b
    while r != 0:
        q = old_r // r
        temp = old_r
        old_r = r
        r = temp - q*r
        temp = old_s
        old_s = s
        s = temp - q*s
        temp = old_t
        old_t = t
        t = temp - q*t
    assert old_r == 1, (
        f"firstBezoutFactor: Expected first argument ({a}) to be "
        + f"a unit mod second argument ({b}).")
    return old_s


# Generate montgomeryDefinitions.cpp

output = ""

for modulus in moduli:
    negModInv = normalize(-firstBezoutFactor(modulus, 2**32), 2**32)
    r2mod = normalize(2**64, modulus)
    output += ("\n"
        +"// NOTE: Assumption: the Montgomery modulus is 1UL<<32\n"
        +"struct MP"+str(modulus)+" {\n"
        +"    // original modulus\n"
        +"    static const scalar_t modulus = "+str(modulus)+"U;\n"
        +"    // R*R mod modulus: used for transforming to Montgomery form\n"
        +"    static const scalar_t R2_mod_modulus = "+str(r2mod)+"U;\n"
        +"    // modulus * neg_modulus_inv === -1 (mod R)\n"
        +"    static const scalar_t neg_modulus_inv = "+str(negModInv)+"U;\n"
        +"};\n\n"
        +"#if !(IGNORE_CUDA)\n"
        +"__device__ inline\n"
        +"mont_t montgomery_mul_device_"+str(modulus)+"(mont_t x, mont_t y)\n"
        +"{\n"
        +"    scalar_t r;\n"
        +"    asm(\"{\\n\\t\"\n"
        +"        \"   .reg .u32  Tl;                  \\n\\t\"\n"
        +"        \"   .reg .u32  Th;                  \\n\\t\"\n"
        +"        \"   .reg .u32  m;                   \\n\\t\"\n"
        +"        \"   .reg .u32  tl;                  \\n\\t\"\n"
        +"        \"   .reg .u32  th;                  \\n\\t\"\n"
        +"        \"   .reg .u32  Mf;                  \\n\\t\"\n"
        +"        \"   .reg .u32  Mo;                  \\n\\t\"\n"
        +"        \"   .reg .pred p;                   \\n\\t\"\n"
        +"        \"   mov.u32         Mo, "+str(modulus)+"U; \\n\\t\"\n"
        +"        \"   mov.u32         Mf, "+str(negModInv)+"U; \\n\\t\"\n"
        +"        \"   mul.lo.u32      Tl, %1, %2;     \\n\\t\"\n"
        +"        \"   mul.hi.u32      Th, %1, %2;     \\n\\t\"\n"
        +"        \"   mul.lo.u32      m, Tl, Mf;      \\n\\t\"\n"
        +"        \"   mad.lo.cc.u32   tl, m, Mo, Tl;  \\n\\t\"\n"
        +"        \"   madc.hi.u32     th, m, Mo, Th;  \\n\\t\"\n"
        +"        \"   setp.ge.u32     p, th, Mo;      \\n\\t\"\n"
        +"        \"@p sub.u32         th, th, Mo;     \\n\\t\"\n"
        +"        \"   mov.u32         %0, th;         \\n\\t\"\n"
        +"        \"}\\n\\t\"\n"
        +"        : // output [and input--output] operands below\n"
        +"        // '=' write-register\n"
        +"        // '+' read-and-write register\n"
        +"        // 'r' .u32 reg\n"
        +"        // 'l' .u64 reg\n"
        +"        \"=r\"(r)\n"
        +"        : // input operands below [omit colon if none]\n"
        +"        \"r\"(x),\n"
        +"        \"r\"(y)\n"
        +"        );\n"
        +"    return r;\n"
        +"}\n"
        +"#endif\n\n"
        +"template <>\n"
        +"__host__ __device__ inline \n"
        +"mont_t montgomery_mul<MP"+str(modulus)+">(mont_t x, mont_t y)\n"
        +"{\n"
        +"#if (defined (__CUDA_ARCH__)) && !(IGNORE_CUDA)\n"
        +"    return montgomery_mul_device_"+str(modulus)+"(x,y);\n"
        +"#else\n"
        +"    return montgomery_mul_host<MP"+str(modulus)+">(x,y);\n"
        +"#endif\n"
        +"}\n")

with open("montgomery_definitions.cpp", "w") as f:
    f.write(output)


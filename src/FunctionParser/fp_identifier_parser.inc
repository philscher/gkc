/* NOTE:
  Do not include this file in your project. The fparser.cc file #includes
this file internally and thus you don't need to do anything (other than keep
this file in the same directory as fparser.cc).

  Part of this file is generated code (by using the make_function_name_parser
utility, found in the development version of this library). It's not intended
to be modified by hand.
*/

        unsigned nameLength = 0;
        const unsigned maximumNameLength = 0x80000000U-8;
        /*
        Due to the manner the identifier lengths are returned from
        the readOpcode() function, the maximum supported length for
        identifiers is 0x7FFFFFFF bytes. We minus 8 here to add some
        buffer, because of the multibyteness of UTF-8.
        Function names are limited to 0xFFFF bytes instead, but because
        function names that long just are not defined, the point is moot.
        */
        const unsigned char* const uptr = (const unsigned char*) input;
        typedef signed char schar;
        while(likely(nameLength < maximumNameLength))
        {
            unsigned char byte = uptr[nameLength+0];
            /* Handle the common case of A-Za-z first */
            if(byte >= 0x40)
            {
                if(byte < 0x80) // 0x40..0x7F - most common case
                {
                    // Valid characters in 40..7F: A-Za-z_
                    // Valid bitmask for 40..5F: 01111111111111111111111111100001
                    // Valid bitmask for 60..7F: 01111111111111111111111111100000
                    if(sizeof(unsigned long) == 8)
                    {
                        const unsigned n = sizeof(unsigned long)*8-32;
                        // ^ avoids compiler warning when not 64-bit
                        unsigned long masklow6bits = 1UL << (byte & 0x3F);
                        if(masklow6bits & ~((1UL << 0) | (0x0FUL << (0x1B  ))
                                          | (1UL << n) | (0x1FUL << (0x1B+n))))
                            { ++nameLength; continue; }
                    }
                    else
                    {
                        unsigned masklow5bits = 1 << (byte & 0x1F);
                        if((masklow5bits & ~(1 | (0x1F << 0x1B))) || byte == '_')
                            { ++nameLength; continue; }
                    }
                    break;
                }
                if(byte < 0xF0)
                {
                    if(byte < 0xE0)
                    {
                        if(byte < 0xC2) break; // 0x80..0xC1
                        if(byte == 0xC2 && uptr[nameLength+1]==0xA0) break; // skip nbsp
                        // C2-DF - next common case when >= 0x40
                        // Valid sequence: C2-DF 80-BF
                        if(schar(uptr[nameLength+1]) > schar(0xBF)) break;
                        nameLength += 2;
                        continue;
                    }
                    if(byte == 0xE0) // E0
                    {
                        // Valid sequence: E0 A0-BF 80-BF
                        if((unsigned char)(uptr[nameLength+1] - 0xA0) > (0xBF-0xA0)) break;
                    }
                    else
                    {
                        if(byte == 0xED) break; // ED is invalid
                        // Valid sequence: E1-EC 80-BF 80-BF
                        //            And: EE-EF 80-BF 80-BF
                        if(byte == 0xE2)
                        {
                            // break on various space characters
                            if(uptr[nameLength+1] == 0x80
                            && (schar(uptr[nameLength+2]) <= schar(0x8B)
                            || (uptr[nameLength+2] == 0xAF))) break;
                            if(uptr[nameLength+1] == 0x81
                            && uptr[nameLength+2] == 0x9F) break;
                        } else
                        if(byte == 0xE3 && uptr[nameLength+1] == 0x80
                        && uptr[nameLength+2] == 0x80) break; // this too

                        if(schar(uptr[nameLength+1]) > schar(0xBF)) break;
                    }
                    if(schar(uptr[nameLength+2]) > schar(0xBF)) break;
                    nameLength += 3;
                    continue;
                }
                if(byte == 0xF0) // F0
                {
                    // Valid sequence: F0 90-BF 80-BF 80-BF
                    if((unsigned char)(uptr[nameLength+1] - 0x90) > (0xBF-0x90)) break;
                }
                else
                {
                    if(byte > 0xF4) break; // F5-FF are invalid
                    if(byte == 0xF4) // F4
                    {
                        // Valid sequence: F4 80-8F
                        if(schar(uptr[nameLength+1]) > schar(0x8F)) break;
                    }
                    else
                    {
                        // F1-F3
                        // Valid sequence: F1-F3 80-BF 80-BF 80-BF
                        if(schar(uptr[nameLength+1]) > schar(0xBF)) break;
                    }
                }
                if(schar(uptr[nameLength+2]) > schar(0xBF)) break;
                if(schar(uptr[nameLength+3]) > schar(0xBF)) break;
                nameLength += 4;
                continue;
            }
            if(nameLength > 0)
            {
                if(sizeof(unsigned long) == 8)
                {
                    // Valid bitmask for 00..1F: 00000000000000000000000000000000
                    // Valid bitmask for 20..3F: 00000000000000001111111111000000
                    const unsigned n = sizeof(unsigned long)*8-32;
                    // ^ avoids compiler warning when not 64-bit
                    unsigned long masklow6bits = 1UL << byte;
                    if(masklow6bits & (((1UL << 10)-1UL) << (16+n)))
                        { ++nameLength; continue; }
                }
                else
                {
                    if(byte >= '0' && byte <= '9')
                        { ++nameLength; continue; }
                }
            }
            break;
        }

        /* This function generated with make_function_name_parser.cc */
#define lR lL 5;}
#define lQ l2 lG
#define lP l3 lG
#define lO if('i'l7
#define lN 'n'l7
#define lM l4 3]={
#define lL return
#define lK )==0)l0
#define lJ l5 3 lK
#define lI 0x80000003U:3;
#define lH lA uptr[
#define lG case
#define lF l6 3;}lG
#define lE std::memcmp(uptr+
#define lD 0x80000005U:5;
#define lC :lH 1]){lG
#define lB 'a'lC
#define lA switch(
#define l9 <<16)|
#define l8 lD lR lG
#define l7 ==uptr[
#define l6 default:lL
#define l5 lE 1,tmp,
#define l4 static const char tmp[
#define l3 lI lL 3;
#define l2 0x80000004U:4;lL 4;
#define l1 .enabled()?(
#define l0 lL Functions[
lA
nameLength){lG
2:lO
0]&&'f'l7
1])l0
cIf]l1
cIf
l9
0x80000002U:2;lL
2;lG
3:lH
0]){lG
lB'b':if('s'l7
2])l0
cAbs]l1
cAbs
l9
lP'r':if('g'l7
2])l0
cArg]l1
cArg
l9
l3
lF'c'lC'o':lH
2]){lG's':l0
cCos]l1
cCos
l9
lI
lG't':l0
cCot]l1
cCot
l9
lI
lF's':if('c'l7
2])l0
cCsc]l1
cCsc
l9
l3
lF'e':if('x'l7
1]&&'p'l7
2])l0
cExp]l1
cExp
l9
lP'i':if(lN
1]&&'t'l7
2])l0
cInt]l1
cInt
l9
lP'l':if('o'l7
1]&&'g'l7
2])l0
cLog]l1
cLog
l9
lP'm'lC'a':if('x'l7
2])l0
cMax]l1
cMax
l9
lP'i':if(lN
2])l0
cMin]l1
cMin
l9
l3
lF'p':if('o'l7
1]&&'w'l7
2])l0
cPow]l1
cPow
l9
lP's'lC'e':if('c'l7
2])l0
cSec]l1
cSec
l9
lP'i':if(lN
2])l0
cSin]l1
cSin
l9
l3
lF't':if('a'l7
1]&&lN
2])l0
cTan]l1
cTan
l9
l3
lF
4:lH
0]){lG
lB'c':if('o'l7
2]&&'s'l7
3])l0
cAcos]l1
cAcos
l9
lQ's':lO
2]&&lN
3])l0
cAsin]l1
cAsin
l9
lQ't':if('a'l7
2]&&lN
3])l0
cAtan]l1
cAtan
l9
l2
l6
4;}
lG'c'lC'b':if('r'l7
2]&&'t'l7
3])l0
cCbrt]l1
cCbrt
l9
lQ'e':lO
2]&&'l'l7
3])l0
cCeil]l1
cCeil
l9
lQ'o':lH
2]){lG'n':if('j'l7
3])l0
cConj]l1
cConj
l9
lQ's':if('h'l7
3])l0
cCosh]l1
cCosh
l9
l2
l6
4;}
l6
4;}
lG'e'lC'v':if('a'l7
2]&&'l'l7
3])l0
cEval]l1
cEval
l9
lQ'x':if('p'l7
2]&&'2'l7
3])l0
cExp2]l1
cExp2
l9
l2
l6
4;}
lG'i':{lM'm','a','g'}
;if(lJ
cImag]l1
cImag
l9
l2}
lG'l':{lM'o','g','2'}
;if(lJ
cLog2]l1
cLog2
l9
l2}
lG'r':{lM'e','a','l'}
;if(lJ
cReal]l1
cReal
l9
l2}
lG's'lC'i':if(lN
2]&&'h'l7
3])l0
cSinh]l1
cSinh
l9
lQ'q':if('r'l7
2]&&'t'l7
3])l0
cSqrt]l1
cSqrt
l9
l2
l6
4;}
lG't':{lM'a','n','h'}
;if(lJ
cTanh]l1
cTanh
l9
l2}
l6
4;}
lG
5:lH
0]){lG
lB'c':{lM'o','s','h'}
;if(lE
2,tmp,3
lK
cAcosh]l1
cAcosh
l9
l8's':{lM'i','n','h'}
;if(lE
2,tmp,3
lK
cAsinh]l1
cAsinh
l9
l8't':if('a'l7
2]){if(lN
3]){lH
4]){lG'2':l0
cAtan2]l1
cAtan2
l9
lD
lG'h':l0
cAtanh]l1
cAtanh
l9
lD
l6
5;}
}
lR
lL
5;l6
5;}
lG'f':{l4
4]={'l','o','o','r'}
;if(l5
4
lK
cFloor]l1
cFloor
l9
l8'h':{l4
4]={'y','p','o','t'}
;if(l5
4
lK
cHypot]l1
cHypot
l9
l8'l':{l4
4]={'o','g','1','0'}
;if(l5
4
lK
cLog10]l1
cLog10
l9
l8'p':{l4
4]={'o','l','a','r'}
;if(l5
4
lK
cPolar]l1
cPolar
l9
l8't':{l4
4]={'r','u','n','c'}
;if(l5
4
lK
cTrunc]l1
cTrunc
l9
lD
lR
l6
5;}
default:break;}
lL
nameLength;
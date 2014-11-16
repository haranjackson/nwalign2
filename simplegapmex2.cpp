//=============================================================================
// Name        : simplegapmex2.cpp
// Version     : 1.0
// Copyright   : 2014 Haran Jackson
//
// Description : A divide and conquer version of MATLAB's simplegapmex function
//                 (found in nwalign.m), based on Hirschberg's Algorithm:
//
//               en.wikipedia.org/wiki/Hirschberg%27s_algorithm
//
//                 simplegapmex2 is used in nwalign2.m.
//
//                 nwalign2 is used in the same way as nwalign.
//=============================================================================


#include <string>
#include <algorithm>
#include <vector>

#include "mex.h"

using namespace std;


struct out          // To allow the return of multiple variables from functions
{
    vector<int> X;
    vector<int> Y;
    string mid;
};


long max(long a, long b, long c)
{
    if (a >= c)
    {
        if (a >= b)  { return a; }
        else         { return b; }
    }
    else
    {
        if (b >= c)  { return b; }
        else         { return c; }
    }
}


vector<long> ScoringMatrix;

int M0, N0;               // Dimensions of scoring matrix

int gapchar;              // Numerical value of a gap character in the alphabet use (amino / dna)

long gap;                 // Gap penalty

bool glocal;              // Whether glocal alignment is performed

bool isAminoAcid;         // True if input sequences are amino acids, false if they are nucleotides


vector<long> NWScore(vector<int> & x, vector<int> & y, bool start)
{
    long xlen = x.size();
    long ylen = y.size();

    vector<long> prev(ylen + 1);
    vector<long> next(ylen + 1);

    // In glocal alignment, no penalties for gaps at the beginning / end of alignment

    if (start)
    {
        for (long j = 0; j < ylen + 1; j++)  { prev[j] = 0; }

        for (long i = 1; i < xlen + 1; i++)
        {
            next[0] = 0;

            for (long j = 1; j < ylen + 1; j++)
            {
                next[j] = max (prev[j-1] + ScoringMatrix[x[i-1] * N0 + y[j-1]], \
                               prev[j]   + gap, \
                               next[j-1] + gap);
            }
            prev = next;
        }
    }
    else
    {
        for (long j = 0; j < ylen + 1; j++)  { prev[j] = j * gap; }

        for (long i = 1; i < xlen + 1; i++)
        {
            next[0] = i * gap;

            for (long j = 1; j < ylen + 1; j++)
            {
                next[j] = max (prev[j-1] + ScoringMatrix[x[i-1] * N0 + y[j-1]], \
                               prev[j]   + gap, \
                               next[j-1] + gap);
            }
            prev = next;
        }
    }

    return next;
}


long PartitionY(vector<long> & scoreL, vector<long> & scoreR, long n)
{
    long max_index = 0;

    long max_sum = scoreL[0] + scoreR[n];

    for (long i = 1; i <= n; i++)
    {
        if (scoreL[i] + scoreR[n - i] > max_sum)
        {
            max_sum = scoreL[i] + scoreR[n - i];
            max_index = i;
        }
    }

    return max_index;
}


out NeedlemanWunsch(vector<int> & x, vector<int> & y, bool start, bool end)
{
    out NWOut;

    long m = x.size() + 1;
    long n = y.size() + 1;

    vector<long> M (m * n);
    M[0] = 0;

    vector<char> Path (m * n);
    Path[0] = 'f';


    if (start)
    {
        for (long i = 1; i < m; i++)
        {
            M[i * n] = 0;
            Path[i * n] = 'u';
        }

        for (long j = 1; j < n; j++)
        {
            M[j] = 0;
            Path[j] = 'l';
        }
    }
    else
    {
        for (long i = 1; i < m; i++)
        {
            M[i * n] = i * gap;
            Path[i * n] = 'u';
        }

        for (long j = 1; j < n; j++)
        {
            M[j] = j * gap;
            Path[j] = 'l';
        }
    }


    long best;

    if (end)
    {
        for (long i = 1; i < m - 1; i++)
        {
            /// i < m - 1,  j < n - 1 ///

            for (long j = 1; j < n - 1; j++) {

                long u = M[(i - 1) * n + j]       + gap;
                long l = M[i * n + (j - 1)]       + gap;
                long d = M[(i - 1) * n + (j - 1)] + ScoringMatrix [x[i-1] * N0 + y[j-1]];

                if (u > l)
                {
                    best = u;
                    Path[i * n + j] = 'u';
                }
                else
                {
                    best = l;
                    Path[i * n + j] = 'l';
                }
                if (d >= best)
                {
                    M[i * n + j] = d;
                    Path[i * n + j] = 'd';
                }
                else
                {
                    M[i * n + j] = best;
                }
            }

            /// i < m - 1,  j = n - 1 ///

            long u = M[(i - 1) * n + (n - 1)];
            long l = M[i * n + (n - 2)]        + gap;
            long d = M[(i - 1) * n + (n - 2)]  + ScoringMatrix [x[i-1] * N0 + y[n-2]];

            if (u > l)
            {
                best = u;
                Path[i * n + (n - 1)] = 'u';
            }
            else
            {
                best = l;
                Path[i * n + (n - 1)] = 'l';
            }
            if (d >= best)
            {
                M[i * n + (n - 1)] = d;
                Path[i * n + (n - 1)] = 'd';
            }
            else
            {
                M[i * n + (n - 1)] = best;
            }
        }

        //// i = m - 1,  j < n - 1 ////

        for (long j = 1; j < n - 1; j++)
        {
            long u = M[(m - 2) * n + j]        + gap;
            long l = M[(m - 1) * n + (j - 1)];
            long d = M[(m - 2) * n + (j - 1)]  + ScoringMatrix[x[m-2] * N0 + y[j-1]];

            if (u > l)
            {
                best = u;
                Path[(m - 1) * n + j] = 'u';
            }
            else
            {
                best = l;
                Path[(m - 1) * n + j] = 'l';
            }
            if (d >= best)
            {
                M[(m - 1) * n + j] = d;
                Path[(m - 1) * n + j] = 'd';
            }
            else
            {
                M[(m - 1) * n + j] = best;
            }
        }

        /// i = m - 1,  j = n - 1 ///

        long u = M[(m - 2) * n + (n - 1)];
        long l = M[(m - 1) * n + (n - 2)];
        long d = M[(m - 2) * n + (n - 2)] + ScoringMatrix [x[m-2] * N0 + y[n-2]];

        if (u > l)
        {
            best = u;
            Path[(m - 1) * n + (n - 1)] = 'u';
        }
        else
        {
            best = l;
            Path[(m - 1) * n + (n - 1)] = 'l';
        }
        if (d >= best)
        {
            M[(m - 1) * n + (n - 1)] = d;
            Path[(m - 1) * n + (n - 1)] = 'd';
        }
        else
        {
            M[(m - 1) * n + (n - 1)] = best;
        }
    }

    else
    {
        for (long i = 1; i < m; i++)
        {
            for (long j = 1; j < n; j++)
            {
                long u = M[(i - 1) * n + j]       + gap;
                long l = M[i * n + (j - 1)]       + gap;
                long d = M[(i - 1) * n + (j - 1)] + ScoringMatrix [x[i-1] * N0 + y[j-1]];

                if (u > l)
                {
                    best = u;
                    Path[i * n + j] = 'u';
                }
                else
                {
                    best = l;
                    Path[i * n + j] = 'l';
                }
                if (d >= best)
                {
                    M[i * n + j] = d;
                    Path[i * n + j] = 'd';
                }
                else
                {
                    M[i * n + j] = best;
                }
            }
        }
    }


    long i = m - 1;
    long j = n - 1;

    while (Path[i * n + j] != 'f')
    {
        if (Path[i * n + j] == 'd')
        {
            NWOut.X.push_back(x[i - 1]);
            NWOut.Y.push_back(y[j - 1]);

            if (x[i - 1] == y[j - 1])  { NWOut.mid.insert(0, "|"); }
            else                       { NWOut.mid.insert(0, " "); }

            i--;
            j--;
        }
        else if (Path[i * n + j] == 'u')
        {
            NWOut.X.push_back(x[i - 1]);
            NWOut.Y.push_back(gapchar);
            NWOut.mid.insert(0, " ");

            i--;
        }
        else
        {
            NWOut.X.push_back(gapchar);
            NWOut.Y.push_back(y[j - 1]);
            NWOut.mid.insert(0, " ");

            j--;
        }
    }

    reverse(NWOut.X.begin(), NWOut.X.end());
    reverse(NWOut.Y.begin(), NWOut.Y.end());

    return (NWOut);
}


out Hirschberg (vector<int> & x, vector<int> & y, bool start, bool end)
{
    out HirschbergOut;

    long xlen = x.size();
    long ylen = y.size();

    if (xlen == 0 || ylen == 0)
    {
        if (xlen == 0)
        {
            vector<int> xvec(ylen, gapchar);
            HirschbergOut.X = xvec;
            HirschbergOut.Y = y;
            HirschbergOut.mid = string(ylen, ' ');
        }
        else
        {
            vector<int> yvec(xlen, gapchar);
            HirschbergOut.X = x;
            HirschbergOut.Y = yvec;
            HirschbergOut.mid = string(xlen, ' ');
        }
    }
    else if (xlen == 1 || ylen == 1)
    {
        HirschbergOut = NeedlemanWunsch(x, y, start, end);
    }
    else
    {
        long xmid = xlen / 2;
        
        vector<int>::const_iterator first1 = x.begin();
        vector<int>::const_iterator last1  = x.begin() + xmid;
        vector<int> xsub1(first1, last1);

        vector<long> scoreL = NWScore(xsub1, y, start);
        
        vector<int>::const_iterator first2 = x.begin() + xmid;
        vector<int>::const_iterator last2  = x.end();
        vector<int> xsub2(first2, last2);

        vector<int> xrev = xsub2;
        vector<int> yrev = y;

        reverse(xrev.begin(), xrev.end());            // Reverses xrev
        reverse(yrev.begin(), yrev.end());            // Reverses yrev

        vector<long> scoreR = NWScore(xrev, yrev, end);

        long ymid = PartitionY(scoreL, scoreR, ylen);
        
        vector<int>::const_iterator first3 = y.begin();
        vector<int>::const_iterator last3  = y.begin() + ymid;
        vector<int> ysub1(first3, last3);

        vector<int>::const_iterator first4 = y.begin() + ymid;
        vector<int>::const_iterator last4  = y.end();
        vector<int> ysub2(first4, last4);

        out LOut, ROut;

        if (ymid == ylen)  { LOut = Hirschberg(xsub1, ysub1, start, end);   }
        else               { LOut = Hirschberg(xsub1, ysub1, start, false); }

        if (ymid == 0)     { ROut = Hirschberg(xsub2, ysub2, start, end);   }
        else               { ROut = Hirschberg(xsub2, ysub2, false, end);   }

        vector<int> XL = LOut.X;
        vector<int> XR = ROut.X;
        vector<int> YL = LOut.Y;
        vector<int> YR = ROut.Y;
        
        XL.insert(XL.end(), XR.begin(), XR.end());
        YL.insert(YL.end(), YR.begin(), YR.end());

        HirschbergOut.X   = XL;
        HirschbergOut.Y   = YL;
        HirschbergOut.mid = LOut.mid + ROut.mid;
    }

    return(HirschbergOut);
}


int CalculateScore(out Output, bool glocal)
{
    long n = Output.X.size();

    int score = 0;

    long startchar = 0;
    long endchar = n;

    if (glocal)
    {
        if (Output.X[0] == gapchar)
        {
            long i = 0;

            while (Output.X[i] == gapchar)
            {
                startchar = i + 1;
                i++;
            }
        }

        if (Output.Y[0] == gapchar)
        {
            long i = 0;

            while (Output.Y[i] == gapchar)
            {
                startchar = i + 1;
                i++;
            }
        }

        if (Output.X[n - 1] == gapchar)
        {
            long i = n - 1;

            while (Output.X[i] == gapchar)
            {
                endchar = i;
                i--;
            }
        }

        if (Output.Y[n - 1] == gapchar)
        {
            long i = n - 1;

            while (Output.Y[i] == gapchar)
            {
                endchar = i;
                i--;
            }
        }
    }


    for (long i = startchar; i < endchar; i++)
    {
        if (Output.X[i] == gapchar || Output.Y[i] == gapchar)
        {
            score += gap;
        }
        else
        {
            score += ScoringMatrix [Output.X[i] * N0 + Output.Y[i]];
        }
    }


    return score;
}


vector<int> nt2int (string & str, long n) {
    
    vector<int> out (n);

    for (long i = 0; i < n; i++)
    {
        switch (str[i])
        {
            case(65): { out[i] = 0;  break; }   // A
            case(67): { out[i] = 1;  break; }   // C
            case(71): { out[i] = 2;  break; }   // G
            case(84): { out[i] = 3;  break; }   // T

            case(78): { out[i] = 14; break; }   // N
            case(45): { out[i] = 15; break; }   // -

            case(85): { out[i] = 3;  break; }   // U
            case(82): { out[i] = 4;  break; }   // R
            case(89): { out[i] = 5;  break; }   // Y

            case(75): { out[i] = 6;  break; }   // K
            case(77): { out[i] = 7;  break; }   // M
            case(83): { out[i] = 8;  break; }   // S
            case(87): { out[i] = 9;  break; }   // W
            case(66): { out[i] = 10; break; }   // B
            case(68): { out[i] = 11; break; }   // D
            case(72): { out[i] = 12; break; }   // H
            case(86): { out[i] = 13; break; }   // V
        }
    }
    return out;
}


vector<int> aa2int (string & str, long n)
{
    vector<int> out (n);

    for (long i = 0; i < n; i++)
    {
        switch (str[i])
        {
            case(65): { out[i] = 0;  break; }   // A
            case(82): { out[i] = 1;  break; }   // R
            case(78): { out[i] = 2;  break; }   // N
            case(68): { out[i] = 3;  break; }   // D
            case(67): { out[i] = 4;  break; }   // C
            case(81): { out[i] = 5;  break; }   // Q
            case(69): { out[i] = 6;  break; }   // E
            case(71): { out[i] = 7;  break; }   // G
            case(72): { out[i] = 8;  break; }   // H
            case(73): { out[i] = 9;  break; }   // I
            case(76): { out[i] = 10; break; }   // L
            case(75): { out[i] = 11; break; }   // K
            case(77): { out[i] = 12; break; }   // M
            case(70): { out[i] = 13; break; }   // F
            case(80): { out[i] = 14; break; }   // P
            case(83): { out[i] = 15; break; }   // S
            case(84): { out[i] = 16; break; }   // T
            case(87): { out[i] = 17; break; }   // W
            case(89): { out[i] = 18; break; }   // Y
            case(86): { out[i] = 19; break; }   // V
            case(66): { out[i] = 20; break; }   // B
            case(90): { out[i] = 21; break; }   // Z
            case(88): { out[i] = 22; break; }   // X
            case(42): { out[i] = 23; break; }   // *
            case(45): { out[i] = 24; break; }   // -
        }
    }
    return out;
}


string int2nt (vector<int> & intseq) {
    
    string strseq;
    
    for (long i = 0; i < intseq.size(); i++) {
        
        switch (intseq[i])
        {
            case (0):  { strseq.append("A"); break; }
            case (1):  { strseq.append("C"); break; }
            case (2):  { strseq.append("G"); break; }
            case (3):  { strseq.append("T"); break; }

            case (14): { strseq.append("N"); break; }
            case (15): { strseq.append("-"); break; }

            case (4):  { strseq.append("R"); break; }
            case (5):  { strseq.append("Y"); break; }
            case (6):  { strseq.append("K"); break; }
            case (7):  { strseq.append("M"); break; }
            case (8):  { strseq.append("S"); break; }
            case (9):  { strseq.append("W"); break; }
            case (10): { strseq.append("B"); break; }
            case (11): { strseq.append("D"); break; }
            case (12): { strseq.append("H"); break; }
            case (13): { strseq.append("V"); break; }
        }
    }
    return strseq;
}


string int2aa(vector<int> & intseq) {

    string strseq;

    for (long i = 0; i < intseq.size(); i++) {

        switch (intseq[i])
        {
            case (24): { strseq.append("-"); break; }
            case (0):  { strseq.append("A"); break; }
            case (1):  { strseq.append("R"); break; }
            case (2):  { strseq.append("N"); break; }
            case (3):  { strseq.append("D"); break; }
            case (4):  { strseq.append("C"); break; }
            case (5):  { strseq.append("Q"); break; }
            case (6):  { strseq.append("E"); break; }
            case (7):  { strseq.append("G"); break; }
            case (8):  { strseq.append("H"); break; }
            case (9):  { strseq.append("I"); break; }
            case (10): { strseq.append("L"); break; }
            case (11): { strseq.append("K"); break; }
            case (12): { strseq.append("M"); break; }
            case (13): { strseq.append("F"); break; }
            case (14): { strseq.append("P"); break; }
            case (15): { strseq.append("S"); break; }
            case (16): { strseq.append("T"); break; }
            case (17): { strseq.append("W"); break; }
            case (18): { strseq.append("Y"); break; }
            case (19): { strseq.append("V"); break; }
            case (20): { strseq.append("B"); break; }
            case (21): { strseq.append("Z"); break; }
            case (22): { strseq.append("X"); break; }
            case (23): { strseq.append("*"); break; }
        }
    }
    return strseq;
}


void StrOut(char *input_buf, size_t buflen, char *output_buf)
{
    mwSize i;

    if (buflen == 0) return;

    for (i = 0; i < buflen-1; i++)  { *(output_buf + i) = *(input_buf + i); }
}
        

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    string seq1 = mxArrayToString(prhs[0]);
    string seq2 = mxArrayToString(prhs[1]);
    
    long n1 = seq1.length();
    long n2 = seq2.length();

    if (mxGetScalar(prhs[5]) == 1)  { isAminoAcid = true;  }
    else                            { isAminoAcid = false; }
    

    vector<int> intseq1, intseq2;

    if (isAminoAcid)
    {
        gapchar = 24;
        intseq1 = aa2int (seq1, n1);
        intseq2 = aa2int (seq2, n2);
    }
    else
    {
        gapchar = 15;
        intseq1 = nt2int (seq1, n1);
        intseq2 = nt2int (seq2, n2);
    }
    
    
    gap = mxGetScalar(prhs[2]);
    
    
    M0 = mxGetM(prhs[3]);
    N0 = mxGetN(prhs[3]);    
    
    ScoringMatrix.resize(M0 * N0);
    
    for (int i = 0; i < M0; i++)
    {
        for (int j = 0; j < N0; j++)
        {
            ScoringMatrix[i * N0 + j] = (mxGetPr(prhs[3]))[i * N0 + j];
        }
    }
    

    if (mxGetScalar(prhs[4]) == 3)  { glocal = true;  }
    else                            { glocal = false; }


    out Output;
    string Xseq, Yseq;


    if (n1 >= n2)
    {
        Output = Hirschberg(intseq1, intseq2, glocal, glocal);

        if (isAminoAcid)
        {
            Xseq = int2aa(Output.X);
            Yseq = int2aa(Output.Y);
        }
        else
        {
            Xseq = int2nt(Output.X);
            Yseq = int2nt(Output.Y);
        }
    }
    else
    {
        Output = Hirschberg(intseq2, intseq1, glocal, glocal);

        if (isAminoAcid)
        {
            Xseq = int2aa(Output.Y);
            Yseq = int2aa(Output.X);
        }
        else
        {
            Xseq = int2nt(Output.Y);
            Yseq = int2nt(Output.X);
        }
    }


    plhs[0] = mxCreateDoubleScalar(CalculateScore(Output, glocal));
    plhs[1] = mxCreateString(Xseq.c_str());
    plhs[2] = mxCreateString(Output.mid.c_str());
    plhs[3] = mxCreateString(Yseq.c_str());
}

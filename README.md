nwalign2
========

Let sequences seq1 and seq2 have lengths n1 and n2, respectively, with  n1<n2.  To align these sequences, MATLAB's nwalign function requires  O(n1*n2)  space and  O(n1*n2)  time. nwalign2 requires  O(n1)  space and  O(n1*n2)  time.

nwalign2 takes the same input as nwalign, and produces the same output. Note, however, that nwalign2 cannot produce the full score and path matrices for the alignment, as this would require  O(n1*n2)  space.


## Improvements

- Incorporate 'ExtendGap' option
- When the 'Glocal' option is chosen, nwalign2 sometimes produces a slightly different result to nwalign. This will be fixed.
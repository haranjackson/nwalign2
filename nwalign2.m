function [score, alignment, startat] = nwalign2(seq1,seq2,varargin)
%
%   NWALIGN2 performs Needleman-Wunsch global alignment of two sequences
%   using only linear space.
%
%   NOTE: Unlike NWALIGN, it is not possible to
%   use the 'showscore' option, or to return the score and pointer
%   matrices (as this would require quadratic space).


%   NWALIGN2(SEQ1, SEQ2) returns the score (in bits) for the optimal
%   alignment. Note: The scale factor used to calculate the score is
%   provided by the scoring matrix info (see below). If this is not
%   defined, then NWALIGN2 returns the raw score.
%
%   [SCORE, ALIGNMENT] = NWALIGN2(SEQ1, SEQ2) returns a string showing an
%   optimal global alignment of amino acid (or nucleotide) sequences SEQ1
%   and SEQ2.
%
%   [SCORE, ALIGNMENT, STARTAT] = NWALIGN2(SEQ1, SEQ2)  returns a 2x1
%   vector with the starting point indices indicating the starting point of
%   the alignment in the two sequences. Note: this output is for
%   consistency with SWALIGN and will always be [1;1] because this is a
%   global alignment.


%   NWALIGN2(..., 'ALPHABET', A) specifies whether the sequences are
%   amino acids ('AA') or nucleotides ('NT'). The default is AA.
%
%   NWALIGN2(..., 'SCORINGMATRIX', matrix) defines the scoring matrix to be
%   used for the alignment. The default is BLOSUM50 for AA or NUC44 for NT.
%
%   NWALIGN2(..., 'SCALE' ,scale) indicates the scale factor of the scoring
%   matrix to return the score using arbitrary units. If the scoring matrix
%   Info also provides a scale factor, then both are used.
%
%   NWALIGN2(..., 'GAPOPEN', penalty) defines the penalty for opening a gap
%   in the alignment. The default gap open penalty is 8.
%
%   NWALIGN2(..., 'EXTENDGAP', penalty) defines the penalty for extending a
%   gap in the alignment. If EXTENDGAP is not specified, then extensions to
%   gaps are scored with the same value as GAPOPEN.
%
%   NWALIGN2(..., 'GLOCAL', true) performs a semi-global alignment. In a
%   semi-global alignment (also known as GLOCAL) gap penalties at the end
%   of the sequences are null.


%   Examples:
%
%       % Return the score in bits and the global alignment using the
%       % default scoring matrix (BLOSUM50).
%       [score, align] = nwalign2('VSPAGMASGYD', 'IPGKASYD')
%
%       % Use user-specified scoring matrix and "gap open" penalty.
%       [score, align] = nwalign2('IGRHRYHIGG', 'SRYIGRG',...
%                               'scoringmatrix', @pam250, 'gapopen',5)
%
%       % Return the score in nat units (nats).
%       [score, align] = nwalign2('HEAGAWGHEE', 'PAWHEAE', 'scale', log(2))


glocal = false;
gapopen = -8;
gapextend = -8;
setGapExtend = false;
isAminoAcid = true;
scale = 1;

% If the input is a structure then extract the Sequence data.
if isstruct(seq1)
    seq1 = bioinfoprivate.seqfromstruct(seq1);
end
if isstruct(seq2)
    seq2 = bioinfoprivate.seqfromstruct(seq2);
end
if nargin > 2
    if rem(nargin,2) == 1
        error(message('bioinfo:nwalign:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'scoringmatrix','gapopen','extendgap','alphabet','scale','glocal'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:nwalign:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:nwalign:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % scoring matrix
                    if isnumeric(pval)
                        ScoringMatrix = pval;
                    else
                        if ischar(pval)
                            pval = lower(pval);
                        end
                        try
                            [ScoringMatrix,ScoringMatrixInfo] = feval(pval);
                        catch allExceptions
                            error(message('bioinfo:nwalign:InvalidScoringMatrix'));
                        end
                    end
                case 2 %gap open penalty
                    gapopen = -pval;
                case 3 %gap extend penalty
                    gapextend = -pval;
                    setGapExtend = true;
                case 4 %if sequence is nucleotide
                    isAminoAcid = bioinfoprivate.optAlphabet(pval,okargs{k}, mfilename);
                case 5 % scale
                    scale=pval;
                case 6 % glocal
                    glocal = bioinfoprivate.opttf(pval, okargs{k}, mfilename);
            end
        end
    end
end


% setting the default scoring matrix
if ~exist('ScoringMatrix','var')
    if isAminoAcid
        [ScoringMatrix,ScoringMatrixInfo] = blosum50;
    else
        [ScoringMatrix,ScoringMatrixInfo] = nuc44;
    end
end


% getting the scale from ScoringMatrixInfo, if it exists
if exist('ScoringMatrixInfo','var') && isfield(ScoringMatrixInfo,'Scale')
    scale=scale*ScoringMatrixInfo.Scale;
end

% handle properly "?" characters typically found in pdb files
if isAminoAcid
    if ischar(seq1)
        seq1 = strrep(seq1,'?','X');
    else
        seq1(seq1 == 26) = 23;
    end
    if ischar(seq2)
        seq2 = strrep(seq2,'?','X');
    else
        seq2(seq2 == 26) = 23;
    end
end

% check input sequences
if isAminoAcid && ~(bioinfoprivate.isaa(seq1) && bioinfoprivate.isaa(seq2))
    error(message('bioinfo:nwalign:InvalidAminoAcidSequences'));
elseif ~isAminoAcid && ~(bioinfoprivate.isnt(seq1) && bioinfoprivate.isnt(seq2))
    error(message('bioinfo:nwalign:InvalidNucleotideSequences'));
end

% use numerical arrays for easy indexing
if ischar(seq1)
    seq1=upper(seq1); %the output alignment will be all uppercase
    if isAminoAcid
        intseq1 = aa2int(seq1);
    else
        intseq1 = nt2int(seq1);
    end
else
    intseq1 = uint8(seq1);
    if isAminoAcid
        seq1 = int2aa(intseq1);
    else
        seq1 = int2nt(intseq1);
    end
end
if ischar(seq2)
    seq2 = upper(seq2); %the output alignment will be all uppercase
    if isAminoAcid
        intseq2 = aa2int(seq2);
    else
        intseq2 = nt2int(seq2);
    end
else
    intseq2 = uint8(seq2);
    if isAminoAcid
        seq2 = int2aa(intseq2);
    else
        seq2 = int2nt(intseq2);
    end
end


m = length(seq1);
n = length(seq2);
if ~n||~m
    error(message('bioinfo:nwalign:InvalidLengthSequences'));
end


scoringMatrixSize = size(ScoringMatrix,1);

if max([intseq1, intseq2]) > scoringMatrixSize
    % if the matrix contains the 'Any' we map to that
    if isAminoAcid
        anyVal = aa2int('X');
    else
        anyVal = nt2int('N');
    end
    if scoringMatrixSize >= anyVal
        if isAminoAcid
            seq1(intseq1>scoringMatrixSize) = 'X';
            seq2(intseq2>scoringMatrixSize) = 'X';
        else
            seq1(intseq1>scoringMatrixSize) = 'N';
            seq2(intseq2>scoringMatrixSize) = 'N';
        end
    else
        error(message('bioinfo:nwalign:InvalidSymbolsInInputSequences'));
    end
end


if glocal
   algorithm = 3;
else
   algorithm = 1;
end


% flip order of input sequences for consistency with older versions
if setGapExtend
    [score, AlignSeq1, Alignmid, AlignSeq2] = affinegapmex2(seq1, seq2, gapopen, gapextend, ScoringMatrix, algorithm, isAminoAcid);
else
    [score, AlignSeq1, Alignmid, AlignSeq2] = simplegapmex2(seq1, seq2, gapopen, ScoringMatrix, algorithm, isAminoAcid);
end

alignment(1,:) = AlignSeq1;
alignment(2,:) = Alignmid;
alignment(3,:) = AlignSeq2;

% re-scaling the output score
score = scale * score;

startat = [1;1];

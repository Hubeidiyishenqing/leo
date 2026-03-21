function decBits = decodeLDPC_helper(llr, ldpcDecoder, numCW, noCodedbits)
%--------------------------------------------------------------------------
%   Decodes multiple LDPC codewords from a concatenated LLR vector.
%
%   llr          numCW*noCodedbits x 1 log-likelihood ratio vector
%   ldpcDecoder  comm.LDPCDecoder system object
%   numCW        Number of codewords
%   noCodedbits  Bits per codeword (e.g. 64800)
%
%   Returns decoded info bits (concatenated).
%--------------------------------------------------------------------------

decBits = [];
buf = llr(:);
for q = 1:numCW
    seg = buf(1:noCodedbits);
    buf(1:noCodedbits) = [];
    decoded = double(ldpcDecoder(seg));
    decBits = [decBits; decoded];
end

end

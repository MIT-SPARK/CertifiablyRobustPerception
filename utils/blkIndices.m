function matrixIndices = blkIndices(blockIndex, blockSize)
% The function returns the indices corresponding to the blockIndex-th block
% for a matrix or vector with blocks of size blockIndex X blockIndex
% 
% Luca Carlone
% Georgia Institute of Technology
% 1/5/2013

numBlocks = length(blockIndex);
matrixIndices = zeros(numBlocks*blockSize,1);

for i=1:numBlocks
  matrixIndices(blockSize*i-(blockSize-1) : blockSize*i) = [blockSize*blockIndex(i)-(blockSize-1) : blockSize*blockIndex(i)];
end

end
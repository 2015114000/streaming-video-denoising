function  [par,parBSM]=ParSet_New(nSig)
% parameters setting for ITS_DeNoising

parBSM.lambda     =   10;
parBSM.rho        =   1.2;
parBSM.mu         =   250;
meanD             =   0.2;
DimUpRank         =   [1,1,1,1];
par.nSig          =   nSig;                                 % Variance of the noise image
par.delta         =   0.1;                                  % Parameter between each iter


if nSig <= 0.1


    par.SigLam        =   0.012*meanD*prod(DimUpRank);  % Noise estimete parameter
    parBSM.maxIter    =   25;                           % max iteration number for ITSReg tensor recovery
    par.deNoisingIter =   2;                            % total iteration numbers
elseif nSig <=0.15
                         
                      
    par.SigLam        =   0.013*meanD*prod(DimUpRank);
    parBSM.maxIter    =   25;
    par.deNoisingIter =   2;
elseif nSig <=0.2
                      
    par.SigLam        =   0.014*meanD*prod(DimUpRank);
    parBSM.maxIter    =   35;
    par.deNoisingIter =   3;
elseif nSig <=0.25
                        
    par.SigLam        =   0.015*meanD*prod(DimUpRank);
    parBSM.maxIter    =   35;
    par.deNoisingIter =   3;
elseif nSig <=0.3
                     
    par.SigLam        =   0.016*meanD*prod(DimUpRank);
    parBSM.maxIter    =   30;
    par.deNoisingIter =   3;
    
elseif nSig >0.3
                     
    par.SigLam        =   0.017*meanD*prod(DimUpRank);
    parBSM.maxIter    =   30;
    par.deNoisingIter =   3;
    
end

parBSM.maxIter    =   parBSM.maxIter-2;
end
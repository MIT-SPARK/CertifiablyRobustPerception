function SDP = rmlindep_sdpt3(SDP)
%% remove linearly dependent constraints
[Atnew,bnew,~,~,~,~,~] = ...
              checkdepconstr(SDP.blk,SDP.At,SDP.b,zeros(length(SDP.b),1),1);        
SDP.At  = Atnew;
SDP.b   = bnew;
end
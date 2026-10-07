function tf = ERW_Should_Apply_TC_Daily_Delta(j,i,NN,Datam)
% Apply each retained T&C daily state change once. The final ERW flush can
% revisit the same j to integrate the remaining sub-daily forcing.

tf = j >= 2 && (i ~= NN || Datam(i,4) == 1);
end

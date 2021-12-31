function R=avalcrs(M,c,r,S,numc,numr,nums)
    R=replace(M,c,numc);
    R=replace(R,r,numr);
    R=replace(R,S,nums);
end
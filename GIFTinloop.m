function result = GIFTinloop(comi,proj)
    global Drug2Sub Protein2Domain Sub2Domain_Recover
    subs = Drug2Sub(comi,:);
    doms = Protein2Domain(proj,:);
    result = 1 - exp(sum(sum(log(1 - Sub2Domain_Recover(subs==1,doms==1)))));
end

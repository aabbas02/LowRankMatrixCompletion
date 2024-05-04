function [tAltMinParfor, tAltMinCntrl, tAltMinPrvt, tAltGD, tAltGDFed, tAltGDMineta1, tAltGDMinCntrl]  = ... 
 findTimes(timeAltMinCntrl, SDAltMinCntrl, timeAltMinParfor, SDAltMinParfor, ...
          timeAltMinPrvt, SDAltMinPrvt, timeAltGD, SDAltGD, timeAltGDFed, SDAltGDFed, ...
          timeAltGDMineta1, SDAltGDMineta1, timeAltGDMinCntrl, SDAltGDMinCntrl)


    lim = 10^-11;
    idx1 = find(SDAltMinParfor < lim,1);
    tAltMinParfor = timeAltMinParfor(idx1);

    idx2 = find(SDAltMinCntrl < lim,1);
    tAltMinCntrl = timeAltMinCntrl(idx2);

    idx3 = find(SDAltMinPrvt < lim,1);
    tAltMinPrvt = timeAltMinPrvt(idx3);

    idx4 =  find(SDAltGD < lim,1);
    tAltGD = timeAltGD(idx4);

    idx5 = find(SDAltGDFed  < lim,1);
    tAltGDFed = timeAltGDFed(idx5);

    idx7 = find(SDAltGDMineta1 < lim,1);
    tAltGDMineta1 =  timeAltGDMineta1(idx7);

    idx8 = find(SDAltGDMinCntrl < lim,1);
    tAltGDMinCntrl =  timeAltGDMinCntrl(idx8);

end

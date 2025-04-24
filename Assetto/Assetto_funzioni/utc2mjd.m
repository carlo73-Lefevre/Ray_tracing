function MJD=utc2mjd(utc)
% n=datenum(utc);
% JD=juliandate(n);
JD = JULIANm(utc(2), utc(3), utc(1))+(utc(4)+(utc(5)+(utc(6)/60))/60)/24;
MJD=JD-2400000.5;
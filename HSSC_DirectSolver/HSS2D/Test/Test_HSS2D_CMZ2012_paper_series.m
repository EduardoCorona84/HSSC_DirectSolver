function Test_HSS2D_CMZ2012_paper_series(fname)

Ng = 7; 
np = [4 8 16 32];% 16];
lnp = length(np); 
EDvC = zeros(lnp,6); ED = EDvC; EC = ED; Ttr = ED; TDI = ED; TCI = ED; 
TDA = ED; TCA = ED; MBD = ED; MBC = ED; 

for i=1:lnp
    [EDvC(i,:),ED(i,:),EC(i,:),Ttr(i,:),TDI(i,:),TCI(i,:),TDA(i,:),TCA(i,:),MBD(i,:),MBC(i,:)] = Test_HSS2D_CMZ2012_paper(Ng,np(i));
    save(fname,'EDvC','ED','EC','Ttr','TDI','TCI','TDA','TCA','MBD','MBC'); 
end

end
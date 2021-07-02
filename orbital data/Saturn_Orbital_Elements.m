function [] = Saturn_Orbital_Elements(inputfile, outputfile)
% read csv or text file with column vector
myt=csvread(inputfile);
[ecco,inco,w_po,capomo,L_po,obl,wo,fo,Mo,ao,Lso,t]=orbelm_swift(myt);
csvwrite(outputfile,[ecco',inco',w_po',capomo',L_po',obl',wo',fo',Mo',ao',Lso',t']);
end


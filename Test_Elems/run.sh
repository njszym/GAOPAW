
for elem in H Li Be B C O F Na Mg Al Si S Cl K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I
do                                                                                                                                     
mkdir $elem                                                                                                                            
cp dakota_input/$elem.dakota.in $elem/dakota.in                                                                                                                     
cp gaopaw_input/$elem.gaopaw.yaml $elem/gaopaw.yaml                                                                                     
cd $elem                                                                                                                               
dakota dakota.in
cd ../                                                                                                                                 
done             


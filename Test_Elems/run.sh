
for elem in H Li Be B C O F Na Mg Al Si S Cl                                                                                                                   
do                                                                                                                                     
mkdir $elem                                                                                                                            
cp dakota_input/$elem.dakota.in $elem/dakota.in                                                                                                                     
cp gaopaw_input/$elem.gaopaw.yaml $elem/gaopaw.yaml                                                                                     
cd $elem                                                                                                                               
dakota dakota.in
cd ../                                                                                                                                 
done             


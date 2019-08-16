
for folder in Cs Ba La Hf Ta W Re Os Ir Pt Au Tl Pb Bi
do
cd $folder
/nfs/working/releases/suite2019-2/run /home/szymansk/scripts/test_PP.py
rm -r *save* *xml */*save*
cd ../
done

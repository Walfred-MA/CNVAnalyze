python summary.py 
python summaryaddinfo.py

cat allpvalues.txt  | grep -v 'Ref' | grep -v 'Dup' | grep -v 'Novel'  > allpvalues_alt.txt
cat allpvalues.txt  | grep 'Dup\|Novel'  > allpvalues_dup.txt 

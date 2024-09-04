cat results/*_anova.txt > all_anova.txt
cat results/*_snpanova.txt > all_snpanova.txt
cat results/*_combineanova.txt > all_combineanova.txt
python summary.py

#!/bin/bash

python collectall.py
python organize.py
python summaryaddinfo.py

cat expr_tissuesummary.txt | awk '{if($5 != $12) print $0}' > expr_tissuesummary.txt_diff.txt


cat expr_tissuesummary.txt_diff.txt  | grep -v 'Ref' | grep -v 'Dup' | grep -v 'Novel'  > expr_tissuesummary.txt_diff.txt_alt.txt
cat expr_tissuesummary.txt_diff.txt  | grep 'Dup\|Novel'  > expr_tissuesummary.txt_diff.txt_dup.txt

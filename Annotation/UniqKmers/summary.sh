ls results/* -1 | xargs -I [] bash -c 'pref=$(echo [] | rev | cut -d/ -f1 | rev | cut -d_ -f1,2 );  awk -v prf="$pref" '\''{print prf, $0}'\'' [] '

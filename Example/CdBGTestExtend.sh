for e in 17 20 23 26 29 32
do
  ./StLiter -M E -r dataset1 -m 10 -L $e -l 3 -p 3 -t 4 -a 0 -i 1
  ./StLiter -M E -r dataset1 -m 10 -L $e -l 3 -p 3 -t 4 -a 0 -i 0
done

for e in 35 45 55 65 75 85 
do
  ./StLiter -M E -r dataset1 -m 10 -L $e -l 10 -p 3 -t 4 -a 0 -i 1
  ./StLiter -M E -r dataset1 -m 10 -L $e -l 10 -p 3 -t 4 -a 0 -i 0
done

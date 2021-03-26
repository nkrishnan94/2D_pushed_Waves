
#!/bin/bash
# A sample Bash script, by Ryan
cnt=0
tot=400



for k in 50000 100000 500000 1000000  ### Inner for loop ###
do
   for b in 0 1 2 3 4 5 6 10
   do
      for (( i = 0; i <= 25; i++ )) ### Outer for loop ###
      do
      	g++ eff_pop_2d_bounds_full.cpp -lgsl
      	./a.out -K $k -B $b 
      	cnt=$[cnt+1]
      	echo $cnt "out of" $tot "\n"
      done 

   done

done


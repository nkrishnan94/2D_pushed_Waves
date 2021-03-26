
#!/bin/bash
# A sample Bash script, by Ryan
cnt=0
tot=450



for j in 1000 5000 10000  ### Inner for loop ###
do
   for b in 0 4 8
   do
      for (( i = 0; i <= 25; i++ )) ### Outer for loop ###
      do
      	g++ eff_pop_2d_bounds.cpp -lgsl
      	./a.out -K $j -B $b 
      	cnt=$[cnt+1]
      	echo $cnt "out of" $tot "\n"
      done 	
   done

done


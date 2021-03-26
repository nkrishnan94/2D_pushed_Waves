
#!/bin/bash
# A sample Bash script, by Ryan
cnt=0
tot=75



for k in 10000  ### Inner for loop ###
do
   for b in  0 3.5 6
   do
 
      for (( i = 0; i < 50; i++ )) ### Outer for loop ###
      do
      	g++ -o inv_out eff_pop_2d_bounds_inv.cpp -lgsl
      	./inv_out -K $k -B $b
      	cnt=$[cnt+1]
      	echo $cnt "out of" $tot "\n"

      done 

   done

done


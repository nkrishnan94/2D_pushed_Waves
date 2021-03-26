
#!/bin/bash
# A sample Bash script, by Ryan
cnt=0
tot=75



for k in 1000000 500000 50000 10000 ### Inner for loop ###
do
   for b in 0 2.5 3.5 6 10
   do

      for (( i = 0; i < 20; i++ )) ### Outer for loop ###
      do
      	g++ -o out_kpz eff_pop_2d_fisher_init_long.cpp -lgsl
      	./out_kpz -K $k -B $b 
      	cnt=$[cnt+1]
      	echo $cnt "out of" $tot "\n"

      done 

   done

done



#!/bin/bash
# A sample Bash script, by Ryan
cnt=0
tot=75



for k in 20000 50000 ### Inner for loop ###
do
   for b in 3.5 6 10
   do
      for (( i = 0; i < 100; i++ )) ### Outer for loop ###
      do
      	g++ -o out_front eff_pop_2d_bounds.cpp -lgsl
      	./out_front -K $k -B $b 
      	cnt=$[cnt+1]
      	echo $cnt "out of" $tot "\n"
      done
      done 

   done

done



#!/bin/bash
# A sample Bash script, by Ryan
cnt=0
tot=75


for (( i = 0; i < 10; i++ )) ### Outer for loop ###
do
   for k in 1000000 ### Inner for loop ###
   do
      for b in 0 2.5 3.5 6 10
      do
         	g++ eff_pop_2d_fisher_init.cpp -lgsl
         	./a.out -K $k -B $b 
         	cnt=$[cnt+1]
         	echo $cnt "out of" $tot "\n"

      done 

   done

done


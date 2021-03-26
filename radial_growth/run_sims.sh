
#!/bin/bash
# A sample Bash script, by Ryan
cnt=0
tot=450


for b in 8
do
   for u in 2 5 10 20 50 80 100 
   do
      for (( i = 0; i <50; i++ )) ### Outer for loop ###
      do
      	g++ eff_pop_2d_radial.cpp -lgsl
      	./a.out -B $b -U $u 
      	cnt=$[cnt+1]
      	echo $cnt "out of" $tot "\n"
      done 	
   done
done




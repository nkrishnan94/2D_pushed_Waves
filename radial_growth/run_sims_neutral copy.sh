
#!/bin/bash
# A sample Bash script, by Ryan
cnt=0
tot=450




for b in 0 2.5 4.5
do
   for u in 250
   do
      for (( i = 0; i < 25 ; i++ )) ### Outer for loop ###
      do
      	g++ -o outfs eff_pop_2d_radial_neutral_smallr.cpp -lgsl
      	./outfs -B $b -U $u
      	cnt=$[cnt+1]
      	echo $cnt "out of" $tot "\n"
      done 	
   done
done




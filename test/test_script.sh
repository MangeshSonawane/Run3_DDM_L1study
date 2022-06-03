i=0
for line in $(cat HTo2LongLivedTo4mu_MH-125_MFF-25_CTau-1500mm.list)
do
	./RunIIIreMCRecipe.sh $line $i
	((i=i+1))
##	if [[ $i -ge 6 ]]; then
#		break
#	fi
done

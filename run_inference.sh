#!/bin/bash
TREECOUNT=32
DO_ASTRAL=true
DO_PAUP=true
TARGETTREES=./QuartetMethods/example/trees_small.txt
FACTORS=(0.5 1.0 2.0 4.0 8.0)

for f in ${FACTORS[@]}; do
    DATASET=./QuartetMethods/example/simulated_data_small-$f
    TREEOUTPUT=./QuartetMethods/inference_outputs-$f/
    ASTRAL_SCOREOUTPUT=$TREEOUTPUT/ASTRAL/allscores.txt # This is for ASTRAL, TODO: change the name so that it reflects this
    PAUP_SCOREOUTPUT=$TREEOUTPUT/PAUP/allscores.txt # This is for ASTRAL, TODO: change the name so that it reflects this
    # initialise tree output space
    mkdir $TREEOUTPUT
    if [[ $DO_ASTRAL ]]; then
        mkdir $TREEOUTPUT/ASTRAL
        mkdir $TREEOUTPUT/ASTRAL/logs
        mkdir $TREEOUTPUT/ASTRAL/trees
        >$ASTRAL_SCOREOUTPUT # set up score output
    fi
    if [[ $DO_ASTRAL ]]; then
        mkdir $TREEOUTPUT/PAUP
        mkdir $TREEOUTPUT/PAUP/trees
        mkdir $TREEOUTPUT/PAUP/scores
        >$PAUP_SCOREOUTPUT # set up score output
    fi
    touch quartet_temp.txt
    touch current_tree.txt
    touch nexus_temp.txt
    SETTINGS=(low mod modhigh high veryhigh)
    for setting in ${SETTINGS[@]}; do
        CSVS=$DATASET/$setting'_noborrowing/no-morph'
        for ((i=1;i<=$TREECOUNT;i++)); do # tree number
            for ((r=1;r<=4;r++)); do # rep number
                pattern="sim_tree$i""_"$r".csv"
                head -"$i" $TARGETTREES | tail -1 > current_tree.txt # write current gold tree to current_tree.txt
                for FILE in $CSVS/*; do
                    if [[ $FILE =~ $pattern ]]; then
                        id=$setting'_'$i'_'$r
                        echo "Factor: $f; ID = $id: target is tree $i"
                        # generate quartets
                        python -c "from QuartetMethods.scripts.getQuartets import *; print_quartets('$FILE')" > quartet_temp.txt
                        echo "✅ Quartet generation" 
                        if $DO_PAUP; then 
                            >nexus_temp.nex
                            Rscript ./QuartetMethods/scripts/commandLineNex.R -f $DATASET/$setting'_noborrowing/sim_tree'$i'_'$r.csv -o nexus_temp.nex -p 3 -m 50.0 > /dev/null 2> /dev/null
                            ./QuartetMethods/scripts/paup -n nexus_temp.nex > tmp_out.txt 2> tmp_run.txt
                            mv paup_out.trees $TREEOUTPUT/PAUP/trees/$id.trees
                            mv paup_out.scores $TREEOUTPUT/PAUP/scores/$id.scores
                            echo "✅ MP4 Tree inference" 

                            echo $FILE >> $PAUP_SCOREOUTPUT
                            Rscript ./QuartetMethods/scripts/QuartetScorer.R -f nexus -r $(<current_tree.txt) -m 1 -p 0 -i $TREEOUTPUT/PAUP/trees/$id.trees >> $PAUP_SCOREOUTPUT
                            echo "✅ PAUP* Tree scoring" 
                        fi
                        if $DO_ASTRAL; then # run ASTRAL 
                            java -jar ./QuartetMethods/ASTRAL/astral.5.7.8.jar -i quartet_temp.txt -o $TREEOUTPUT/ASTRAL/trees/$id.tre -x > /dev/null 2> $TREEOUTPUT/ASTRAL/logs/$id.log # Run ASTRAL in exact mode
                            echo "✅ ASTRAL Tree inference" 

                            echo $FILE >> $ASTRAL_SCOREOUTPUT
                            Rscript ./QuartetMethods/scripts/QuartetScorer.R -f newick -r $(<current_tree.txt) -m 1 -p 0 -i $TREEOUTPUT/ASTRAL/trees/$id.tre >> $ASTRAL_SCOREOUTPUT
                            echo "✅ ASTRAL Tree scoring" 
                        fi
                    fi
                done
            done
        done
    done
    rm quartet_temp.txt
    rm current_tree.txt
    rm nexus_temp.txt
done

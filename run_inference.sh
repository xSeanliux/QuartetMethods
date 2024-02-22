#!/bin/bash
TREECOUNT=32
OS_TYPE='RedHat' # RedHat / OSX
DO_ASTRAL=true
DO_MP4=true
USE_MORPH=true
TARGETTREES=./QuartetMethods/example/trees_small.txt
FACTORS=(0.5 1.0 2.0 4.0 8.0)

if [[ $OS_TYPE = "RedHat" ]]; then 
    PAUP_PATH=./QuartetMethods/scripts/paup4a168_centos64
    chmod a+x $PAUP_PATH
elif [[ $OS_TYPE = "OSX" ]]; then 
    PAUP_PATH=./QuartetMethods/scripts/paup
else 
    echo "PAUP_PATH could not be set (OSTYPE ="$OSTYPE" may be invalid )"
fi

for f in ${FACTORS[@]}; do
    DATASET=./QuartetMethods/example/simulated_data_small-$f
    if $USE_MORPH; then 
        TREEOUTPUT=./QuartetMethods/outputs_with_morphology/inference_outputs-$f/
    else
        TREEOUTPUT=./QuartetMethods/outputs_no_morphology/inference_outputs-$f/
    fi
    ASTRAL_SCOREOUTPUT=$TREEOUTPUT/ASTRAL/allscores.txt # This is for ASTRAL, TODO: change the name so that it reflects this
    MP4_SCOREOUTPUT=$TREEOUTPUT/MP4/allscores.txt # This is for ASTRAL, TODO: change the name so that it reflects this
    # initialise tree output space
    mkdir -p $TREEOUTPUT
    if [[ $DO_ASTRAL ]]; then
        mkdir $TREEOUTPUT/ASTRAL
        mkdir -p $TREEOUTPUT/ASTRAL/logs
        mkdir -p $TREEOUTPUT/ASTRAL/trees
        >$ASTRAL_SCOREOUTPUT # set up score output
    fi
    if [[ $DO_MP4 ]]; then
        mkdir $TREEOUTPUT/MP4
        mkdir -p $TREEOUTPUT/MP4/trees
        mkdir -p $TREEOUTPUT/MP4/scores
        >$MP4_SCOREOUTPUT # set up score output
    fi
    touch quartet_temp.txt
    touch current_tree.txt
    touch nexus_temp.txt
    SETTINGS=(low mod modhigh high veryhigh)
    for setting in ${SETTINGS[@]}; do
        CSVS=$DATASET/$setting'_noborrowing/no-morph'
        if $USE_MORPH; then 
            CSVS=$DATASET/$setting'_noborrowing/'
        else
            CSVS=$DATASET/$setting'_noborrowing/no-morph'
        fi
        for ((i=1;i<=$TREECOUNT;i++)); do # tree number
            for ((r=1;r<=4;r++)); do # rep number
                pattern="sim_tree$i""_"$r".csv"
                head -"$i" $TARGETTREES | tail -1 > current_tree.txt # write current gold tree to current_tree.txt
                for FILE in $CSVS/*; do
                    if [[ $FILE =~ $pattern ]]; then
                        id=$setting'_'$i'_'$r
                        echo "Factor: $f; ID = $id: target is tree $i"
                        # generate quartets
                        
                        if $DO_MP4; then 
                            >nexus_temp.nex
                            Rscript ./QuartetMethods/scripts/commandLineNex.R -f $CSVS/sim_tree$i'_'$r.csv -o nexus_temp.nex -p 3 -m 1.0 > /dev/null 2> /dev/null
                            echo "✅ MP4 nexus files"
                            $PAUP_PATH -n nexus_temp.nex > tmp_out.txt 2> tmp_run.txt
                            mv paup_out.trees $TREEOUTPUT/MP4/trees/$id.trees
                            mv paup_out.scores $TREEOUTPUT/MP4/scores/$id.scores
                            echo "✅ MP4 tree inference" 

                            echo $FILE >> $MP4_SCOREOUTPUT
                            Rscript ./QuartetMethods/scripts/QuartetScorer.R -f nexus -r $(<current_tree.txt) -m 1 -p 0 -i $TREEOUTPUT/MP4/trees/$id.trees >> $MP4_SCOREOUTPUT
                            echo "✅ MP4 tree scoring" 
                        fi
                        if $DO_ASTRAL; then # run ASTRAL 
                            python -c "from QuartetMethods.scripts.getQuartets import *; print_quartets('$FILE')" > quartet_temp.txt
                            echo "✅ ASTRAL quartet generation" 
                            java -jar ./QuartetMethods/ASTRAL/astral.5.7.8.jar -i quartet_temp.txt -o $TREEOUTPUT/ASTRAL/trees/$id.tre -x > /dev/null 2> $TREEOUTPUT/ASTRAL/logs/$id.log # Run ASTRAL in exact mode
                            echo "✅ ASTRAL tree inference" 

                            echo $FILE >> $ASTRAL_SCOREOUTPUT
                            Rscript ./QuartetMethods/scripts/QuartetScorer.R -f newick -r $(<current_tree.txt) -m 1 -p 0 -i $TREEOUTPUT/ASTRAL/trees/$id.tre >> $ASTRAL_SCOREOUTPUT
                            echo "✅ ASTRAL tree scoring" 
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

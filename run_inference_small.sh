#!/bin/bash
shopt -s expand_aliases # expand alias so that mb works
TREECOUNT=32
REPLICA_COUNT=4  # a number from 1 to 4
OS_TYPE='RedHat' # RedHat / OSX
DO_ASTRAL=false  # a(stral)
DO_MP4=false     # p(arsimony)
DO_GA=false      # g(ray & atkinson)
USE_MORPH=false  # m(orphology)

# ASTRAL modes
QT_MODE=1
BP_MODE=1

TARGETTREES=./QuartetMethods/example/trees_small.txt
FACTORS=(0.5 1.0 2.0 4.0 8.0) # f
SETTINGS=(low mod modhigh high veryhigh) # s 
RUNID=$(tr -dc A-Za-z0-9 </dev/urandom | head -c 13; echo) # random string so that runs don't use the same name (i.e. for temp files)
MB_EXEC=/projects/tallis/zxliu2/bin/bin/mb
while getopts 'apgq:b:ms:f:' o; do 
    echo $o' '$OPTARG
    case $o in 
        a) DO_ASTRAL=true;;
        p) DO_MP4=true;;
        g) DO_GA=true;;
        q) QT_MODE=$OPTARG;;
        b) BP_MODE=$OPTARG;;
        m) USE_MORPH=true;;
        f) 
        if [[ ${FACTORS[@]} =~ $OPTARG ]]; then  
            echo "FOUND FACTOR "$OPTARG 
            FACTORS=($OPTARG)
        else 
            echo "FACTOR NOT FOUND. USING ALL"
        fi;;
        s)
        if [[ ${SETTINGS[@]} =~ $OPTARG ]]; then  
            echo "FOUND SETTING "$OPTARG 
            SETTINGS=($OPTARG)
        else 
            echo "SETTING NOT FOUND. USING ALL"
        fi;;
        *) echo "Unknown argument: "$o
    esac
done 

echo "ASTRAL: $DO_ASTRAL"
if $DO_ASTRAL; then 
    echo "QT MODE: $QT_MODE, BP MODE: $BP_MODE"
fi
echo "MP4: $DO_MP4"
echo "GA: $DO_GA"
echo "Settings:"${SETTINGS[@]}
echo "Factors:"${FACTORS[@]}

if [[ $OS_TYPE = "RedHat" ]]; then 
    PAUP_PATH=./QuartetMethods/scripts/paup4a168_centos64
    chmod a+x $PAUP_PATH
elif [[ $OS_TYPE = "OSX" ]]; then 
    PAUP_PATH=./QuartetMethods/scripts/paup
else 
    echo "PAUP_PATH could not be set (OSTYPE ="$OSTYPE" may be invalid )"
fi

for f in ${FACTORS[@]}; do
    echo "RUNNING "$f
    DATASET=./QuartetMethods/example/simulated_data_small-$f
    for setting in ${SETTINGS[@]}; do
        if $USE_MORPH; then 
            TREEOUTPUT=./QuartetMethods/outputs_small/outputs_with_morph/inference_outputs-$f/$setting'_noborrowing'
            # TREEOUTPUT=./QuartetMethods/outputs_with_morph/inference_outputs-$f/
        else
            TREEOUTPUT=./QuartetMethods/outputs_small/outputs_no_morph/inference_outputs-$f/$setting'_noborrowing'
            # TREEOUTPUT=./QuartetMethods/outputs_no_morph/inference_outputs-$f/
        fi
        ASTRAL_VARIANT=ASTRAL\($QT_MODE,$BP_MODE\)
        ASTRAL_SCOREOUTPUT=$TREEOUTPUT/$ASTRAL_VARIANT/allscores.txt # This is for ASTRAL, TODO: change the name so that it reflects this
        MP4_SCOREOUTPUT=$TREEOUTPUT/MP4/allscores.txt # This is for ASTRAL, TODO: change the name so that it reflects this
        GA_SCOREOUTPUT=$TREEOUTPUT/GA/allscores.txt # This is for ASTRAL, TODO: change the name so that it reflects this
        # initialise tree output space
        mkdir -p $TREEOUTPUT
        if $DO_ASTRAL; then
            mkdir -p $TREEOUTPUT/$ASTRAL_VARIANT/logs
            mkdir -p $TREEOUTPUT/$ASTRAL_VARIANT/trees
            >$ASTRAL_SCOREOUTPUT # set up score output
        fi
        if $DO_MP4; then
            mkdir -p $TREEOUTPUT/MP4/trees
            mkdir -p $TREEOUTPUT/MP4/scores
            >$MP4_SCOREOUTPUT # set up score output
        fi
        if $DO_GA; then
            mkdir -p $TREEOUTPUT/GA/trees
            mkdir -p $TREEOUTPUT/GA/trees1
            mkdir -p $TREEOUTPUT/GA/scores
            >$GA_SCOREOUTPUT # set up score output
        fi
        touch scratch/tmp_quartet_$RUNID.txt
            CSVS=""
        if $USE_MORPH; then 
            CSVS=$DATASET/$setting'_noborrowing/'
        else
            CSVS=$DATASET/$setting'_noborrowing/no-morph'
        fi
        for ((i=1;i<=$TREECOUNT;i++)); do # tree number
            for ((r=1;r<=$REPLICA_COUNT;r++)); do # rep number
                pattern="sim_tree$i""_"$r".csv"
                CURRENT_TREE=`head -"$i" $TARGETTREES | tail -1`
                for FILE in $CSVS/*; do
                    if [[ $FILE =~ $pattern ]]; then
                        id=$setting'_'$i'_'$r
                        echo "Factor: $f; ID = $id: target is tree $i"
                        # generate quartets
                        if $DO_GA; then 
                            if ! test -f $TREEOUTPUT/GA/trees1/$id.trees; then 
                                >scratch/tmp_mb_$RUNID.nex
                                Rscript ./QuartetMethods/scripts/commandLineNex.R -H $RUNID -f $FILE -o scratch/tmp_mb_$RUNID.nex --resolve-poly 4 --morph-weight 1.0  > /dev/null 2> /dev/null
                                echo "✅ GA nexus files"
                                $MB_EXEC scratch/tmp_mb_$RUNID.nex > /dev/null 2> /dev/null # > tmp_mb_out_$RUNID.txt 2> tmp_mb_log_$RUNID.txt
                                echo "✅ GA sampling"
                                mv Bayes_out_$RUNID.t $TREEOUTPUT/GA/trees1/$id.trees 
                                mv Bayes_out_$RUNID.con.tre $TREEOUTPUT/GA/trees/$id.trees
                                rm Bayes_out_$RUNID.* # tmp_mb*
                            fi
                        fi
                        if $DO_MP4; then 
                            if ! grep -q $FILE $MP4_SCOREOUTPUT; then # You can run this multiple times to continue where you left off
                                >scratch/tmp_mp4_$RUNID.nex
                                Rscript ./QuartetMethods/scripts/commandLineNex.R -H $RUNID -f $FILE -o scratch/tmp_mp4_$RUNID.nex -p 3 -m 1.0 # > /dev/null 2> /dev/null
                                echo "✅ MP4 nexus files"
                                $PAUP_PATH -n scratch/tmp_mp4_$RUNID.nex > /dev/null 2> /dev/null
                                mv scratch/paup_out_$RUNID.trees $TREEOUTPUT/MP4/trees/$id.trees # If we run lots of instances of this script in parallel, paup_out might be overwritten so we can't have that
                                mv scratch/paup_out_$RUNID.scores $TREEOUTPUT/MP4/scores/$id.scores
                                echo "✅ MP4 tree inference" 

                                echo $FILE >> $MP4_SCOREOUTPUT
                                Rscript ./QuartetMethods/scripts/QuartetScorer.R -f nexus -r $CURRENT_TREE -m 1 -p 0 -i $TREEOUTPUT/MP4/trees/$id.trees >> $MP4_SCOREOUTPUT
                                echo "✅ MP4 tree scoring" 
                            else 
                                echo "Skipping "$id
                            fi
                        fi
                        if $DO_ASTRAL; then # run ASTRAL 
                            if ! grep -q $FILE $ASTRAL_SCOREOUTPUT; then
                                python -c "from QuartetMethods.scripts.getQuartets import *; print_quartets('$FILE', mode=$QT_MODE)" > scratch/tmp_quartet_$RUNID.txt
                                echo "✅ ASTRAL quartet generation, $(wc -l scratch/tmp_quartet_$RUNID.txt | awk '{ print $1 }') quartets" 

                                python -c "from QuartetMethods.scripts.getBipartitions import *; get_bipartitions(morph='$USE_MORPH', f='$f', poly='$setting', treeidx=$i, replica=$r, mode=$BP_MODE)" > scratch/tmp_bipartitions_$RUNID.bootstrap.trees
                                echo "✅ Heuristic ASTRAL Get bipartitions" 
                                java -jar ./QuartetMethods/ASTRAL/astral.5.7.8.jar -i scratch/tmp_quartet_$RUNID.txt -o $TREEOUTPUT/$ASTRAL_VARIANT/trees/$id.tre -x > /dev/null 2> $TREEOUTPUT/$ASTRAL_VARIANT/logs/$id.log # Run ASTRAL in default mode
                                echo "✅ Heuristic ASTRAL tree inference" 

                                echo $FILE >> $ASTRAL_SCOREOUTPUT
                                Rscript ./QuartetMethods/scripts/QuartetScorer.R --prune r -f newick -r $CURRENT_TREE -m 1 -p 0 -i $TREEOUTPUT/$ASTRAL_VARIANT/trees/$id.tre >> $ASTRAL_SCOREOUTPUT
                                echo "✅ Heuristic ASTRAL tree scoring" 
                                rm scratch/tmp_quartet_$RUNID.txt
                                rm scratch/tmp_bipartitions_$RUNID.bootstrap.trees
                            else 
                                echo "Heuristic ASTRAL: Skipping "$id
                            fi
                        fi
                    fi
                done
            done
        done
    done
done

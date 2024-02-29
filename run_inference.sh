#!/bin/bash
shopt -s expand_aliases # expand alias so that mb works
TREECOUNT=32
OS_TYPE='RedHat' # RedHat / OSX
DO_ASTRAL=false  # a(stral)
DO_MP4=false     # p(arsimony)
DO_GA=false      # g(ray & atkinson)
DO_HEURISTIC_ASTRAL=false # h(euristic)
USE_MORPH=false  # m(orphology)
TARGETTREES=./QuartetMethods/example/trees_small.txt
FACTORS=(0.5 1.0 2.0 4.0 8.0) # f
SETTINGS=(low mod modhigh high veryhigh) # s 
RUNID=$(tr -dc A-Za-z0-9 </dev/urandom | head -c 13; echo) # random string so that runs don't use the same name (i.e. for temp files)
MB_EXEC=/projects/tallis/zxliu2/bin/bin/mb
while getopts 'apghmf:' o; do 
    echo $o' '$OPTARG
    case $o in 
        a) DO_ASTRAL=true;;
        p) DO_MP4=true;;
        g) DO_GA=true;;
        m) USE_MORPH=true;;
        h) DO_HEURISTIC_ASTRAL=true;;
        f) 
        if [[ ${FACTORS[@]} =~ $OPTARG ]]; then  
            echo "FOUND FACTOR "$OPTARG 
            FACTORS=($OPTARG)
        else 
            echo "FACTOR NOT FOUND. USING ALL"
        fi;;
    esac
done 

echo "ASTRAL_HEURISTIC $DO_HEURISTIC_ASTRAL"
echo "ASTRAL: $DO_ASTRAL"
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
    if $USE_MORPH; then 
        TREEOUTPUT=./QuartetMethods/outputs_with_morph/inference_outputs-$f/
    else
        TREEOUTPUT=./QuartetMethods/outputs_no_morph/inference_outputs-$f/
    fi
    ASTRAL_SCOREOUTPUT=$TREEOUTPUT/ASTRAL/allscores.txt # This is for ASTRAL, TODO: change the name so that it reflects this
    ASTRAL_H_SCOREOUTPUT=$TREEOUTPUT/ASTRAL-H/allscores.txt # This is for ASTRAL, TODO: change the name so that it reflects this
    MP4_SCOREOUTPUT=$TREEOUTPUT/MP4/allscores.txt # This is for ASTRAL, TODO: change the name so that it reflects this
    GA_SCOREOUTPUT=$TREEOUTPUT/GA/allscores.txt # This is for ASTRAL, TODO: change the name so that it reflects this
    # initialise tree output space
    mkdir -p $TREEOUTPUT
    if $DO_HEURISTIC_ASTRAL; then
        mkdir -p $TREEOUTPUT/ASTRAL-H/logs
        mkdir -p $TREEOUTPUT/ASTRAL-H/trees
        >$ASTRAL_H_SCOREOUTPUT # set up score output
    fi
    if $DO_ASTRAL; then
        mkdir -p $TREEOUTPUT/ASTRAL/logs
        mkdir -p $TREEOUTPUT/ASTRAL/trees
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
    touch tmp_quartet_$RUNID.txt
    touch tmp_nexus_$RUNID.txt
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
                CURRENT_TREE=`head -"$i" $TARGETTREES | tail -1`
                for FILE in $CSVS/*; do
                    if [[ $FILE =~ $pattern ]]; then
                        id=$setting'_'$i'_'$r
                        echo "Factor: $f; ID = $id: target is tree $i"
                        # generate quartets
                        if $DO_GA; then 
                            >tmp_mb_$RUNID.nex
                            Rscript ./QuartetMethods/scripts/commandLineNex.R -h $RUNID -f $FILE -o tmp_mb_$RUNID.nex --resolve-poly 4 --morph-weight 1.0 > /dev/null 2> /dev/null
                            echo "✅ GA nexus files"
                            $MB_EXEC tmp_mb_$RUNID.nex > tmp_mb_out_$RUNID.txt 2> tmp_mb_log_$RUNID.txt
                            echo "✅ GA sampling"
                            mv Bayes_out_$RUNID.con.tre $TREEOUTPUT/GA/trees/$id.trees
                            mv Bayes_out_$RUNID.t $TREEOUTPUT/GA/trees1/$id.trees 
                            rm Bayes_out_$RUNID.* mb_tmp*
                        fi
                        if $DO_MP4; then 
                            >tmp_mp4_$RUNID.nex
                            Rscript ./QuartetMethods/scripts/commandLineNex.R -f $FILE -o tmp_mp4_$RUNID.nex -p 3 -m 1.0 > /dev/null 2> /dev/null
                            echo "✅ MP4 nexus files"
                            $PAUP_PATH -n tmp_mp4_$RUNID.nex > /dev/null 2> /dev/null
                            mv paup_out.trees $TREEOUTPUT/MP4/trees/$id.trees # If we run lots of instances of this script in parallel, paup_out might be overwritten so we can't have that
                            mv paup_out.scores $TREEOUTPUT/MP4/scores/$id.scores
                            echo "✅ MP4 tree inference" 

                            echo $FILE >> $MP4_SCOREOUTPUT
                            Rscript ./QuartetMethods/scripts/QuartetScorer.R -f nexus -r $CURRENT_TREE -m 1 -p 0 -i $TREEOUTPUT/MP4/trees/$id.trees >> $MP4_SCOREOUTPUT
                            echo "✅ MP4 tree scoring" 
                            rm mp4_nexus_temp.nex
                        fi
                        if $DO_ASTRAL || $DO_HEURISTIC_ASTRAL; then # run ASTRAL 
                            python -c "from QuartetMethods.scripts.getQuartets import *; print_quartets('$FILE')" > tmp_quartet_$RUNID.txt
                            echo "✅ ASTRAL quartet generation" 

                            if $DO_ASTRAL; then 
                                java -jar ./QuartetMethods/ASTRAL/astral.5.7.8.jar -i tmp_quartet_$RUNID.txt -o $TREEOUTPUT/ASTRAL/trees/$id.tre -x > /dev/null 2> $TREEOUTPUT/ASTRAL/logs/$id.log # Run ASTRAL in exact mode
                                echo "✅ ASTRAL tree inference" 

                                echo $FILE >> $ASTRAL_SCOREOUTPUT
                                Rscript ./QuartetMethods/scripts/QuartetScorer.R -f newick -r $CURRENT_TREE -m 1 -p 0 -i $TREEOUTPUT/ASTRAL/trees/$id.tre >> $ASTRAL_SCOREOUTPUT
                                echo "✅ ASTRAL tree scoring" 
                                rm tmp_quartet_$RUNID.txt
                            fi

                            if $DO_HEURISTIC_ASTRAL; then 
                                if ! [[ -f $TREEOUTPUT/GA/trees1/$id.trees ]]; then 
                                    echo "❌ Heuristic ASTRAL NEXUS file not found: "$TREEOUTPUT/GA/trees1/$id.trees
                                    exit 0
                                fi
                                python -c "from QuartetMethods.scripts.getQuartets import *; nexus_to_newick('$TREEOUTPUT/GA/trees1/$id.trees')" > tmp_bipartitions_$RUNID.bootstrap.trees
                                echo "✅ Heuristic ASTRAL NEXUS to Newick" 
                                java -jar ./QuartetMethods/ASTRAL/astral.5.7.8.jar -i tmp_quartet_$RUNID.txt -o $TREEOUTPUT/ASTRAL-H/trees/$id.tre -e tmp_bipartitions_$RUNID.bootstrap.trees > /dev/null 2> $TREEOUTPUT/ASTRAL-H/logs/$id.log # Run ASTRAL in exact mode
                                echo "✅ Heuristic ASTRAL tree inference" 

                                echo $FILE >> $ASTRAL_H_SCOREOUTPUT
                                Rscript ./QuartetMethods/scripts/QuartetScorer.R -f newick -r $CURRENT_TREE -m 1 -p 0 -i $TREEOUTPUT/ASTRAL-H/trees/$id.tre >> $ASTRAL_H_SCOREOUTPUT
                                echo "✅ Heuristic ASTRAL tree scoring" 
                                rm tmp_quartet_$RUNID.txt
                                rm tmp_bipartitions_$RUNID.bootstrap.trees
                            fi
                        fi
                    fi
                done
            done
        done
    done
    rm tmp_quartet_$RUNID.txt
    rm tmp_nexus_$RUNID.txt
done

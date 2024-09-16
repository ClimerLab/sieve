ORIG_GML=
OUTPUT_DIR=
RUNTAG=
CLEANUP_FLAG=true

# -----------------------------------------------
# Seperate graph into components
# -----------------------------------------------
# Check for dir, if not found create it using mkdir #
if [ ! -d $OUTPUT_DIR ]; then
	mkdir -p $OUTPUT_DIR
fi

# Remove 'NUM_COMP.txt' if it exists
if [ -f $OUTPUT_DIR"NUM_COMP.txt" ]; then
	rm $OUTPUT_DIR"NUM_COMP.txt"
fi

# Run 'Seperate Components program'
./SepComp $ORIG_GML $OUTPUT_DIR $RUNTAG

# Source values from 'Num_COMP.txt'
if [ -f $OUTPUT_DIR"NUM_COMP.txt" ]; then
	source $OUTPUT_DIR"NUM_COMP.txt"
	rm $OUTPUT_DIR"NUM_COMP.txt"
else
	NUM_COMP=-1
	TOTAL_NODES=0
	TOTAL_EDGES=0
fi

if [ $NUM_COMP -eq -1 ]; then
	echo "Error seperating graph...exiting"
	exit -1
fi

# -----------------------------------------------
# Optimize Sieve value for each component
# -----------------------------------------------
# Loop through all components found by "SepComp"
for ((CUR_COMP=1; CUR_COMP<=$NUM_COMP; CUR_COMP++))
do
  echo "Processing component $CUR_COMP"

  COMP_GML="$OUTPUT_DIR"/$RUNTAG"_comp"$CUR_COMP".gml"
  ./S_MIP $COMP_GML $OUTPUT_DIR $RUNTAG $TOTAL_NODES $CUR_COMP
done

# -----------------------------------------------
# Combine results from each component
# -----------------------------------------------
./CombClustOut $ORIG_GML $OUTPUT_DIR $RUNTAG $NUM_COMP

# Clean up temp files
if [ "$CLEANUP_FLAG" = true ]; then
  if [ -f $OUTPUT_DIR$RUNTAG"_comp0.clust" ]; then
	  rm $OUTPUT_DIR$RUNTAG"_comp0.clust"
  fi
  if [ -f $OUTPUT_DIR$RUNTAG"_comp0.nn" ]; then
	  rm $OUTPUT_DIR$RUNTAG"_comp0.nn"
  fi

  for ((CUR_COMP=1; CUR_COMP<=$NUM_COMP; CUR_COMP++))
  do
    if [ -f $OUTPUT_DIR$RUNTAG"_comp1_S.out" ]; then
	    rm $OUTPUT_DIR$RUNTAG"_comp1_S.out"
    fi
    if [ -f $OUTPUT_DIR$RUNTAG"_comp1.nn" ]; then
	    rm $OUTPUT_DIR$RUNTAG"_comp1.nn"
    fi
    if [ -f $OUTPUT_DIR$RUNTAG"_comp1.gml" ]; then
	    rm $OUTPUT_DIR$RUNTAG"_comp1.gml"
    fi
  done
fi

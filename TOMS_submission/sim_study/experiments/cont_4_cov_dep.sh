cd ../

# save dir
save_dir="./experiments/z4"

# experiment name
exp=cont_4_cov_dep

# number of trials
ntrials=50

dims=( 10 25 50 100 )
n=225

for (( i=0; i<${#dims[@]}; i++ )); do

  now=$(date '+%Y-%m-%d_%H-%M')
  p=${dims[i]}
  out_name=$save_dir/$exp\_ntrials$ntrials\_p$p\_n\_$n\_$now.Rout
  printf "\n\nRunning $out_name...\n"
  SECONDS=0
  R CMD BATCH --no-save --no-restore "--args \
  save_dir='$save_dir' \
  experiment='$exp' \
  p=$p \
  n=$n \
  n_trials=$ntrials" \
  main.R $out_name
  duration=$SECONDS
  echo "Done. $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"

done

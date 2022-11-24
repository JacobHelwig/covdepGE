cd ../../

# save dir
save_dir="./experiments/discrete"

# experiment name
exp=disc_cov_dep

# sample size and signal strength
n1=80
n2=20
lambda=15

# number of trials
ntrials=50

dims=( 51 31 11 )

for p in ${dims[@]}; do

  now=$(date '+%Y-%m-%d_%H:%M')
  out_name=$save_dir/$exp\_ntrials$ntrials\_p$p\_n1\_$n1\_n2\_$n2\_lambda$lambda\_$now.Rout
  printf "\n\nRunning $out_name...\n"
  SECONDS=0
  R CMD BATCH --no-save --no-restore "--args \
  save_dir='$save_dir' \
  experiment='$exp' \
  p=$p \
  n1=$n1 \
  n2=$n2 \
  lambda=$lambda \
  n_trials=$ntrials" \
  main.R $out_name
  duration=$SECONDS
  echo "Done. $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"

done

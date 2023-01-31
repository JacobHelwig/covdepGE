cd ../../

# save dir
save_dir="./experiments/discrete"

# experiment name
exp=disc_cov_free

# dimension
p=11

# number of trials
ntrials=50

# sample sizes and signal strengths
n1s=( 50 50 80 80 )
n2s=( 50 50 20 20 )
lambdas=( 15 3 15 3 )

for (( i=0; i<${#lambdas[@]}; i++ )); do

  now=$(date '+%Y-%m-%d_%H:%M')
  n1=${n1s[i]}
  n2=${n2s[i]}
  lambda=${lambdas[i]}
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
